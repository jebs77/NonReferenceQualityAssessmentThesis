import numpy as np
import pandas as pd
from skimage import color
from nss_functions import * 
from pyntcloud import PyntCloud
import os 
import cupy as cp
import pandas as pd
from scipy.stats import entropy
from scipy.special import gamma

def my_scale(vec):
    vec = (vec - cp.mean(vec)) / cp.std(vec, ddof=1)
    return vec

def get_color_nss_param(vec):
    return [estimate_basic_param(vec)]

def get_geometry_nss_param(vec):
    return [estimate_basic_param(vec), estimate_ggd_param(vec),
            estimate_aggd_param(my_scale(vec)), estimate_gamma_param(vec)]

def Entropy(labels):
    # Ensure labels is a CuPy array
    labels_cp = cp.asarray(labels)
    
    # Convert CuPy array to NumPy array before using in Pandas Series
    labels_np = labels_cp.get()
    probs = pd.Series(labels_np).value_counts(bins=2000) / len(labels_np)
    en = entropy(probs)
    return en

def estimate_basic_param(vec):
    result = [cp.mean(vec), cp.std(vec, ddof=1), Entropy(vec)]
    return result   

def estimate_ggd_param(vec):
    gam = cp.arange(0.2, 10 + 0.001, 0.001)
    r_gam = (gamma(1.0 / gam) * gamma(3.0 / gam) / (gamma(2.0 / gam) ** 2))

    sigma_sq = cp.mean(vec ** 2)
    sigma = cp.sqrt(sigma_sq)
    E = cp.mean(cp.abs(vec))
    rho = sigma_sq / E ** 2

    differences = cp.abs(rho - r_gam)
    array_position = cp.argmin(differences)
    gamparam = gam[array_position]
    result = [gamparam, sigma]
    return result

def estimate_aggd_param(vec):
    gam = cp.arange(0.2, 10 + 0.001, 0.001)
    r_gam = ((gamma(2.0 / gam)) ** 2) / (
                gamma(1.0 / gam) * gamma(3.0 / gam))

    left_std = cp.sqrt(cp.mean((vec[vec < 0]) ** 2))
    right_std = cp.sqrt(cp.mean((vec[vec > 0]) ** 2))
    gamma_hat = left_std / right_std
    rhat = (cp.mean(cp.abs(vec))) ** 2 / cp.mean((vec) ** 2)
    rhat_norm = (rhat * (gamma_hat ** 3 + 1) * (gamma_hat + 1)) / (
            (gamma_hat ** 2 + 1) ** 2)

    differences = (r_gam - rhat_norm) ** 2
    array_position = cp.argmin(differences)
    alpha = gam[array_position]
    const = cp.sqrt(gamma(1 / alpha)) / cp.sqrt(gamma(3 / alpha))
    mean_param = (right_std - left_std) * (
                    gamma(2 / alpha) / gamma(1 / alpha)) * const
    result = [alpha, mean_param, left_std, right_std]
    return result

def estimate_gamma_param(vec):
    mean = cp.mean(vec)
    std = cp.std(vec)
    shape = (mean / std) ** 2
    scale = (std ** 2) / mean
    result = [shape, scale]
    return result

def get_feature_vector(objpath):  
  #load colored point cloud
  print("Begin loading point cloud")
  cloud = PyntCloud.from_file(objpath)
  
  #begin geometry projection
  print("Begin geometry feature extraction.")
  k_neighbors = cloud.get_neighbors(k=10)
  ev = cloud.add_scalar_field("eigen_values", k_neighbors=k_neighbors)
  cloud.add_scalar_field("curvature", ev=ev)
  cloud.add_scalar_field("anisotropy",ev=ev)
  cloud.add_scalar_field("linearity",ev=ev)
  cloud.add_scalar_field("planarity",ev=ev)
  cloud.add_scalar_field("sphericity",ev=ev)
  curvature = cloud.points['curvature(11)'].to_numpy()
  anisotropy = cloud.points['anisotropy(11)'].to_numpy()
  linearity = cloud.points['linearity(11)'].to_numpy()
  planarity = cloud.points['planarity(11)'].to_numpy()
  sphericity = cloud.points['sphericity(11)'].to_numpy()


  #begin color projection
  print("Begin color feature extraction.")
  rgb_color = cloud.points[['red','green','blue']].to_numpy()/255
  lab_color = color.rgb2lab(rgb_color)
  l = lab_color[:,0]
  a = lab_color[:,1]
  b = lab_color[:,2]
  
  
  print("Begin NSS parameters estimation.")
  # compute nss parameters
  nss_params = []
  # compute color nss features
  for tmp in [l,a,b]:
      tmp_cupy = cp.array(tmp)
      params = get_color_nss_param(tmp_cupy)
      #flatten the feature vector
      nss_params = nss_params + [i for item in params for i in item]
  # compute geomerty nss features
  for tmp in [curvature,anisotropy,linearity,planarity,sphericity]:
      tmp_cupy = cp.array(tmp)
      params = get_geometry_nss_param(tmp_cupy)
      #flatten the feature vector
      nss_params = nss_params + [i for item in params for i in item]
  return nss_params

