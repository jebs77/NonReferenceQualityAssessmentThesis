import feature_extract as fe
import feature_functions as ff
import time
import csv
import pandas as pd
import os

# get all point clouds in the folder
objpath = "models/Longdress/"
point_clouds = []
for files in os.listdir(objpath):
    if files.endswith(".ply"):
        point_clouds.append("Longdress/"+files)
print(point_clouds)

all_pcs = []
for pc in point_clouds:
    objpath = "models/" + pc
    all_pcs.append(objpath)
print(all_pcs)

first_row = ["name","l_mean,l_std","l_entropy","a_mean,a_std","a_entropy,b_mean","b_std","b_entropy","curvature_mean","curvature_std","curvature_entropy","curvature_ggd1","curvature_ggd2","curvature_aggd1","curvature_aggd2","curvature_aggd3","curvature_aggd4","curvature_gamma1","curvature_gamma2","anisotropy_mean","anisotropy_std","anisotropy_entropy","anisotropy_ggd1,anisotropy_ggd2","anisotropy_aggd1","anisotropy_aggd2","anisotropy_aggd3","anisotropy_aggd4","anisotropy_gamma1","anisotropy_gamma2","linearity_mean","linearity_std","linearity_entropy","linearity_ggd1","linearity_ggd2","linearity_aggd1","linearity_aggd2","linearity_aggd3","linearity_aggd4","linearity_gamma1","linearity_gamma2","planarity_mean","planarity_std","planarity_entropy","planarity_ggd1","planarity_ggd2","planarity_aggd1","planarity_aggd2","planarity_aggd3","planarity_aggd4","planarity_gamma1","planarity_gamma2","sphericity_mean","sphericity_std","sphericity_entropy","sphericity_ggd1","sphericity_ggd2","sphericity_aggd1","sphericity_aggd2","sphericity_aggd3","sphericity_aggd4","sphericity_gamma1","sphericity_gamma2"]
ff.create_csv("house.csv",first_row)
df = pd.read_csv("house.csv")

start = time.time()

for pc in all_pcs:
    features = fe.get_feature_vector(pc)
    ff.print_info(features)
    features.insert(0, pc)
    ff.append_csv("house.csv", features)
df = pd.read_csv("house.csv")
all_names = df["name"].tolist()

print(all_names)