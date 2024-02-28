import feature_extract as fe
import feature_functions as ff
import time
import csv
import pandas as pd
import os
import glob

# get all point clouds in the folder

base_folder = "VSENSE"
ply_pattern = f"{base_folder}/**/*.ply"

point_clouds = glob.glob(ply_pattern, recursive=True)
print(point_clouds)
# Print the list of .ply files
# for file in point_clouds:
#     print(file)

# first_row = ["name","l_mean,l_std","l_entropy","a_mean,a_std","a_entropy,b_mean","b_std","b_entropy","curvature_mean","curvature_std","curvature_entropy","curvature_ggd1","curvature_ggd2","curvature_aggd1","curvature_aggd2","curvature_aggd3","curvature_aggd4","curvature_gamma1","curvature_gamma2","anisotropy_mean","anisotropy_std","anisotropy_entropy","anisotropy_ggd1,anisotropy_ggd2","anisotropy_aggd1","anisotropy_aggd2","anisotropy_aggd3","anisotropy_aggd4","anisotropy_gamma1","anisotropy_gamma2","linearity_mean","linearity_std","linearity_entropy","linearity_ggd1","linearity_ggd2","linearity_aggd1","linearity_aggd2","linearity_aggd3","linearity_aggd4","linearity_gamma1","linearity_gamma2","planarity_mean","planarity_std","planarity_entropy","planarity_ggd1","planarity_ggd2","planarity_aggd1","planarity_aggd2","planarity_aggd3","planarity_aggd4","planarity_gamma1","planarity_gamma2","sphericity_mean","sphericity_std","sphericity_entropy","sphericity_ggd1","sphericity_ggd2","sphericity_aggd1","sphericity_aggd2","sphericity_aggd3","sphericity_aggd4","sphericity_gamma1","sphericity_gamma2"]
# ff.create_csv("house.csv",first_row)
# df = pd.read_csv("house.csv")

start = time.time()
i = 0
start_time = time.time()
average_time = 0
for pc in point_clouds:
    i += 1
    if i%10 == 0:
        new_time = time.time()
        average_time += new_time
        average_time /= i
    print("average time :" + str(average_time))
    print(pc)
    features = fe.get_feature_vector(pc)
    ff.print_info(features)
    features.insert(0, pc)
    ff.append_csv("test.csv", features)
df = pd.read_csv("test.csv")
# all_names = df["name"].tolist()
end_time = time.time()-start_time
print("final time is"+ str(end_time))
# print(all_names)