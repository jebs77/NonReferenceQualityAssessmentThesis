import feature_extract as fe
import feature_functions as ff
import time
import csv
import pandas as pd
import os
import glob

# get all point clouds in the folder
# the base folder for the feature extraction
base_folder = "statistical_features/WPC"
ply_pattern = f"{base_folder}/**/*.ply"

point_clouds = glob.glob(ply_pattern, recursive=True)
print(point_clouds)


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
    ff.append_csv("WPC_NSS.csv", features)
df = pd.read_csv("WPC_NSS.csv")
# all_names = df["name"].tolist()
end_time = time.time()-start_time
print("final time is"+ str(end_time))
# print(all_names)