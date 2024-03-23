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

# print(point_clouds)


start = time.time()
i = 0
start_time = time.time()
average_time = 0
for pc in point_clouds:
    
    print(pc)
    features = fe.get_feature_vector(pc)
    ff.print_info(features)
    features.insert(0, pc)
    ff.append_csv("test.csv", features)
df = pd.read_csv("test.csv")
