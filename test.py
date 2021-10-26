import os
import numpy as np

dataset_path = "dataset/"
result_path = "results/"
graph_1 = "mz_test"
graph_2 = "mz_test"
sim = "sim"
ini = "edit"
label_constraint="N"
opt= "null"
num_thus=32
w_o = 0.4
w_i = 0.4
w_l = 0.2
theta = 0
w_l=0.2

out_root = "./edit_results/"
command = "/usr/bin/time -v timeout 30h ./calsim "+ dataset_path + " " + result_path + " " + graph_1 + " " + graph_2 + " " + ini + " " + opt + " " + sim + " " + label_constraint + " " + str(w_o) + " " + str(w_i) + " " + str(w_l) + " " + str(theta) + " " + str(num_thus) + " > " + out_root + graph_1 + "_" + ini + "_" + opt + "_" + sim + "_" + label_constraint + "_" + str(w_o) + "_" + str(w_l) +"_theta" + str(theta) + "_" + str(num_thus) + "_out 2>&1"
# print(command)
os.system(command)

# 2>&1 redirects std_err to std_out