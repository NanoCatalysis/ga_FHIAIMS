#this will create the python file to run 
#Ga generated BY Rodrigo Espinola 
#created on  2020-08-31 14:13:03.484281 
import atompacking_functions as af 
import subprocess

N= 6
atom = "Au"
path = "pools_au6/Au6_6"
R_min = 2.0
Num_decimals =4
Dist_min= 2
generation =0

N = int (N)
pool_size = af.Pool_size(N)
dist = af.Cluster_size(N)
energies =[]
directories =[]
#gen_path = path + "Gen" + str(generation) +"/"
#preff = atom + str(N)
#gen_dir = create_directory(preff, gen_path , add =0) + "/"
#text =["path =", path]
af.print_wami()
for x in range(pool_size):
	directory= af.create_files(Size = N, Atom = atom, Path = path, r_min = R_min ,r_max = dist ,num_decimals =Num_decimals ,dist_min =Dist_min , dist_max =dist)
	directories.append(directory)

print("\n directories : ", directories)
for x in directories:
	try:
		run_raw = "./" +x +"/shforrunning.sh"
		subprocess.call(run_raw,universal_newlines = True, shell = True)
	except :
		print("Error running")
