import random 
import math 
import numpy as np
import subprocess 
import shlex
import time 
import os 
import datetime
import re

###### Geometric functions

def distance_1(atom1,atom2):
    dist = math.sqrt(math.pow(atom1[0] - atom2[0], 2) +
                     math.pow(atom1[1] - atom2[1], 2) +
                     math.pow(atom1[2] - atom2[2], 2))
    return dist


#This function generate coordinates using distribution over spheric coordinates and the making them Euclidean

def generate_coordinates(r_min=0,r_max=1):
    theta = random.uniform(0, 2 * math.pi)
    phi = random.uniform(0, 2 * math.pi)
    r = random.uniform(r_min,r_max)
    x = (r * math.cos(theta) * math.sin(phi))
    y = (r * math.sin(theta) * math.sin(phi))
    z = (r * math.cos(phi))
    vector = [x,y,z]
    return vector

#This Function sorts atoms 
def select_atom(atoms =[]):
    at = np.array(atoms)
    shape = np.shape(at)
    dim = shape[1]
    number =shape[0]
    if number >=2 :
        rm = random.randint(0,number-1)
        selection = atoms[rm]
    
        
    else :    
        selection = atoms[0]
    return selection

def add_coordinates(vector1, vector2):
    result =[]
    if vector1 ==[]:
        result =vector2
        
    else:
        for i in range(len(vector2)):
            result = [vector1[i] + vector2[i] for i in range(len(vector2))];
    return result
        
def proof_distance(atom1, atom2, r_min=0,r_max=10):
    distance = distance_1(atom1, atom2)
    condition = None
    if distance <= r_min:
        condition = False
    elif distance >= r_max:
        condition = False
    else:
        condition = True
    return condition       






def generate_atom (atoms=[],r_min = 2.0,r_max = 7,num_decimals =4,dist_min =2, dist_max =7):
    condition =0
    #print("r_max ", r_max, " dist_max", dist_max)
    if atoms == [] :
        at1 =[0,0,0]
    else :
        temporal_atom = None 
        num_atoms = len(atoms)
        while condition ==0:
            random_coordinates =generate_coordinates(r_min,r_max)
            base_atom = select_atom(atoms)
            temporal_atom = add_coordinates(base_atom,random_coordinates)
            at = np.array(atoms)
            number = at.ndim
            if number !=1 :
               
                for atom in atoms:
                    condition =proof_distance(atom, temporal_atom,dist_min ,dist_max)
                    if condition == 0:
                        break
                    else:
                        condition = 1
            else:
                atom = atoms
                condition =proof_distance(atom, temporal_atom,dist_min ,dist_max)
                if condition == 0:
                    break
                else:
                    condition = 1
          
                    
        at1 = temporal_atom;             

                       

    return at1

def print_xyz(size , matrix, atom , path=""):
	
	# print("shape: "+ str(shape[0]))
	Path_xyz = path 
	number= size
	file_name = Path_xyz+"/" + atom + str(number)
	print(Path_xyz)
	print(file_name)
	
	lines_of_text =[]
	for x in matrix:
		temp_string = "atom\t" +str(x[0])+"\t"+ str(x[1])+"\t"+str(x[2])+"\t"+ atom+ "\n"
		lines_of_text.append(temp_string)
	with  open("%s.xyz"%file_name, "w") as fh :
		fh.write(str(number)+ "\n")
		fh.write("\n")
		fh.writelines(lines_of_text)
		fh.close()
	return("Done")	 

def print_matrix_geometryin(matrix =[[],[]],name = "", path= ""):

	#name_folder=path.split('/')
	#original_path= str([string += x for x in name_folder[:-1]])
	at = np.array(matrix)
	shape = np.shape(at)
	number =shape[0]
	create_folder(name="", path=path)
	file_name = path +  "/" + name 
	lines_of_text =[]
	for x in matrix:
		temp_string = "atom   " + str(x[0]) + "   " + str(x[1]) + "   " + str(x[2]) + "   " + str(x[3])+ "\n"
		lines_of_text.append(temp_string)
	with  open(file_name, "w") as fh :
		#fh.write(str(number)+ "\n")
		#fh.write("\n")
		fh.writelines(lines_of_text)
		fh.close()




def print_geometryin(matrix,atom ="Au",path=""):
	file_name = path + "/geometry.in"
	
	lines_of_text=[]
	for x in matrix:
		temp_string = "atom"+"  "+str(x[0]) +"   "+str(x[1])+"   "+str(x[2])+"   " +atom + "\n"
		lines_of_text.append(temp_string)
	with open(file_name, "w") as fh :
		fh.writelines(lines_of_text)
		fh.close()
	return "Done"

def create_cluster(size =55, atom="Au",path ="", R_min = 2.0,R_max = 7,Num_decimals =4,Dist_min =2, Dist_max =7,cores =16):
	cluster =[]
	Size = size
	Atom = atom
	Path_cluster= path
	print("size : ", Size , "atom : ",Atom,"r_max",R_max ,"dist = ", Dist_max )
	for i in range(int(size)):
 		at=generate_atom(cluster,r_min = R_min,r_max = R_max , dist_max = Dist_max )
 		cluster.append(at)

	print_geometryin(cluster,Atom ,Path_cluster)
	print_xyz(size,cluster, Atom, Path_cluster)
	#create_shforrunning(size,Atom,Path_cluster,cores =cores)
	create_shforrunning_name(name= Atom+str(size),path = path, cores=cores)
	create_control_in(path =Path_cluster)

	return "Done"



#########################################################
# qsub is for miztli
#########################################################
def create_qsub(size =55, atom ="Au", path =""):
	file_name_out =  atom + 	str(size) +".out"
	file_name_sh = path + "/qsub_fhi.sh"
	print("Creating :" + file_name_sh)
	
	text = ["!/bin/bash \n", 
	"#BSUB -q  q_residual \n",
	"#BSUB -oo fhi-aims.%J.o \n",
	"#BSUB -eo fhi-aims.%J.e \n",
	# num cores 
	"#BSUB -n  16 \n",
	#nodos 
	'#BSUB -m "g1" \n',
	"module purge \n",
	"module load use.own\n",
	"module load fhi-aims/1\n",
	"mpirun aims.171221_1.scalapack.mpi.x < control.in > " + file_name_out]
	
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)
		
	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
	return "qsub_fhi.sh"
###########################################################
# nueva forma de calculo en archivo init_ga.py
#
#def create_qsub_init(size =55, atom ="Au", path ="", cores ="16", node= "g1"):
#	file_name_out =  atom + 	str(size) +".out"
#	file_name_sh = path + "/qsub_fhi.sh"
#	print("Creating :" + file_name_sh)
#	
#	text = ["!/bin/bash \n", 
#	"#BSUB -q  q_residual \n",
#	"#BSUB -oo fhi-aims.%J.o \n",
#	"#BSUB -eo fhi-aims.%J.e \n",
#	# num cores 
#	"#BSUB -n  {} \n".format(cores),
#	#nodos 
#	'#BSUB -m "{}"\n'.format(node),
#	"module purge \n",
#	"module load use.own\n",
#	"module load fhi-aims/1\n",
#	"module load python/3.7.6",
#	"python3 "]
#	#"mpirun aims.171221_1.scalapack.mpi.x < control.in > " + file_name_out]
#	#print(text)
#	with open(file_name_sh, "w") as fh: 
#		fh.writelines(text)
#		
#	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
#	return "qsub_fhi.sh"
#
#
#

#########################################################
# runsh is for ela 
#########################################################
def create_runsh(size =55, atom ="Au", path =""):
        file_name_out =  atom + str(size) +".out"
        file_name_sh = path + "/run.sh"
        print("Creating :" + file_name_sh)
        
        text = ["nohup mpiexec -np 24 -f /home/raet/machinefiles.txt /opt/fhi-aims/bin/aims.171221_1.scalapack.mpi.x < control.in"]
        
        with  open(file_name_sh, "w") as fh :

        	fh.writelines(text)
        
        subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
        #fh.close;
        return "run.sh"

##########################################################
def create_shforrunning(size =55, atom ="Au", path ="", cores = "16"):
	file_name_out =  atom + 	str(size) +".out"
	file_name_sh = path + "/shforrunning.sh"
	print("Creating :" + file_name_sh)
	root_dir =  os.path.dirname(os.path.abspath(__file__))
	complete_dir = os.path.join(root_dir, path)
	text = ["#!/bin/bash \n", 
	"mpirun -np " + str(cores)+ " aims.171221_1.scalapack.mpi.x < "+complete_dir +"/control.in > "+ complete_dir + "/" +file_name_out+"\n"]
	
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)
		fh.close()	
	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
	return "shforrunning.sh"
		
def create_shforrunning_name(name="Name.out", path ="", cores = "16"):
	file_name_out =  name +".out"
	file_name_sh = path + "/shforrunning.sh"
	print("Creating :" + file_name_sh)
	root_dir =  os.path.dirname(os.path.abspath(__file__))
	complete_dir = os.path.join(root_dir, path)
	text = ["#!/bin/bash \n", 
	"mpirun -np " + str(cores)+ " aims.171221_1.scalapack.mpi.x < "+complete_dir +"/control.in > "+ complete_dir + "/" +file_name_out+"\n"]
	
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)
		fh.close()	
	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
	return "shforrunning.sh"		



#Not needed 
#def create_directory(size =55, atom ="Au", path ="", add =0):
#	count = add
#	original_size = size
#	original_atom = atom 
#	original_path = path 
#	directory_name =original_atom + str(original_size)
#	directory_path = original_path+"/"+ directory_name	
#	answer = "not changing answer"
#	if count != 0:
#		 directory_path = original_path+ "/" + directory_name + "_" + str(count)
#		 answer = directory_path
#	if not os.path.exists(directory_path):
#		os.mkdir(directory_path)
#		print("creating folder '{}' ".format(directory_path))
#		answer =  directory_path	
#	else:
#		print("folder {} already exists".format(directory_path))
#		directory_path=create_directory(original_size, original_atom, original_path, count +1)
#		answer = directory_path		
#	
#	return answer	
#


def check_and_rename(file, add=0):
    original_file = file
    if add != 0:
        split = file.split("_")
        part_1 = split[0] + "_" + str(add)
        file = "_".join([part_1, split[1]])
    if not os.path.exists(file):
    	os.mkdir(file)
        # save here
    else:
        check_and_rename(original_file, add= add + 1)




def create_control_in(path =""):
	text = [' # DFT details\n',
	'xc                     pbe\n',
	'spin                   collinear            # non-collinear spin\n',
	'relativistic           atomic_zora scalar  # basis set (used zora for single-point, atomic_zora for opt)\n',
	'charge                 0.\n',
	'default_initial_moment 1\n',
	'\n',
	'# SCF CONVERGENCE\n',
	'occupation_type        gaussian 0.01         # this is required for metals\n',
	'charge_mix_param       0.2\n',
	'sc_accuracy_rho        1E-4\n',
	'sc_accuracy_eev        5E-3\n',
	'sc_accuracy_etot       5E-4\n',
	'sc_iter_limit          1000\n',
	'\n',
	'#  Relaxation\n',
	'relax_geometry               bfgs 5.e-2\n',
	'hessian_to_restart_geometry  .false.\n',
	'write_restart_geometry       .true.\n',
	'\n',
	'\n',
	'################################################################################\n',
	'#\n',
	'#  FHI-aims code project\n',
	'#  VB, Fritz-Haber Institut, 2009\n',
	'#\n',
	'#  Suggested "light" defaults for Au atom (to be pasted into control.in file)\n',
	'#  Be sure to double-check any results obtained with these settings for post-processing,\n',
	'#  e.g., with the "tight" defaults and larger basis sets.\n',
	'#\n',
	'################################################################################\n',
	'  species        Au\n',
	'#     global species definitions\n',
	'    nucleus             79\n',
	'    mass                196.966569\n',
	'#\n',
	'    l_hartree           4\n',
	'#\n',
	'    cut_pot             3.5  1.5  1.0\n',
	'    basis_dep_cutoff    1e-4\n',
	'#\n',
	'    radial_base         73 5.0\n',
	'    radial_multiplier   1\n',
	'    angular_grids specified\n',
	'      division   0.5066   50\n',
	'      division   0.9861  110\n',
	'      division   1.2821  194\n',
	'      division   1.5344  302\n',
	'#      division   2.0427  434\n',
	'#      division   2.1690  590\n',
	'#      division   2.2710  770\n',
	'#      division   2.3066  974\n',
	'#      division   2.7597 1202\n',
	'#      outer_grid 974\n',
	'      outer_grid 302\n',
	'################################################################################\n',
	'#\n',
	'#  Definition of "minimal" basis\n',
	'#\n',
	'################################################################################\n',
	'#     valence basis states\n',
	'    valence      6  s   1.\n',
	'    valence      5  p   6.\n',
	'    valence      5  d  10.\n',
	'    valence      4  f  14.\n',
	'#     ion occupancy\n',
	'    ion_occ     6  s   0.\n',
	'    ion_occ     5  p   6.\n',
	'    ion_occ     5  d   9.\n',
	'    ion_occ     4  f   14.\n',
	'################################################################################\n',
	'#\n',
	'#  Suggested additional basis functions. For production calculations, \n',
	'#  uncomment them one after another (the most important basis functions are\n',
	'#  listed first).\n',
	'#\n',
	'#  Constructed for dimers: 2.10, 2.45, 3.00, 4.00 AA\n',
	'#\n',
	'################################################################################\n',
	'#  "First tier" - max. impr. -161.60  meV, min. impr. -4.53 meV\n',
	'     ionic 6 p auto\n',
	'     hydro 4 f 7.4\n',
	'     ionic 6 s auto\n',
	'#     hydro 5 g 10\n',
	'#     hydro 6 h 12.8\n',
	'     hydro 3 d 2.5\n',
	'#  "Second tier" - max. impr. -2.46  meV, min. impr. -0.28 meV\n',
	'#     hydro 5 f 14.8\n',
	'#     hydro 4 d 3.9\n',
	'#     hydro 3 p 3.3\n',
	'#     hydro 1 s 0.45\n',
	'#     hydro 5 g 16.4\n',
	'#     hydro 6 h 13.6\n',
	'#  "Third tier" - max. impr. -0.49  meV, min. impr. -0.09 meV\n',
	'#     hydro 4 f 5.2\n',
	'#     hydro 4 d 5\n',
	'#     hydro 5 g 8\n',
	'#     hydro 5 p 8.2\n',
	'#     hydro 6 d 12.4\n',
	'#     hydro 6 s 14.8\n',
	'#  Further basis functions: -0.08 meV and below\n',
	'#     hydro 5 f 18.8\n',
	'#     hydro 5 g 20\n',
	'#    hydro 5 g 15.2']

	file_controlin = path + "/control.in"
	with open(file_controlin, "w") as fh:
		fh.writelines(text)
		fh.close()

	subprocess.call(["chmod", "754",file_controlin], universal_newlines=True)
	
	return "Done"

def get_hostname():
	host = subprocess.check_output(["hostname"], universal_newlines=True)
	print(" You're currently in %s"%host)
	return host

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def check_if_files_exists(path="", filename=""):
	if os.path.isfile(filename):
	    print ("File exist")
	else:
	    print ("File not exist")

def Cluster_size(N=55, R_ws=1.44):
	##Revisar bibliografia de esto 
	dist_max = round(2 * R_ws* math.pow(int(N) , 1/3), 4)
	print("For ", N , "atoms distance is : ", dist_max)
	return dist_max

def Pool_size(N= 55):
	##Revisar bibliografia de esto 
	pool = int(math.pow(int(N),2/3))
	return pool
	
def Number_ofGenerations(N=55):
	##Revisar bibliografia de esto 
	generations = int(math.pow(int(N),3/2))
	return generations

def Convergence(path):
	file = path+ "/geometry.in.next_step"
	if not os.path.exists(file):
		last_dir = str(path.split("/")[-1])
		command = "rm -r " + last_dir
		run_command = shlex.split(command)
		subprocess.call(run_command, universal_newlines = True, shell = True)
		return False
	else:
		return True

def create_files(Size = 55, Atom = "Au", Path ="", r_min = 2.0,r_max = 7,num_decimals =4,dist_min =2, dist_max =7,cores =16):
	

	atom = Atom 
	size = Size
	path = Path


	dist_max = Cluster_size(Size, R_ws= 1.44)
	#directory_name =create_directory(size, atom, path)
	directory_name = create_folder(name=Atom+str(Size), path=path)
	create_cluster(size, atom, path = directory_name, R_min = r_min,R_max = r_max,Num_decimals =num_decimals,Dist_min =dist_min, Dist_max =dist_max, cores=cores)
	
	#host =get_hostname()

	#if host == "basie":
	#	run_raw = "./" + create_runsh(size, atom, path= directory_name)
	#else:
	#	run_raw = "bsub < " + create_qsub(size,atom, path = directory_name)
	#
#
	#run_ready = shlex.split(run_raw)
	#create_control_in(path =directory_name)
#
	#with cd(directory_name):
	#	print(run_raw)
	#	subprocess.call(run_raw,universal_newlines = True, shell = True)
	#	#grep_cmd =shlex.split('grep  \"| Total energy of the DFT / Hartree-Fock s.c.f. calculation      :\" ' +  directory_name + "/nohup.out")
		#print(grep_cmd)
		#while subprocess.call(grep_cmd,universal_newlines =True, shell = True) == 1:
		#print("Listo")
        

	return directory_name
	
def check_convergence(filename="",path =""):
	converged = False 
	energy = 0 
	path.replace("\n", "")
	try :		
		grep_cmd ='grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation"  {}/{}'.format(path.replace("\n", ""), filename)	
		output =subprocess.check_output(grep_cmd,shell=True)
		output_string=str(output)
		vec1=output_string.split("   :   ")
		energy = float(vec1[1].split("eV")[0])
		converged = True
	
	except :
		print("Cluster didn't converged")
		command = "rm -r " + path
		print("Run: " , command)
	return converged, energy




def Proof_convergence(atom, size,  complete_path):

	converged = False 
	energy = 0 
	complete_path.replace("\n", "")
	try :		
		#with cd(complete_path):
			#grep_cmd =shlex.split('grep " Total energy of the DFT / Hartree-Fock s.c.f. calculation"      {}/nohup.out'.format(directory_name))
		grep_cmd ='grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation"  {}/{}{}.out'.format(complete_path.replace("\n", ""), atom, str(size))	
		#process =subprocess.run(grep_cmd,shell=True, check=True, universal_newlines=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
		#output = process.stdout
		#print(output)
		#ster = process.stderr
		#print(ster)
		output =subprocess.check_output(grep_cmd,shell=True)
		output_string=str(output)
		vec1=output_string.split("   :   ")
		energy = float(vec1[1].split("eV")[0])
		converged = True
	
	except :
		print("Cluster didn't converged")
		#last_dir = str(complete_path.split("/")[-1])
		command = "rm -r " + complete_path
		print("Run: " , command)
		#run_command = shlex.split(command)
		#subprocess.call(command, universal_newlines = True, shell = True)
	return converged, energy
	 
	 

##### Energia Raw
#grep_cmd ='''grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation"  /tmpu/lopb_g/raet_a/FHIaims/light_coarse/code/ga_FHIAIMS/pools_au6/Au10/Au10/Au10.out'''
#
#print("Command = " , grep_cmd, "\n")
##process =subprocess.call(grep_cmd,shell=True,  universal_newlines=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
#output =subprocess.check_output(grep_cmd,shell=True)
#
#print(output)
#output_string=str(output)
#print("string : ", output_string)
##ster = process.stderr
#vec1=output_string.split("   :   ")
##for x in vec1:
##       print(x,"\n")
#energy = float(vec1[1].split("eV")[0])
#print("Energy =" , energy)
#











def print_wami():
	process= subprocess.run(["pwd"], check=True, stdout=subprocess.PIPE, universal_newlines=True) 
	output = process.stdout
	print("I am here :", output)
	a_string = output.rstrip("\n")
	return a_string

def create_pool(N= 55, atom = "Au", path = "", R_min = 2.0, Num_decimals =4, Dist_min= 2 ,generation =0, cores = 16):
	pool_size = Pool_size(N)
	dist = Cluster_size(N)
	directories =[]
	#gen_path = path + "Gen" + str(generation) +"/"
	#preff = atom + str(N)
	#gen_dir = create_directory(preff, gen_path , add =0) + "/"
	#text =["path =", path]
	for x in range(pool_size):
		directory= create_files(Size = N, Atom = atom, Path = path, r_min = R_min ,r_max = dist ,num_decimals =Num_decimals ,dist_min =Dist_min , dist_max =dist)
		directories.append(directory)

	THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
	
	print("\n directories : ")

	#for x in directories:
	#    print(x , "\n ")

	#time.sleep(20)	
	used_directories =[]
	for x in directories:
		path_dum = os.path.join(THIS_FOLDER, x)
		dir_list = os.listdir(path_dum)  
		#print("Files and directories in '", path_dum, "' :", dir_list)
		file_running = path_dum +'/shforrunning.sh'
		if os.path.exists(file_running) == True: 
			#print("sh created")
			used_directories.append(path_dum)
	return used_directories		


	#		try:
	#			#run_raw = "." + file_running
	#			subprocess.call(file_running,universal_newlines = True, shell = True)
	#		except :
	#			print("Error running")
	#		converged, energy = Proof_convergence(directory_name=x, path = path)
	#		energies.append(energy)
	#		directories.append(x)
	#		text_line = ["directory= ", x ,"  energy= ", energy]
	#		text.append(text_line)
	#	except:
	#		energy=0
	#		directories.append(x)
	#		text_line = ["directory= ", x ,"  energy= ", energy]
	#		text.append(text_line)
#
	#with open(path + "Energy.txt", "w") as fh:
	#	fh.writelines(text)

#####

#
#
#
#
#Normalized_energies=Normalize_energies(Energies)
#fitnessed_energies= calculate_fitnesscalculate_fitness(Normalized_energies,func = "tanh")
#probabilities = probability_i(fitnessed_energies)

#########################################################Normalization			
def Normalization(E_i, E_min, E_max):
	p_i = (E_i - E_min)/(E_max - E_min)
	return round(p_i,5) 

def Normalize_energies(Energies):
	E_max = max(Energies)
	E_min = min(Energies)
	Normalized_energies = [round((E_i - E_min)/(E_max - E_min),5) for E_i in Energies]
	return Normalized_energies




####################################################### Fitness
def f_exp(alpha, p_i):
	f_i =math.pow(math.e, -alpha* p_i)	
	return round(f_i,5)

def f_tanh( p_i):
	f_i = (1/2)*(1-math.tanh(2*p_i -1))	
	return round(f_i, 5)

def f_lin(p_i,b =0.7 ):
	f_i = 1- b*p_i
	return round(f_i,5) 

def calculate_fitness(Normalized_energies,func = "tanh", alpha = 1):
	fitnessed_energies =[]
	if func == "tanh":
		fitnessed_energies = [ f_tanh(e_i) for e_i in Normalized_energies ] 
	elif func == "exp":
		fitnessed_energies	 = [ f_exp(alpha,e_i) for e_i in Normalized_energies ] 
	elif func == "lin":
		fitnessed_energies = [ f_lin(e_i, alpha) for e_i in Normalized_energies ] 
	else:
		print("Undefined fitness function using tanh")
		fitnessed_energies = [ f_tanh(e_i) for e_i in Normalized_energies ]
	return fitnessed_energies	

##################################################### Probability 
#### p_i = f_i /sum(f_i)
###1 == sum(p_i)
##### random
def probability_i(fitnessed):
	sum_fit = sum(fitnessed)
	p =[round((p_i/sum_fit),5) for p_i in fitnessed]
	return p

def selection_energy(Energies, fitnessed_energies):
	energies_rand= random.choices(Energies,weights =fitnessed_energies)
	return energies_rand



########################################################## code for Mutate

def kick_atom(atom,r_min=0, r_max=1):
    
    x_i = atom[0]
    y_i = atom[1]
    z_i= atom[2]
    theta = random.uniform(0, 2 * math.pi)
    phi = random.uniform(0, 2 * math.pi)
    r = random.uniform(r_min,r_max)
    x =round((r * math.cos(theta) * math.sin(phi)) +x_i,5)
    y =round((r * math.sin(theta) * math.sin(phi)) + y_i,5)
    z =round((r * math.cos(phi)) + z_i,5)
    vector = [x,y,z]
    return  vector



def read_geometry_nextstep(filename=""):
    matriz=[]
    lineas =[]
    with open(filename, "r") as fh:
        lineas = fh.readlines()
        fh.close()
    lineas_utiles= lineas[5:]
    for i in range(len(lineas_utiles)):
        atomo=lineas_utiles[i].replace("atom", "").replace("\n","").split(" ")
        str_list = list(filter(None, atomo))
        matriz.append(str_list)
    return matriz


def kick(filename = ""):
	lineas= read_geometry_nextstep(filename=filename)
	coord_txt = [linea[:-1] for linea in  lineas]
	atom_list = [linea[-1] for linea in  lineas]
	coordenadas= [[float(texto) for texto in linea] for linea in coord_txt]
	nuevas_coordenadas =[kick_atom(coordenada) for coordenada in coordenadas]
	[nuevas_coordenadas[i].append(atom_list[i]) for i in range(len(atom_list))]
	
	return nuevas_coordenadas

def kick_mutation(filename_mutated = "geometry.in", path="", original_file=""):	
	print("path:",path,"file",original_file, "filename_mutated", filename_mutated)
	matrix_used= kick(filename = original_file)
	print_matrix_geometryin(matrix= matrix_used, name=filename_mutated, path=path)
	
	
def choose_mutation(mutation ="",filename_mutated = "geometry.in", path="", original_file=""):
	if mutation =="kick":
		kick_mutation(filename_mutated = filename_mutated, path=path, original_file=original_file)

def create_files_mutation(size=52, atom="",path ="",cores =16):
	create_control_in(path =path)	
	create_shforrunning_name(name=atom+ str(size), path=path ,cores =cores)


#################################################### code for rotations
def rotate_over_z(point,angle):
    if angle == 0:
        new_point = point
    else:
        psi = angle
        qx =  math.cos(psi) * point[0] - math.sin(psi) * point[1]
        qy =  math.sin(psi) * point[0] + math.cos(psi) * point[1]
        qz = point[2] 
        new_point= [qx,qy,qz]
    return new_point

def rotate_over_x(point,angle):
    if angle == 0:
        new_point = point
    else:
        phi= angle
        qx = point[0]
        qy =  math.cos(phi) * point[1] - math.sin(phi) * point[2]
        qz =  math.sin(phi) * point[1] + math.cos(phi) * point[2]
        new_point= [qx,qy,qz]
    return new_point

def rotate_over_y(point,angle):
    if angle == 0:
        new_point = point
    else:
        theta =angle
        qx =   math.cos(theta) * point[0] + math.sin(theta) * point[2]
        qz =  -math.sin(theta) * point[0] + math.cos(theta) * point[2]
        qy = point[1]
        new_point= [qx,qy,qz]
    return new_point

def choose_rotations(number_of_rotations =2):
    list_of_rotations= list(range(3))
    choices_of_rotations = random.sample(list_of_rotations,number_of_rotations)
    #print(choices_of_rotations)
    #angle 1 is over x , angle 2 is over y and angle3 is over z
    if 0 in choices_of_rotations:
        angle1 = random.uniform(0, 2 * math.pi)
    else:
        angle1 = 0
    if 1 in choices_of_rotations:
        angle2 = random.uniform(0, 2 * math.pi)
    else:
        angle2 =0
    if 2 in choices_of_rotations:
        angle3= random.uniform(0, 2 * math.pi)
    else:
        angle3 =0
    angles=[angle1, angle2, angle3]
    return angles

def rotate_overxyz(point ,angles):
    
    point_x = rotate_over_x(point, angles[0])
    point_y= rotate_over_y(point_x, angles[1])
    point_z = rotate_over_z(point_y,angles[2])
    return point_z

def radius(point):
    r = math.sqrt(point[0]**2 + point[1]**2 + point[2]**2)
    return r

def rotate_matrix_over_xyz(matrix):
    angles_for_rotation = choose_rotations(number_of_rotations=2)
    print(angles_for_rotation)
    matrix_rotated=[]
    for point_str in matrix:
        point = []
        for i in range(3):
            point.append(float(point_str[i]))

        vector = rotate_overxyz(point, angles=angles_for_rotation)
        #print(round(radius(point),4),"=", round(radius(vector),4))
        vector_str=[ str(round(x,4)) for x in vector]
        vector_str.append(point_str[3])
        matrix_rotated.append(vector_str)

    return matrix_rotated

	
#### 3 angle rotation must be applied AzAyAx to preserve orientation 

#################################################################  code for Mating 

def order_matrix(matrix):
    order_z = [x[2] for x in matrix]
    order_z.sort(reverse = True)
    ordered_matrix=[]
    for order_i in order_z:
        #print(order_i)
        for x in matrix:
            if order_i == x[2]:
                ordered_matrix.append(x)
    return ordered_matrix   

def cut_n_splice(matrix1, matrix2, percentage_cut):
    if len(matrix1) == len(matrix2) and len(matrix1[0])== len(matrix2[0]):
        size_of_cut_int_1= int(round(len(matrix1)*(percentage_cut/100),0))
        print("size1",size_of_cut_int_1)
        size_of_cut_int_2 = len(matrix)-size_of_cut_int_1
        print("size2",size_of_cut_int_2)
        new_matrix =[]
        for i in range(size_of_cut_int_1):
            new_matrix.append(matrix1[i])
        for j in range(size_of_cut_int_2):
            new_matrix.append(matrix2[j+ size_of_cut_int_1])
        return new_matrix;  




########################################################################
def create_py(size=55, atom="Au", path="", cores =16):
	today = datetime.datetime.now()
	THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
	file_name_out =  path +"/"+ "run_"+atom + str(size) +".py"
	text=["#this will create the python file to run \n",
	"#Ga generated BY Rodrigo Espinola \n",
	"#created on  {} \n".format(today),
	"import sys\n",
	"sys.path.append('{}')\n".format(THIS_FOLDER),	
	'import atompacking_functions as af \n',
	'print("running python for run dirs ") \n'
	'af.run_dirs("{}/{}") \n'.format(THIS_FOLDER, path),
	'af.complete_cicle_mutation(file_dirs ="{}/{}/file_dirs.txt", Atom = "{}", Size = {}, path = "{}", cores ={}) \n'.format(THIS_FOLDER, path,atom,size,path,cores )
	#'af.Mutate(data_last_step="{}/{}/data_last_step.txt", path = "{}",  cores ={}, file_energies="{}/{}/energies.txt", Atom ="{}", Size={})\n'.format(THIS_FOLDER, path,path,cores,THIS_FOLDER, path,atom,size )
	#'af.Cicle_mutation(data_last_step="{}/{}/data_last_step.txt", path = "{}", name="", cores ={}, file_energies="{}", Atom ="{}", Size={})\n'.format(THIS_FOLDER, path,path,cores,THIS_FOLDER, path,atom,size )
	]
	with open(file_name_out, "w") as fh: 
		fh.writelines(text)
		fh.close()
	subprocess.call(["chmod", "754",file_name_out], universal_newlines=True)

def create_qsub_init(size =55, atom ="Au", path ="", cores ="16", node= "g1"):
	
	file_name_sh =  path+"/" + "qsub_fhi.sh"
	THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
	complete_path =  os.path.join(THIS_FOLDER,path)
	
	text = ["#!/bin/bash \n",
	#'''#BSUB -R "same[model] span[ptile='!',Intel_EM64T:16,Intel_a:20,Intel_b:20,Intel_c:32,Intel_d:32]" \n''',  
	"#BSUB -q  q_residual \n",
	"#BSUB -oo fhi-aims.%J.o \n",
	"#BSUB -eo fhi-aims.%J.e \n",
	# num cores 
	"#BSUB -n  {} \n".format(cores),
	#nodos 
	'#BSUB -m "{}"\n'.format(node),
	"module purge \n",
	"module load use.own\n",
	"module load fhi-aims/1\n",
	"module load mpi/intel-2017_update3 \n",		
	"module load python/3.7.6 \n",
	"python3 {}/run_{}.py \n".format(complete_path,atom + str(size))]
	#"mpirun aims.171221_1.scalapack.mpi.x < control.in > " + file_name_out]
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)

	subprocess.call(["chmod", "777",file_name_sh], universal_newlines=True)
	return file_name_sh

def run_calc(filename):    
	host =get_hostname()
	if host == "basie":
		run_raw = "./" + filename
	else:
		run_raw = "bsub < " + filename
	#subprocess.call(run_raw,universal_newlines = True, shell = True)
	process= subprocess.run(run_raw, check=True, stdout=subprocess.PIPE, universal_newlines=True) 
	output = process.stdout
	print("Output", output)	


def init_calc(Size =55, Atom ="Au", Path ="", Cores ="16", Node= "g1"):
	Path = create_folder(name=Atom+str(Size), path= Path)
	create_py(size=Size, atom=Atom, path=Path, cores =int(Cores))
	file_bsub = create_qsub_init(size=Size, atom=Atom,path=Path,cores=Cores, node=Node)
	print("file bsub:", file_bsub)
	run_calc(file_bsub)


def create_folder( name ="Au_6", path ="", add =0):
	count = add
	original_name = name 
	original_path = path
	directory_name =original_name
	directory_path = original_path +"/"+ directory_name	
	answer = "not changing answer"
	if count != 0:
		 directory_path = original_path+"/" + directory_name + "_" + str(count)
		 answer = directory_path
	if not os.path.exists(directory_path):
		os.mkdir(directory_path)
		answer =  directory_path	
	else:
		directory_path=create_folder( original_name, original_path, count +1)
		answer = directory_path		
	
	return answer	 			



def create_all_files(Size =55, Atom ="Au", Path ="", Cores ="16", Node= "g1"):
	path_master = create_folder(name=Atom+str(Size), path= Path)
	print("Path :", path_master)
	dirs = create_pool(N= Size, atom = Atom, path =path_master, R_min = 2.0, Num_decimals =4, Dist_min= 2 ,generation =0, cores = int(Cores))
	print("Creating : file of directories" )
	file_dirs= path_master + "/file_dirs.txt"
	with open(file_dirs, "w") as fh:
		for x in dirs:
			fh.write(x + "\n")
		fh.close()
	
	create_py(size=Size, atom=Atom, path=path_master, cores =int(Cores))
	file_bsub = create_qsub_init(size=Size, atom=Atom,path=path_master,cores=Cores, node=Node)	
	print("For running write: bsub < ", file_bsub)

	return file_dirs

def append_file(text_to_append="", file_to_append=""):
	with open(file_to_append, "a") as fh:
		fh.write(text_to_append + "\n")
		fh.close()



def read_files(file_dirs):
	with open(file_dirs, "r") as fh:
		directories = fh.readlines()
		fh.close()
	return directories


def check_convergence_pool_first_step( file_dirs ="", Atom = "Au", Size = 52, path ="",cores= 16 ):
	name = Atom + str(Size) 
	directories = read_files(file_dirs)
	Energies =[]
	Converged = []
	Normalized_energies =[]
	

	for x in directories:
		#print("currently in ", x)
		converged = False
		energy =0
		converged , energy_not_rounded = Proof_convergence(atom = Atom, size = Size,  complete_path = x)
		energy = round(energy_not_rounded, 5)
		Converged.append(converged)
		Energies.append(energy)

	E_max = max(Energies)
	E_min = min(Energies)
	E_min_2 = min([x for x in a if x != min(Energies)]) 
	Normalized_energies=Normalize_energies(Energies)
	fitnessed_energies= calculate_fitness(Normalized_energies,func = "tanh")
	probabilities = probability_i(fitnessed_energies)
	
	print("energies", Energies)
	print("Normalized ", Normalized_energies)
	print("fitness," , fitnessed_energies)
	print("probabilities", probabilities)
	
	#data_last_step=[Energies, Normalized_energies, fitnessed_energies, probabilities, directories]
	Number_ofsteps=Number_ofGenerations(Size)
	file_energies= print_energies(filename="energies.txt",path=path, Energies=Energies, Normalized_energies=Normalized_energies, fitnessed_energies=fitnessed_energies, probabilities=probabilities, directories=directories, text_option="a")
	data_last_step= print_energies(filename="data_last_step.txt",path=path, Energies=Energies, Normalized_energies=Normalized_energies, fitnessed_energies=fitnessed_energies, probabilities=probabilities, directories=directories, text_option="w",step=0, Number_ofGenerations=Number_ofsteps)
	return file_energies, data_last_step

def check_convergence_pool( file_dirs ="", Atom = "Au", Size = 52, path ="",cores= 16 ):
	name = Atom + str(Size) 
	directories = read_files(file_dirs)
	Energies =[]
	Converged = []
	Normalized_energies =[]
	

	for x in directories:
		#print("currently in ", x)
		converged = False
		energy =0
		converged , energy_not_rounded = Proof_convergence(atom = Atom, size = Size,  complete_path = x)
		energy = round(energy_not_rounded, 5)
		Converged.append(converged)
		Energies.append(energy)

	E_max = max(Energies)
	E_min = min(Energies)
	E_min_2 = min([x for x in a if x != min(Energies)])

	Normalized_energies=Normalize_energies(Energies)
	fitnessed_energies= calculate_fitness(Normalized_energies,func = "tanh")
	probabilities = probability_i(fitnessed_energies)
	
	print("energies", Energies)
	print("Normalized ", Normalized_energies)
	print("fitness," , fitnessed_energies)
	print("probabilities", probabilities)
	
	#data_last_step=[Energies, Normalized_energies, fitnessed_energies, probabilities, directories]
	Number_ofsteps=Number_ofGenerations(Size)
	file_energies= print_energies(filename="energies.txt",path=path, Energies=Energies, Normalized_energies=Normalized_energies, fitnessed_energies=fitnessed_energies, probabilities=probabilities, directories=directories, text_option="a")
	data_last_step= print_energies(filename="data_last_step.txt",path=path, Energies=Energies, Normalized_energies=Normalized_energies, fitnessed_energies=fitnessed_energies, probabilities=probabilities, directories=directories, text_option="a",step=0, Number_ofGenerations=Number_ofsteps)
	########################selection for mutation
	selected_energy = selection_energy(Energies, fitnessed_energies)
	print("Selected Energy: ", selected_energy)
	index_selected = Energies.index(selected_energy[0])

	text_selecting ="Selected Energy: "+ str(selected_energy[0]) + " ,Index of Energy:" + str(index_selected) + " ,directory : " + str(directories[index_selected])
	print(text_selecting)
	mating_text = "Best suited energies for mating =" + str(E_min) +"with index"+str(Energies.index(E_min))+" and "+ str(E_min_2) +"with index"+str(Energies.index(E_min_2))

	with open(file_energies, "a") as fa:
		fa.write(text_selecting)
		fa.write(mating_text)
		fa.close()

	############################  mutation

	path_mutated = create_folder(name="{}_mutated".format(name), path= path)
	#append_file(text_to_append=path_mutated, file_to_append=file_dirs)
	kick_mutation(filename_mutated = "geometry.in", path= path_mutated, original_file=str(directories[index_selected]).replace("\n", "")+"/geometry.in.next_step")
	create_files_mutation(size=Size, atom=Atom,path =path_mutated,cores =cores)
	run_file(path=path_mutated, filename= './shforrunning.sh')


	converged, mutated_energy = check_convergence(filename=Atom+ str(Size)+".out",path =path_mutated)
	Normalized_mutated_energy =Normalization(mutated_energy, E_min, E_max)
	fitnessed_mutated = f_tanh( Normalized_mutated_energy)
	
	with open(file_energies, "a") as fh:
		fh.write("After mutation:\n")
		fh.write("Energies,\t  Normalized_energies,\t fitnessed_energies,\t prob,\t dir \n")	
		fh.write(str( mutated_energy)+",\t"+ str(Normalized_mutated_energy) + ",\t"+ str(fitnessed_mutated) + ",\t" + path_mutated+"\n")	
		fh.close()

	#path_mating = create_folder(name="{}_mating".format(name), path= path)
	#kick_mutation(filename_mutated = "geometry.in", path= path_mutated, original_file=str(directories[index_selected]).replace("\n", "")+"/geometry.in.next_step")
	#create_files_mutation(size=Size, atom=Atom,path =path_mutated,cores =cores)
	#run_file(path=path_mutated, filename= './shforrunning.sh')
#
#
	#converged, mutated_energy = check_convergence(filename=Atom+ str(Size)+".out",path =path_mutated)
	#Normalized_mutated_energy =Normalization(mutated_energy, E_min, E_max)
	#fitnessed_mutated = f_tanh( Normalized_mutated_energy)
	#
	#with open(file_energies, "a") as fh:
	#	fh.write("After mutation:\n")
	#	fh.write("Energies,\t  Normalized_energies,\t fitnessed_energies,\t prob,\t dir \n")	
	#	fh.write(str( mutated_energy)+",\t"+ str(Normalized_mutated_energy) + ",\t"+ str(fitnessed_mutated) + ",\t" + path_mutated+"\n")	
	#	fh.close()
#text option is "a" for append , "w" for write
def print_energies(filename="",path="./", Energies=[], Normalized_energies=[], fitnessed_energies=[], probabilities=[], directories=[],text_option="a",step=None,Number_ofGenerations=32): 
	file_energies = path+"/" + filename
	with open(file_energies, text_option) as fh:
		if step is not None:
			fh.write(str(step)+"/"+str(Number_ofGenerations)+ "\n")
		fh.write("Energies,\t  Normalized_energies,\t fitnessed_energies,\t prob,\t dir \n")	
		for i in range(len(Energies)):
			fh.write(str( Energies[i])+",\t"+ str(Normalized_energies[i]) + ",\t"+ str(fitnessed_energies[i]) + ",\t"+ str(probabilities[i])+",\t" + directories[i]+"\n")
		fh.close()	
	return file_energies

def remove_file(filename=""):
	if os.path.exists(filename):
  		os.remove(filename)
	else:
  		print("The file does not exist")

def read_data(filename="",path="./"):
    vector_headers = ['Energies', '  Normalized_energies', ' fitnessed_energies', ' prob', ' dir ']
    with open(path+filename, "r") as f:
        lines_read=f.readlines()
        f.close()
    string_steps = lines_read[0]
    vector_steps = string_steps.split("/")
    step= int(vector_steps[0])
    Number_ofGenerations= int(vector_steps[1])
    lines=lines_read
    print(lines)
    Energies_1=[]
    Normalized_energies_1=[]
    fitnessed_energies_1=[]
    probabilities_1=[]
    directories_1=[]
    print("Lines:")
    #x = re.search("^The.*Spain$", txt)
    [print(x) for x in lines ]
    for line in lines:
        vector_line=line.replace("\t","").replace("\n","").split(",")
        print("index: ",lines.index(line),"Vector line")
        [print(x) for x in vector_line ]
        print(vector_line)
        len_vector = len(vector_line)
        if vector_line != [''] and vector_line!= vector_headers and len_vector > 1:
            Energies_1.append(float(vector_line[0]))
            Normalized_energies_1.append(float(vector_line[1]))
            fitnessed_energies_1.append(float(vector_line[2]))
            probabilities_1.append(float(vector_line[3]))
            directories_1.append(vector_line[4])
    return Energies_1, Normalized_energies_1, fitnessed_energies_1, probabilities_1, directories_1, step, Number_ofGenerations

def Cicle_mutation(data_last_step= "", path = "", name="", cores =16, file_energies="", Atom ="Au", Size=52):
	try:	
		#check_if_files_exists()
		with open(data_last_step, "r") as fa:
			text=fa.readline()
			fa.close()
		text=text.replace("/n","")
		vector_1 = text.split("/")
		vector = [str(x) for x in vector_1]
		for x in range (int(vector[0]), int(vector[1])):
			print("Step:",vector[0]," of", vector[0])
			Mutate(data_last_step= data_last_step, path = path, name=name, cores =cores, file_energies=file_energies, Atom =Atom, Size=Size)
	except IOError:
		print("Data last step accesible")

def Mutate(data_last_step= "", path = "", name="", cores =16, file_energies="", Atom ="Au", Size=52):
	pool_size = Pool_size(N=Size)
	print("selecting energy from: ", data_last_step)
	Energies, Normalized_energies, fitnessed_energies, probabilities, directories, step, Number_ofGenerations = read_data(filename=data_last_step,path="") 
	selected_energy = selection_energy(Energies, fitnessed_energies)
	#print("Selected Energy: ", selected_energy)
	index_selected = Energies.index(selected_energy[0])
	text_selecting ="Selected Energy: "+ str(selected_energy[0]) + " ,Index of Energy:" + str(index_selected) + " ,directory : " + str(directories[index_selected]+"\n")
	print(text_selecting)
	with open(file_energies, "a") as fa:
		fa.write(text_selecting)
		fa.close()

	############################mutation

	path_mutated = create_folder(name="{}_mutated".format(Atom+str(Size) ), path= path)
	kick_mutation(filename_mutated = "geometry.in", path= path_mutated, original_file=str(directories[index_selected]).replace("\n", "")+"/geometry.in.next_step")
	create_files_mutation(size=Size, atom=Atom,path =path_mutated,cores =cores)
	run_file(path=path_mutated, filename= './shforrunning.sh')

	#E_min= min(Energies)
	#E_max = max(Energies)
	#index_min = Energies.index(E_min)

	converged, mutated_energy = check_convergence(filename=Atom+ str(Size)+".out",path =path_mutated)
	#Normalized_mutated_energy =Normalization(mutated_energy, E_min, E_max)
	#fitnessed_mutated = f_tanh( Normalized_mutated_energy)
	if converged == True:
		new_energies = Energies
		new_energies.append(mutated_energy)
		new_E_max = max(new_energies)
		new_E_max_index= Energies.index(new_E_max)
		new_directories= directories
		new_directories.append(path_mutated)

		new_energies.remove(new_E_max)
		new_directories.remove(new_directories[new_E_max_index])
		#new_index = new_energies.index(mutated_energy)

		Normalized_new_energies=Normalize_energies(new_energies)
		fitnessed_new_energies= calculate_fitness(Normalized_new_energies,func = "tanh")
		new_probabilities = probability_i(fitnessed_new_energies)

		#new_directories[index_min] = path_mutated
		
	

		with open(file_energies, "a") as fh:
			fh.write("After mutation:\n")
			fh.write("Energies,\t  Normalized_energies,\t fitnessed_energies,\t prob,\t dir \n")	
			for i in range(len(new_energies)):
				fh.write(str( new_energies[i])+",\t"+ str(Normalized_new_energies[i]) + ",\t"+ str(fitnessed_new_energies[i]) + ",\t"+ str(new_probabilities[i])+",\t" + new_directories[i]+"\n")
			fh.close()

		remove_file(filename=path +"/"+ "data_last_step.txt")
		data_last_step= print_energies(filename="data_last_step.txt",path=path, Energies=new_energies, Normalized_energies=Normalized_new_energies, fitnessed_energies=fitnessed_new_energies, probabilities=new_probabilities, directories=new_directories, text_option="a",step=step +1 , Number_ofGenerations=Number_ofGenerations)


		#if float(mutated_energy) < float(E_min):
		#	print("Mutation worked :")
		#	
		#	remove_file(filename=path +"/"+ "data_last_step.txt")
		#	data_last_step= print_energies(filename="data_last_step.txt",path=path, Energies=new_energies, Normalized_energies=Normalized_new_energies, fitnessed_energies=fitnessed_new_energies, probabilities=new_probabilities, directories=directories, text_option="a",step=step +1 , Number_ofGenerations=Number_ofGenerations)
		#else :
		#	print("Mutation didn't work")
	else:
		print("Mutation didn't work")
		remove_file(filename=path +"/"+ "data_last_step.txt")
		#Energies, Normalized_energies, fitnessed_energies, probabilities, directories, step, Number_ofGenerations
		data_last_step= print_energies(filename="data_last_step.txt",path=path, Energies=Energies, Normalized_energies=Normalized_energies, fitnessed_energies=fitnessed_energies, probabilities=probabilities, directories=directories, text_option="a",step=step +1 , Number_ofGenerations=Number_ofGenerations)

 	



def run_file(path="", filename= './shforrunning.sh'):
	with cd(path):
		subprocess.call(filename,universal_newlines = True, shell = True)



def complete_cicle_mutation(file_dirs ="", Atom = "Au", Size = 52, path ="",cores= 16 ):
	file_energies, data_last_step = check_convergence_pool_first_step( file_dirs =file_dirs, Atom = Atom, Size = Size, path =path,cores=cores )
	Mutate(data_last_step= data_last_step, path = path, name="{}+{}".format(Atom,str(Size)), cores =cores, file_energies=file_energies, Atom =Atom, Size=Size)
	Cicle_mutation(data_last_step= data_last_step, path = path, name="{}+{}".format(Atom,str(Size)), cores =16, file_energies=file_energies, Atom =Atom, Size=Size)




def run_dirs(path =""):
	with  open(path + '/file_dirs.txt', 'r') as file1:
		directories = file1.readlines()
		file1.close()

	for x in directories:
		s= x.replace('\n', '')

		with cd(s):
			subprocess.call('./shforrunning.sh',universal_newlines = True, shell = True)
			#process =subprocess.run('./shforrunning.sh', shell=True, universal_newlines=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
			#output = process.stdout
			#print("output: ", output)
			#ster = process.stderr
			#print("ster: ", ster)
			#process =subprocess.call(my_file, universal_newlines=True)
		
##rerankear el pool y conservar los n mejores
def minN(elements, n):
    return sorted(elements, reverse=False)[:n] 
## rotar sobre angulos pphi y rho y cortar en el nuevo z 
