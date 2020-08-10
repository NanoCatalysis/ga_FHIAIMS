import atompacking_functions as af
import datetime

def create_qsub_init(size =55, atom ="Au", path ="", cores ="16", node= "g1"):
    today = datetime.datetime.now()
	file_name_out =  atom + 	str(size) +".out"
	file_name_sh = path + "/qsub_fhi.sh"
	print("Creating :" + file_name_sh)
	
	text = ["!/bin/bash \n", 
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
	"module load python/3.7.6",
	"python3 "]
	#"mpirun aims.171221_1.scalapack.mpi.x < control.in > " + file_name_out]
	#print(text)
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)
		
	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
	return "qsub_fhi.sh"

def run_calc(size =55, atom ="Au",directory_name=""):    
host =af.get_hostname()
if host == "basie":
	run_raw = "./" + create_runsh(size, atom, path= directory_name)
else:
	run_raw = "bsub < " + create_qsub(size,atom, path = directory_name)
