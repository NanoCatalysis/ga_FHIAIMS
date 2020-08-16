import atompacking_functions as af
import datetime
import subprocess
def create_py(size =55, atom ="Au", path =""):
	today = datetime.datetime.now()
	file_name_out =  "run_"+atom + str(size) +".py"
	text=["#this will create the python file to run",
	"#Ga generated BY Rodrigo Espinola",
	"#created on  {}".format(today),
	'import atompacking_functions as af',
	'af.create_pool("{}", "{}","{}")'.format(size,atom,path)
	]
	with open(file_name_out, "w") as fh: 
		fh.writelines(text)
		
	subprocess.call(["chmod", "754",file_name_out], universal_newlines=True)

def create_qsub_init(size =55, atom ="Au", path ="", cores ="16", node= "g1"):
	#today = datetime.datetime.now()
	#file_name_out =  atom + 	str(size) +".out"
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
	"python3 run_"+atom + str(size) +".py" ]
	#"mpirun aims.171221_1.scalapack.mpi.x < control.in > " + file_name_out]
	#print(text)
	with open(file_name_sh, "w") as fh: 
		fh.writelines(text)

	subprocess.call(["chmod", "754",file_name_sh], universal_newlines=True)
	return file_name_sh

def run_calc(filename):    
	host =af.get_hostname()
	if host == "basie":
		run_raw = "./" + filename
	else:
		run_raw = "bsub < " + filename
	subprocess.call(run_raw,universal_newlines = True, shell = True)	


def init_calc(Size =55, Atom ="Au", Path ="", Cores ="16", Node= "g1"):
	create_py(size=Size, atom=Atom, path=Path)
	file_bsub = create_qsub_init(size=Size, atom=Atom,path=Path,cores=Cores, node=Node)
	run_calc(file_bsub)
