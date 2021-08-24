import atompacking_functions as af
for i in range(3,5):
	file_dirs= af.create_all_files(Size= i, Atom="Au", Path="pools_test",Cores=16,Node= "g1", queue = "q_residual", vdw = True)