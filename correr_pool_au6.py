import atompacking_functions as af
#af.create_pool(6, "Au","pools_au6/")

#import init_ga as init 
#af.init_calc(6, "Au", "pools_au6/",6,"g1")

file_dirs= af.create_all_files(Size= 6, Atom="Au", Path="pools_au6",Cores= 16,Node= "g1")
#af.check_convergence_pool( file_dirs =file_dirs, Atom = "Au", Size = 6, path = "pools_au6")