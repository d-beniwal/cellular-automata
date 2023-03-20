from helpers_nuclei_growth import *

print("Enter parameters to create initial state of the system...")
sys_size_input = input("System size [height width] (e.g. 100 150): ")
sys_size = tuple(int(i) for i in sys_size_input.split())
nuclei_pos = input("Nuclei position ('c'-centre, 'r'-random): ")
n_nuclei = int(input("No. of nuclei to initiate: "))
nuclei_shape = input("Nuclei shape ('c'-circle, 's'-square): ")
nuclei_size = int(input("Nuclei size: "))
min_nuclei_spacing = float(input("Minimum spacing between nuclei: "))

print("\nEnter parameters to evolve the system...")
nb_size = int(input("Size of neighborhood (odd integer;e.g. 3, 5, 7): "))
nb_type = input("Neighborhood type ('m'-Moore, 'vn'-Von Newmann): ")
BC_type = input("Boundary condition ('p'-periodic): ")
rule = int(input("RULE-Minimum neighbors need to change state: "))
t_steps = int(input("Time steps: "))
savename_gif = input("Savename for gif of time evolution ('n'-don't save): ")


print("\n| Initializing system... ", end="")
sys_init_state, sys_init_grainID = f_sys_initialize(sys_size, nuclei_pos, n_nuclei,
                                                    nuclei_shape, nuclei_size, min_nuclei_spacing)

print("DONE.")

print("-- Time evolution initiated --")
sys_store_state, sys_store_grainID = f_evolve_sys(sys_init_state, sys_init_grainID,
                                                  nb_size, nb_type, BC_type,
                                                  rule, t_steps)

print("\nFinished.")

if savename_gif.lower() not in ["n", "no", "false"]:
    print("\n| Saving gif...", end="")
    f_2Darray_list_to_gif(sys_store_state, f"{savename_gif}-state")
    f_2Darray_list_to_gif(sys_store_grainID, f"{savename_gif}-grainID")
    print("DONE.")