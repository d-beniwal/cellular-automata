from helpers_clustering_of_states import *

print("Enter parameters to create initial state of the system...")
sys_size_input = input("System size [height width] (e.g. 100 150): ")
sys_size = tuple(int(i) for i in sys_size_input.split())
n_states = int(input("No. of states in system (integer, minimum value= 2): "))
state_fractions_input = input("Phase fraction of states except last state (e.g.: if two states: 0.45, if three states: 0.3 0.4): ")
state_fractions = [float(i) for i in state_fractions_input.split(" ")]

print("\nEnter parameters to evolve the system...")
nb_size = int(input("Size of neighborhood (odd integer;e.g. 3, 5, 7): "))
nb_type = input("Neighborhood type ('m'-Moore, 'vn'-Von Newmann): ")
BC_type = input("Boundary condition ('p'-periodic): ")
t_steps = int(input("Time steps: "))
savename_gif = input("Savename for gif of time evolution ('n'-don't save): ")


print("\n| Initializing system... ", end="")
sys_init_state = f_sys_initialize(sys_size, n_states, state_fractions)
print("DONE.")

print("-- Time evolution initiated --")
sys_store_state = f_evolve_sys(sys_init_state, nb_size, nb_type, BC_type, t_steps)

print("\nFinished.")

if savename_gif.lower() not in ["n", "no", "false"]:
    print("\n| Saving gif...", end="")
    f_2Darray_list_to_gif(sys_store_state, f"{savename_gif}-state")
    print("DONE.")