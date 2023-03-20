from core_functions import f_evolve_WolframCA, f_space_time_plot

# -----------------------------
# Take inputs
sys_size = int(input("system size (integer): "))
init_type = input("Type of initialization ('r-random' / 'c-centre'): ")

init_rand_state = 0
if init_type.lower() in ["random", "r"]:
    init_rand_state = int(input("Random state for initialization (integer): "))

BC_type = input("Boundary condition ('p-periodic' / 'fix-L-R'-fixed [e.g. 'fix-1-0']): ")
c_map = "summer"
rule_number = int(input("Rule to use (integer: [0, 255]): "))
time_steps = int(input("Time steps (integer): "))
plot_type = "space-time"

# -----------------------------
# Fixed parameters
nb_size = 3
print("Neighborhood size = 3 (FIXED)")
    
n_states = 2
print("Number of possible cell states = 2 (FIXED)")

# -----------------------------
# Evolve system
sys_store_list = f_evolve_WolframCA(sys_size, init_type, init_rand_state, nb_size, n_states, BC_type, rule_number, time_steps)

# -----------------------------
# Create title of time-space plot
if init_type.lower() in ["random", "r"]:
    plt_title = f"Rule: {rule_number} || Initialize: '{init_type}-{init_rand_state}', BC: '{BC_type}'"
else:
    plt_title = f"Rule: {rule_number} || Initialize: '{init_type}', BC: '{BC_type}'"

# -----------------------------
# Create space time plot
f_space_time_plot(sys_store_list, c_map, plt_title)