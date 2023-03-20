import numpy as np
from matplotlib import pyplot as plt


# -----------------------------
def f_sys_initialize(sys_size, init_type, init_rand_state, n_states=2):
    
    """generates initial state of the system
    
    Arguments:
        sys_size (int): size of system
        init_type (str): type of initialization
            - "random" or "r": systems elements randomly initiated as 1
            - "centre" or "c": single element at centre initiated as 1
        init_rand_state (int): random state to use (used only when 'init_type' is 'random'; otherwise ignored)
        n_states (int): number of states
        
    Returns:
        sys_init (numpy array): initialized system
    """
    
    if init_type.lower() in ["random", "r"]:
        rng_init = np.random.default_rng(seed=init_rand_state) # random generator for initial state
        sys_init = rng_init.integers(n_states, size=sys_size) # create initial state

    if init_type.lower() in ["centre", "center", "c"]:
        sys_init = np.zeros(shape=(sys_size, )) # create initial state
        sys_init[int(sys_size/2)] = 1
        
    return (sys_init)


# -----------------------------
def f_WolframCA_rule_lookup(rule_number, nb_size, n_states):
    
    """generates lookup for a give wolfram CA rule number

    Arguments:
        rule_number (int): rule to use
        nb_size (int): size of neighborhood to use
        n_states (int): possible states of a cell
    
    Returns:
        input_pattern (numpy array): array where each row represents possible local neighborhood configuration
        output_pattern (list of integers): output corresponding to each neighborhood (row) in input_pattern
    """
    
    n_nb_configs = n_states ** nb_size # no. of unique configurations of neighborhood
    
    # list storing N-bit representation of rule number where N=no of possible neighborhood configurations
    # e.g.: rule 110 becomes [0 1 1 0 1 1 1 0] when N=8 (i.e. n_states=2 & nb_size=3)
    output_pattern = [int(x) for x in np.binary_repr(rule_number, width=n_nb_configs)]

    # creating array with (n_nb_configs) rows and (nb_size) columns; each row will represent a possible neighborhood
    # During evolution, if local neighborhood is ith row in 'input_pattern', new cell state will be ith element in 'output_pattern'
    input_pattern = np.zeros([n_nb_configs, nb_size])
    for i in range(n_nb_configs):
        input_pattern[i, :] = [int(x) for x in np.binary_repr(n_nb_configs-1-i, width=nb_size)]
        
    return (input_pattern, output_pattern)


# -----------------------------
def f_evolve_WolframCA(sys_size, init_type, init_rand_state, nb_size, n_states, BC_type, rule_number, time_steps):
    
    """initialize and evolve system using wolfram CA rules

    Args:
        sys_size (int): size of system
        init_type (str): type of initialization
            - "random" or "r": systems elements randomly initiated as 1
            - "centre" or "c": single element at centre initiated as 1
        init_rand_state (int): random state to use (used only when 'init_type' is 'random'; otherwise ignored)
        nb_size (int): size of neighborhood to use
        n_states (int): possible states of a cell
        BC_type (str): Boundary condition
            - "periodic" or "p": implements periodic boundary condition
            - "fix-L-R": fixed boundary condition where L and R are fixed cell states for left and right boundary (example: "fix-1-1")
        rule_number (int): wolfram rule to use (integer between 0 and 255)
        time_steps (int): number of time steps

    Returns:
        sys_store_list (list of arrays): list that stores each time state of system
    """
    
    n_nb_configs = n_states ** nb_size # no. of unique configurations of neighborhood
    nb_order = int((nb_size - 1)/2) # order of neighborhood (nearest neighbour order)

    sys_init = f_sys_initialize(sys_size, init_type, init_rand_state, n_states)
    sys_store_list = [] # list to store system configration after each time step
    sys_store_list.append(sys_init) # add to system store list

    input_pattern, output_pattern = f_WolframCA_rule_lookup(rule_number, nb_size, n_states) # create input and output arrays for looking up rules

    for t in range(1, time_steps + 1):
        sys_state_old = sys_store_list[t-1]

        if BC_type.lower() in ["periodic", "p"]:
            sys_state_old_bc = np.concatenate((sys_state_old[-nb_order:],
                                               sys_state_old,
                                               sys_state_old[0:nb_order]))
        
        if ("fix" in BC_type) or ("Fix" in BC_type):
            bc_L, bc_R = int(BC_type.split("-")[1]), int(BC_type.split("-")[2])
            sys_state_old_bc = np.concatenate((np.full((nb_order,), bc_L),
                                               sys_state_old,
                                               np.full((nb_order,), bc_R)))

        sys_state_new_bc = np.zeros(shape=sys_state_old_bc.shape, dtype=int)

        for i in range(nb_order, sys_size+nb_order):
            nb = sys_state_old_bc[i-nb_order:i+nb_order+1]

            for k in range(0, n_nb_configs):
                # if neighborhood matches; assign new state based on 'output_pattern'
                if np.array_equal(input_pattern[k], nb):
                    sys_state_new_bc[i] = output_pattern[k]

        sys_state_new = sys_state_new_bc[nb_order:-nb_order]
        sys_store_list.append(sys_state_new) # add to store
    
    return (sys_store_list)


# -----------------------------
def f_space_time_plot(sys_store_list, c_map, plt_title):
    """create space time plot

    Args:
        sys_store_list (list of arrays): list that stores each time state of system
        c_map (str): colormap to use (default="summer"; any matplotlib compatible cmap can be used)
        plt_title (str): title of plot

    Returns:
        None
    """
    plt.figure(figsize=(8,10))
    fig = plt.imshow(np.array(sys_store_list), cmap=c_map)
    plt.title(plt_title)
    plt.ylabel("Time steps", fontsize=12)
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.spines['left'].set_visible(False)
    fig.axes.spines['top'].set_visible(False)
    fig.axes.spines['right'].set_visible(False)
    fig.axes.spines['bottom'].set_visible(False)
    plt.show()
    
    return None