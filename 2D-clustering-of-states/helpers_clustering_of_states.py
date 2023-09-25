import numpy as np
    
    
def f_sys_initialize(sys_size, n_states, state_fractions):
    """creates initial state of system containing all states
    
    Args:
        sys_size (tuple): (height, width) of system
        n_states (int): number of different states present in system
        state_fractions (list): Fraction of different states. Contains a value for each state except
                                the last which is calculated as 1 - (sum of all other state fractions)
        
    Returns:
        sys_init_state (numpy array): initial state of system
    """
    
    state_ids = np.arange(1, n_states+1, 1) # ids to assign each state (1, 2, 3, ....)
    sys_n_cells = sys_size[0] * sys_size[1] # total cells in system
    sys_flattened = np.zeros(shape=sys_n_cells) # 1D array with number of cells
    
    n_unassigned_cells = len(np.where(sys_flattened==0)[0]) # number of unassigned cells (i.e. with zero value)

    # Run untill a state is assigned to each cell
    while n_unassigned_cells > 0:
        n_unassigned_cells = len(np.where(sys_flattened==0)[0]) # number of unassigned cells (i.e. with zero value)
        state_n_cells = [int(round(state_frac*n_unassigned_cells)) for state_frac in state_fractions] # number of cells to be assigned for each state
        state_n_cells.append(n_unassigned_cells - np.sum(state_n_cells)) # number of cells for last state: total_cells - number of other state cells
        
        # Randomly select cells based on fraction of each state and assign state values
        for (current_state_id, current_state_n_cells) in zip(state_ids, state_n_cells):
            current_state_locs = np.random.choice(np.where(sys_flattened==0)[0], size=current_state_n_cells)
            sys_flattened[current_state_locs] = current_state_id
            
    
    sys_2D_init_state = sys_flattened.reshape(sys_size)
    
    
    return (sys_2D_init_state)
    
    

def f_expand_array_for_bc(sys_array, BC_type, nb_order):
    """expands array beyond boundaries based on boundary condition used

    Args:
        sys_array (numpy array): actual system
        BC_type (str): type of boundary condition. Accepted values:
            - 'p'-periodic
        nb_order (int): neighbor interaction order; 1 means nearest neighbor, 2 mean next-nearest also

    Returns:
        sys_array_bc (numpy array): expanded array based on boundary condition
    """
    
    sys_size = sys_array.shape
    sys_bc_size = (sys_size[0] + 2*nb_order, sys_size[1] + 2*nb_order)
    
    if BC_type.lower() in ["periodic", "p"]:
        sys_array_bc = np.zeros(shape=sys_bc_size) # initiate as all zero
        sys_array_bc[nb_order:-nb_order, nb_order:-nb_order] = sys_array # main
        sys_array_bc[0:nb_order, 0:nb_order] = sys_array[-nb_order:, -nb_order:] # top left 1
        sys_array_bc[0:nb_order, nb_order:-nb_order] = sys_array[-nb_order:, :] # top 2
        sys_array_bc[0:nb_order, -nb_order:] = sys_array[-nb_order:, 0:nb_order] # top right 3
        sys_array_bc[nb_order:-nb_order, -nb_order:] = sys_array[:, -nb_order:] # right 4
        sys_array_bc[-nb_order:, -nb_order:] = sys_array[0:nb_order, 0:nb_order] # right bottom 5
        sys_array_bc[-nb_order:, nb_order:-nb_order] = sys_array[0:nb_order, :] # bottom 6
        sys_array_bc[-nb_order:, 0:nb_order] = sys_array[0:nb_order, -nb_order:] # bottom left 7
        sys_array_bc[nb_order:-nb_order, 0:nb_order] = sys_array[:,  -nb_order:] # left 8

        
    return (sys_array_bc)



def f_2Darray_list_to_gif(array_list, gif_savename, remap_values=True, apply_cmap=True):
    """create a gif from an array list

    Args:
        array_list (list): list of numpy arrays; each array is treated as a frame
        gif_savename (str): filename for output gif (don't specify any extension)
        remap_values (bool): default=True; rescale array values by factor of (255/maxValue)
        apply_cmap (bool): default=True; applies a colormap to images

    Returns:
        None
    """
    
    import cv2
    from PIL import Image

    heatmap_used = cv2.COLORMAP_JET
    max_cell_value = np.max(array_list)
    images = []

    for i in range(0, len(array_list)):
        
        if remap_values == True:
            frame = (array_list[i] / max_cell_value) * 255 # changes to intensity levels corresponding to 8-bit image
        else:
            frame = array_list[i]
        
        frame = Image.fromarray(frame).convert('P') # converts frame from array to image format
        
        # get 3 channel color image by replicating same frame to 3 channels
        frame = cv2.cvtColor(np.asarray(frame), cv2.COLOR_RGB2BGR)
        
        if apply_cmap == True:
            frame = cv2.applyColorMap(frame, heatmap_used) # apply colormap to image
        
        frame = Image.fromarray(frame).convert('P') # converts 3 channel frame to color image
        images.append(frame) # add frame to images list
    
    # create gif
    images[0].save(f"{gif_savename}.gif", save_all=True, append_images=images[1:], optimize=True, duration=100)
    
    
    return None
    
    

def f_evolve_sys(sys_init_state, nb_size, nb_type, BC_type, t_steps):
    """creates initial state of system containing nuclei

    Args:
        sys_init_state (numpy array): initial state of the system
        nb_type (str): type of neighborhood. Accepted values:
            - 'm'-Moore
            - 'vn'-Von Newmann
        BC_type (str): type of boundary condition. Accepted values:
            - 'p'-periodic
        t_steps (int): number of time steps

    Returns:
        sys_init_state (numpy array): initial state of system
        sys_init_grainID (numpy array): initial grainID map of system
    """
    

    nb_order = int((nb_size - 1)/2) # order of neighborhood (nearest neighbour order)
    sys_size = sys_init_state.shape
    state_ids = np.unique(sys_init_state) # unique state ids present in initial system

    sys_store_state = [sys_init_state]

    print("Evolving system: Time step ", end="", flush=True)
    for t in range(1, t_steps + 1):

        print(t, end=" ", flush=True)

        sys_old_state = sys_store_state[t-1]
        if BC_type.lower() in ["periodic", "p"]:
            sys_old_state_bc = f_expand_array_for_bc(sys_old_state, BC_type, nb_order)

        sys_new_state_bc = np.zeros(shape=sys_old_state_bc.shape, dtype=int)

        if nb_type.lower() in ["m", "moore"]:
            for Ri in range(nb_order, sys_size[0] + nb_order):
                for Ci in range(nb_order, sys_size[1] + nb_order):

                    #----------------
                    # get neighborhood (state)
                    nb_old_state = sys_old_state_bc[Ri-nb_order:Ri+nb_order+1, Ci-nb_order:Ci+nb_order+1]

                    nb_cell_state_counts = [np.count_nonzero(nb_old_state == i) for i in state_ids]

                    # Assign the state ID that has maximum count in neighborhood
                    sys_new_state_bc[Ri, Ci] = state_ids[np.argmax(nb_cell_state_counts)]

            sys_new_state = sys_new_state_bc[nb_order:-nb_order, nb_order:-nb_order]

            sys_store_state.append(sys_new_state)
            
            
    return (sys_store_state)