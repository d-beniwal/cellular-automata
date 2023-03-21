import numpy as np

def f_sys_initialize(sys_size, pos, n_seed, shape, size, min_spacing=1.5):
    """creates initial state of system containing nuclei

    Args:
        sys_size (tuple): (height, width) of system
        pos (str): position where nuclei have to be initiated. Possible values:
            - 'r'-random; nuclei generated randomly
            - 'c'-centre; single nuclei at centre of system
        n_seed (int): number of nuclei (fixed to 1 if pos='c')
        shape (str): shape of nuclei. Possible values:
            - 'c'-circular
            - 's'-square
        size (int): size (radius/side length) of nuclei in pixel
        min_spacing (float, optional): minimum spacing between centres of nuclei.
            It represents multiple of size i.e. '1' means nuclei can just touch. Defaults to 1.5.

    Returns:
        sys_init_state (numpy array): initial state of system
        sys_init_grainID (numpy array): initial grainID map of system
    """
    
    if pos == "c":
        n_seed = 1
        
    sys_init_state = np.zeros(shape=sys_size) # cell state: 0-liq; 1-solid
    sys_init_grainID = np.zeros(shape=sys_size) # grain id: 0-liq; [1, n_seed]-solid
    
    (H, W) = sys_size
    grain_ids = np.arange(1, n_seed+1, step=1)
    
    if pos.lower() in ["c", "centre", "center"]:
        loc_R, loc_C = [int(H/2)], [int(W/2)]
    
    if pos.lower() in ["r", "random"]:

        # D_threshold controls that minimum spacing maintained between all seeds
        D_threshold = min_spacing*(2*size)
        D = D_threshold - 1 # initialize D value
        
        while D < D_threshold:
            # generate (row, column) locations for seeds
            loc_R = np.random.randint(low=size+1, high=H-size-1, size=(n_seed))
            loc_C = np.random.randint(low=size+1, high=W-size-1, size=(n_seed))

            D_list = []
            # Calculate distance between all seeds and store in D_list
            for i in range(0, n_seed-1):
                for j in range(i+1, n_seed):
                    Dij = ((loc_R[i] - loc_R[j])**2 + (loc_C[i] - loc_C[j])**2)**(0.5)
                    D_list.append(Dij)
            
            D = np.min(np.array(D_list))
            
    if shape.lower() in ["s", "sq", "square"]:
        for i in range(0, n_seed):
            # Change grain ID of cells corresponding to seed
            sys_init_grainID[(loc_R[i]-size):(loc_R[i]+size+1) ,
                            (loc_C[i]-size):(loc_C[i]+size+1)] = grain_ids[i]
        
            # Change state of cell from 0 to 1
            sys_init_state[(loc_R[i]-size):(loc_R[i]+size+1) ,
                            (loc_C[i]-size):(loc_C[i]+size+1)] = 1
    
    if shape.lower() in ["c", "circle"]:
        for i in range(0, n_seed):
            
            # Defining range to search for creating seed
            R_min, R_max = loc_R[i]-size, loc_R[i]+size+1
            C_min, C_max = loc_C[i]-size, loc_C[i]+size+1

            for R in range(R_min, R_max):
                for C in range(C_min, C_max):
                    dist = ((R - loc_R[i])**2 + (C - loc_C[i])**2)**(0.5)
                    if (dist) < size:
                        sys_init_grainID[R, C] = grain_ids[i]
                        sys_init_state[R, C] = 1
                
            
    return (sys_init_state, sys_init_grainID)



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



def f_evolve_sys(sys_init_state, sys_init_grainID, nb_size, nb_type, BC_type, rule, t_steps):
    """evolves the system over time

    Args:
        sys_init_state (numpy array): _description_
        sys_init_grainID (numpy array): _description_
        nb_size (int): size of neighborhood (odd integer;e.g. 3, 5, 7)
        nb_type (str): type of neighborhood. Accepted values:
            - 'm'-Moore
            - 'vn'-Von Newmann
        BC_type (str): type of boundary condition. Accepted values:
            - 'p'-periodic
        rule (int): minimum neighbors needed to change state
        t_steps (int): number of time steps

    Returns:
        sys_store_state (list): list of system state arrays at each time
        sys_store_grainID (list): list of system grainID arrays at each time
    """
    
    nb_order = int((nb_size - 1)/2) # order of neighborhood (nearest neighbour order)
    sys_size = sys_init_state.shape
    
    sys_store_state = [sys_init_state]
    sys_store_grainID = [sys_init_grainID]
    
    print("Evolving system: Time step ", end="")
    for t in range(1, t_steps + 1):
        
        print(t, end=" ")
        
        sys_old_state = sys_store_state[t-1]
        sys_old_grainID = sys_store_grainID[t-1]
        
        if BC_type.lower() in ["periodic", "p"]:
            sys_old_state_bc = f_expand_array_for_bc(sys_old_state, BC_type, nb_order)
            sys_old_grainID_bc = f_expand_array_for_bc(sys_old_grainID, BC_type, nb_order)


        sys_new_state_bc = np.zeros(shape=sys_old_state_bc.shape, dtype=int)
        sys_new_grainID_bc = np.zeros(shape=sys_old_grainID_bc.shape, dtype=int)


        if nb_type.lower() in ["m", "moore"]:
            for Ri in range(nb_order, sys_size[0] + nb_order):
                for Ci in range(nb_order, sys_size[1] + nb_order):
                    
                    #----------------
                    # get neighborhood (state and grainID)
                    nb_old_state = sys_old_state_bc[Ri-nb_order:Ri+nb_order+1, Ci-nb_order:Ci+nb_order+1]
                    nb_old_grainID = sys_old_grainID_bc[Ri-nb_order:Ri+nb_order+1, Ci-nb_order:Ci+nb_order+1]
                    
                    #----------------
                    # check if rule follows
                    if np.sum(nb_old_state) >= rule:
                        
                        sys_new_state_bc[Ri, Ci] = 1 # update cell state
                        
                        # get unique grain IDs and their counts in neighborhood
                        gID_unique, gID_counts = np.unique(nb_old_grainID[~np.isin(nb_old_grainID, [0])], return_counts=True)
                        
                        # get list of grain ID(s) that are most frequent in neighborhood
                        gID_possibe = list(gID_unique[np.where(gID_counts==max(gID_counts))[0]])
                        
                        # assign new grain ID randomly from gID_possible
                        sys_new_grainID_bc[Ri, Ci] = np.random.choice(gID_possibe)
                        
                    else:
                        sys_new_state_bc[Ri, Ci] = 0
                        sys_new_grainID_bc[Ri, Ci] = 0
                        
                        
            sys_new_state = sys_new_state_bc[nb_order:-nb_order, nb_order:-nb_order]
            sys_new_grainID = sys_new_grainID_bc[nb_order:-nb_order, nb_order:-nb_order]

            sys_store_state.append(sys_new_state)
            sys_store_grainID.append(sys_new_grainID)
    
    
    return (sys_store_state, sys_store_grainID)
