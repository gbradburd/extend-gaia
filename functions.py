import msprime, tskit
import numpy as np
import gaiapy as gp
import pandas as pd
from scipy.spatial.distance import pdist, squareform

#function that takes the tree sequence and returns a array with the columns 
# node id | wether its a sample or ancestor | x coordinate | y coordinate 
def locations(ts):
    nodes = ts.nodes()
    locs_array = []
    for node in nodes:
        if node.individual != -1:
            ind = ts.individual(node.individual)
            x = ind.location[0]
            y = ind.location[1]
            is_sample = node.is_sample()
            locs_array.append(node.id)
            locs_array.append(is_sample)
            locs_array.append(x)
            locs_array.append(y)
        
    locs_array = np.array(locs_array)
    locs = locs_array.reshape(-1, 4)
    return locs


#function that takes the tree sequence and finds the unary nodes in the tree sequence
# returns a list of all the unary indices by node id 
def findUnary(ts):
    unary_nodes = np.zeros(ts.num_nodes) # binary vector specifying if a node is unary or not anywhere on the tree sequence
    for tree in ts.trees():
        num_children = tree.num_children_array[:ts.num_nodes]
        is_unary = num_children == 1
        for i, condition in enumerate(is_unary):
            if is_unary[i] == True:
                unary_nodes[i] = 1
    mask = unary_nodes == 1
    #mask, unary_nodes[mask]
    unary_list = np.where(mask)
    unary_indices = unary_list[0]
    return unary_indices


#function for testing small tree sequences that dont have associated location data
# cannot use the location function, so takes in made up location data from samples, ancestors, everything
def small_getAccOut(ts, unary_indices, samples, ancestors, everything, outPrefix):
    # tqdm - package to test timing of things 

    sample_locations = samples # comment this out when not usign test data set 
    ancestor_locations = ancestors

    sample_centroid = np.mean(sample_locations[:, 1:], axis=0)
    # print(sample_centroid)
    mpr = gp.quadratic_mpr(ts, sample_locations)
    map_x = gp.quadratic_mpr_minimize(mpr)
    
    is_unary = np.isin(everything[:, 0], unary_indices) # is_unary = np.isin(locs[:, 0], unary_indices)
    is_ancestor = np.isin(everything[:, 0], ancestor_locations[:, 0])
    target_nodes = is_unary & is_ancestor 
   

    target_locations = everything[target_nodes][:, [0, 1, 2]]

    e = np.sqrt(np.sum((map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    
    dist_from_sample_centroid0 = np.sqrt(np.sum((target_locations[:, 1:3] - sample_centroid)**2, axis=1))
    dist_from_sample_centroid = np.sqrt(np.sum((map_x[target_nodes] - sample_centroid)**2, axis=1))
    
    target_node_ids = everything[target_nodes, 0]
    target_node_times = ts.nodes_time[target_nodes]

    ancestor_df = pd.DataFrame({
        'node_id': target_node_ids,
        'node_time': target_node_times,
        'error': e,
        'dist_from_sample_centroid0': dist_from_sample_centroid0,
        'dist_from_sample_centroid': dist_from_sample_centroid
    })

    ancestor_df.to_csv(f"{outPrefix}_ancestor-estimates.csv",
                      mode='a', header=False, index=False)
    
    return ancestor_df



# the full version of getAccOut that takes in large tree sequences with associated location data
# so locations can be used and made up location data arrays dont need to be passed in 
# also includes the sigma and rep in the returned data frame so the csv file can be 
# used to generate heatmaps from the R files 
def big_getAccOut(ts, unary_indices, outPrefix, sigma, rep):
    # tqdm - package to test timing of things 

    locs = locations(ts)
    sample_locations = locs[locs[:, 1] == 1][:, [0, 2, 3]]
    ancestor_locations = locs[locs[:, 1] != 1][:, [0, 2, 3]]


    sample_centroid = np.mean(sample_locations[:, 1:], axis=0)
    # print(sample_centroid)
    mpr = gp.quadratic_mpr(ts, sample_locations)
    map_x = gp.quadratic_mpr_minimize(mpr)
    
    is_unary = np.isin(locs[:, 0], unary_indices)
    is_ancestor = np.isin(locs[:, 0], unary_indices)
    target_nodes = is_unary & is_ancestor 
   
    target_locations = locs[target_nodes][:, [0, 2, 3]]

    e = np.sqrt(np.sum((map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    
    dist_from_sample_centroid0 = np.sqrt(np.sum((target_locations[:, 1:3] - sample_centroid)**2, axis=1))
    dist_from_sample_centroid = np.sqrt(np.sum((map_x[target_nodes] - sample_centroid)**2, axis=1))
    
    target_node_ids = locs[target_nodes, 0]
    target_node_times = ts.nodes_time[target_nodes]

    ancestor_df = pd.DataFrame({
        'sigma': np.repeat(sigma, target_node_ids.shape[0]),
        'rep': np.repeat(rep, target_node_ids.shape[0]),
        'node_id': target_node_ids,
        'node_time': target_node_times,
        'error': e,
        'dist_from_sample_centroid0': dist_from_sample_centroid0,
        'dist_from_sample_centroid': dist_from_sample_centroid
    })

    ancestor_df.to_csv(f"{outPrefix}_ancestor-estimates.csv",
                      mode='a', header=False, index=False)
    
    return ancestor_df