import msprime, tskit
import numpy as np
import gaiapy as gp
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

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

def node_spans(ts, include_missing=False):
    """
    Returns the array of "node spans", i.e., the `j`th entry gives
    the total span over which node `j` is in the tree sequence.
    Sample nodes that are isolated are "missing data"; inclusion
    of these spans are controlled by `include_missing`. (If
    `include_missing` is `True` then the span of each sample is
    always equal to the sequence length.)

    :param bool include_missing: Whether to include spans of nodes
        on which they have missing data.
    """
    child_spans = np.bincount(
        ts.edges_child,
        weights=ts.edges_right - ts.edges_left,
        minlength=ts.num_nodes,
    )
    for t in ts.trees():
        span = t.span
        for r in t.roots:
            # do this check to exempt 'missing data'
            if include_missing or (t.num_children(r) > 0):
                child_spans[r] += span
    return child_spans

# here ts is the simplified tree 
def get_span_stats(ts, ets):
    added_span = np.zeros(ets.num_nodes)
    wrong_added_span = np.zeros(ets.num_nodes)
    
    ts_node_ids = set(n.id for n in ts.nodes())  # all valid node IDs in original ts

    for interval, t, et in ts.coiterate(ets):
        interval_length = interval[1] - interval[0]
        t_nodes = set(t.nodes())  # node IDs present in this specific tree
        for n in et.nodes():
            if et.num_children(n) == 1:
                added_span[n] += interval_length
            if n not in t_nodes:
                assert et.num_children(n) == 1
                wrong_added_span[n] += interval_length

    return added_span, wrong_added_span

#here ts is the unsimplified tree
def get_node_stats(ts):
    #nodes = list(ts.nodes())
    nodes = np.array(list(ts.nodes()))
    node_ids = np.array([n.id for n in nodes])
    # Find the samples
    #is_sample = np.asarray(np.isin(nodes, ts.samples()), dtype=int)
    is_sample = np.asarray(np.isin(node_ids, ts.samples()), dtype=int)
    # Find all the other things (this requires checking tree by tree)
    tree = ts.first()
    start = tree.interval[0] 
    end = ts.sequence_length
    num_children = np.zeros(nodes.shape[0])
    num_parents = np.zeros(nodes.shape[0])
    distinct_children, distinct_parents, distinct_populations = list(), list(), list()
    is_root = np.zeros(nodes.shape[0])
    #for i,node in tqdm(enumerate(nodes)):
    for i, node_id in tqdm(enumerate(node_ids)): 
        tree.seek(start)
        children = list()
        parents = list()
        is_root[i] = tree.is_root(node_id)
        node_children = list(tree.children(node_id))

        children.extend(node_children)
        parents.append(tree.parent(node_id))

        w = (None, tree.interval[0])
        while w[1] < end and tree.next():
            w = (w[1], min(tree.interval[1], end))
            is_root[i] = tree.is_root(node_id)
            node_children = list(tree.children(node_id))

            children.extend(node_children)
            parents.append(tree.parent(node_id))

        distinct_parents.append(np.unique(parents))
        distinct_children.append(np.unique(children))
        num_children[i] = np.unique(children).shape[0]
        num_parents[i] = np.unique(parents).shape[0]
    
    data_dict = {
        'id': node_ids,
        'num_children': num_children,
        'distinct_children': distinct_children,
        'distinct_parents': distinct_parents,
        'num_parents': num_parents,
        'is_sample': is_sample,
        'is_root': is_root,
    }
    return pd.DataFrame(data_dict)


def getAccOut(ts, outPrefix, sigma, rep):
    # tqdm - package to test timing of things 
    simp = ts.simplify()
    ets = simp.extend_haplotypes()
    unary_indices = findUnary(ets)

    # use simplified tree
    total_added_span, wrongly_added_span = get_span_stats(simp, ets)
    # use original tree 
    node_stats_df = get_node_stats(ts)
    
    simp_spans = node_spans(simp)
    ets_spans = node_spans(ets)

    simp_num_trees = simp.num_trees
    ets_num_trees = ets.num_trees
    simp_num_edges = simp.num_edges
    ets_num_edges = ets.num_edges

    locs = locations(ets)
    sample_locations = locs[locs[:, 1] == 1][:, [0, 2, 3]]
    ancestor_locations = locs[locs[:, 1] != 1][:, [0, 2, 3]]


    sample_centroid = np.mean(sample_locations[:, 1:], axis=0)
    # print(sample_centroid)

    simp_mpr = gp.quadratic_mpr(simp, sample_locations)
    simp_map_x = gp.quadratic_mpr_minimize(simp_mpr)

    ets_mpr = gp.quadratic_mpr(ets, sample_locations)
    ets_map_x = gp.quadratic_mpr_minimize(ets_mpr)
    
    is_unary = np.isin(locs[:, 0], unary_indices)
    is_ancestor = np.isin(locs[:, 0], unary_indices)
    target_nodes = is_unary & is_ancestor 
    target_locations = locs[target_nodes][:, [0, 2, 3]]

    simp_e = np.sqrt(np.sum((simp_map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    ets_e = np.sqrt(np.sum((ets_map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    
    dist_from_sample_centroid0 = np.sqrt(np.sum((target_locations[:, 1:3] - sample_centroid)**2, axis=1))

    simp_dist_from_sample_centroid = np.sqrt(np.sum((simp_map_x[target_nodes] - sample_centroid)**2, axis=1))
    ets_dist_from_sample_centroid = np.sqrt(np.sum((ets_map_x[target_nodes] - sample_centroid)**2, axis=1))
    
    target_node_ids = locs[target_nodes, 0]
    target_node_times = ets.nodes_time[target_nodes]

    target_simp_spans = simp_spans[target_nodes]
    target_ets_spans = ets_spans[target_nodes]

    target_total_added_span = total_added_span[target_nodes]
    target_wrong_added_span = wrongly_added_span[target_nodes]
   
    

    ancestor_df = pd.DataFrame({
        'node_id': target_node_ids,
        'node_time': target_node_times,
        'simp_error': simp_e,
        'ets_error': ets_e,
        'simp_span': target_simp_spans,
        'ets_span': target_ets_spans,
        'added_span': target_total_added_span,
        'wrongly_added_span': target_wrong_added_span,
        'dist_from_sample_centroid0': dist_from_sample_centroid0,
        'simp_dist_from_sample_centroid': simp_dist_from_sample_centroid,
        'ets_dist_from_sample_centroid': ets_dist_from_sample_centroid
    })

    static_df = pd.DataFrame({
        'sigma': [sigma],
        'rep': [rep],
        'simp_num_trees': [simp_num_trees],
        'ets_num_trees': [ets_num_trees], 
        'simp_num_edges': [simp_num_edges],
        'ets_num_edges': [ets_num_edges]
    })



    ancestor_df.to_csv(f"{outPrefix}_results.csv",
                      mode='a', header=True, index=False)
    
    static_df.to_csv(f"{outPrefix}_static_info.csv",
                      mode='a', header=True, index=False)
    
    node_stats_df.to_csv(f"{outPrefix}_node_stats.csv",
                      mode='a', header=True, index=False)

    
    return ancestor_df, static_df, node_stats_df

ts = tskit.load("/home/islar/bradburdlab/tree_project/unary_project/tree-S0.5-R3.trees")
sigma = 0.5
rep = 3
getAccOut(ts, "tree-S0.5-R3_new", sigma, rep)

# could edit locs so it only creates the location array for nodes that are unary?






# /*
# more nodes now than ancient nodes 
# start with linear disection of nodes - take like max time / however many bins want 

# some percentage of sewquence length - liked iveide by ten 

# */


# figure out pllotly !! look at plotly!! 




# Look at direct descendants and ancestors of nodes about unary nodes before and after extension 
# (look at parent and child errors of unary NotADirectoryError 
                                                                                                
# for marginal trees nodes could have differnet parents children 
# overtrees 

# preprocessing step ? for ev

# have column in csv thats immediate descendents and immediate ancestors 
# - those entries arent gonna be a single nodes - like a list of nodes

# prepost immediate ancestors and descendents - checking correctness - tscompare 

# tie percentage of wrong span addded to every node add to csv - use code halley sent 



#function for testing small tree sequences that dont have associated location data
# cannot use the location function, so takes in made up location data from samples, ancestors, everything
def small_getAccOut(ts, samples, ancestors, everything):
    simp = ts.simplify()
    ext = ts.extend_haplotypes()
    unary_indices = findUnary(ext)
    print("1")

    total_added_span, wrongly_added_span = get_span_stats(ts, ext)

    print("1.5")

    simp_spans = node_spans(simp)
    ext_spans = node_spans(ext)

    print("2")

    simp_num_trees = simp.num_trees
    ext_num_trees = ext.num_trees
    simp_num_edges = simp.num_edges
    ext_num_edges = ext.num_edges

    print("3")
    # locs = locations(ext)
    sample_locations = samples # comment this out when not usign test data set 
    ancestor_locations = ancestors

    sample_centroid = np.mean(sample_locations[:, 1:], axis=0)
    # print(sample_centroid)

    print("4")

    simp_mpr = gp.quadratic_mpr(simp, sample_locations)
    simp_map_x = gp.quadratic_mpr_minimize(simp_mpr)

    print("5")

    ext_mpr = gp.quadratic_mpr(ext, sample_locations)
    ext_map_x = gp.quadratic_mpr_minimize(ext_mpr)

    print("6")

    is_unary = np.isin(everything[:, 0], unary_indices)
    is_ancestor = np.isin(everything[:, 0], unary_indices)
    target_nodes = is_unary & is_ancestor 
    target_locations = everything[target_nodes][:, [0, 1, 2]]

    print("7")

    simp_e = np.sqrt(np.sum((simp_map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    ext_e = np.sqrt(np.sum((ext_map_x[target_nodes] - target_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))

    print("8")

    dist_from_sample_centroid0 = np.sqrt(np.sum((target_locations[:, 1:3] - sample_centroid)**2, axis=1))

    simp_dist_from_sample_centroid = np.sqrt(np.sum((simp_map_x[target_nodes] - sample_centroid)**2, axis=1))
    ext_dist_from_sample_centroid = np.sqrt(np.sum((ext_map_x[target_nodes] - sample_centroid)**2, axis=1))

    print("9")

    target_node_ids = everything[target_nodes, 0]
    target_node_times = ext.nodes_time[target_nodes]

    print("10")

    target_simp_spans = simp_spans[target_nodes]
    target_ext_spans = ext_spans[target_nodes]

    print("11")

    target_total_added_span = total_added_span[target_nodes]
    target_wrong_added_span = wrongly_added_span[target_nodes]

    print("12")


    ancestor_df = pd.DataFrame({
        'node_id': target_node_ids,
        'node_time': target_node_times,
        'simp_error': simp_e,
        'ext_error': ext_e,
        'simp_span': target_simp_spans,
        'ext_span': target_ext_spans,
        'added_span': target_total_added_span,
        'wrongly_added_span': target_wrong_added_span,
        'dist_from_sample_centroid0': dist_from_sample_centroid0,
        'simp_dist_from_sample_centroid': simp_dist_from_sample_centroid,
        'ext_dist_from_sample_centroid': ext_dist_from_sample_centroid
    })

    print("13")

    static_df = pd.DataFrame({
        'simp_num_trees': [simp_num_trees],
        'ext_num_trees': [ext_num_trees], 
        'simp_num_edges': [simp_num_edges],
        'ext_num_edges': [ext_num_edges]
    })
    print("14")

    # ancestor_df.to_csv(f"shits.csv",
    #                   mode='a', header=True, index=False)
    
    # static_df.to_csv(f"shits_static_info.csv",
    #                   mode='a', header=True, index=False)
    
    return ancestor_df, static_df



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



