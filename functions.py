
import msprime, tskit
import numpy as np
import gaiapy as gp
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import time

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

# finds the distinct children of unary nodes and returns their indices 
def find_children(unary_indices, node_stats_df):
    children = set()
    for node_id in unary_indices:
        n = node_stats_df.loc[node_id, 'distinct_children']
        if n is None or len(n) == 0:
            continue
        n = n[n != -1]
        children.update(n)
    return np.array(list(children))


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

def get_span_stats(ts, ets):
    node_map = {}
    added_span = np.zeros(ets.num_nodes)
    wrong_added_span = np.zeros(ets.num_nodes)
    for n in ts.nodes():
        slim_id = n.metadata["slim_id"]
        assert slim_id not in node_map
        node_map[slim_id] = n.id
    for interval, t, et in ts.coiterate(ets):
        interval_length = interval[1] - interval[0]
        t_nodes = list(t.nodes())
        for n in et.nodes():
            if et.num_children(n) == 1:
                added_span[n] += interval_length
            node = ets.node(n)
            on = node_map[node.metadata["slim_id"]]
            if on not in t_nodes:
                assert et.num_children(n) == 1, print(interval, n, et.num_children(n), et.time(n))
                wrong_added_span[n] += interval_length
    return added_span, wrong_added_span


# try deletign the time part - parse on node id 

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

    # use original tree
    total_added_span, wrongly_added_span = get_span_stats(ts, ets)
    # use original tree 
    node_stats_df = get_node_stats(ts)

    # the indices of unary nodes
    unary_indices = findUnary(ets)
    # the indices of the children of unary nodes 
    children = find_children(unary_indices, node_stats_df)

   
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
    simp_quadratic_mpr_start = time.time()
    simp_mpr = gp.quadratic_mpr(simp, sample_locations)
    simp_quadratic_mpr_end = time.time()

    simp_quadratic_mpr_minimize_start = time.time()
    simp_map_x = gp.quadratic_mpr_minimize(simp_mpr)
    simp_quadratic_mpr_minimize_end = time.time()

    ets_quadratic_mpr_start = time.time()
    ets_mpr = gp.quadratic_mpr(ets, sample_locations)
    ets_quadratic_mpr_end = time.time()

    ets_quadratic_mpr_minimize_start = time.time()
    ets_map_x = gp.quadratic_mpr_minimize(ets_mpr)
    ets_quadratic_mpr_minimize_end = time.time()

    simp_quadratic_mpr_runtime = simp_quadratic_mpr_end - simp_quadratic_mpr_start
    simp_quadratic_mpr_minimize_runtime = simp_quadratic_mpr_minimize_end - simp_quadratic_mpr_minimize_start

    ets_quadratic_mpr_runtime = ets_quadratic_mpr_end - ets_quadratic_mpr_start
    ets_quadratic_mpr_minimize_runtime = ets_quadratic_mpr_minimize_end - ets_quadratic_mpr_minimize_start

    
    is_unary = np.isin(locs[:, 0], unary_indices)
    is_ancestor = np.isin(locs[:, 0], unary_indices)
    is_child = np.isin(locs[:,0], children)
    unary_nodes = is_unary & is_ancestor 
    child_nodes = is_child

    
    unary_locations = locs[unary_nodes][:, [0, 2, 3]]
    child_locations = locs[child_nodes][:, [0, 2, 3]]

    unary_simp_e = np.sqrt(np.sum((simp_map_x[unary_nodes] - unary_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    unary_ets_e = np.sqrt(np.sum((ets_map_x[unary_nodes] - unary_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))

    child_simp_e = np.sqrt(np.sum((simp_map_x[child_nodes] - child_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    child_ets_e = np.sqrt(np.sum((ets_map_x[child_nodes] - child_locations[:, 1:2])**2, axis=1)) / np.max(pdist(sample_locations[:, 1:2]))
    
    unary_dist_from_sample_centroid0 = np.sqrt(np.sum((unary_locations[:, 1:3] - sample_centroid)**2, axis=1))
    child_dist_from_sample_centroid0 = np.sqrt(np.sum((child_locations[:, 1:3] - sample_centroid)**2, axis=1))

    unary_simp_dist_from_sample_centroid = np.sqrt(np.sum((simp_map_x[unary_nodes] - sample_centroid)**2, axis=1))
    unary_ets_dist_from_sample_centroid = np.sqrt(np.sum((ets_map_x[unary_nodes] - sample_centroid)**2, axis=1))

    child_simp_dist_from_sample_centroid = np.sqrt(np.sum((simp_map_x[child_nodes] - sample_centroid)**2, axis=1))
    child_ets_dist_from_sample_centroid = np.sqrt(np.sum((ets_map_x[child_nodes] - sample_centroid)**2, axis=1))
    
    unary_node_ids = locs[unary_nodes, 0]
    unary_node_times = ets.nodes_time[unary_nodes]

    child_node_ids = locs[child_nodes, 0]
    child_node_times = ets.nodes_times[child_nodes]


    simp_spans = simp_spans[unary_nodes]
    ets_spans = ets_spans[unary_nodes]

    total_added_span = total_added_span[unary_nodes]
    wrong_added_span = wrongly_added_span[unary_nodes]
   


    ancestor_df = pd.DataFrame({
        'unary_node_id': unary_node_ids,
        'unary_node_time': unary_node_times,
        'child_node_id': child_node_ids,
        'child_node_time': child_node_times,
        'unary_simp_error': unary_simp_e,
        'unary_ets_error': unary_ets_e,
        'child_simp_error': child_simp_e,
        'child_ets_error': child_ets_e,
        'simp_span': simp_spans,
        'ets_span': ets_spans,
        'added_span': total_added_span,
        'wrongly_added_span': wrong_added_span,
        'unary_dist_from_sample_centroid0': unary_dist_from_sample_centroid0,
        'unary_simp_dist_from_sample_centroid': unary_simp_dist_from_sample_centroid,
        'unary_ets_dist_from_sample_centroid': unary_ets_dist_from_sample_centroid,
        'child_dist_from_sample_centroid0': child_dist_from_sample_centroid0,
        'child_simp_dist_from_sample_centroid': child_simp_dist_from_sample_centroid,
        'child_ets_dist_from_sample_centroid': child_ets_dist_from_sample_centroid
    })

    static_df = pd.DataFrame({
        'sigma': [sigma],
        'rep': [rep],
        'simp_num_trees': [simp_num_trees],
        'ets_num_trees': [ets_num_trees], 
        'simp_num_edges': [simp_num_edges],
        'ets_num_edges': [ets_num_edges]
    })

    time_df = pd.DataFrame({
        'simp_quad_mpr': [simp_quadratic_mpr_runtime], 
        'simp_minimize': [simp_quadratic_mpr_minimize_runtime],
        'ets_quad_mpr': [ets_quadratic_mpr_runtime],
        'ets_minimize': [ets_quadratic_mpr_minimize_runtime]
    })


    ancestor_df.to_csv(f"{outPrefix}_results.csv",
                      mode='a', header=True, index=False)
    
    static_df.to_csv(f"{outPrefix}_static_info.csv",
                      mode='a', header=True, index=False)
    
    time_df.to_csv(f"{outPrefix}_time.csv",
                   mode='a', header=True, index=False)
    
    node_stats_df.to_csv(f"{outPrefix}_node_stats.csv",
                      mode='a', header=True, index=False)

    return ancestor_df, static_df, time_df, node_stats_df