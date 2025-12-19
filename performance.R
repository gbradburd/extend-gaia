library(gaia)
library(reticulate)
#py_install("tskit")
tskit <- reticulate::import("tskit")

locations = function(ts)
{
  ind = treeseq_individuals(ts)
  nodes = treeseq_nodes(ts)
  locs = t(apply(ind, 1, function(i) {
    c(i$individual_id, i$location[1:2])
  })) #locs
  locs = cbind(
    nodes$node_id
    , nodes$is_sample
    , locs[match(nodes$individual_id, locs[, 1]), -1]
  )
  colnames(locs) = c("node_id","is_sample","x","y")
  locs
} #locations

getAccOut <- function(ts,outPrefix,SIGMA=0.2,REP=0, only_unary=FALSE, original_ts=NULL){
  # og was SIGMA=0.025,REP= 1
  nodes = treeseq_nodes(ts)
  edges = treeseq_edges(ts)
  inds = treeseq_individuals(ts)
  internal_nodes = which(nodes$is_sample != 1L)
  locs = locations(ts)
  sample_locations = locs[locs[,2] == 1, c(1,3,4)]
  ancestor_locations = locs[locs[,2] != 1,c(1,3,4)]
  sample_centroid = colMeans(sample_locations[,2:3])
  mpr = treeseq_quadratic_mpr(ts, sample_locations, TRUE)
  mean_map_rate = sqrt(mpr[[1]] / 2)
  effective_rate = sqrt(mean(apply(data.matrix(edges), 1L, function(e) {
    parent = nodes[e[4]+1, 5]
    child = nodes[e[5]+1, 5]
    len = nodes[e[4]+1, 3] - nodes[e[5]+1, 3]
    from = inds[parent+1, ]$location[1:2]
    to = inds[child+1, ]$location[1:2]
    sum((from - to)^2) / len
  })) / 2)
  ans = data.frame(
    as.numeric(SIGMA),
    as.integer(REP),
    mean_map_rate,
    effective_rate
  )
  write.table(
    ans,
    file=paste0(outPrefix,"_rate-estimates.csv",sep=""),
    sep=",",
    col.names=FALSE,
    row.names=FALSE,
    append=TRUE
  )
  map_x = treeseq_quadratic_mpr_minimize(mpr)
  
  is_ancestor = locs[, "is_sample"] != 1 
  #find the rows that are ancestors 
  
  is_unary = locs[, "node_id"] %in% unary_indices
  
  if (sum(is_unary) != 0) {
    target_nodes = is_unary & is_ancestor 
  }
  else {
    target_nodes = is_ancestor
  }
  
  ancestor_locations = locs[target_nodes, c(1,3,4)] # new ancestor locations  
  
  e = sqrt(rowSums((map_x[target_nodes, ] - ancestor_locations[,2:3])^2)) / max(dist(sample_locations[,2:3]))

  dist_from_sample_centroid0 = sqrt(rowSums(sweep(ancestor_locations[, 2:3], 2, sample_centroid)^2))

  dist_from_sample_centroid = sqrt(rowSums(sweep(map_x[target_nodes, ], 2, sample_centroid)^2))
  
  target_node_ids = nodes$node_id[target_nodes] - 1L
  target_node_times = nodes$time[target_nodes]
  
  ans = data.frame(
    as.numeric(SIGMA),
    as.integer(REP),
    target_node_ids, # internal_nodes-1L,
    target_node_times, # nodes[internal_nodes, "time"], 
    e,
    dist_from_sample_centroid0,
    dist_from_sample_centroid
  )

  write.table(
    ans,
    file=paste0(outPrefix,"_ancestor-estimates.csv",sep=""),
    sep=",",
    col.names=FALSE,
    row.names=FALSE,
    append=TRUE
  )
  
  
}


SIGMA = 0.2
REP = 0
#TREESEQ <- file.path("C:/Users/islar/OneDrive/Documents/bradburd_lab", sprintf("tree-S%s-R%s.trees", SIGMA, REP))

EXTTREESEQ <- file.path("C:/Users/islar/OneDrive/Documents/bradburd_lab/extended-tree-S0.2-R0.trees")
SIMPTREESEQ <- file.path("C:/Users/islar/OneDrive/Documents/bradburd_lab/simp-tree-S0.2-R0.trees")
unary_indices <- read.csv('unary_indices.csv', header=FALSE)$V1

#TTREESEQ = file.path(getwd(), "../simulations/trees", sprintf("tree-S%s-R%s.trees", SIGMA, REP))
#hard coded the name of the tree file im looking at rn 

ts = treeseq_load(SIMPTREESEQ)
tsex = treeseq_load(EXTTREESEQ)

nodes = treeseq_nodes(ts) ## is this correct????


idx = match(sample(unique(nodes$individual_id[nodes$is_sample == 1L]), 100L),
            nodes$individual_id)


# for the extended tree:
getAccOut(tsex,outPrefix="ext")
# tsex for non-simplified version
# ts2ex for simplified after extending 


# for the non-extended tree:
getAccOut(ts,outPrefix="normie") 
#  ts2 = simplified
# ts = non-simplified version 


