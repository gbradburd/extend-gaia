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

getAccOut <- function(ts,outPrefix,SIGMA=0.8,REP=9){
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
  
  
  
  e = sqrt(rowSums((map_x[is_ancestor, ] - ancestor_locations[,2:3])^2)) / max(dist(sample_locations[,2:3]))
  #e = sqrt(rowSums((map_x[-(1:200),] - ancestor_locations[,2:3])^2)) / max(dist(sample_locations[,2:3]))

  dist_from_sample_centroid0 = sqrt(
    rowSums(sweep(ancestor_locations[, 2:3], 2, sample_centroid)^2))

  dist_from_sample_centroid = sqrt(
    #rowSums(sweep(map_x[-(1:200), ], 2, sample_centroid)^2))
    rowSums(sweep(map_x[is_ancestor, ], 2, sample_centroid)^2))

  ans = data.frame(
    as.numeric(SIGMA),
    as.integer(REP),
    internal_nodes-1L,
    nodes[internal_nodes, "time"],
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

SIGMA = 0.8
REP = 9
TREESEQ <- file.path("C:/Users/islar/OneDrive/Documents/bradburd_lab", sprintf("tree-S%s-R%s.trees", SIGMA, REP))

#TTREESEQ = file.path(getwd(), "../simulations/trees", sprintf("tree-S%s-R%s.trees", SIGMA, REP))
#hard coded the name of the tree file im looking at rn 

ts = treeseq_load(TREESEQ)
nodes = treeseq_nodes(ts)

idx = match(sample(unique(nodes$individual_id[nodes$is_sample == 1L]), 100L),
            nodes$individual_id)

#ts2 = treeseq_simplify(ts, nodes$node_id[c(rbind(idx, idx+1L))])
#this is where it simplifies 


extend = "yes" #just so i can control when this runs
if (extend == "yes"){
  treeseq_write(ts, "temp.trees") #bc cant extend on gaia object
  
  # changing ts2 to ts to run extension on non-simplified tree 
  
  temp <- tskit$load("temp.trees") # loading the temporary simplified tree file
  
  extended_temp   <- temp$extend_haplotypes() # extending the temp simplified tree
  
  extended_temp$dump("temp_extended.trees") # putting extened back into temp 
  tsex <- treeseq_load("temp_extended.trees")
  
  #ts2ex = treeseq_simplify(tsex, nodes$node_id[c(rbind(idx, idx+1L))])
  # trying simplifying after extending 
 
}

# for the extended tree:
getAccOut(tsex,outPrefix="ext")
# tsex for non-simplified version
# ts2ex for simplified after extending 


# for the non-extended tree:
getAccOut(ts,outPrefix="normie") 
#  ts2 = simplified
# ts = non-simplified version 


