library(reticulate)
library(gaia)
use_virtualenv("~/home/islar/bradburdlab/tree_project/extend-gaia/bioenv", required = TRUE)
tskit <- import("tskit")

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

timing = function(file_path)
{
    path = file.path(file_path)
    ts = treeseq_load(path)

    # all this bc cannot extend on a  gaia tree object so must write to disk to extend then back to gaia tree 
    treeseq_write(ts, "temp.trees")
    temp <- tskit$load("temp.trees")
    simplified_temp <- temp$simplify()
    extended_temp <- simplified_temp$extend_haplotypes()
    simplified_temp$dump("temp_simplified.trees")
    extended_temp$dump("temp_extended.trees")
    ets <- treeseq_load("temp_extended.trees")
    sts <- treeseq_load("temp_simplified.trees")


    locs = locations(ets)
    sample_locations = locs[locs[,2] == 1, c(1,3,4)]

    sts_start = Sys.time()
    mpr = treeseq_quadratic_mpr(sts, sample_locations, TRUE)
    sts_stop = Sys.time()
    sts_elapsed = unclass(sts_stop - sts_start)[1]


    ets_start = Sys.time()
    mpr = treeseq_quadratic_mpr(ets, sample_locations, TRUE)
    ets_stop = Sys.time()
    ets_elapsed = unclass(ets_stop - ets_start)[1]

    R_timings_




}# timing 

path1 = "tree-files/tree-S0.2-R0.trees"
timing(path1)






ets_start = Sys.time()
mpr = treeseq_quadratic_mpr(ets_S02_R0, sample_locations, TRUE)
ets_stop = Sys.time()
ets_elapsed = unclass(stop - start)[1]



# sts_S02_R0 = treeseq_simplify(S02_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S05_R0 = treeseq_simplify(S05_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S08_R0 = treeseq_simplify(S08_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S11_R0 = treeseq_simplify(S11_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S14_R0 = treeseq_simplify(S14_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S17_R0 = treeseq_simplify(S17_R0, nodes$node_id[c(rbind(idx, idx+1L))])
# sts_S20_R0 = treeseq_simplify(S20_R0, nodes$node_id[c(rbind(idx, idx+1L))])




S02_R0_path = file.path("tree-files/tree-S0.2-R0.trees")
# S05_R0_path = file.path("tree-files/tree-S0.5-R0.trees")
# S08_R0_path = file.path("tree-files/tree-S0.8-R0.trees")
# S11_R0_path = file.path("tree-files/tree-S1.1-R0.trees")
# S14_R0_path = file.path("tree-files/tree-S1.4-R0.trees")
# S17_R0_path = file.path("tree-files/tree-S1.7-R0.trees")
# S20_R0_path = file.path("tree-files/tree-S2.0-R0.trees")

S02_R0 = treeseq_load(S02_R0_path)
# S05_R0 = treeseq_load(S05_R0_path)
# S08_R0 = treeseq_load(S08_R0_path)
# S11_R0 = treeseq_load(S11_R0_path)
# S14_R0 = treeseq_load(S14_R0_path)
# S17_R0 = treeseq_load(S17_R0_path)
# S20_R0 = treeseq_load(S20_R0_path)


treeseq_write(S02_R0, "temp.trees")
temp <- tskit$load("temp.trees")
simplified_temp <- temp$simplify()
extended_temp <- simplified_temp$extend_haplotypes()
simplified_temp$dump("temp_simplified.trees")
extended_temp$dump("temp_extended.trees")
ets_S02_R0 <- treeseq_load("temp_extended.trees")
sts_S02_R0 <- treeseq_load("temp_simplified.trees")
