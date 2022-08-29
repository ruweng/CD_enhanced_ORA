######################

# Enrichment Graphical Model
# topology_small_high

######################

library(igraph)
library(tidyverse)
library(graphite)
library(pheatmap)
library(tidygraph)
library(focus)
library(foreach)
library(doParallel)

#########

sharp::Graph
load("../data/small_topology.RData")

sim_small_top_high <- small_topology[[3]]

data_list = NULL
theta_list = NULL

for (i in 1:length(sim_small_top_high)){
  data_list[[i]] = sim_small_top_high[[i]]$data
  theta_list[[i]] = sim_small_top_high[[i]]$theta
}


###

parallel:::setDefaultClusterOptions(setup_strategy = "parallel")
n.cores <- parallel::detectCores()

my.cluster <- parallel::makeCluster(
  30, 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)

clusterEvalQ(my.cluster, {
  library(focus)
  library(igraph)
})

func <- function(x) {
  stability <- GraphicalModel(xdata = x)
  return(stability)}
clusterExport(cl = my.cluster, c("data_list", "theta_list", "func"))


stab <- parLapply(cl = my.cluster, data_list, fun = func)
adj_graphical <- parLapply(cl = my.cluster, stab, fun = Adjacency)
network_graphical <- parLapply(cl = my.cluster, adj_graphical, fun = Graph)

output <- list ("data" = data_list,
                "theta" = theta_list,
                "stability_selection" = stab,
                "adj_graphical" = adj_graphical,
                "network_graphical" = network_graphical
)


save(output, file = "../data/Graphical_top_small_high.RData")