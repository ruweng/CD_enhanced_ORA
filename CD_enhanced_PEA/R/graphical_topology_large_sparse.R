######################

# Enrichment Graphical Model
# topology_large_sparse

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

load("../data/large_topology.RData")

sim_large_top_sparse <- large_topology[[2]]

data_list = NULL
theta_list = NULL

for (i in 1:length(sim_large_top_sparse)){
  data_list[[i]] = sim_large_top_sparse[[i]]$data
  theta_list[[i]] = sim_large_top_sparse[[i]]$theta
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


save(output, file = "../data/Graphical_top_large_sparse.RData")