######################

# Network Simulation

######################

# install.packages("~/Downloads/sharp-1.0-beta.tar", repos = NULL, type="source")
# install.packages("sharp")
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("graphite")

library(clusterProfiler)
library(igraph)
library(tidyverse)
library(graphite)
library(sharp)
library(pheatmap)
library(tidygraph)
library(focus)

######################

#' A. [Network Topolgy Bias study] Simulation of synthetic networks with random pathway topology 
#' varying in the size of pathways and connectivity between pathways

pk1 = c(20,20,20,20) # m = 4 small pathways 
pk2 = c(50,50,50,50) # m = 4 large pathways

nu_1 = 0 # disconnected pathways
nu_2 = 0.01 # sparsely connected pathways
nu_3 = 0.05 # highly connected pathways

SimulateTopology <- function(k, n, pk, nu){
  
    
    sim_network_list = NULL
    
    for (i in 1:k){
      
      # simulate network
      
      sim_network <- sharp::SimulateGraphical(n = n, pk = pk, topology = "random", nu_within = 0.1, nu_between = nu)    
      sim_network_list[[i]] <- sim_network
      
    }
    
    return(sim_network_list)
}

set.seed(2222)
small_dis <- SimulateTopology(100, 1000, pk1, nu_1)
set.seed(2222)
small_sparse <- SimulateTopology(100, 1000, pk1, nu_2)
set.seed(2222)
small_high <- SimulateTopology(100, 1000, pk1, nu_3)

set.seed(2222)
large_dis <- SimulateTopology(100, 1000, pk2, nu_1)
set.seed(2222)
large_sparse <- SimulateTopology(100, 1000, pk2, nu_2)
set.seed(2222)
large_high <- SimulateTopology(100, 1000, pk2, nu_3)


small_topology <- list(small_dis, small_sparse, small_high)
save(small_topology, file = "../data/small_topology.RData")
large_topology <- list(large_dis,large_sparse, large_high)
save(large_topology, file = "../data/large_topology.RData")







