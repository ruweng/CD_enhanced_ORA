######################

# Enrichment Graphical Model

######################

# install.packages("~/Downloads/sharp-1.0-beta.tar", repos = NULL, type="source")
# install.packages("sharp")
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("graphite")
# library(clusterProfiler)
# library(sharp)

library(igraph)
library(tidyverse)
library(graphite)
library(pheatmap)
library(tidygraph)
library(focus)
library(foreach)
library(doParallel)

####################


EnrichGraphicalModel <- function(x){
  
  ## Function to create networks using stability enhanced Graphical LASSO 
  # from EnrichmentSimulation data
  
  n.cores <- parallel::detectCores()
    
    my.cluster <- parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
  
    doParallel::registerDoParallel(cl = my.cluster)
    
    clusterEvalQ(my.cluster, {
      library(focus)
      library(igraph)
    })
    
  # initiate items
  
  stab = NULL # List of stability selection Graphical LASSO objects
  adj_graphical = NULL #  List of adjacency matrices from stability selection Graphical LASSO
  network_graphical = NULL # List of Graphs of Graphical LASSO networks
  performance = NULL # List of Graphical LASSO performance metrics
  performance_graph = NULL # List of Performance Graphs for Graphical LASSO
  
  # iterate over simulations
  
  Graphical_network <- foreach (i = 1:length(sim_enrich_small$simulations_complete)) %dopar% {
    stab[[i]] <- focus::GraphicalModel(sim_enrich_small$simulations_complete[[i]]$data)
    adj_graphical[[i]] <- Adjacency(stab[[i]])
    network_graphical[[i]] <- Graph(adj_graphical[[i]])
    performance[[i]] <- SelectionPerformance(theta = adj_graphical[[i]], theta_star = x$simulations_complete[[i]]$theta)
    performance_graph[[i]] <- SelectionPerformanceGraph(
      theta = adj_graphical[[i]],
      theta_star = x$simulations_complete[[i]]$theta,
    )
  }
  
  Graphical_network
  
  ## Returns
  
  parallel::stopCluster(cl = my.cluster)
  
  EnrichGraphicalLASSO = list("graphical_stab" = stab,
                              "adjacency_graphical" = adj_graphical,
                              "network_graphical" = network_graphical,
                              "performance" = performance,
                              "performance_graph" = performance_graph)
  
  return(EnrichGraphicalLASSO)
  
}

