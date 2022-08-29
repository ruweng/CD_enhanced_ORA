######################

# Network Simulation

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

######################

#' B. [Enrichment Detection Study] Data-base-derived Simulation
#' simulate network by randomly selecting items from pathway database
#' 1) randomly select q pathways from KEGG database 
#' 2) generate adjacency matrix of network of selected pathway
#' 3) include detection call (remove adjacency) and random noise (both annotated and unannotated variables)
#' 4) simulate data with sharp (multivariate normal distribution) using adjacency matrix as true underlying structure
#' 
#' 
#' ---
#' 
#' @Parameters
#' k number of iterations
#' n observations in simulation
#' m pathways selected 
#' dc detection call (percentage of (randomly selected) true underlying structure detected per pathway) # interesting for clustering applications as it breaks up structure
#' f amount of random noise from KEGG annotated variables
#' g amount of random noise from un-annotated variables
#' 
#' 
#' ---
#' 
#' @Returns
#' pathways_list - list of selected / enriched pathways list per iteration
#' 
#' topology_true_list - list of symmetric pathway adjacency matrices lists
#' adj_network_true_list - list of adjacency matrices of true networks
#' sim_network_true_list - simulated data based on adjacency matrices of true networks
#' 
#' topology_detected_list - list of detected pathway adjacency matrices lists 
#' adj_network_detected_list - list of adjacency matrices of detected networks
#' sim_network_detected_list - list of simulated data based on adjacency matrices of detected networks
#' 
#' topology_complete_list - list of complete pathway adjacency matrices lists
#' adj_network_complete_list - list of adjacency matrices of complete networks
#' sim_network_complete_list - list of simulated data based on adjacency matrices of complete networks
#' 
#' 
#' for each iteration k:
#' pathways - list of selected / enriched pathways 
#' pathway_names_list - list of pathway node names
#' path_list - list of directed graphs of enriched pathways
#' theta_list - list of directed adjacency matrices of enriched pathways
#' det_theta_list - list of detected directed adjacency matrices of enriched pathways
#' graph_list - list of undirected graphs of enriched pathways
#' detected_graph_list - list of detected undirected graphs of enriched pathways
#' adj_list - list of undirected / symmetric adjacency matrices of enriched pathways
#' detected_adj_list - list of detected symmetric adjacency matrices of enriched pathways
#' ---


SimulateEnrichment <- function(k, n, m, dc, f, g){
  
  ### 
  
  # KEGG Simulation
  
  # select full KEGG network
  KEGG <- pathways("hsapiens", "kegg")
  
  # iterate sampling over k iterations
  
  pathways_list = NULL
  topology_list_true = NULL
  topology_list_detected = NULL
  adj_network_true_list = NULL
  adj_network_detected_list = NULL
  adj_network_complete_list = NULL
  sim_network_true_list = NULL
  sim_network_detected_list = NULL
  sim_network_complete_list = NULL
  
  for (j in 1:k){
    
    # randomly select m amount of pathways enriched
    pathways <- KEGG[sample(1:length(KEGG), m)]
    
    # initiate iteration over i = number of pathways
    pathway_nodes_list = NULL # list of pathway node names
    path_list = NULL # directed enriched pathways
    theta_list = NULL # directed adjacency matrices of enriched pathways
    det_theta_list = NULL # detected directed adjacency matrices of enriched pathways
    
    adj_list = NULL # undirected adjacency matrices of enriched pathways
    detected_adj_list = NULL # undirected adjacency matrices of detected enriched pathways
    graph_list = NULL # undirected enriched pathways 
    detected_graph_list = NULL # detected undirected pathways
    
    network_true = graph_from_literal() # true network of enriched pathways
    network_detected = graph_from_literal() # detected network of enriched pathways
    
    for (i in 1:length(pathways)){
      
      ## Directed pathways
      
      # create list of pathway nodes
      
      pathway_nodes <- nodes(pathways[[i]], which = c("proteins"))
      pathway_nodes_list <- c(pathway_nodes_list,pathway_nodes)
      
      # create directed graph for each pathway
      path <- pathwayGraph(pathways[[i]])
      path_list[[i]] <- path
      
      # create directed adjacency matrix for each pathway
      theta <- as(path, "matrix")
      theta_list[[i]] <- theta
      
      ## Undirected pathways
      
      # create undirected graph for each pathway
      graph <- graph_from_adjacency_matrix(theta_list[[i]], mode = "undirected")
      graph_list[[i]] <- graph
      
      # create symmetric adjacency matrix for each pathway
      adj <- as.matrix(as_adjacency_matrix(graph_list[[i]]))
      adj_list[[i]] <- adj
      
      ## Network
      
      # create true (full) network of enriched pathways
      network_true <- network_true + graph
      
      
      ### Detected Topology and Network
      
      ## Detected Pathway
      
      # randomly sample (1-dc) percent of pathway adjacency matrix and set to 0
      
      det_theta <- theta # detected adjacency matrix
      ind <- which(det_theta == 1) # Gives you the indices of all 1 in the adjacency matrix
      ind_to_change <- sample(ind, floor(length(ind)*(1-dc))) # randomly sample dc% of the indices
      det_theta[ind_to_change] <- 0 # set the sample indices to 0
      det_theta_list[[i]] <- det_theta
      
      ## graph detected adjacency matrices
      
      det_graph <- graph_from_adjacency_matrix(det_theta_list[[i]], mode = "undirected")
      detected_graph_list[[i]] <- det_graph
      det_adj <- as.matrix(as_adjacency_matrix(det_graph))
      detected_adj_list[[i]] <- det_adj
      
      ## Detected Network
      
      # create detected network of enriched pathways
      network_detected <- network_detected + det_graph
      
      
    }
    
    ### Simulate Networks
    
    ## True Network
    
    # create adjacency matrix 
    adj_network_true <- as.matrix(as_adjacency_matrix(network_true))
    
    # simulate data based on network 
    sim_network_true <- SimulateGraphical(n = n, theta = adj_network_true)
    
    ## Detected Network
    
    # create adjacency matrix of detected network
    adj_network_detected <- as.matrix(as_adjacency_matrix(network_detected))
    
    # simulate data based on detected network
    sim_network_detected <- SimulateGraphical(n = n, theta = adj_network_detected)
    
    ## Complete network with random noise
    
    # add random noise of f KEGG annotated variables
    
    # create list of unique KEGG IDs
    
    KEGG_list <- as.list(KEGG)
    
    KEGG_nodes_list = NULL
    
    for (r in 1:length(KEGG_list)) {
      KEGG_nodes <- nodes(KEGG[[r]], which = c("proteins"))
      KEGG_nodes_list <- c(KEGG_nodes_list, KEGG_nodes)
    }
    
    KEGG_IDs <- unique(KEGG_nodes_list)
    
    # create list of unique pathway IDs
    
    pathway_IDs <- unique(pathway_nodes_list)
    
    # create list of unique KEGG IDs not in pathway
    
    non_path_IDs <- KEGG_IDs[KEGG_IDs %in% pathway_IDs == FALSE]
    
    # randomly sample f variables from KEGG outside of pathways
    
    f_IDs <- sample(non_path_IDs, f)
    
    # add f variables to detected network
    
    annot_graph <- make_empty_graph(n = f, directed = FALSE)
    
    V(annot_graph)$name <- f_IDs
    
    network_complete <- network_detected + annot_graph
    
    # add random noise of g un-annotated variables
    
    unannot_graph <- make_empty_graph(n = g, directed = FALSE)
    
    V(unannot_graph)$name <- paste0("Var",1:g)
    
    network_complete <- network_complete + unannot_graph
    
    # create adjacency matrix 
    adj_network_complete <- as.matrix(as_adjacency_matrix(network_complete))
    
    # simulate data based on network 
    sim_network_complete <- SimulateGraphical(n = n, theta = adj_network_complete)
    colnames(sim_network_complete$data) <- colnames(adj_network_complete)
    rownames(sim_network_complete$theta) <- rownames(adj_network_complete)
    colnames(sim_network_complete$theta) <- colnames(adj_network_complete)
    
    ### Returns
    
    # list of pathway list
    pathways_list[[j]] <- pathways
    
    # list of adjacency matrix list
    topology_list_true[[j]] <- adj_list
    topology_list_detected[[j]] <- detected_adj_list
    
    # list of adjacency matrices of true networks
    adj_network_true_list[[j]] <- adj_network_true
    adj_network_detected_list[[j]] <- adj_network_detected
    adj_network_complete_list[[j]] <- adj_network_complete
    
    # list of simulated data based on adjacency matrices of true networks
    sim_network_true_list[[j]] <- sim_network_true
    sim_network_detected_list[[j]] <- sim_network_detected
    sim_network_complete_list[[j]] <- sim_network_complete
    
  }
  
  EnrichmentSimulation <- list("pathways" = pathways_list,
                               "pathway_topologies_true" = topology_list_true,
                               "pathway_topologies_detected" = topology_list_detected,
                               "network_adjacencies_true" = adj_network_true_list,
                               "network_adjacencies_detected" = adj_network_detected_list,
                               "network_adjacencies_complete" = adj_network_complete_list,
                               "simulations_true" = sim_network_true_list,
                               "simulations_detected" = sim_network_detected_list,
                               "simulations_complete" = sim_network_complete_list)
  return(EnrichmentSimulation)
}

### detection call

set.seed(2222)
sim_enrichDC_small <- SimulateEnrichment(100, 1000, 4, 1.0, 100, 100)
save(sim_enrichDC_small, file = "../data/sim_enrichDC_small.RData")

set.seed(2222)
sim_enrichDC_medium <- SimulateEnrichment(100, 1000, 4, 0.7, 100, 100)
save(sim_enrichDC_medium, file = "../data/sim_enrichDC_medium.RData")

set.seed(2222)
sim_enrichDC_large <- SimulateEnrichment(100, 1000, 4, 0.5, 100, 100)
save(sim_enrichDC_large, file = "../data/sim_enrichDC_large.RData")

### pathway number

set.seed(2222)
sim_enrichP_small <- SimulateEnrichment(100, 1000, 2, 0.7, 100, 100)
save(sim_enrichP_small, file = "../data/sim_enrichP_small.RData")

set.seed(2222)
sim_enrichP_medium <- SimulateEnrichment(100, 1000, 4, 0.7, 100, 100)
save(sim_enrichP_medium, file = "../data/sim_enrichP_medium.RData")

set.seed(2222)
sim_enrichP_large <- SimulateEnrichment(100, 1000, 8, 0.7, 100, 100)
save(sim_enrichP_large, file = "../data/sim_enrichP_large.RData")

