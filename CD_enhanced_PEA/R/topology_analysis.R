######################

# Topology Analysis 

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
library(caret)
library(enrichplot)

######################


######

#' BenchTopologyAnalysis 
#' CD enhanced ORA of true networks
#'
#' @Parameters
#' x = list of true network adjacency matrices
#' pz = size of pathways
#' fz = size of full networks
#' 
#' @returns 
#' detection = dataframe of variables per iteration, true pathways, community membership and true_pathways
#' prop_true = proportion of variables accurately detected 


BenchTopologyAnalysis <- function(x, pz, fz){
  
  network_list = NULL
  graph_list = NULL
  cluster_list = NULL
  communities_list = NULL
  pathways = NULL
  variables = NULL
  true_mapping = NULL
  membership = NULL
  community_membership = NULL
  pathway_members = NULL
  community_mapping = NULL
  count_table = NULL
  pathways_detected = NULL
  communities_detected = NULL
  detection = NULL
  comparison = NULL
  comparison_logical = NULL
  prop = NULL
  
  for (i in 1:length(x)){
    
    # perform community detection for each network graph
    network_list[[i]] <- x[[i]]$theta
    graph_list[[i]] <- Graph(network_list[[i]])
    cluster_list[[i]] <- cluster_edge_betweenness(graph_list[[i]])
    communities_list[[i]] <- communities(cluster_list[[i]])

    # construct df of variables names belonging to each true pathway
    pathways[[i]] <- as.data.frame( unlist(list(rep(1,pz), rep(2,pz), rep(3,pz), rep(4,pz))))
    colnames(pathways[[i]]) <- c("true_pathway")
    variables[[i]] <- as.data.frame(rownames(network_list[[i]]))
    colnames(variables[[i]]) <- c("variable")
    true_mapping[[i]] <- cbind(variables[[i]],pathways[[i]])
    
    # construct df of variable names belonging to each community
    membership[[i]] <- as.matrix(membership(cluster_list[[i]]))
    community_membership[[i]] <- cbind(rownames(membership[[i]]), data.frame(membership[[i]], row.names=NULL))
    colnames(community_membership[[i]]) <- c("variable", "community")
  
    # construct df of variable names, pathways and community memberships
    pathway_members[[i]] <- merge(true_mapping[[i]], community_membership[[i]], all = TRUE)
    community_mapping[[i]] <- pathway_members[[i]][,-1]
    rownames(community_mapping[[i]]) <- pathway_members[[i]][,1]
  
    # count number of variables in each community per pathway
    count_table[[i]] <- as.data.frame(table(community_mapping[[i]]))
    
    ## perform fishers exact test for each community and each pathway
    
    community_pval = NULL
    detected_pathway = NULL
    community_table = NULL
    
    for (j in 1:length(summary(count_table[[i]]$community))){
      
      community_table[[j]] <- subset(count_table[[i]], community == j)
      
      pathway_pval = NULL
      c_p = NULL
      c_np = NULL
      nc_p = NULL
      nc_np = NULL
      
      for (p in 1:length(summary(count_table[[i]]$true_pathway))){
        
        # overlap between community and pathway of interest
        
        c_p[[p]] <- community_table[[j]]$Freq[community_table[[j]]$true_pathway == p]
          
        # overlap between community and outside pathway of interest
        
        c_np[[p]] <- sum(community_table[[j]]$Freq[community_table[[j]]$true_pathway != p])
          
        # overlap between pathway of interest and outside community
          
        nc_p[[p]] <- (pz - c_p[[p]])
        
        # overlap between outside pathway and outside community
        
        nc_np[[p]] <- (fz - pz)
        
        # build matrix
        
        matrix <- matrix(c(c_p[[p]] , c_np[[p]], nc_p[[p]], nc_np[[p]]), nrow = 2, dimnames = list(c("pathway","non-pathway"), c("community","non-community")))
        
        
        # perform fisher exact test
        
        p_value <- fisher.test(matrix, alternative = "greater")[["p.value"]]
        
        pathway_pval[[p]] <- p_value
        
      }
      
      community_pval[[j]] <- pathway_pval
      
      # for each variable in each community select pathway detected
      
      detected_pathway[[j]] <- which.min(community_pval[[j]])
      
      
    }
    
    pathways_detected[[i]] <- t(as.data.frame((detected_pathway)))
    communities_detected[[i]] <- as.data.frame(seq(1:length(summary(count_table[[i]]$community))))
    detection[[i]] <- cbind(communities_detected[[i]], pathways_detected[[i]])
    colnames(detection[[i]]) <- c("community", "detected_pathway")
    
    # calculate proportion of true pathways detected via community enhanced ORA
    comparison[[i]] <- merge(community_mapping[[i]], detection[[i]], by = "community", all = TRUE)
    comparison_logical[[i]] <- (comparison[[i]]$detected_pathway == comparison[[i]]$true_pathway)
    comparison_logical[[i]][is.na(comparison_logical[[i]])] <- FALSE 
    prop[[i]] <- sum(comparison_logical[[i]]) / fz
    
   
  }
  
    output = list("detection" = comparison,
                  "prop_true" = prop)
  
    return(output)
  
}

#####


load("../data/small_topology.RData")
load("../data/large_topology.RData")

bench_small_dis <- small_topology[[1]]
bench_small_sparse <- small_topology[[2]]
bench_small_high <- small_topology[[3]]

b_s_d <- BenchTopologyAnalysis(bench_small_dis, 20, 80)
b_s_s <- BenchTopologyAnalysis(bench_small_sparse, 20, 80)
b_s_h <- BenchTopologyAnalysis(bench_small_high, 20, 80)

bench_large_dis <- large_topology[[1]]
bench_large_sparse <- large_topology[[2]]
bench_large_high <- large_topology[[3]]

b_l_d <- BenchTopologyAnalysis(bench_large_dis, 50, 200)
b_l_s <- BenchTopologyAnalysis(bench_large_sparse, 50, 200)
b_l_h <- BenchTopologyAnalysis(bench_large_high, 50, 200)

######

#' GraphicalTopologyAnalysis 
#' CD enhanced ORA of Graphical networks models
#'
#' @Parameters
#' x = list of true network adjacency matrices
#' pz = size of pathways
#' fz = size of full networks
#' 
#' @returns 
#' detection = dataframe of variables per iteration, true pathways, community membership and true_pathways
#' prop_true = proportion of variables accurately detected 


GraphicalTopologyAnalysis <- function(x, pz, fz){
  
  network_list = NULL
  graph_list = NULL
  cluster_list = NULL
  communities_list = NULL
  pathways = NULL
  variables = NULL
  true_mapping = NULL
  membership = NULL
  community_membership = NULL
  pathway_members = NULL
  community_mapping = NULL
  count_table = NULL
  pathways_detected = NULL
  communities_detected = NULL
  detection = NULL
  comparison = NULL
  comparison_logical = NULL
  prop = NULL
  
  for (i in 1:length(x$adj_graphical)){
    
    # perform community detection for each network graph
    network_list[[i]] <- x$adj_graphical[[i]]
    graph_list[[i]] <- Graph(network_list[[i]])
    cluster_list[[i]] <- cluster_edge_betweenness(graph_list[[i]])
    communities_list[[i]] <- communities(cluster_list[[i]])
    
    # construct df of variables names belonging to each true pathway
    pathways[[i]] <- as.data.frame( unlist(list(rep(1,pz), rep(2,pz), rep(3,pz), rep(4,pz))))
    colnames(pathways[[i]]) <- c("true_pathway")
    variables[[i]] <- as.data.frame(rownames(network_list[[i]]))
    colnames(variables[[i]]) <- c("variable")
    true_mapping[[i]] <- cbind(variables[[i]],pathways[[i]])
    
    # construct df of variable names belonging to each community
    membership[[i]] <- as.matrix(membership(cluster_list[[i]]))
    community_membership[[i]] <- cbind(rownames(membership[[i]]), data.frame(membership[[i]], row.names=NULL))
    colnames(community_membership[[i]]) <- c("variable", "community")
    
    # construct df of variable names, pathways and community memberships
    pathway_members[[i]] <- merge(true_mapping[[i]], community_membership[[i]], all = TRUE)
    community_mapping[[i]] <- pathway_members[[i]][,-1]
    rownames(community_mapping[[i]]) <- pathway_members[[i]][,1]
    
    # count number of variables in each community per pathway
    count_table[[i]] <- as.data.frame(table(community_mapping[[i]]))
    
    ## perform fishers exact test for each community and each pathway
    
    community_pval = NULL
    detected_pathway = NULL
    community_table = NULL
    
    for (j in 1:length(summary(count_table[[i]]$community))){
      
      community_table[[j]] <- subset(count_table[[i]], community == j)
      
      pathway_pval = NULL
      c_p = NULL
      c_np = NULL
      nc_p = NULL
      nc_np = NULL
      
      for (p in 1:length(summary(count_table[[i]]$true_pathway))){
        
        # overlap between community and pathway of interest
        
        c_p[[p]] <- community_table[[j]]$Freq[community_table[[j]]$true_pathway == p]
        
        # overlap between community and outside pathway of interest
        
        c_np[[p]] <- sum(community_table[[j]]$Freq[community_table[[j]]$true_pathway != p])
        
        # overlap between pathway of interest and outside community
        
        nc_p[[p]] <- (pz - c_p[[p]])
        
        # overlap between outside pathway and outside community
        
        nc_np[[p]] <- (fz - pz)
        
        # build matrix
        
        matrix <- matrix(c(c_p[[p]] , c_np[[p]], nc_p[[p]], nc_np[[p]]), nrow = 2, dimnames = list(c("pathway","non-pathway"), c("community","non-community")))
        
        
        # perform fisher exact test
        
        p_value <- fisher.test(matrix, alternative = "greater")[["p.value"]]
        
        pathway_pval[[p]] <- p_value
        
      }
      
      community_pval[[j]] <- pathway_pval
      
      # for each variable in each community select pathway detected
      
      detected_pathway[[j]] <- which.min(community_pval[[j]])
      
      
    }
    
    pathways_detected[[i]] <- t(as.data.frame((detected_pathway)))
    communities_detected[[i]] <- as.data.frame(seq(1:length(summary(count_table[[i]]$community))))
    detection[[i]] <- cbind(communities_detected[[i]], pathways_detected[[i]])
    colnames(detection[[i]]) <- c("community", "detected_pathway")
    
    # calculate proportion of true pathways detected via community enhanced ORA
    comparison[[i]] <- merge(community_mapping[[i]], detection[[i]], by = "community", all = TRUE)
    comparison_logical[[i]] <- (comparison[[i]]$detected_pathway == comparison[[i]]$true_pathway)
    comparison_logical[[i]][is.na(comparison_logical[[i]])] <- FALSE 
    prop[[i]] <- sum(comparison_logical[[i]]) / fz
    
    
  }
  
  output = list("detection" = comparison,
                "prop_true" = prop)
  
  return(output)
  
}

######

load("../data/Graphical_top_small_dis.RData")
Graphical_small_dis <- output

load("../data/Graphical_top_small_sparse.RData")
Graphical_small_sparse <- output

load("../data/Graphical_top_small_high.RData")
Graphical_small_high <- output


load("../data/Graphical_top_large_dis.RData")
Graphical_large_dis <- output

load("../data/Graphical_top_large_sparse.RData")
Graphical_large_sparse <- output

load("../data/Graphical_top_large_high.RData")
Graphical_large_high <- output

g_s_d <- GraphicalTopologyAnalysis(Graphical_small_dis, 20, 80)

g_s_s <- GraphicalTopologyAnalysis(Graphical_small_sparse, 20, 80)

g_s_h <- GraphicalTopologyAnalysis(Graphical_small_high, 20, 80)

g_l_d <- GraphicalTopologyAnalysis(Graphical_large_dis, 50, 200)
  
g_l_s <- GraphicalTopologyAnalysis(Graphical_large_sparse, 50, 200)

g_l_h <- GraphicalTopologyAnalysis(Graphical_large_high, 50, 200)

######

## df small

bench_small_dis = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_small_dis) <- c("model", "pathway_size", "nu", "prop_true")
bench_small_dis <- bench_small_dis %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0,
    prop_true = b_s_d$prop_true
  )

bench_small_sparse = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_small_sparse) <- c("model", "pathway_size", "nu", "prop_true")
bench_small_sparse <- bench_small_sparse %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0.01,
    prop_true = b_s_s$prop_true
  )

bench_small_high = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_small_high) <- c("model", "pathway_size", "nu", "prop_true")
bench_small_high <- bench_small_high %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0.05,
    prop_true = b_s_h$prop_true
  )

graph_small_dis <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_small_dis) <- c("model", "pathway_size", "nu", "prop_true")
graph_small_dis <- graph_small_dis %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0,
    prop_true = g_s_d$prop_true
  )

graph_small_sparse <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_small_sparse) <- c("model", "pathway_size", "nu", "prop_true")
graph_small_sparse <- graph_small_sparse %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0.01,
    prop_true = g_s_s$prop_true
  )

graph_small_high <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_small_high) <- c("model", "pathway_size", "nu", "prop_true")
graph_small_high <- graph_small_high %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0.05,
    prop_true = g_s_h$prop_true
  )

topological_small <- rbind(
  bench_small_dis,
  bench_small_sparse,
  bench_small_high,
  graph_small_dis,
  graph_small_sparse,
  graph_small_high
)

topological_small <- topological_small %>%
  mutate(
    model = as.factor(model),
    pathway_size = as.factor(pathway_size),
    nu = as.factor(nu),
    prop_true = as.numeric(prop_true)
  )

str(topological_small)

g1 <- ggplot(data = topological_small, aes( y = prop_true, fill = model)) + geom_boxplot() + facet_wrap(topological_small$nu)
plot(g1)

## df large

bench_large_dis = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_large_dis) <- c("model", "pathway_size", "nu", "prop_true")
bench_large_dis <- bench_large_dis %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0,
    prop_true = b_l_d$prop_true
  )

bench_large_sparse = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_large_sparse) <- c("model", "pathway_size", "nu", "prop_true")
bench_large_sparse <- bench_large_sparse %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0.01,
    prop_true = b_l_s$prop_true
  )

bench_large_high = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(bench_large_high) <- c("model", "pathway_size", "nu", "prop_true")
bench_large_high <- bench_large_high %>%
  mutate(
    model = "bench",
    pathway_size = 20,
    nu = 0.05,
    prop_true = b_l_h$prop_true
  )

graph_large_dis <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_large_dis) <- c("model", "pathway_size", "nu", "prop_true")
graph_large_dis <- graph_large_dis %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0,
    prop_true = g_l_d$prop_true
  )

graph_large_sparse <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_large_sparse) <- c("model", "pathway_size", "nu", "prop_true")
graph_large_sparse <- graph_large_sparse %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0.01,
    prop_true = g_l_s$prop_true
  )

graph_large_high <- as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(graph_large_high) <- c("model", "pathway_size", "nu", "prop_true")
graph_large_high <- graph_large_high %>%
  mutate(
    model = "graphical",
    pathway_size = 20,
    nu = 0.05,
    prop_true = g_l_h$prop_true
  )

topological_large <- rbind(
  bench_large_dis,
  bench_large_sparse,
  bench_large_high,
  graph_large_dis,
  graph_large_sparse,
  graph_large_high
)

topological_large <- topological_large %>%
  mutate(
    model = as.factor(model),
    pathway_size = as.factor(pathway_size),
    nu = as.factor(nu),
    prop_true = as.numeric(prop_true)
  )

str(topological_large)

g2 <- ggplot(data = topological_large, aes( y = prop_true, fill = model)) + geom_boxplot() + facet_wrap(topological_large$nu)
plot(g2)

save(topological_small, file = "../data/final_df_topsmall.RData")
save(topological_large, file= "../data/final_df_toplarge.RData")
