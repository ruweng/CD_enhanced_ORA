######################

# Enrichment Analysis (ORA)

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

# load("../data/sim_enrich_small.RData")
# load("../data/sim_enrich_medium.RData")
# load("../data/sim_enrich_large.RData")


### Community detection (potentially stability enhanced edge-betweenness clustering?)

bench <- function(x,m){
  
  ## Community detection and ORA
  
  true_adj = NULL
  true_network = NULL
  true_network_clean = NULL
  cluster = NULL
  communities = NULL
  communities_ORA_list = NULL
  communities_ORA_list_clean = NULL
  network_names = NULL
  network_names_clean = NULL
  network_ORA_list = NULL
  
  for (i in 1:length(x$simulations_complete)){
  
  # use true underlying structure as adjacency matrix
  
  true_adj[[i]] <- x$simulations_complete[[i]]$theta
  
  # graph true network
  
  true_network[[i]] <- graph_from_adjacency_matrix(true_adj[[i]], mode = "undirected")
  
  true_network_clean[[i]] <- delete.vertices(true_network[[i]], V(true_network[[i]])[ degree(true_network[[i]]) == 0])
  
  
  # community detection
  
  cluster[[i]] <- cluster_edge_betweenness(true_network_clean[[i]]) 
  
  communities[[i]] <- communities(cluster[[i]])
  
  # perform ORA for communities 
  
  community = NULL
  community_clean = NULL
  community_ORA = NULL
  
    for (j in 1:length(communities[[i]])){
      community[[j]] <- as.vector(communities[[i]][[j]])
      community_clean[[j]] <- gsub("ENTREZID:","", community[[j]])
      community_clean[[j]] <- community_clean[[j]][!(community_clean[[j]] = str_detect(community_clean[[j]], "Var"))]                
      community_ORA[[j]] <- enrichKEGG(community_clean[[j]], 
                             keyType = "kegg", 
                             pAdjustMethod = "bonferroni",
                             pvalueCutoff = 0.05)
    }
  communities_ORA_list[[i]] <- community_ORA
  communities_ORA_list_clean[[i]] <- compact(communities_ORA_list[[i]])
  
  # perform ORA for networks
  
  network_names[[i]] <- as.vector(V(true_network_clean[[i]])$name)
  network_names_clean[[i]] <- gsub("ENTREZID:","", network_names[[i]])
  # network_names_clean[[i]] <-  network_names[[i]][!(network_names[[i]] = str_detect(network_names[[i]], "Var"))]      
  network_ORA_list[[i]] <- enrichKEGG(network_names_clean[[i]], 
                                      keyType = "kegg", 
                                      pAdjustMethod = "bonferroni",
                                      pvalueCutoff = 0.05)
  
  }
  
  ### results confusion matrix
  
    ## return list of KEGG pathway titles
    
    KEGG <- pathways("hsapiens", "kegg")
    
    KEGG_titles = NULL # list of pathways in KEGG
    
    for (i in 1:length(KEGG@entries))(
      KEGG_titles[[i]] <- KEGG@entries[[i]]@title
    )
    
    ## return list of true pathway titles per simulation
    
    true_pathways_list = NULL # list of true pathways per simulation
    true_pathways_titles_list = NULL # list of pathway titles per simulation
    
    
    for (i in 1:length(x$simulations_complete)){
     
      
      true_pathways_list[[i]] <- x$pathways[[i]]@entries
      
      true_pathways_titles = NULL # titles of pathways in simulation
      
      for (j in 1:length(true_pathways_list[[i]])){
      true_pathways_titles[[j]] <- true_pathways_list[[i]][[j]]@title
      }
      
      true_pathways_titles_list[[i]] <- true_pathways_titles
    }
   
    ## return list of detected pathway titles per simulation
      # for each community take first title 
    

    detected_pathways_list = NULL

    for (i in 1:length(communities_ORA_list_clean)){
      
        significant_pathway_titles = NULL
          
        for (j in 1:length(communities_ORA_list_clean[[i]])){
          significant_pathway_titles[[j]] <- communities_ORA_list_clean[[i]][[j]]@result$Description[1]
        }
        detected_pathways_list[[i]] <- significant_pathway_titles
        
    }
    
    # for each network take all below p-val threshold
      
    network_ORA_pathway_titles = NULL
    m_network_ORA_pathway_titles = NULL
    
      for (i in 1:length(network_ORA_list)){
        network_ORA_pathway_titles[[i]] <- subset(network_ORA_list[[i]]@result$Description, network_ORA_list[[i]]@result$p.adjust < 0.05)
        m_network_ORA_pathway_titles[[i]] <- network_ORA_list[[i]]@result$Description[1:m]
      }
    
    ## confusion matrix CD ORA
    
    KEGG_negatives = NULL
    true_positives = NULL
    false_positives = NULL
    true_negatives = NULL
    false_negatives = NULL
    
    for (i in 1:length(detected_pathways_list)){
      true_positives[[i]] <- sum(table(detected_pathways_list[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
      false_positives[[i]] <- sum(table(detected_pathways_list[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
      
      KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
      true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% detected_pathways_list[[i]])["FALSE"])
      false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% detected_pathways_list[[i]])["TRUE"])
    }
    
    # confusion matrix of frequencies
    true_positives <- do.call(rbind.data.frame, true_positives)
    false_positives <- do.call(rbind.data.frame, false_positives)
    true_negatives <- do.call(rbind.data.frame, true_negatives)
    false_negatives <- do.call(rbind.data.frame, false_negatives)
    CD_confusion <- cbind(true_positives, false_positives, true_negatives, false_negatives)
    CD_confusion[is.na(CD_confusion)] <- 0
    colnames(CD_confusion) <- c("TP", "FP", "TN", "FN")
    
    # calculate proportion of correct / incorrect detection
    CD_confusion <- CD_confusion %>%
      mutate(TP_prop = TP / (FP + TP),
             FP_prop = FP / (FP + TP),
             TN_prop = TN / (FN + TN),
             FN_prop = FN / (FN + TN),
             specificity = TN / (TN + FP),
             sensitivity = TP / (TP + FN)
             )
    
    ## confusion matrix Network ORA
    
    network_true_positives = NULL
    network_false_positives = NULL
    network_true_negatives = NULL
    network_false_negatives = NULL
    
    for (i in 1:length(network_ORA_pathway_titles)){
      network_true_positives[[i]] <- sum(table(network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
      network_false_positives[[i]] <- sum(table(network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
      
      KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
      network_true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% network_ORA_pathway_titles[[i]])["FALSE"])
      network_false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% network_ORA_pathway_titles[[i]])["TRUE"])
    }
    
    # confusion matrix of frequencies
    network_true_positives <- do.call(rbind.data.frame, network_true_positives)
    network_false_positives <- do.call(rbind.data.frame, network_false_positives)
    network_true_negatives <- do.call(rbind.data.frame, network_true_negatives)
    network_false_negatives <- do.call(rbind.data.frame, network_false_negatives)
    network_confusion <- cbind(network_true_positives, network_false_positives, network_true_negatives, network_false_negatives)
    network_confusion[is.na(network_confusion)] <- 0
    colnames(network_confusion) <- c("TP", "FP", "TN", "FN")
    
    # calculate proportion of correct / incorrect detection
    network_confusion <- network_confusion %>%
      mutate(TP_prop = TP / (FP + TP),
             FP_prop = FP / (FP + TP),
             TN_prop = TN / (FN + TN),
             FN_prop = FN / (FN + TN),
             specificity = TN / (TN + FP),
             sensitivity = TP / (TP + FN)
      )
    
    ## confusion matrix m Network ORA
    
    m_network_true_positives = NULL
    m_network_false_positives = NULL
    m_network_true_negatives = NULL
    m_network_false_negatives = NULL
    
    for (i in 1:length(m_network_ORA_pathway_titles)){
      m_network_true_positives[[i]] <- sum(table(m_network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
      m_network_false_positives[[i]] <- sum(table(m_network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
      
      KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
      m_network_true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% m_network_ORA_pathway_titles[[i]])["FALSE"])
      m_network_false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% m_network_ORA_pathway_titles[[i]])["TRUE"])
    }
    
    # confusion matrix of frequencies
    m_network_true_positives <- do.call(rbind.data.frame, m_network_true_positives)
    m_network_false_positives <- do.call(rbind.data.frame, m_network_false_positives)
    m_network_true_negatives <- do.call(rbind.data.frame, m_network_true_negatives)
    m_network_false_negatives <- do.call(rbind.data.frame, m_network_false_negatives)
    m_network_confusion <- cbind(m_network_true_positives, m_network_false_positives, m_network_true_negatives, m_network_false_negatives)
    m_network_confusion[is.na(m_network_confusion)] <- 0
    colnames(m_network_confusion) <- c("TP", "FP", "TN", "FN")
    
    # calculate proportion of correct / incorrect detection
    m_network_confusion <- m_network_confusion %>%
      mutate(TP_prop = TP / (FP + TP),
             FP_prop = FP / (FP + TP),
             TN_prop = TN / (FN + TN),
             FN_prop = FN / (FN + TN),
             specificity = TN / (TN + FP),
             sensitivity = TP / (TP + FN)
      )
    
    
    output = list( "Community_detection_confusion_matrix" = CD_confusion,
                   "Network_confusion_matrix" = network_confusion,
                   "m_network_confusion_matrix" = m_network_confusion,
                   "m_network_titles" = m_network_ORA_pathway_titles,
                   "detected_pathways" = detected_pathways_list,
                   "true_pathways" = true_pathways_titles_list,
                   "network_titles" = network_ORA_pathway_titles,
                   "network_ORA" = network_ORA_list,
                   "network_names" = network_names,
                   "network_names_clean" = network_names_clean
    )
    
    return(output)
}

# limitations 
# not accounting for overlaps outside of enriched pathways
# not accounting for bad community structure
## how many pathways to choose ?

#####

load("../data/sim_enrichDC_small.RData")
load("../data/sim_enrichDC_medium.RData")
load("../data/sim_enrichDC_large.RData")

load("../data/sim_enrichP_small.RData")
load("../data/sim_enrichP_medium.RData")
load("../data/sim_enrichP_large.RData")

bench_DC_small <- bench(sim_enrichDC_small, 4)
bench_DC_medium <- bench(sim_enrichDC_medium, 4)
bench_DC_large <- bench(sim_enrichDC_large, 4)

bench_P_small <- bench(sim_enrichP_small, 2)
bench_P_medium <- bench(sim_enrichP_medium, 4)
bench_P_large <- bench(sim_enrichP_large, 8)

bench_DC_summary <- list("bench_DC_small" = bench_DC_small,
                         "bench_DC_medium" = bench_DC_medium,
                         "bench_DC_large" = bench_DC_large)

bench_P_summary <- list("bench_P_small" = bench_P_small,
                        "bench_P_medium" = bench_P_medium,
                        "bench_P_large" = bench_P_large)

save(bench_DC_summary, file = "../data/bench_DC.Rdata")
save(bench_P_summary, file = "../data/bench_P.Rdata")

load("../data/Graphical_enrichDC_small.RData")
load("../data/sim_enrichDC_small.Rdata")

#######

experiment <- function(x, g, m){
  
  ## Community detection and ORA
  

  true_network_clean = NULL
  cluster = NULL
  communities = NULL
  communities_ORA_list = NULL
  communities_ORA_list_clean = NULL
  network_names = NULL
  network_names_clean = NULL
  network_ORA_list = NULL
  
  for (i in 1:length(x$simulations_complete)){
    
    true_network_clean[[i]] <- g$network_graphical[[i]]
    
    # community detection
    
    cluster[[i]] <- cluster_edge_betweenness(true_network_clean[[i]]) 
    
    communities[[i]] <- communities(cluster[[i]])
    
    # perform ORA for communities 
    
    community = NULL
    community_clean = NULL
    community_ORA = NULL
    
    for (j in 1:length(communities[[i]])){
      community[[j]] <- as.vector(communities[[i]][[j]])
      community_clean[[j]] <- gsub("ENTREZID:","", community[[j]])
      community_clean[[j]] <- community_clean[[j]][!(community_clean[[j]] = str_detect(community_clean[[j]], "Var"))]                
      community_ORA[[j]] <- enrichKEGG(community_clean[[j]], 
                                       keyType = "kegg", 
                                       pAdjustMethod = "bonferroni",
                                       pvalueCutoff = 0.05)
    }
    communities_ORA_list[[i]] <- community_ORA
    communities_ORA_list_clean[[i]] <- compact(communities_ORA_list[[i]])
    
    # perform ORA for networks
    
    network_names[[i]] <- as.vector(V(true_network_clean[[i]])$name)
    network_names_clean[[i]] <- gsub("ENTREZID:","", network_names[[i]])
    # network_names_clean[[i]] <-  network_names[[i]][!(network_names[[i]] = str_detect(network_names[[i]], "Var"))]      
    network_ORA_list[[i]] <- enrichKEGG(network_names_clean[[i]], 
                                        keyType = "kegg", 
                                        pAdjustMethod = "bonferroni",
                                        pvalueCutoff = 0.05)
  }
  
  ### results confusion matrix
  
  ## return list of KEGG pathway titles
  
  KEGG <- pathways("hsapiens", "kegg")
  
  KEGG_titles = NULL # list of pathways in KEGG
  
  for (i in 1:length(KEGG@entries))(
    KEGG_titles[[i]] <- KEGG@entries[[i]]@title
  )
  
  ## return list of true pathway titles per simulation
  
  true_pathways_list = NULL # list of true pathways per simulation
  true_pathways_titles_list = NULL # list of pathway titles per simulation
  
  
  for (i in 1:length(x$simulations_complete)){
    
    
    true_pathways_list[[i]] <- x$pathways[[i]]@entries
    
    true_pathways_titles = NULL # titles of pathways in simulation
    
    for (j in 1:length(true_pathways_list[[i]])){
      true_pathways_titles[[j]] <- true_pathways_list[[i]][[j]]@title
    }
    
    true_pathways_titles_list[[i]] <- true_pathways_titles
  }
  
  ## return list of detected pathway titles per simulation
  # for each community take first title 
  
  
  detected_pathways_list = NULL
  
  for (i in 1:length(communities_ORA_list_clean)){
    
    significant_pathway_titles = NULL
    
    for (j in 1:length(communities_ORA_list_clean[[i]])){
      significant_pathway_titles[[j]] <- communities_ORA_list_clean[[i]][[j]]@result$Description[1]
    }
    detected_pathways_list[[i]] <- significant_pathway_titles
    
  }
  
  # for each network take all below p-val threshold
  
  network_ORA_pathway_titles = NULL
  m_network_ORA_pathway_titles = NULL
  
  for (i in 1:length(network_ORA_list)){
    network_ORA_pathway_titles[[i]] <- subset(network_ORA_list[[i]]@result$Description, network_ORA_list[[i]]@result$p.adjust < 0.05)
    m_network_ORA_pathway_titles[[i]] <- network_ORA_list[[i]]@result$Description[1:m]
  }
  
  ## confusion matrix CD ORA
  
  KEGG_negatives = NULL
  true_positives = NULL
  false_positives = NULL
  true_negatives = NULL
  false_negatives = NULL
  
  for (i in 1:length(detected_pathways_list)){
    true_positives[[i]] <- sum(table(detected_pathways_list[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
    false_positives[[i]] <- sum(table(detected_pathways_list[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
    
    KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
    true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% detected_pathways_list[[i]])["FALSE"])
    false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% detected_pathways_list[[i]])["TRUE"])
  }
  
  # confusion matrix of frequencies
  true_positives <- do.call(rbind.data.frame, true_positives)
  false_positives <- do.call(rbind.data.frame, false_positives)
  true_negatives <- do.call(rbind.data.frame, true_negatives)
  false_negatives <- do.call(rbind.data.frame, false_negatives)
  CD_confusion <- cbind(true_positives, false_positives, true_negatives, false_negatives)
  CD_confusion[is.na(CD_confusion)] <- 0
  colnames(CD_confusion) <- c("TP", "FP", "TN", "FN")
  
  # calculate proportion of correct / incorrect detection
  CD_confusion <- CD_confusion %>%
    mutate(TP_prop = TP / (FP + TP),
           FP_prop = FP / (FP + TP),
           TN_prop = TN / (FN + TN),
           FN_prop = FN / (FN + TN),
           specificity = TN / (TN + FP),
           sensitivity = TP / (TP + FN))
  
  ## confusion matrix Network ORA
  
  network_true_positives = NULL
  network_false_positives = NULL
  network_true_negatives = NULL
  network_false_negatives = NULL
  
  for (i in 1:length(network_ORA_pathway_titles)){
    network_true_positives[[i]] <- sum(table(network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
    network_false_positives[[i]] <- sum(table(network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
    
    KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
    network_true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% network_ORA_pathway_titles[[i]])["FALSE"])
    network_false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% network_ORA_pathway_titles[[i]])["TRUE"])
  }
  
  # confusion matrix of frequencies
  network_true_positives <- do.call(rbind.data.frame, network_true_positives)
  network_false_positives <- do.call(rbind.data.frame, network_false_positives)
  network_true_negatives <- do.call(rbind.data.frame, network_true_negatives)
  network_false_negatives <- do.call(rbind.data.frame, network_false_negatives)
  network_confusion <- cbind(network_true_positives, network_false_positives, network_true_negatives, network_false_negatives)
  network_confusion[is.na(network_confusion)] <- 0
  colnames(network_confusion) <- c("TP", "FP", "TN", "FN")
  
  # calculate proportion of correct / incorrect detection
  network_confusion <- network_confusion %>%
    mutate(TP_prop = TP / (FP + TP),
           FP_prop = FP / (FP + TP),
           TN_prop = TN / (FN + TN),
           FN_prop = FN / (FN + TN),
           specificity = TN / (TN + FP),
           sensitivity = TP / (TP + FN)
           )
  
  ## confusion matrix m Network ORA
  
  m_network_true_positives = NULL
  m_network_false_positives = NULL
  m_network_true_negatives = NULL
  m_network_false_negatives = NULL
  
  for (i in 1:length(m_network_ORA_pathway_titles)){
    m_network_true_positives[[i]] <- sum(table(m_network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["TRUE"])
    m_network_false_positives[[i]] <- sum(table(m_network_ORA_pathway_titles[[i]] %in% true_pathways_titles_list[[i]])["FALSE"])
    
    KEGG_negatives[[i]] <- subset((KEGG_titles), !KEGG_titles %in% true_pathways_titles_list[[i]])
    m_network_true_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% m_network_ORA_pathway_titles[[i]])["FALSE"])
    m_network_false_negatives[[i]] <- sum(table(KEGG_negatives[[i]] %in% m_network_ORA_pathway_titles[[i]])["TRUE"])
  }
  
  # confusion matrix of frequencies
  m_network_true_positives <- do.call(rbind.data.frame, m_network_true_positives)
  m_network_false_positives <- do.call(rbind.data.frame, m_network_false_positives)
  m_network_true_negatives <- do.call(rbind.data.frame, m_network_true_negatives)
  m_network_false_negatives <- do.call(rbind.data.frame, m_network_false_negatives)
  m_network_confusion <- cbind(m_network_true_positives, m_network_false_positives, m_network_true_negatives, m_network_false_negatives)
  m_network_confusion[is.na(m_network_confusion)] <- 0
  colnames(m_network_confusion) <- c("TP", "FP", "TN", "FN")
  
  # calculate proportion of correct / incorrect detection
  m_network_confusion <- m_network_confusion %>%
    mutate(TP_prop = TP / (FP + TP),
           FP_prop = FP / (FP + TP),
           TN_prop = TN / (FN + TN),
           FN_prop = FN / (FN + TN),
           specificity = TN / (TN + FP),
           sensitivity = TP / (TP + FN)
           )
  
  
  output = list( "Community_detection_confusion_matrix" = CD_confusion,
                 "Network_confusion_matrix" = network_confusion,
                 "M_Network_confusion_matrix" = m_network_confusion,
                 "network_ORA_pathway_titles" = network_ORA_pathway_titles,
                 "m_network_ORA_pathway_titles" = m_network_ORA_pathway_titles,
                 "detected_pathways" = detected_pathways_list,
                 "true_pathways" = true_pathways_titles_list
  )
  
  return(output)
}

#####

load("../data/Graphical_enrichDC_small.RData")
Graphical_enrichDC_small <- output

load("../data/Graphical_enrichDC_medium.RData")
Graphical_enrichDC_medium <- output

load("../data/Graphical_enrichDC_large.RData")
Graphical_enrichDC_large <- output

exp_enrichDC_small <- experiment(sim_enrichDC_small, Graphical_enrichDC_small, 4)
exp_enrichDC_medium <- experiment(sim_enrichDC_medium, Graphical_enrichDC_medium, 4)
exp_enrichDC_large <- experiment(sim_enrichDC_large, Graphical_enrichDC_large, 4)

exp_enrichDC_summary <- list("exp_enrichDC_small" = exp_enrichDC_small,
                             "exp_enrichDC_medium" = exp_enrichDC_medium,
                             "exp_enrichDC_large" = exp_enrichDC_large)

save(exp_enrichDC_summary, file = "../data/exp_enrichDC.RData")

load("../data/Graphical_enrichP_small.RData")
Graphical_enrichP_small <- output

load("../data/Graphical_enrichP_medium.RData")
Graphical_enrichP_medium <- output

load("../data/Graphical_enrichP_large.RData")
Graphical_enrichP_large <- output


exp_enrichP_small <- experiment(sim_enrichP_small, Graphical_enrichP_small, 2)
exp_enrichP_medium <- experiment(sim_enrichP_medium, Graphical_enrichP_medium, 4)
exp_enrichP_large <- experiment(sim_enrichP_large, Graphical_enrichP_large, 8)

exp_enrichP_summary <- list("exp_enrichP_small" = exp_enrichP_small,
                             "exp_enrichP_medium" = exp_enrichP_medium,
                             "exp_enrichP_large" = exp_enrichP_large)

save(exp_enrichP_summary, file = "../data/exp_enrichP.RData")


#####

# Figures


load("../data/exp_enrichDC.Rdata")
load("../data/bench_DC.Rdata")

# DC

  ## specificitiy 
  
    ### bench 
    
      #### exp
    
        DC1_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_bench_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_bench_exp_specificity <- DC1_bench_exp_specificity %>%
          mutate(DC = 1, 
                 model = "bench",
                 arm = "exp",
                 specificity = bench_DC_summary[["bench_DC_small"]][["Community_detection_confusion_matrix"]][["specificity"]])
        
        DC0.7_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_bench_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_bench_exp_specificity <- DC0.7_bench_exp_specificity %>%
          mutate(DC = 0.7, 
                 model = "bench",
                 arm = "exp",
                 specificity = bench_DC_summary[["bench_DC_medium"]][["Community_detection_confusion_matrix"]][["specificity"]])
        
        DC0.5_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_bench_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_bench_exp_specificity <- DC0.5_bench_exp_specificity %>%
          mutate(DC = 0.5, 
                 model = "bench",
                 arm = "exp",
                 specificity = bench_DC_summary[["bench_DC_large"]][["Community_detection_confusion_matrix"]][["specificity"]])
    
      #### CtrlALL
        
        DC1_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_bench_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_bench_ctrlall_specificity <- DC1_bench_ctrlall_specificity %>%
          mutate(DC = 1,
                 model = "bench",
                 arm = "ctrl_all",
                 specificity = bench_DC_summary[["bench_DC_small"]][["Network_confusion_matrix"]][["specificity"]])
        
        DC0.7_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_bench_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_bench_ctrlall_specificity <- DC0.7_bench_ctrlall_specificity %>%
          mutate(DC = 0.7,
                 model = "bench",
                 arm = "ctrl_all",
                 specificity = bench_DC_summary[["bench_DC_medium"]][["Network_confusion_matrix"]][["specificity"]])
        
        DC0.5_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_bench_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_bench_ctrlall_specificity <- DC0.5_bench_ctrlall_specificity %>%
          mutate(DC = 0.5,
                 model = "bench",
                 arm = "ctrl_all",
                 specificity = bench_DC_summary[["bench_DC_large"]][["Network_confusion_matrix"]][["specificity"]])
        
      #### CtrlM
        
        DC1_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_bench_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_bench_ctrlm_specificity <- DC1_bench_ctrlm_specificity %>%
          mutate(DC = 1,
                 model = "bench",
                 arm = "ctrl_m",
                 specificity = bench_DC_summary[["bench_DC_small"]][["m_network_confusion_matrix"]][["specificity"]])
        
        DC0.7_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_bench_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_bench_ctrlm_specificity <- DC0.7_bench_ctrlm_specificity %>%
          mutate(DC = 0.7,
                 model = "bench",
                 arm = "ctrl_m",
                 specificity = bench_DC_summary[["bench_DC_medium"]][["m_network_confusion_matrix"]][["specificity"]])
        
        DC0.5_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_bench_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_bench_ctrlm_specificity <- DC0.5_bench_ctrlm_specificity %>%
          mutate(DC = 0.5,
                 model = "bench",
                 arm = "ctrl_m",
                 specificity = bench_DC_summary[["bench_DC_large"]][["m_network_confusion_matrix"]][["specificity"]])
        
    ### Graphical
        
      #### Exp
        
        DC1_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_graphical_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_graphical_exp_specificity <- DC1_graphical_exp_specificity %>%
          mutate(DC = 1,
                 model = "model",
                 arm = "exp",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_small"]][["Community_detection_confusion_matrix"]][["specificity"]])
        
        DC0.7_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_graphical_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_graphical_exp_specificity <- DC0.7_graphical_exp_specificity %>%
          mutate(DC = 0.7,
                 model = "model",
                 arm = "exp",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["Community_detection_confusion_matrix"]][["specificity"]])
        
        DC0.5_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_graphical_exp_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_graphical_exp_specificity <- DC0.5_graphical_exp_specificity %>%
          mutate(DC = 0.5,
                 model = "model",
                 arm = "exp",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_large"]][["Community_detection_confusion_matrix"]][["specificity"]])
        
      #### CtrlALL
        
        DC1_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_graphical_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_graphical_ctrlall_specificity <- DC1_graphical_ctrlall_specificity %>%
          mutate(DC = 1,
                 model = "model",
                 arm = "ctrl_all",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_small"]][["Network_confusion_matrix"]][["specificity"]])
        
        DC0.7_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_graphical_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_graphical_ctrlall_specificity <- DC0.7_graphical_ctrlall_specificity %>%
          mutate(DC = 0.7,
                 model = "model",
                 arm = "ctrl_all",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["Network_confusion_matrix"]][["specificity"]])
        
        DC0.5_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_graphical_ctrlall_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_graphical_ctrlall_specificity <- DC0.5_graphical_ctrlall_specificity %>%
          mutate(DC = 0.5,
                 model = "model",
                 arm = "ctrl_all",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_large"]][["Network_confusion_matrix"]][["specificity"]])
        
      #### CtrlM
        
        DC1_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC1_graphical_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC1_graphical_ctrlm_specificity <- DC1_graphical_ctrlm_specificity %>%
          mutate(DC = 1,
                 model = "model",
                 arm = "ctrl_m",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_small"]][["M_Network_confusion_matrix"]][["specificity"]])
        
        DC0.7_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.7_graphical_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.7_graphical_ctrlm_specificity <- DC0.7_graphical_ctrlm_specificity %>%
          mutate(DC = 0.7,
                 model = "model",
                 arm = "ctrl_m",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["M_Network_confusion_matrix"]][["specificity"]])
        
        DC0.5_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
        colnames(DC0.5_graphical_ctrlm_specificity) <- c("DC", "model", "arm", "specificity")
        DC0.5_graphical_ctrlm_specificity <- DC0.5_graphical_ctrlm_specificity %>%
          mutate(DC = 0.5,
                 model = "model",
                 arm = "ctrl_m",
                 specificity = exp_enrichDC_summary[["exp_enrichDC_large"]][["M_Network_confusion_matrix"]][["specificity"]])
        
    ##
        
  df_DC_specificity <- as.data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_DC) <- c("DC", "model", "arm", "specificity")
  
  df_DC_specificity <- df_DC_specificity %>% 
    rbind(DC1_bench_exp_specificity,
          DC0.7_bench_exp_specificity,
          DC0.5_bench_exp_specificity,
          DC1_bench_ctrlall_specificity,
          DC0.7_bench_ctrlall_specificity,
          DC0.5_bench_ctrlall_specificity,
          DC1_bench_ctrlm_specificity,
          DC0.7_bench_ctrlm_specificity,
          DC0.5_bench_ctrlm_specificity,
          DC1_graphical_exp_specificity,
          DC0.7_graphical_exp_specificity,
          DC0.5_graphical_exp_specificity,
          DC1_graphical_ctrlall_specificity,
          DC0.7_graphical_ctrlall_specificity,
          DC0.5_graphical_ctrlall_specificity,
          DC1_graphical_ctrlm_specificity,
          DC0.7_graphical_ctrlm_specificity,
          DC0.5_graphical_ctrlm_specificity)
  
  df_DC_specificity <- df_DC_specificity %>%
    mutate(DC = as.factor(DC),
           model = as.factor(model),
           arm = as.factor(arm))
        
  ## sensitivity
        
    ### bench
        
      #### exp
  
  DC1_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_bench_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_bench_exp_sensitivity <- DC1_bench_exp_sensitivity %>%
    mutate(DC = 1, 
           model = "bench",
           arm = "exp",
           sensitivity = bench_DC_summary[["bench_DC_small"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_bench_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_bench_exp_sensitivity <- DC0.7_bench_exp_sensitivity %>%
    mutate(DC = 0.7, 
           model = "bench",
           arm = "exp",
           sensitivity = bench_DC_summary[["bench_DC_medium"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_bench_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_bench_exp_sensitivity <- DC0.5_bench_exp_sensitivity %>%
    mutate(DC = 0.5, 
           model = "bench",
           arm = "exp",
           sensitivity = bench_DC_summary[["bench_DC_large"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
  
      #### CtrlALL
  
  DC1_bench_ctrlall_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_bench_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_bench_ctrlall_sensitivity <- DC1_bench_ctrlall_sensitivity %>%
    mutate(DC = 1,
           model = "bench",
           arm = "ctrl_all",
           sensitivity = bench_DC_summary[["bench_DC_small"]][["Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_bench_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_bench_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_bench_ctrlall_sensitivity <- DC0.7_bench_ctrlall_sensitivity %>%
    mutate(DC = 0.7,
           model = "bench",
           arm = "ctrl_all",
           sensitivity = bench_DC_summary[["bench_DC_medium"]][["Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_bench_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_bench_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_bench_ctrlall_sensitivity <- DC0.5_bench_ctrlall_sensitivity %>%
    mutate(DC = 0.5,
           model = "bench",
           arm = "ctrl_all",
           sensitivity = bench_DC_summary[["bench_DC_large"]][["Network_confusion_matrix"]][["sensitivity"]])
  
      #### CtrlM
  
  DC1_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_bench_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_bench_ctrlm_sensitivity <- DC1_bench_ctrlm_sensitivity %>%
    mutate(DC = 1,
           model = "bench",
           arm = "ctrl_m",
           sensitivity = bench_DC_summary[["bench_DC_small"]][["m_network_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_bench_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_bench_ctrlm_sensitivity <- DC0.7_bench_ctrlm_sensitivity %>%
    mutate(DC = 0.7,
           model = "bench",
           arm = "ctrl_m",
           sensitivity = bench_DC_summary[["bench_DC_medium"]][["m_network_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_bench_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_bench_ctrlm_sensitivity <- DC0.5_bench_ctrlm_sensitivity %>%
    mutate(DC = 0.5,
           model = "bench",
           arm = "ctrl_m",
           sensitivity = bench_DC_summary[["bench_DC_large"]][["m_network_confusion_matrix"]][["sensitivity"]])
  
    ### graphical
  
      #### exp
  
  DC1_graphical_exp_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_graphical_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_graphical_exp_sensitivity <- DC1_graphical_exp_sensitivity %>%
    mutate(DC = 1,
           model = "model",
           arm = "exp",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_small"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_graphical_exp_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_graphical_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_graphical_exp_sensitivity <- DC0.7_graphical_exp_sensitivity %>%
    mutate(DC = 0.7,
           model = "model",
           arm = "exp",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_graphical_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_graphical_exp_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_graphical_exp_sensitivity <- DC0.5_graphical_exp_sensitivity %>%
    mutate(DC = 0.5,
           model = "model",
           arm = "exp",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_large"]][["Community_detection_confusion_matrix"]][["sensitivity"]])
  
      #### CtrlALL
  
  DC1_graphical_ctrlall_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_graphical_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_graphical_ctrlall_sensitivity <- DC1_graphical_ctrlall_sensitivity %>%
    mutate(DC = 1,
           model = "model",
           arm = "ctrl_all",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_small"]][["Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_graphical_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_graphical_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_graphical_ctrlall_sensitivity <- DC0.7_graphical_ctrlall_sensitivity %>%
    mutate(DC = 0.7,
           model = "model",
           arm = "ctrl_all",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_graphical_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_graphical_ctrlall_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_graphical_ctrlall_sensitivity <- DC0.5_graphical_ctrlall_sensitivity %>%
    mutate(DC = 0.5,
           model = "model",
           arm = "ctrl_all",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_large"]][["Network_confusion_matrix"]][["sensitivity"]])
  
      #### CtrlM
  
  DC1_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC1_graphical_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC1_graphical_ctrlm_sensitivity <- DC1_graphical_ctrlm_sensitivity %>%
    mutate(DC = 1,
           model = "model",
           arm = "ctrl_m",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_small"]][["M_Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.7_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.7_graphical_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.7_graphical_ctrlm_sensitivity <- DC0.7_graphical_ctrlm_sensitivity %>%
    mutate(DC = 0.7,
           model = "model",
           arm = "ctrl_m",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_medium"]][["M_Network_confusion_matrix"]][["sensitivity"]])
  
  DC0.5_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
  colnames(DC0.5_graphical_ctrlm_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  DC0.5_graphical_ctrlm_sensitivity <- DC0.5_graphical_ctrlm_sensitivity %>%
    mutate(DC = 0.5,
           model = "model",
           arm = "ctrl_m",
           sensitivity = exp_enrichDC_summary[["exp_enrichDC_large"]][["M_Network_confusion_matrix"]][["sensitivity"]])
  
  
  ##
  
  df_DC_sensitivity <- as.data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_DC_sensitivity) <- c("DC", "model", "arm", "sensitivity")
  
  df_DC_sensitivity <- df_DC_sensitivity %>% 
    rbind(DC1_bench_exp_sensitivity,
          DC0.7_bench_exp_sensitivity,
          DC0.5_bench_exp_sensitivity,
          DC1_bench_ctrlall_sensitivity,
          DC0.7_bench_ctrlall_sensitivity,
          DC0.5_bench_ctrlall_sensitivity,
          DC1_bench_ctrlm_sensitivity,
          DC0.7_bench_ctrlm_sensitivity,
          DC0.5_bench_ctrlm_sensitivity,
          DC1_graphical_exp_sensitivity,
          DC0.7_graphical_exp_sensitivity,
          DC0.5_graphical_exp_sensitivity,
          DC1_graphical_ctrlall_sensitivity,
          DC0.7_graphical_ctrlall_sensitivity,
          DC0.5_graphical_ctrlall_sensitivity,
          DC1_graphical_ctrlm_sensitivity,
          DC0.7_graphical_ctrlm_sensitivity,
          DC0.5_graphical_ctrlm_sensitivity)

  df_DC_sensitivity <- df_DC_sensitivity %>%
    mutate(DC = as.factor(DC),
           model = as.factor(model),
           arm = as.factor(arm))
  
## DC dataset
  
  df_DC <- as.data.frame(matrix(ncol = 5, nrow = 1800))
  colnames(df_DC) <- c("DC", "model", "arm", "specificity", "sensitivity")
  
 df_DC <- df_DC %>%
   mutate(DC = df_DC_specificity$DC,
          model = df_DC_specificity$model,
          arm = df_DC_specificity$arm,
          specificity = df_DC_specificity$specificity,
          sensitivity = df_DC_sensitivity$sensitivity)
 
save(df_DC, file = "../data/final_df_enrichDC.RData")



#####

load("../data/exp_enrichP.RData")
load("../data/bench_P.RData")

# P

## specificitiy 

### bench 

#### exp

P2_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_exp_specificity) <- c("P", "model", "arm", "specificity")
P2_bench_exp_specificity <- P2_bench_exp_specificity %>%
  mutate(P = 2, 
         model = "bench",
         arm = "exp",
         specificity = bench_P_summary[["bench_P_small"]][["Community_detection_confusion_matrix"]][["specificity"]])

P4_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_exp_specificity) <- c("P", "model", "arm", "specificity")
P4_bench_exp_specificity <- P4_bench_exp_specificity %>%
  mutate(P = 4, 
         model = "bench",
         arm = "exp",
         specificity = bench_P_summary[["bench_P_medium"]][["Community_detection_confusion_matrix"]][["specificity"]])

P8_bench_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_exp_specificity) <- c("P", "model", "arm", "specificity")
P8_bench_exp_specificity <- P8_bench_exp_specificity %>%
  mutate(P = 8, 
         model = "bench",
         arm = "exp",
         specificity = bench_P_summary[["bench_P_large"]][["Community_detection_confusion_matrix"]][["specificity"]])

#### CtrlALL

P2_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P2_bench_ctrlall_specificity <- P2_bench_ctrlall_specificity %>%
  mutate(P = 2,
         model = "bench",
         arm = "ctrl_all",
         specificity = bench_P_summary[["bench_P_small"]][["Network_confusion_matrix"]][["specificity"]])

P4_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P4_bench_ctrlall_specificity <- P4_bench_ctrlall_specificity %>%
  mutate(P = 4,
         model = "bench",
         arm = "ctrl_all",
         specificity = bench_P_summary[["bench_P_medium"]][["Network_confusion_matrix"]][["specificity"]])

P8_bench_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P8_bench_ctrlall_specificity <- P8_bench_ctrlall_specificity %>%
  mutate(P = 8,
         model = "bench",
         arm = "ctrl_all",
         specificity = bench_P_summary[["bench_P_large"]][["Network_confusion_matrix"]][["specificity"]])

#### CtrlM

P2_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P2_bench_ctrlm_specificity <- P2_bench_ctrlm_specificity %>%
  mutate(P = 2,
         model = "bench",
         arm = "ctrl_m",
         specificity = bench_P_summary[["bench_P_small"]][["m_network_confusion_matrix"]][["specificity"]])

P4_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P4_bench_ctrlm_specificity <- P4_bench_ctrlm_specificity %>%
  mutate(P = 4,
         model = "bench",
         arm = "ctrl_m",
         specificity = bench_P_summary[["bench_P_medium"]][["m_network_confusion_matrix"]][["specificity"]])

P8_bench_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P8_bench_ctrlm_specificity <- P8_bench_ctrlm_specificity %>%
  mutate(P = 8,
         model = "bench",
         arm = "ctrl_m",
         specificity = bench_P_summary[["bench_P_large"]][["m_network_confusion_matrix"]][["specificity"]])

### Graphical

#### Exp

P2_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_exp_specificity) <- c("P", "model", "arm", "specificity")
P2_graphical_exp_specificity <- P2_graphical_exp_specificity %>%
  mutate(P = 2,
         model = "model",
         arm = "exp",
         specificity = exp_enrichP_summary[["exp_enrichP_small"]][["Community_detection_confusion_matrix"]][["specificity"]])

P4_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_exp_specificity) <- c("P", "model", "arm", "specificity")
P4_graphical_exp_specificity <- P4_graphical_exp_specificity %>%
  mutate(P = 4,
         model = "model",
         arm = "exp",
         specificity = exp_enrichP_summary[["exp_enrichP_medium"]][["Community_detection_confusion_matrix"]][["specificity"]])

P8_graphical_exp_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_exp_specificity) <- c("P", "model", "arm", "specificity")
P8_graphical_exp_specificity <- P8_graphical_exp_specificity %>%
  mutate(P = 8,
         model = "model",
         arm = "exp",
         specificity = exp_enrichP_summary[["exp_enrichP_large"]][["Community_detection_confusion_matrix"]][["specificity"]])

#### CtrlALL

P2_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P2_graphical_ctrlall_specificity <- P2_graphical_ctrlall_specificity %>%
  mutate(P = 2,
         model = "model",
         arm = "ctrl_all",
         specificity = exp_enrichP_summary[["exp_enrichP_small"]][["Network_confusion_matrix"]][["specificity"]])

P4_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P4_graphical_ctrlall_specificity <- P4_graphical_ctrlall_specificity %>%
  mutate(P = 4,
         model = "model",
         arm = "ctrl_all",
         specificity = exp_enrichP_summary[["exp_enrichP_medium"]][["Network_confusion_matrix"]][["specificity"]])

P8_graphical_ctrlall_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_ctrlall_specificity) <- c("P", "model", "arm", "specificity")
P8_graphical_ctrlall_specificity <- P8_graphical_ctrlall_specificity %>%
  mutate(P = 8,
         model = "model",
         arm = "ctrl_all",
         specificity = exp_enrichP_summary[["exp_enrichP_large"]][["Network_confusion_matrix"]][["specificity"]])

#### CtrlM

P2_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P2_graphical_ctrlm_specificity <- P2_graphical_ctrlm_specificity %>%
  mutate(P = 2,
         model = "model",
         arm = "ctrl_m",
         specificity = exp_enrichP_summary[["exp_enrichP_small"]][["M_Network_confusion_matrix"]][["specificity"]])

P4_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P4_graphical_ctrlm_specificity <- P4_graphical_ctrlm_specificity %>%
  mutate(P = 4,
         model = "model",
         arm = "ctrl_m",
         specificity = exp_enrichP_summary[["exp_enrichP_medium"]][["M_Network_confusion_matrix"]][["specificity"]])

P8_graphical_ctrlm_specificity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_ctrlm_specificity) <- c("P", "model", "arm", "specificity")
P8_graphical_ctrlm_specificity <- P8_graphical_ctrlm_specificity %>%
  mutate(P = 8,
         model = "model",
         arm = "ctrl_m",
         specificity = exp_enrichP_summary[["exp_enrichP_large"]][["M_Network_confusion_matrix"]][["specificity"]])

##

df_P_specificity <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_P_specificity) <- c("P", "model", "arm", "specificity")

df_P_specificity <- df_P_specificity %>% 
  rbind(P2_bench_exp_specificity,
        P4_bench_exp_specificity,
        P8_bench_exp_specificity,
        P2_bench_ctrlall_specificity,
        P4_bench_ctrlall_specificity,
        P8_bench_ctrlall_specificity,
        P2_bench_ctrlm_specificity,
        P4_bench_ctrlm_specificity,
        P8_bench_ctrlm_specificity,
        P2_graphical_exp_specificity,
        P4_graphical_exp_specificity,
        P8_graphical_exp_specificity,
        P2_graphical_ctrlall_specificity,
        P4_graphical_ctrlall_specificity,
        P8_graphical_ctrlall_specificity,
        P2_graphical_ctrlm_specificity,
        P4_graphical_ctrlm_specificity,
        P8_graphical_ctrlm_specificity)

df_P_specificity <- df_P_specificity %>%
  mutate(P = as.factor(P),
         model = as.factor(model),
         arm = as.factor(arm))

  
## sensitivity

### bench

#### exp

P2_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_bench_exp_sensitivity <- P2_bench_exp_sensitivity %>%
  mutate(P = 2, 
         model = "bench",
         arm = "exp",
         sensitivity = bench_P_summary[["bench_P_small"]][["Community_detection_confusion_matrix"]][["sensitivity"]])

P4_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_bench_exp_sensitivity <- P4_bench_exp_sensitivity %>%
  mutate(P = 4, 
         model = "bench",
         arm = "exp",
         sensitivity = bench_P_summary[["bench_P_medium"]][["Community_detection_confusion_matrix"]][["sensitivity"]])

P8_bench_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_bench_exp_sensitivity <- P8_bench_exp_sensitivity %>%
  mutate(P = 8, 
         model = "bench",
         arm = "exp",
         sensitivity = bench_P_summary[["bench_P_large"]][["Community_detection_confusion_matrix"]][["sensitivity"]])


#### CtrlALL

P2_bench_ctrlall_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_bench_ctrlall_sensitivity <- P2_bench_ctrlall_sensitivity %>%
  mutate(P = 2,
         model = "bench",
         arm = "ctrl_all",
         sensitivity = bench_P_summary[["bench_P_small"]][["Network_confusion_matrix"]][["sensitivity"]])

P4_bench_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_bench_ctrlall_sensitivity <- P4_bench_ctrlall_sensitivity %>%
  mutate(P = 4,
         model = "bench",
         arm = "ctrl_all",
         sensitivity = bench_P_summary[["bench_P_medium"]][["Network_confusion_matrix"]][["sensitivity"]])

P8_bench_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_bench_ctrlall_sensitivity <- P8_bench_ctrlall_sensitivity %>%
  mutate(P = 8,
         model = "bench",
         arm = "ctrl_all",
         sensitivity = bench_P_summary[["bench_P_large"]][["Network_confusion_matrix"]][["sensitivity"]])

#### CtrlM

P2_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_bench_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_bench_ctrlm_sensitivity <- P2_bench_ctrlm_sensitivity %>%
  mutate(P = 2,
         model = "bench",
         arm = "ctrl_m",
         sensitivity = bench_P_summary[["bench_P_small"]][["m_network_confusion_matrix"]][["sensitivity"]])

P4_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_bench_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_bench_ctrlm_sensitivity <- P4_bench_ctrlm_sensitivity %>%
  mutate(P = 4,
         model = "bench",
         arm = "ctrl_m",
         sensitivity = bench_P_summary[["bench_P_medium"]][["m_network_confusion_matrix"]][["sensitivity"]])

P8_bench_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_bench_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_bench_ctrlm_sensitivity <- P8_bench_ctrlm_sensitivity %>%
  mutate(P = 8,
         model = "bench",
         arm = "ctrl_m",
         sensitivity = bench_P_summary[["bench_P_large"]][["m_network_confusion_matrix"]][["sensitivity"]])

### graphical

#### exp

P2_graphical_exp_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_graphical_exp_sensitivity <- P2_graphical_exp_sensitivity %>%
  mutate(P = 2,
         model = "model",
         arm = "exp",
         sensitivity = exp_enrichP_summary[["exp_enrichP_small"]][["Community_detection_confusion_matrix"]][["sensitivity"]])

P4_graphical_exp_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_graphical_exp_sensitivity <- P4_graphical_exp_sensitivity %>%
  mutate(P = 4,
         model = "model",
         arm = "exp",
         sensitivity = exp_enrichP_summary[["exp_enrichP_medium"]][["Community_detection_confusion_matrix"]][["sensitivity"]])

P8_graphical_exp_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_exp_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_graphical_exp_sensitivity <- P8_graphical_exp_sensitivity %>%
  mutate(P = 8,
         model = "model",
         arm = "exp",
         sensitivity = exp_enrichP_summary[["exp_enrichP_large"]][["Community_detection_confusion_matrix"]][["sensitivity"]])

#### CtrlALL

P2_graphical_ctrlall_sensitivity= as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_graphical_ctrlall_sensitivity <- P2_graphical_ctrlall_sensitivity %>%
  mutate(P = 2,
         model = "model",
         arm = "ctrl_all",
         sensitivity = exp_enrichP_summary[["exp_enrichP_small"]][["Network_confusion_matrix"]][["sensitivity"]])

P4_graphical_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_graphical_ctrlall_sensitivity <- P4_graphical_ctrlall_sensitivity %>%
  mutate(P = 4,
         model = "model",
         arm = "ctrl_all",
         sensitivity = exp_enrichP_summary[["exp_enrichP_medium"]][["Network_confusion_matrix"]][["sensitivity"]])

P8_graphical_ctrlall_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_ctrlall_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_graphical_ctrlall_sensitivity <- P8_graphical_ctrlall_sensitivity %>%
  mutate(P = 8,
         model = "model",
         arm = "ctrl_all",
         sensitivity = exp_enrichP_summary[["exp_enrichP_large"]][["Network_confusion_matrix"]][["sensitivity"]])

#### CtrlM

P2_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P2_graphical_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P2_graphical_ctrlm_sensitivity <- P2_graphical_ctrlm_sensitivity %>%
  mutate(P = 2,
         model = "model",
         arm = "ctrl_m",
         sensitivity = exp_enrichP_summary[["exp_enrichP_small"]][["M_Network_confusion_matrix"]][["sensitivity"]])

P4_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P4_graphical_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P4_graphical_ctrlm_sensitivity <- P4_graphical_ctrlm_sensitivity %>%
  mutate(P = 4,
         model = "model",
         arm = "ctrl_m",
         sensitivity = exp_enrichP_summary[["exp_enrichP_medium"]][["M_Network_confusion_matrix"]][["sensitivity"]])

P8_graphical_ctrlm_sensitivity = as.data.frame(matrix(ncol = 4, nrow = 100))
colnames(P8_graphical_ctrlm_sensitivity) <- c("P", "model", "arm", "sensitivity")
P8_graphical_ctrlm_sensitivity <- P8_graphical_ctrlm_sensitivity %>%
  mutate(P = 8,
         model = "model",
         arm = "ctrl_m",
         sensitivity = exp_enrichP_summary[["exp_enrichP_large"]][["M_Network_confusion_matrix"]][["sensitivity"]])


##

df_P_sensitivity <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_P_sensitivity) <- c("P", "model", "arm", "sensitivity")

df_P_sensitivity <- df_P_sensitivity %>% 
  rbind(P2_bench_exp_sensitivity,
        P4_bench_exp_sensitivity,
        P8_bench_exp_sensitivity,
        P2_bench_ctrlall_sensitivity,
        P4_bench_ctrlall_sensitivity,
        P8_bench_ctrlall_sensitivity,
        P2_bench_ctrlm_sensitivity,
        P4_bench_ctrlm_sensitivity,
        P8_bench_ctrlm_sensitivity,
        P2_graphical_exp_sensitivity,
        P4_graphical_exp_sensitivity,
        P8_graphical_exp_sensitivity,
        P2_graphical_ctrlall_sensitivity,
        P4_graphical_ctrlall_sensitivity,
        P8_graphical_ctrlall_sensitivity,
        P2_graphical_ctrlm_sensitivity,
        P4_graphical_ctrlm_sensitivity,
        P8_graphical_ctrlm_sensitivity)

df_P_sensitivity <- df_P_sensitivity %>%
  mutate(P = as.factor(P),
         model = as.factor(model),
         arm = as.factor(arm))

## P dataset

df_P <- as.data.frame(matrix(ncol = 5, nrow = 1800))
colnames(df_P) <- c("P", "model", "arm", "specificity", "sensitivity")

df_P <- df_P %>%
  mutate(P = df_P_specificity$P,
         model = df_P_specificity$model,
         arm = df_P_specificity$arm,
         specificity = df_P_specificity$specificity,
         sensitivity = df_P_sensitivity$sensitivity)

save(df_P, file = "../data/final_df_enrichP.RData")


#####

Figures 
load("../data/final_df_enrichDC.Rdata")
load("../data/final_df_enrichP.Rdata")

g1 <- ggplot(data = df_DC, aes(x = arm, y = specificity, fill = model)) + geom_boxplot() + facet_wrap(df_DC$DC)
plot(g1)

g2 <- ggplot(data = df_DC, aes(x = arm, y = sensitivity, fill = model)) + geom_boxplot() + facet_wrap(df_DC$DC)
plot(g2)

g3 <- ggplot(data = df_P, aes(x = arm, y = specificity, fill = model)) + geom_boxplot() + facet_wrap(df_P$P)
plot(g3)

g4 <- ggplot(data = df_P, aes(x = arm, y = sensitivity, fill = model)) + geom_boxplot() + facet_wrap(df_P$P)
plot(g4)
