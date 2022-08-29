library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(tidyverse)
library(table1)
library(sharp)


### Figures topology analysis

load("../data/final_df_topsmall.RData")
load("../data/final_df_toplarge.RData")

## 

levels(topological_small$nu) <- c("disconnected (nu = 0)", "sparsely connected (nu = 0.01)", "highly connected (nu = 0.05)")

g1 <- ggplot(data = topological_small, aes(y = prop_true, x = model)) + 
  geom_boxplot(data = topological_small, aes(fill = model), width = 0.5) + 
  facet_wrap(topological_small$nu) + 
  geom_hline(aes(yintercept=0.25, linetype = str_wrap("Probability of detecting the true pathway membership of a given variable by chance alone", 20)), show.legend = TRUE) +
  scale_linetype_manual(name = "", values = c(2, 2)) +
  stat_compare_means( aes(label = ..p.signif..),comparisons = list(c("bench", "graphical")), method = "t.test") +
  theme_gray() +
  scale_fill_npg() +
  ggtitle("A) Proportion of accurately detected pathway memberships for m = 4 small sized pathways with size pk = 20") +
  xlab("") +
  ylab("Proportion of accurately detected pathway memberships") +
  ylim(0.2,1) +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
 
plot(g1)

label(topological_small$prop_true) <- str_wrap("Proportion of accurately detected pathway memberships",20)

top_small_table <- table1(~ prop_true |nu + model, data = topological_small, overall = F)
top_small_table

## large 

g2 <- ggplot(data = topological_large, aes(y = prop_true, x = model)) + 
  geom_boxplot(data = topological_large, aes(fill = model), width = 0.5) + 
  facet_wrap(topological_small$nu) + 
  geom_hline(aes(yintercept=0.25, linetype = str_wrap("Probability of detecting the true pathway membership of a given variable by chance alone", 20)), show.legend = TRUE) +
  scale_linetype_manual(name = "", values = c(2, 2)) +
  stat_compare_means( aes(label = ..p.signif..),comparisons = list(c("bench", "graphical")), method = "t.test") +
  theme_gray() +
  scale_fill_npg() +
  ggtitle("B) Proportion of accurately detected pathway memberships for m = 4 large sized pathways with size pk = 50") +
  xlab("") +
  ylab("Proportion of accurately detected pathway memberships") +
  ylim(0.2,1.05) +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

plot(g2)

top_large_table <- table1(~ prop_true | nu + model , data = topological_large, overall = F)
top_large_table

# Figures enrichment analysis
load("../data/final_df_enrichDC.Rdata")
load("../data/final_df_enrichP.Rdata")

summary(df_DC$model)

levels(df_DC$model) <- c("bench", "graphical")

df_DC_1 <- df_DC[which(df_DC$DC == 1),]
df_DC_0.7<- df_DC[which(df_DC$DC == 0.7),]
df_DC_0.5<- df_DC[which(df_DC$DC == 0.5),]

DC_1 <- table1(~ specificity + sensitivity | arm + model, data = df_DC_1, overall = F)
DC_0.7 <- table1(~ specificity + sensitivity | arm + model, data = df_DC_0.7, overall = F)
DC_0.5 <- table1(~ specificity + sensitivity | arm + model, data = df_DC_0.5, overall = F)

summary(df_P$model)
levels(df_P$model) <- c("bench", "graphical")

df_P_2 <- df_P[which(df_P$P == 2),]
df_P_4 <- df_P[which(df_P$P == 4),]
df_P_8 <- df_P[which(df_P$P == 8),]

P_2 <- table1(~ specificity + sensitivity | arm + model, data = df_P_2, overall = F)
P_4 <- table1(~ specificity + sensitivity | arm + model, data = df_P_4, overall = F)
P_8 <- table1(~ specificity + sensitivity | arm + model, data = df_P_8, overall = F)
P_8

### graphics enrichment

library(ggstatsplot)

str(df_P)

levels(df_DC$DC) <- c("DC = 0.5", "DC = 0.7", "DC = 1.0")
levels(df_P$P) <- c("m = 2", "m = 4", "m = 8")



g3 <- ggplot(data = df_DC, aes(x = arm, y = specificity, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ df_DC$DC)+
  theme_gray() +
  scale_fill_npg() +
  ggtitle("A) Specificity of Pathway Enrichment Analysis methods dependent on Detection Call (DC)") +
  xlab("") +
  ylab("Specificity") +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

plot(g3)
  
g4 <- ggplot(data = df_DC, aes(x = arm, y = sensitivity, fill = model)) +
  geom_boxplot() +
  facet_wrap(df_DC$DC) +
  theme_gray() +
  scale_fill_npg() +
  ggtitle("B) Sensitivity of Pathway Enrichment Analysis methods dependent on Detection Call (DC)") +
  xlab("") +
  ylab("Sensitvity") +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

plot(g4)

g5 <- ggplot(data = df_P, aes(x = arm, y = specificity, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ df_P$P)+
  theme_gray() +
  scale_fill_npg() +
  ggtitle("A) Specificity of Pathway Enrichment Analysis methods dependent on the number of enriched pathways m") +
  xlab("") +
  ylab("Specificity") +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

plot(g5)

g6 <- ggplot(data = df_P, aes(x = arm, y = sensitivity, fill = model)) +
  geom_boxplot() +
  facet_wrap(df_P$P) +
  theme_gray() +
  scale_fill_npg() +
  ggtitle("B) Sensitivity of Pathway Enrichment Analysis methods dependent on the number of enriched pathways m") +
  xlab("") +
  ylab("Sensitvity") +
  theme(plot.title = element_text(size = 14, vjust = 3),
        axis.title.y = element_text(vjust= 3),
        axis.text=element_text(size=14), 
        axis.text.x=element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

plot(g6)


