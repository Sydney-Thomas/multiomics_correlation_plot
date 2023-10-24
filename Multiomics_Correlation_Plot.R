library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(randomForest)
library(pROC)


setwd('C:/Users/sthomas/OneDrive - University of California, San Diego Health/Reference Files/Example R Scripts/multiomics_correlation')

## Read in metabolomics data 
## Data has been preprocessed using rclr and missing data imputation
metabolites_US <- read.csv("./Example_Data/Metabolites_normalized_imputed.csv")

## Read in microbiome data
## Data is preprocessed using rclr and annotated at genus level
microbiome_US <- read.csv("./Example_Data/Microbiome_table_SEED_genus_rclr.csv")

## Read in HMO data
## HMO data is scaled and centered
HMO <- read.csv("./Example_Data/HMO_scaled_centered.csv")

## For this graph, we just want to look at how metabolites compare between the HMOs 6DDL and microbes Rothia/Veillonella
## So we'll filter the microbiome data and HMO data to reduce dataset size
## If you want to look at everything, you can skip this step
microbiome_US <- microbiome_US %>% dplyr::select(Sample, contains("othia"), contains("eillon"))
HMO <- HMO %>% dplyr::select(Sample, z6DDL)

## Merge datasets
US <- inner_join(metabolites_US, microbiome_US, by = "Sample")
US <- inner_join(HMO, US, by = "Sample")

## Define groups 
## These are the variables that you want to compare everything against
## For this experiment, these groups are the microbes Rothia and Veillonella, and the HMOs 6DDL
## Anything should work here as long as it's a continuous variable
G1 <- "z6DDL"
G2 <- "Rothia"
G3 <- "Veillonella"

## Do spearman correlations
US <- US %>% dplyr::select(Sample, !!sym(G1), !!sym(G2), !!sym(G3), everything())
US <- US %>% pivot_longer(5:ncol(US), names_to = "To", values_to = "To_Abundance")
US <- US %>% pivot_longer(2:4, names_to = "From", values_to = "From_Abundance")
US_corr <- US %>% group_by(To, From) %>% rstatix::cor_test(vars=c("To_Abundance"), vars2=c("From_Abundance"), method="spearman")

## BH adjust p values
US_corr$padj <- p.adjust(US_corr$p, method = "BH")

## Find metabolites that are significantly correlated with Rothia, Veillonella, and 6DDL
## This filters out metabolites that are significantly associated with all three 
US_corr <- US_corr %>%
  group_by(To) %>%
  filter(all(padj < 0.05)) %>%
  ungroup()

## Make table of library and canopus IDs
## These dataframes come from iimn collapsed outputs from GNPS and Sirius
## For more information, see this script: https://github.com/Sydney-Thomas/iimn_clean_normalize
lib_US <- US_corr %>% dplyr::select(To, From, cor)
lib_US <- lib_US %>% pivot_wider(names_from = "From", values_from = "cor")
lib <- read.csv("./Example_Data/library_matches.csv")
lib$Metabolite <- paste0("X", lib$Metabolite)
can <- read.csv("./Example_Data/canopus_matches.csv")
can$Metabolite <- paste0("X", can$Metabolite)
lib_US <- right_join(lib, lib_US, by = c("Metabolite" = "To"))
lib_US <- left_join(lib_US, can, by = "Metabolite")
write_csv(lib_US, "Significant_Metabolites.csv")
rm(lib, can)

## Make graph ###############
## Generate network table
network_table <- US_corr %>% ungroup() %>% dplyr::select(From, To, cor, everything())
network_table$direction <- ifelse(network_table$cor < 0, "neg", "pos")
network_table$weight <- abs(network_table$cor)
network_table$sig <- ifelse(network_table$padj<0.001, 1, ifelse(network_table$padj<0.01, 2, ifelse(network_table$padj<0.05, 3, 4)))
network_table <- network_table %>% unique()

nodes <- network_table %>% arrange(From, cor) %>% dplyr::select(From) %>% unique()
nodes_2 <- network_table %>% arrange(From, cor) %>% filter(!To %in% nodes$From) %>% dplyr::select(To) %>% unique() %>% rename(From = To)
nodes <- rbind(nodes, nodes_2)
names(nodes)[1] <- "node"
rm(nodes_2)
nodes <- nodes %>% unique()

## Add in metabolite colors 
nodes <- left_join(nodes, lib_US, by = c("node" = "Metabolite"), multiple = "first")
## Here you can specify which canopus category (class, subclass) you want to visualize
nodes <- nodes %>% dplyr::select(node, ClassyFire.level.5)

## Add in x/y coordinates (for manual layouts)
## For this manual layout, we are putting Rothia and Veillonella on one side and 6DDL on the other
## You'll need to change this manually depending on which factors you want in each column 
nodes$x <- ifelse(str_detect(nodes$node, "othia|eillon"), 1, ifelse(str_detect(nodes$node, "^X"), 2, 3))
nodes <- nodes %>% group_by(x) %>% mutate(y = row_number())
nodes$y <- ifelse(str_detect(nodes$node, "othia"), (nrow(nodes)/3), 
           ifelse(str_detect(nodes$node, "eillon"), (2*nrow(nodes)/3), 
           ifelse(str_detect(nodes$node, "z6DDL"), (nrow(nodes)/2), nodes$y)))

## Generate network
network_seed <- graph_from_data_frame(d = network_table, vertices = nodes, directed = FALSE)


# Create circular correlation plot
SEED_Network <- ggraph(network_seed, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(color = cor, edge_linetype = as.factor(sig)), edge_width = 1) +
  scale_edge_color_gradientn(colors = brewer.pal(n=8, name="RdBu")) +
  geom_node_point(aes(color = ClassyFire.level.5)) +
  theme_void() +
  #theme(legend.position = "none") +
  geom_node_text(aes(label = name), size = 3, repel = TRUE, max.overlaps = Inf) +
  labs(edge_color = "Spearman")
SEED_Network


# Create column correlation plot
V(network_seed)$x <- nodes$x
V(network_seed)$y <- nodes$y
layout_manual <- as.matrix(nodes[,c("x","y")])

SEED_Network <- ggraph(network_seed, layout = layout_manual) +
  geom_edge_link(aes(color = cor, edge_linetype = as.factor(sig)), edge_width = 1) +
  scale_edge_color_gradientn(colors = brewer.pal(n=8, name="RdBu")) +
  geom_node_point(aes(x = x, y = y, color = ClassyFire.level.5)) +
  theme_void() +
  geom_node_text(aes(x = x, y = y, label = name), size = 3, repel = TRUE, max.overlaps = Inf) +
  labs(edge_color = "Spearman")
SEED_Network
