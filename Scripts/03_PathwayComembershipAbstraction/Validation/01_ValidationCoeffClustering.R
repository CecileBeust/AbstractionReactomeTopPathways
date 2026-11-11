library(igraph)
library(ggplot2)
library(dplyr)
library(purrr)

##############################################
# Load Reactome weighted comembership clique #
##############################################

graph_data <- read.csv("../../../Results/PathwayComembership/06_ScalingReactome/Reactome96WeightedComembershipClique.csv")
head(graph_data)
#g <- graph_from_data_frame(
#  graph_data[, c("prot1", "prot2", "comembershipScore")],
#  directed = FALSE
#)
g <- graph_from_data_frame(
  d = graph_data,
  directed = FALSE
)
g <- set_edge_attr(g, "weight", value = graph_data$comembershipScore)
is_weighted(g)

#g <- set_edge_attr(g, "provenance", value = graph_data$provenanceScore)
E(g)$weight
#E(g)$provenance 

str(g)
#############################################
########### Load cancer hallmarks ###########
#############################################

hallmarks <- list(
  TISSUE_INVASION_METASTASIS = read.csv("../../../Data/ValidationData/CancerHallmarks/TissueInvasionAndMetastasis-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  REPROGRAMMING_ENERGY_METABOLISM = read.csv("../../../Data/ValidationData/CancerHallmarks/ReprogrammingEnergyMetabolism-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  REPLICATIVE_IMMORTALITY = read.csv("../../../Data/ValidationData/CancerHallmarks/ReplicativeImmortality-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  GENOME_INSTABILITY = read.csv("../../../Data/ValidationData/CancerHallmarks/GenomeInstability-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  RESISTING_CELL_DEATH = read.csv("../../../Data/ValidationData/CancerHallmarks/ResistingCellDeath-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  EVADING_GROWTH_SUPPRESSORS = read.csv("../../../Data/ValidationData/CancerHallmarks/EvadingGrowthSuppressors-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  EVADING_IMMUNE_DESTRUCTION = read.csv("../../../Data/ValidationData/CancerHallmarks/EvadingImmuneDestruction-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  SUSTAINED_ANGIOGENESIS = read.csv("../../../Data/ValidationData/CancerHallmarks/SustainedAngiogenesis-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  SUSTAINING_PROLIFERATIVE_SIGNALING = read.csv("../../../Data/ValidationData/CancerHallmarks/SustainingProliferativeSignaling-UniProt.tsv", sep = "\t", header = TRUE)$Entry,
  TUMOR_PROMOTING_INFLAMMATION = read.csv("../../../Data/ValidationData/CancerHallmarks/TumorPromotingInflammation-UniProt.tsv", sep = "\t", header = TRUE)$Entry
)

#################################################################
# Filter hallmarks that are in the Reactome comembership clique #
#################################################################

all_proteins_in_graph <- V(g)$name
length(all_proteins_in_graph)
hallmarks_filtered <- lapply(hallmarks, function(hallmark) {
  intersect(hallmark, all_proteins_in_graph)
})

#################################
# Create subgraphs of hallmarks #
#################################

subgraphs <- lapply(hallmarks_filtered, function(hallmark_proteins) {
  induced_subgraph(g, V(g)[hallmark_proteins], impl = 'create_from_scratch')
})
names(subgraphs) <- names(hallmarks_filtered)
subgraphs$TISSUE_INVASION_METASTASIS

#########################################################
# Compute Barrat weighted global clustering coefficient #
#########################################################

#E(g)$weight
#plot(subgraphs$TISSUE_INVASION_METASTASIS)
#length(subgraphs$TISSUE_INVASION_METASTASIS)
#g2 <- subgraphs$TISSUE_INVASION_METASTASIS
#e2 <- E(g2)$weight
#e2

#transitivity(graph = g2, type="barrat", weights = e2)
#transitivity(graph = g2, type="weighted", weights = e2)
#transitivity(subgraphs$TISSUE_INVASION_METASTASIS, type="globalundirected")

#barrat_coefficients <- sapply(names(subgraphs), function(hm_name) {
#  subgraph <- subgraphs[[hm_name]]
#  edge_weights <- E(subgraph)$weight
#  edge_weights
#  mean(transitivity(subgraph, type = "barrat", weights = edge_weights), na.rm = TRUE)
#})
#barrat_coefficients

##########################################
# COMPUTE DENSITY IN HALLMARKS SUBGRAPHS #
##########################################

density <- sapply(names(subgraphs), function(hm_name) {
  subgraph <- subgraphs[[hm_name]]
  edge_weights <- E(subgraph)$weight
  2 * sum(edge_weights) / (vcount(subgraph) * (vcount(subgraph) - 1))
})
density

###########################################

generate_random_subgraph_weighted <- function(graph, nodes, n_iter = 1000) {
  # Extraire le sous-graphe induit
  subgraph <- induced_subgraph(graph, V(graph)[nodes])
  
  # Si le sous-graphe a moins de 3 nœuds, retourner NA (densité indéfinie)
  if (vcount(subgraph) < 3) {
    return(rep(NA, n_iter))
  }
  
  # Extraire les poids et la distribution de degrés
  edge_weights <- E(subgraph)$weight
  degrees <- degree(subgraph, mode = "all")
  
  # Calculer la densité pondérée observée
  observed_density <- 2 * sum(edge_weights) / (vcount(subgraph) * (vcount(subgraph) - 1))
  
  # Générer des sous-graphes aléatoires avec même distribution de poids
  random_densities <- map_dbl(1:n_iter, ~ {
    # Créer un graphe aléatoire avec la même distribution de degrés (configuration model)
    random_graph <- sample_degseq(degrees, method = "vl")
    
    # Assigner les poids aléatoirement (mais conserver la distribution)
    new_edge_weights <- sample(edge_weights, ecount(random_graph), replace = FALSE)
    set_edge_attr(random_graph, "weight", value = new_edge_weights)
    
    # Calculer la densité pondérée
    2 * sum(new_edge_weights) / (vcount(random_graph) * (vcount(random_graph) - 1))
  })
  
  return(random_densities)
}

# Calculer les densités pondérées observées
density_observed <- sapply(names(subgraphs), function(hm_name) {
  subgraph <- subgraphs[[hm_name]]
  edge_weights <- E(subgraph)$weight
  2 * sum(edge_weights) / (vcount(subgraph) * (vcount(subgraph) - 1))
})

# Générer des densités aléatoires pour chaque hallmark
set.seed(123)  # Pour la reproductibilité
random_densities_list <- map(names(subgraphs), function(hm_name) {
  tryCatch(
    {
      generate_random_subgraph_weighted(g, hallmarks_filtered[[hm_name]], n_iter = 1000)
    },
    error = function(e) {
      message("Erreur pour ", hm_name, ": ", e$message)
      return(rep(NA, 1000))  # Retourne des NA si échec
    }
  )
})

# Calculer les p-valeurs (proportion de densités aléatoires ≥ densité observée)
p_values <- map2_dbl(
  names(subgraphs),
  random_densities_list,
  ~ {
    observed <- density_observed[.x]
    random_densities <- .y
    valid_random <- random_densities[!is.na(random_densities)]
    if (length(valid_random) == 0) {
      return(NA)  # Pas de comparaison possible
    }
    mean(valid_random >= observed, na.rm = TRUE)
  }
)

# Préparer les données pour la visualisation
results_df <- data.frame(
  Hallmark = rep(names(subgraphs), each = 1001),  # 1000 aléatoires + 1 observée
  Density = c(unlist(random_densities_list), density_observed),
  Type = rep(c(rep("Random", 1000), "Observed"), length(names(subgraphs)))
)

# Visualisation avec ggplot2
ggplot(results_df, aes(x = Density, fill = Type)) +
  geom_histogram(data = filter(results_df, Type == "Random"), bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(aes(xintercept = Density[Type == "Observed"], color = "Observed"), linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Hallmark, scales = "free") +
  labs(title = "Comparaison des densités pondérées observées vs aléatoires (distribution de poids équivalente)",
       x = "Densité pondérée",
       y = "Fréquence") +
  theme_minimal() +
  theme(legend.position = "top")

# Afficher les p-valeurs
print(p_values)
















