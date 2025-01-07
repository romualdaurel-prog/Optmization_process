source("Rscript/0_load_library.R")
source("Rscript/0_load_functions.R")

run_traceplot <- FALSE
run_WAIC <- TRUE
model_name <- "base_model"

modelfile <- "model/model_scorff_LCM_v1.R"
datafile <- "input_data/data_scorff_LCM_v2.rds"
Constfile <- "input_data/Const_scorff_LCM_v2.rds"
configfile <- "Rscript/1_configuration_base.R"
initfile <- "Rscript/0_generate_inits_base.R"
Data_nimble <- readRDS(datafile)
Const_nimble <- readRDS(Constfile)

n.iter <- 500e3 # total number of model iterations
n.burnin  <- 0 # burnin define in total number of model iteration
n.chains  <- 3 # number of parallel chains
n.thin  <- 250 # thin applied to keep a reasonable number of posteriors

burnin_output <- 200 # burnin defined relative to number of posterior used in outputs only

# Should chains be run in parallel
parallel_run <- TRUE

# Create log file -----------------------------------------------------
if (!dir.exists("log")) {
  dir.create("log")
}

# Names files to store results ----------------------------------------


# 1 - Configuration & Model compilation -------------------------------
source(configfile)

# 2 - Run Model -------------------------------------------------------
# Model will not be run if name_file_mcmc exists
source("Rscript/2_run_model.R")


# 3 - Convergence Analysis & plot--------------------------------------
source("Rscript/3_convergence.R")

# 4 - extract median & plot parameters  -------------------------------
source("Rscript/4_parameters_summary_and_plot.R")

# 5 - Fit assessment  -------------------------------------------------
source("Rscript/5_fit_assessment.R")

# 6 - WAIC ------------------------------------------------------------
if(run_WAIC){
  source("Rscript/6_WAIC.R")
}

# 7 - Traceplot -------------------------------------------------------
if(run_traceplot){
  source("Rscript/7_traceplot.R")
}

# Get sampler times for each node -------------------------------------
#sampler_times <- compiled.MCMC.model$getTimes()  # Get the time for each node update
#print(sampler_times)
#sampler_times <- MCMC_samplertime
#print(sampler_times)
compiled.MCMC.model$run(1000,time=TRUE)
sampler_times <-compiled.MCMC.model$getTimes()

###Visualisation des temps d'execution des samplers
# Tracer un histogramme des temps
times <- sampler_times
# Scatter plot des temps avec lignes et courbes
plot(times, type = "p", main = "Temps d'exécution par sampler", 
     xlab = "Index du sampler", ylab = "Temps (secondes)", 
     col = "darkgreen", pch = 19)

# Ajouter une ligne reliant les points (ligne brute)
lines(times, col = "blue", lwd = 2)

# Ajouter une courbe lissée (loess)
lines(smooth.spline(times), col = "red", lwd = 2, lty = 2)

# Ajouter une légende
legend("topright", 
       legend = c("Temps par sampler", "Ligne brute", "Courbe lissée (Loess)"),
       col = c("darkgreen", "blue", "red"), 
       pch = c(19, NA, NA), 
       lty = c(NA, 1, 2), 
       lwd = c(NA, 2, 2))

# Calcul automatique du nombre de lignes
nrow_auto <- ceiling(length(times) / 10)  # Ajustez la valeur 10 selon le besoin

# Créer une matrice compatible
heatmap_data <- matrix(times, nrow = nrow_auto, byrow = TRUE)

# Tracer le heatmap
heatmap(heatmap_data, 
        main = "Heatmap des temps d'exécution des samplers", 
        xlab = "Samplers", 
        ylab = "Groupes", 
        col = colorRampPalette(c("yellow", "orange", "red"))(50), 
        margins = c(5, 5))

# Ajouter une légende
legend("topright", 
       legend = c("Temps faible", "Temps élevé"), 
       fill = c("yellow", "red"), 
       title = "Temps d'exécution")

# 3 - Load MCMC results and check format ------------------------------
results <- readRDS("outputs_hindcast/mcmc_results_base_model.rds")
str(results)

# Convert MCMC results to a 3D array format for diagnostics
if (is.list(results) && length(dim(results[[1]])) == 2) {
  iterations <- nrow(results[[1]])
  parameters <- ncol(results[[1]])
  chains <- length(results)
  mcmc_array <- array(NA, dim = c(iterations, chains, parameters))
  for (chain in 1:chains) {
    mcmc_array[, chain, ] <- as.matrix(results[[chain]])
  }
  parameter_names <- colnames(results[[1]])
  dimnames(mcmc_array) <- list(NULL, NULL, parameter_names)
} else {
  stop("Error: 'results' is not in the expected format for conversion to an array.")
}

# Run convergence diagnostics using effectiveSize from coda
mcmc_list <- lapply(results, as.mcmc)

##SUPPRESSION DES NA
# Fonction pour nettoyer les lignes contenant des NA dans toutes les chaînes
cleaned_mcmc_list <- function(mcmc_list) {
  # Vérifier que toutes les chaînes ont initialement les mêmes dimensions
  initial_dims <- lapply(mcmc_list, dim)
  if (!all(sapply(initial_dims, function(x) all(x == initial_dims[[1]])))) {
    stop("Les dimensions initiales des chaînes ne sont pas identiques.")
  }
  
  # Identifier les lignes à conserver (aucun NA dans toutes les chaînes)
  rows_to_keep <- Reduce(`&`, lapply(mcmc_list, function(chain) {
    complete.cases(chain)  # TRUE si une ligne ne contient aucun NA
  }))
  
  # Supprimer les lignes contenant des NA dans toutes les chaînes
  mcmc_list_cleaned <- lapply(mcmc_list, function(chain) {
    chain[rows_to_keep, , drop = FALSE]  # Garder uniquement les lignes valides
  })
  
  # Vérifier que toutes les chaînes nettoyées ont les mêmes dimensions
  cleaned_dims <- lapply(mcmc_list_cleaned, dim)
  if (!all(sapply(cleaned_dims, function(x) all(x == cleaned_dims[[1]])))) {
    stop("Les dimensions des chaînes ne sont pas identiques après le nettoyage.")
  }
  
  return(mcmc_list_cleaned)
}


# Appliquer le nettoyage
mcmc_list_cleaned <- cleaned_mcmc_list(mcmc_list)

# Vérification des résultats
print("Chaînes nettoyées :")
print(mcmc_list_cleaned)

###VERIFICATION
# Fonction de vérification
verify_cleaned_mcmc_list <- function(mcmc_list) {
  # Vérifier qu'il n'y a plus de NA dans aucune chaîne
  no_na <- all(sapply(mcmc_list, function(chain) {
    all(!is.na(chain))
  }))
  
  if (!no_na) {
    stop("Erreur : Des NA sont encore présents dans l'une ou plusieurs des chaînes.")
  }
  
  # Vérifier que toutes les chaînes ont les mêmes dimensions
  dims <- lapply(mcmc_list, dim)
  consistent_dims <- all(sapply(dims, function(x) all(x == dims[[1]])))
  
  if (!consistent_dims) {
    stop("Erreur : Les chaînes n'ont pas toutes les mêmes dimensions.")
  }
  
  # Retourner un message de succès si tout est conforme
  return("Vérification réussie : Pas de NA et toutes les chaînes ont les mêmes dimensions.")
}
# Vérifier
result <- verify_cleaned_mcmc_list(mcmc_list_cleaned)
print(result)

##DETECTION DES PARAMETRES CONSTANTS
# Extraire toutes les colonnes de chaque chaîne
#all_chains_matrix <- lapply(mcmc_list_cleaned, function(chain) as.matrix(chain))

# 2. Combiner toutes les chaînes dans une seule matrice
# Pour combiner les chaînes, on va les mettre en une seule matrice de manière appropriée.
#combined_matrix <- do.call(cbind, all_chains_matrix)

# 3. Calculer la variance de chaque paramètre dans la matrice combinée
#param_variance <- apply(combined_matrix, 2, var)

# 4. Identifier les paramètres constants (variance proche de zéro)
#constant_params <- names(param_variance)[param_variance < 1e-5]  # Ajuster le seuil si nécessaire

# 5. Afficher les paramètres constants
#cat("Les paramètres constants sont :\n")
#print(constant_params)

## Filtrer les paramètres constants
#params_to_remove <- c("theta5", "theta8", "logit_theta5", "logit_theta8", "p_smolt1_migr", "h6_hw","h6_hw")

# Supprimer ces paramètres de chaque chaîne dans mcmc_list
#mcmc_list_cleaned <- lapply(mcmc_list, function(chain) {
  # Garder uniquement les colonnes dont le nom ne correspond pas aux paramètres à supprimer
  #valid_columns <- !grepl(paste(params_to_remove, collapse = "|"), colnames(chain))
  #chain[, valid_columns, drop = FALSE]
#})

#####UN AUTRE OPTION
# Étape 1 : Détecter les paramètres constants
# On va vérifier si la variance des paramètres est nulle (ou proche de zéro) sur toutes les chaînes

# Fonction pour vérifier la variance d'un paramètre dans chaque chaîne
detect_constant_params <- function(mcmc_list_cleaned) {
  constant_params <- NULL
  
  # Vérification de la variance pour chaque paramètre
  for (param_name in colnames(mcmc_list_cleaned[[1]])) {
    param_values <- do.call(cbind, lapply(mcmc_list_cleaned, function(chain) chain[, param_name]))  # Rassembler les valeurs de chaque chaîne pour ce paramètre
    param_variance <- apply(param_values, 1, var)  # Calculer la variance sur les chaînes pour chaque échantillon
    
    # Si la variance est proche de zéro, c'est un paramètre constant
    if (all(abs(param_variance) < 1e-6)) {  # Tolérance de 1e-6 pour les valeurs proches de zéro
      constant_params <- c(constant_params, param_name)
    }
  }
  
  return(constant_params)
}

# Étape 2 : Supprimer les paramètres constants
constant_params <- detect_constant_params(mcmc_list_cleaned)

# Vérification des paramètres constants trouvés
print("Paramètres constants détectés :")
print(constant_params)

# Étape 3 : Supprimer ces paramètres constants de chaque chaîne dans mcmc_list
mcmc_list_cleaned2 <- lapply(mcmc_list_cleaned, function(chain) {
  # Garder uniquement les colonnes dont le nom ne correspond pas aux paramètres constants
  valid_columns <- !colnames(chain) %in% constant_params
  chain[, valid_columns, drop = FALSE]
})

# Vérification des résultats
print("Paramètres après suppression des constants :")
print(colnames(mcmc_list_cleaned2[[1]]))  # Exemple pour la première chaîne



###calcul des criteres

ess_values <- sapply(mcmc_list_cleaned2, effectiveSize)
#ess_values<-totdf$neff
# Extract minimum effective sample size (ESS) across all parameters and chains
eff_min <- min(ess_values, na.rm = TRUE)

# Calculate time to reach ESS = 1000 for each parameter
total_time <- as.numeric(time_mcmc)  # in seconds

# Replace ESS values equal to zero by NA to avoid division by zero
#ess_values[ess_values == 0] <- NA

# Calculate time to reach ESS = 1000 for each parameter
#time_to_reach_ess_1000 <- (1000 * total_time) / ess_values

# Calculate computational efficiency (ESS per second) for each parameter/Algorithmic efficiency
Algorithmic_efficiency <- ess_values / 1800
Computational_efficiency <- (1000 * total_time) / ess_values
#
# Sort the nodes by lgorithmic_efficiency in descending order to identify the worst cases
bottleneck_nodes_by_Algorithmic_efficiency <- order(Algorithmic_efficiency, decreasing = TRUE, na.last = NA)

# Sort the nodes by computational efficiency in ascending order to identify the worst cases
bottleneck_nodes_by_Computational_efficiency <- order(Computational_efficiency, decreasing = FALSE, na.last = NA)

# Sort the nodes by ess_values in descending order to identify the worst cases
bottleneck_nodes_by_ess <- order(ess_values, decreasing = TRUE, na.last = NA)
# Extraction des noms des nœuds
all_param_names <- unique(rownames(ess_values))

# Vérification
cat("Noms des nœuds extraits :\n")
print(head(all_param_names))  # Afficher les premiers noms pour vérifier

# Filtrer les indices valides
valid_bottleneck_nodes_by_Algorithmic_efficiency <- bottleneck_nodes_by_Algorithmic_efficiency[bottleneck_nodes_by_Algorithmic_efficiency <= length(all_param_names)]
valid_bottleneck_nodes_by_Computational_efficiency <- bottleneck_nodes_by_Computational_efficiency[bottleneck_nodes_by_Computational_efficiency <= length(all_param_names)]
valid_bottleneck_nodes_by_ess <- bottleneck_nodes_by_ess[bottleneck_nodes_by_ess <= length(all_param_names)]
# Extraire les noms de nœuds avec des indices valides
Algorithmic_efficiency_bottleneck_node_names <- all_param_names[valid_bottleneck_nodes_by_Algorithmic_efficiency]
Computational_efficiency_bottleneck_node_names <- all_param_names[valid_bottleneck_nodes_by_Computational_efficiency]
Ess_bottleneck_node_names <- all_param_names[valid_bottleneck_nodes_by_ess]
#Algorithmic_efficiency_bottleneck_node_names <- all_param_names[bottleneck_nodes_by_Algorithmic_efficiency]
#Computational_efficiency_bottleneck_node_names <- all_param_names[bottleneck_nodes_by_Computational_efficiency]
#ESS_nodes <- all_param_names[ess_values]
# Extract the 100 nodes with the lowest Algorithmic Efficiency (last 100 elements)
worst_100_algorithmic_efficiency_nodes <- tail(Algorithmic_efficiency_bottleneck_node_names, 100)

# Extract the 100 nodes with the highest Computational Efficiency (first 100 elements)
worst_100_computational_efficiency_nodes <- head(Computational_efficiency_bottleneck_node_names, 100)
#ESS_nodes
# Extract the 100 nodes with the lowest ESS (last 100 elements)
ESS_100_nodes <- tail(Ess_bottleneck_node_names, 100)

# Print the worst bottlenecks based on Algorithmic Efficiency
cat("Nodes with the lowest Algorithmic Efficiency (worst bottlenecks):\n")
print(worst_100_algorithmic_efficiency_nodes)

# Print the best bottlenecks based on Computational Efficiency
cat("Nodes with the highest Computational Efficiency (best nodes):\n")
print(worst_100_computational_efficiency_nodes)

# Print the best bottlenecks based on ESS
cat("ESS (worst nodes):\n")
print(ESS_100_nodes)
# Vérifiez que les longueurs de `ess_values` et `node_names` correspondent
#if (length(ess_values) != length(node_names)) {
  #stop("Error: ESS values and node names must have the same length.")
#}
####VISUALISATION

# Fonction pour extraire les familles automatiquement en fonction des préfixes des noms de nœuds
get_family_from_name <- function(node_name) {
  # Retirer les indices entre crochets pour n'utiliser que le préfixe du nom
  node_base_name <- sub("\\[.*\\]", "", node_name)
  
  # Retourner le préfixe comme famille
  return(node_base_name)
}

# Fusionner les deux listes de nœuds
all_nodes <- c(worst_100_algorithmic_efficiency_nodes, worst_100_computational_efficiency_nodes,ESS_100_nodes)

# Appliquer la fonction pour créer une liste de familles pour tous les nœuds
families <- unique(sapply(all_nodes, get_family_from_name))

# Afficher le nombre total de familles et les familles uniques
length(families)
families
# Appliquer la fonction pour créer une liste de familles pour les nœuds bottleneck
node_families_algorithmic <- unique(sapply(worst_100_algorithmic_efficiency_nodes, get_family_from_name))
node_families_computational <- unique(sapply(worst_100_computational_efficiency_nodes, get_family_from_name))
node_families_ESS <- unique(sapply(ESS_100_nodes, get_family_from_name))

# Créer un DataFrame avec les nœuds, les familles et les métriques
# Étape 1 : Nettoyer les noms des nœuds dans totdf (si nécessaire)
totdf$clean_name <- sub("\\[.*\\]", "", totdf$fullname)  # Retirer les indices entre crochets

# Étape 2 : Aligner les valeurs d'Algorithmic_Efficiency avec les nœuds de node_families_algorithmic
# Utilisation de match() pour obtenir l'alignement entre les nœuds dans node_families_algorithmic et totdf$clean_name
Algorithmic_efficiency_aligned <- Algorithmic_efficiency[match(node_families_algorithmic, totdf$clean_name)]

# Vérification que les longueurs des vecteurs sont correctes
if (length(Algorithmic_efficiency_aligned) == length(node_families_algorithmic)) {
  # Étape 3 : Créer le DataFrame avec les familles et les valeurs d'efficacité algorithmique
  metrics_df_algorithmic <- data.frame(
    Family = node_families_algorithmic,
    Algorithmic_Efficiency = Algorithmic_efficiency_aligned
  )
  
  # Étape 4 : Éliminer les lignes contenant des NA
  metrics_df_algorithmic_clean <- metrics_df_algorithmic %>%
    filter(!is.na(Algorithmic_Efficiency))  # Filtrer les lignes où Algorithmic_Efficiency est NA
  
  # Vérification du DataFrame après nettoyage des NA
  head(metrics_df_algorithmic_clean)
}


# Étape 1 : Nettoyer les noms des nœuds dans totdf (si nécessaire)
totdf$clean_name <- sub("\\[.*\\]", "", totdf$fullname)  # Retirer les indices entre crochets

# Étape 2 : Calculer l'efficacité computationnelle
Computational_efficiency <- (1000 * total_time) / ess_values  # Calcul de Computational_Efficiency

# Étape 3 : Aligner les valeurs de Computational_Efficiency avec les nœuds de node_families_computational
Computational_efficiency_aligned <- Computational_efficiency[match(node_families_computational, totdf$clean_name)]

# Vérification que les longueurs des vecteurs sont correctes
if (length(Computational_efficiency_aligned) == length(node_families_computational)) {
  # Étape 4 : Créer le DataFrame avec les familles et les valeurs d'efficacité computationnelle
  metrics_df_computational <- data.frame(
    Family = node_families_computational,
    Computational_Efficiency = Computational_efficiency_aligned
  )
  
  # Étape 5 : Éliminer les lignes contenant des NA
  metrics_df_computational_clean <- metrics_df_computational %>%
    filter(!is.na(Computational_Efficiency))  # Filtrer les lignes où Computational_Efficiency est NA
  
  # Vérification du DataFrame après nettoyage des NA
  head(metrics_df_computational_clean)
} else {
  print("Les tailles des données ne correspondent pas. Vérifiez les correspondances des nœuds.")
}

# Étape 1 : Nettoyer les noms des nœuds dans totdf (si nécessaire)
totdf$clean_name <- sub("\\[.*\\]", "", totdf$fullname)  # Retirer les indices entre crochets

# Étape 2 : Extraire ess_values à partir de totdf$neff
ess_values <- totdf$neff  # Assurez-vous que totdf$neff contient les valeurs correctes

# Étape 3 : Aligner les valeurs d'ess_values avec les nœuds de node_families_algorithmic
ess_values_aligned <- ess_values[match(node_families_algorithmic, totdf$clean_name)]

# Vérification que les longueurs des vecteurs sont correctes
if (length(ess_values_aligned) == length(node_families_algorithmic)) {
  # Étape 4 : Créer le DataFrame avec les familles et les valeurs d'ESS
  metrics_df_ESS <- data.frame(
    Family = node_families_algorithmic,
    ESS_values = ess_values_aligned
  )
  
  # Étape 5 : Éliminer les lignes contenant des NA
  metrics_df_ESS_clean <- metrics_df_ESS %>%
    filter(!is.na(ESS_values))  # Filtrer les lignes où ESS_values est NA
  
  # Vérification du DataFrame après nettoyage des NA
  head(metrics_df_ESS_clean)
} else {
  print("Les tailles des données ne correspondent pas. Vérifiez les correspondances des nœuds.")
}





library(dplyr)
library(ggplot2)

# Calcul de la médiane pour chaque famille dans chaque DataFrame
metrics_by_family_algorithmic <- metrics_df_algorithmic_clean %>%
  group_by(Family) %>%
  summarise(
    Median_Algorithmic_Efficiency = median(Algorithmic_Efficiency, na.rm = TRUE)
  )

metrics_by_family_computational <- metrics_df_computational_clean %>%
  group_by(Family) %>%
  summarise(
    Median_Computational_Efficiency = median(Computational_Efficiency, na.rm = TRUE)
  )

metrics_by_family_ESS <- metrics_df_ESS_clean %>%
  group_by(Family) %>%
  summarise(
    Median_ESS = median(ESS_values, na.rm = TRUE)
  )

# Fonction pour tracer les graphiques par famille (sans afficher les messages d'erreur)
plot_bottlenecks_by_family <- function(metrics_by_family, metric_name, title, fill_color) {
  # Sélection de la colonne de la médiane en fonction de la métrique
  if (metric_name == "Algorithmic_Efficiency") {
    median_column <- "Median_Algorithmic_Efficiency"
  } else if (metric_name == "Computational_Efficiency") {
    median_column <- "Median_Computational_Efficiency"
  } else if (metric_name == "ESS") {
    median_column <- "Median_ESS"
  }
  
  # Vérification si la colonne existe avant de tracer, sans message d'erreur
  if (median_column %in% colnames(metrics_by_family)) {
    ggplot(metrics_by_family, aes(x = reorder(Family, .data[[median_column]]), 
                                  y = .data[[median_column]])) +
      geom_bar(stat = "identity", fill = fill_color) +
      coord_flip() +
      labs(
        title = title,
        x = "Node families",
        y = "Median"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
      )
  }
}

# Liste des métriques
metrics_list <- list(
  Algorithmic_Efficiency = "Algorithmic_Efficiency",
  Computational_Efficiency = "Computational_Efficiency",
  ESS = "ESS"
)

# Titres pour les graphiques
titles <- list(
  Algorithmic_Efficiency = "Algorithmic_Efficiency",
  Computational_Efficiency = "Computational_Efficiency",
  ESS = "ESS median by Family"
)

# Couleurs pour les graphiques
colors <- list(
  Algorithmic_Efficiency = "darkorange",
  Computational_Efficiency = "steelblue",
  ESS = "firebrick"
)

# Générer les graphiques pour les trois DataFrames (algorithmic, computational et ESS)
metrics_dfs <- list(metrics_by_family_algorithmic, metrics_by_family_computational, metrics_by_family_ESS)

# Tracer les graphiques pour chaque métrique (médiane)
for (metrics_df in metrics_dfs) {
  for (metric in names(metrics_list)) {
    print(
      plot_bottlenecks_by_family(
        metrics_by_family = metrics_df,
        metric_name = metrics_list[[metric]],
        title = titles[[metric]],
        fill_color = colors[[metric]]
      )
    )
  }
}










