suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

# Argumentos
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Ruta al archivo TSV de la red (aristas: from, to, weight)"),
  make_option(c("-m", "--modo"), type = "character", default = "completo",
              help = "Modo de ejecución: 'rapido' o 'completo'"),
  make_option(c("-n", "--n_perm"), type = "integer", default = 1000,
              help = "Número de permutaciones para modelos nulos [defecto: 1000]"),
  make_option(c("-s", "--semilla"), type = "integer", default = 123,
              help = "Semilla aleatoria para reproducibilidad [defecto: 123]"),
  make_option(c("-o", "--output"), type = "character", default = "../results/topology",
              help = "Carpeta de salida"),
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "Archivo con lista de genes semilla (un gen por línea)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Error: Debe especificar un archivo de entrada con -i o --input")
}

set.seed(opt$semilla) # Establece la semilla para sample_gnm y rewire
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE) # Directorio de salida

# 2. Cargar red e instanciar grafo igraph
edges_df <- read.csv2(opt$input,stringsAsFactors = FALSE)
g <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)

# 3. Métricas Globales (Ambos modos)
# El modo rapido calcula sólo las métricas globales
# el modo ocmpleto calcula añade métricas por nodo, distribucion, modelos nulos y robustez
comp_info <- igraph::components(g)
max_comp_size <- max(comp_info$csize)

metricas_globales <- data.frame(
  nodos = igraph::gorder(g),
  aristas = igraph::gsize(g),
  densidad = igraph::edge_density(g),
  grado_medio = mean(igraph::degree(g)),
  num_componentes = comp_info$no,
  tam_componente_mayor = max_comp_size,
  prop_componente_mayor = max_comp_size / igraph::gorder(g),
  diametro = igraph::diameter(g, directed = FALSE),
  distancia_media = igraph::mean_distance(g, directed = FALSE),
  transitividad_global = igraph::transitivity(g, type = "global"),
  transitividad_media = igraph::transitivity(g, type = "average"),
  asortatividad_grado = igraph::assortativity_degree(g, directed = FALSE)
)

write.table(metricas_globales, file = file.path(opt$output, "metricas_globales.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Métricas globales guardadas exitosamente.\n")

if (opt$modo == "rapido") {
  cat("Modo rápido finalizado sin ejecutar modelos nulos ni robustez.\n")
  q("no")
}

# 4. Métricas por Nodo (solo Modo Completo)
# Closeness se calcula solo sobre la componente mayor (LCC)
main_comp_vids <- V(g)[comp_info$membership == which.max(comp_info$csize)]
g_main <- igraph::induced_subgraph(g, main_comp_vids)

closeness_main <- igraph::closeness(g_main)
closeness_all <- setNames(rep(NA, igraph::gorder(g)), V(g)$name)
closeness_all[names(closeness_main)] <- closeness_main

metricas_nodos <- data.frame(
  gene = V(g)$name,
  grado = igraph::degree(g),
  betweenness = igraph::betweenness(g, directed = FALSE),
  harmonic_centrality = igraph::harmonic_centrality(g),
  closeness_lcc = closeness_all
)

write.table(metricas_nodos, file = file.path(opt$output, "metricas_nodos.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5. Distribución de Grado
degree_vals <- igraph::degree(g)
fit_power <- igraph::fit_power_law(degree_vals)

pdf(file.path(opt$output, "distribucion_grado.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2))
hist(degree_vals, main = "Distribución de Grado", xlab = "Grado", col = "lightblue", border = "white")
deg_freq <- table(degree_vals)
plot(as.numeric(names(deg_freq)), as.numeric(deg_freq), log = "xy",
     main = paste("Log-Log (alpha =", round(fit_power$alpha, 2), ")"),
     xlab = "Grado (log)", ylab = "Frecuencia (log)", pch = 19, col = "darkblue")
dev.off()

# 6. Modelos Nulos: de dos formas
# Con sample_gnm
n_v <- igraph::gorder(g)
n_e <- igraph::gsize(g)
obs_trans <- igraph::transitivity(g, type = "global")
obs_dist <- igraph::mean_distance(g, directed = FALSE)
obs_mod <- igraph::modularity(igraph::cluster_louvain(g))

null_trans <- numeric(opt$n_perm)
null_dist <- numeric(opt$n_perm)
null_mod <- numeric(opt$n_perm)

for (i in 1:opt$n_perm) {
  g_null <- igraph::sample_gnm(n_v, n_e, directed = FALSE)
  null_trans[i] <- igraph::transitivity(g_null, type = "global")
  null_dist[i] <- igraph::mean_distance(g_null, directed = FALSE)
  null_mod[i] <- igraph::modularity(igraph::cluster_louvain(g_null))
}

p_trans <- mean(null_trans >= obs_trans)
p_dist <- mean(null_dist <= obs_dist)
p_mod <- mean(null_mod >= obs_mod)

pdf(file.path(opt$output, "modelos_nulos_histogramas.pdf"), width = 9, height = 3)
par(mfrow = c(1, 3))
hist(null_trans, main = paste("Transitividad (p =", round(p_trans, 4), ")"), xlab = "Transitividad Nula", col = "gray90")
abline(v = obs_trans, col = "red", lwd = 2)

hist(null_dist, main = paste("Distancia Media (p =", round(p_dist, 4), ")"), xlab = "Distancia Nula", col = "gray90")
abline(v = obs_dist, col = "red", lwd = 2)

hist(null_mod, main = paste("Modularidad (p =", round(p_mod, 4), ")"), xlab = "Modularidad Nula", col = "gray90")
abline(v = obs_mod, col = "red", lwd = 2)
dev.off()

# Y preservando secuencia de grados (rewire)

null_trans_rewire <- numeric(opt$n_perm)
null_dist_rewire  <- numeric(opt$n_perm)
null_mod_rewire   <- numeric(opt$n_perm)

for (i in seq_len(opt$n_perm)) {
  g_rewired <- igraph::rewire(g, with = igraph::keeping_degseq( niter = 10 * igraph::gsize(g)))

  null_trans_rewire[i] <- igraph::transitivity(g_rewired, type = "global")
  null_dist_rewire[i] <-igraph::mean_distance(g_rewired, directed = FALSE)
  null_mod_rewire[i] <-igraph::modularity(igraph::cluster_louvain(g_rewired))
}

p_trans_rewire <- mean(null_trans_rewire >= obs_trans)
p_dist_rewire <- mean(null_dist_rewire <= obs_dist)
p_mod_rewire <- mean(null_mod_rewire >= obs_mod)

resultados_nulos <- data.frame(
  metrica = c("transitividad", "distancia_media", "modularidad"),
  observado = c(obs_trans, obs_dist, obs_mod),
  p_gnm = c(p_trans, p_dist, p_mod),
  p_rewire = c(p_trans_rewire, p_dist_rewire, p_mod_rewire)
)

write.table(
  resultados_nulos, file = file.path(opt$output, "resultados_modelos_nulos.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 7. Prueba del módulo de enfermedad (opcional, dado por el argumento -g)
#Comprueba si los genes del fenotipo están más conectados entre sí en STRING que conjuntos aleatorios del mismo tamaño

seed_genes <- NULL

if (!is.null(opt$genes)) {
  seed_genes <- readLines(opt$genes)
}

if (!is.null(seed_genes)) {
  seed_genes_in_net <- intersect(seed_genes, V(g)$name)

  if (length(seed_genes_in_net) > 1) {
    g_seed_sub <-igraph::induced_subgraph(g, seed_genes_in_net)
    obs_seed_lcc <- max(igraph::components(g_seed_sub)$csize)
    background_genes <- V(g)$name
    n_seeds <- length(seed_genes_in_net)

    null_seed_lcc <- replicate(opt$n_perm, {
        rand_nodes <- sample(background_genes, n_seeds)
        g_rand_sub <- igraph::induced_subgraph(g, rand_nodes)
        max(igraph::components(g_rand_sub)$csize)
      }
    )

    p_val_disease_module <- mean(null_seed_lcc >= obs_seed_lcc)

    res_disease_module <-
      data.frame(genes_en_red = n_seeds, obs_lcc_size = obs_seed_lcc,
        mean_null_lcc = mean(null_seed_lcc), p_value = p_val_disease_module)

    write.table(res_disease_module,
      file = file.path(opt$output, "disease_module_test.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

    pdf(file.path(opt$output, "disease_module_test.pdf"))
    hist(null_seed_lcc,
      main = paste( "Disease Module Test (p =", round(p_val_disease_module, 4), ")"),
      xlab = "LCC aleatoria", col = "gray90")

    abline(v = obs_seed_lcc, col = "red", lwd = 2)
    dev.off()
  }
}

# 8. Análisis de Robustez (Ataque a Hubs vs. Azar)
simular_robustez <- function(graph, modo = "hub") {
  g_temp <- graph
  sizes <- numeric(igraph::gorder(graph))
  sizes[1] <- max(igraph::components(g_temp)$csize)
  
  for (step in 2:igraph::gorder(graph)) {
    if (modo == "hub") {
      node_to_remove <- names(which.max(igraph::degree(g_temp)))
    } else {
      node_to_remove <- sample(V(g_temp)$name, 1)
    }
    g_temp <- igraph::delete_vertices(g_temp, node_to_remove)
    sizes[step] <- if (igraph::gorder(g_temp) > 0) max(igraph::components(g_temp)$csize) else 0
  }
  return(sizes / igraph::gorder(graph))
}

# Función para simular robustez y promediar resultados
simular_robustez_random_promedio <- function(graph, n_sim = 100) {
  n_nodes <- igraph::gorder(graph)

  mat_results <-matrix(0, nrow = n_sim, ncol = n_nodes)

  for (s in seq_len(n_sim)) {
    g_temp <- graph
    order_random <- sample(V(graph)$name)
    mat_results[s, 1] <- max(igraph::components(g_temp)$csize)

    for (step in 2:n_nodes) {
      g_temp <- igraph::delete_vertices(g_temp, order_random[step - 1])
      mat_results[s, step] <- if (igraph::gorder(g_temp) > 0) {
          max(igraph::components(g_temp)$csize)
        } else {0}
    }
  }

  colMeans(mat_results) / n_nodes
}

# Ejecución de simulaciones de robustez
frac_hubs <- simular_robustez(g, modo = "hub")
frac_random <- simular_robustez_random_promedio(g, n_sim = 100)

pdf(file.path(opt$output, "robustez_red.pdf"), width = 6, height = 5)
plot(seq(0, 1, length.out = length(frac_hubs)), frac_hubs, type = "l", col = "red", lwd = 2,
     xlab = "Proporción de nodos retirados", ylab = "Fracción Componente Mayor",
     main = "Análisis de Robustez")
lines(seq(0, 1, length.out = length(frac_random)), frac_random, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Ataque a Hubs", "Remoción Aleatoria"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
dev.off()

cat("Script 03_topologia.R completado con éxito.\n")