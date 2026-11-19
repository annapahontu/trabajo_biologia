#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(STRINGdb)
    library(igraph)
})

# 1. Argumentos
option_list <- list(
    make_option(c("-p", "--particion"), type = "character", default = NULL,
                help = "Partición final del paso 04: results/clusters/umbral_700/particion_final.tsv"),
    make_option(c("-g", "--genes"), type = "character", default = NULL,
                help = "Genes semilla del paso 01: results/data/hpo_genes.tsv"),
    make_option(c("-u", "--umbral"), type = "integer", default = 700,
                help = "Umbral de STRING [defecto: 700]"),
    make_option(c("--string_dir"), type = "character", default = NULL,
                help = "Carpeta con los archivos de STRING del paso 02: results/network/string_data"),
    make_option(c("-d", "--max_distancia"), type = "integer", default = 3,
                help = "Distancia máxima para buscar genes intermedios [defecto: 3]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Carpeta de salida: results/clusters/umbral_700/distancias")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$particion) || is.null(opt$genes) || is.null(opt$string_dir) || is.null(opt$output)) {
    stop("Uso: Rscript 07_distancias_clusters.R -p particion_final.tsv -g hpo_genes.tsv --string_dir string_data -o salida")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# 2. Partición del paso 04 y genes semilla
particion <- read.delim(opt$particion, colClasses = "character")
genes <- read.delim(opt$genes, colClasses = "character", quote = "")

# 3. Red completa de STRING con el umbral elegido
# Usamos la carpeta del paso 02 para no volver a descargar nada
cat("Cargando la red completa de STRING, puede tardar unos minutos\n")
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = opt$umbral, input_directory = opt$string_dir)
mapeo <- string_db$map(genes, "gene_symbol", removeUnmappedRows = TRUE)
cat("\n")
g <- simplify(string_db$get_graph())
proteinas <- string_db$get_proteins()
cat("Red de STRING (umbral", opt$umbral, "):", vcount(g), "proteínas y", ecount(g), "interacciones\n")

# Las semillas sin ninguna interacción en STRING no tienen distancia a nada
mapeo <- mapeo[mapeo$STRING_id %in% V(g)$name, ]

# 4. Grupos: cada cluster es un grupo, y cada semilla que no entró en ningún
# cluster (gen aislado en la red semilla) es un grupo de un solo gen
mapeo$grupo <- mapeo$gene_symbol
en_cluster <- mapeo$gene_symbol %in% particion$gene
mapeo$grupo[en_cluster] <- paste0("cluster_",
                                  particion$cluster[match(mapeo$gene_symbol[en_cluster], particion$gene)])
grupos <- unique(mapeo$grupo)
n <- length(grupos)
cat("Grupos:", sum(startsWith(grupos, "cluster_")), "clusters y",
    sum(!startsWith(grupos, "cluster_")), "genes aislados\n")

# 5. Distancias
# Primero la distancia entre todas las semillas, y para cada pareja de grupos
# nos quedamos con la menor entre un gen de uno y un gen del otro
d_genes <- distances(g, v = mapeo$STRING_id, to = mapeo$STRING_id)

dist_grupos <- matrix(0, n, n, dimnames = list(grupos, grupos))
for (i in 1:n) {
    for (j in 1:n) {
        dist_grupos[i, j] <- min(d_genes[mapeo$grupo == grupos[i], mapeo$grupo == grupos[j]])
    }
}
write.table(data.frame(grupo = grupos, dist_grupos, check.names = FALSE),
            file.path(opt$output, "distancias_grupos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 6. Camino más corto entre cada pareja de grupos cercanos
# Los genes que quedan dentro del camino son los que faltan para unir los dos grupos
cat("\nBuscando los caminos más cortos entre", n * (n - 1) / 2, "parejas de grupos\n")
barra <- txtProgressBar(min = 0, max = n * (n - 1) / 2, style = 3)
pareja <- 0
caminos <- data.frame()
for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
        pareja <- pareja + 1
        setTxtProgressBar(barra, pareja)
        d <- dist_grupos[i, j]
        if (d < 2 || d > opt$max_distancia) {
            next
        }
        # pareja de genes, uno de cada grupo, que está a la distancia mínima
        genes_a <- mapeo$gene_symbol[mapeo$grupo == grupos[i]]
        genes_b <- mapeo$gene_symbol[mapeo$grupo == grupos[j]]
        sub <- d_genes[mapeo$grupo == grupos[i], mapeo$grupo == grupos[j], drop = FALSE]
        pos <- which(sub == d, arr.ind = TRUE)[1, ]

        camino <- shortest_paths(g, from = rownames(sub)[pos[1]], to = colnames(sub)[pos[2]])$vpath[[1]]
        intermedios <- names(camino)[2:d]
        simbolos <- proteinas$preferred_name[match(intermedios, proteinas$protein_external_id)]

        caminos <- rbind(caminos, data.frame(
            grupo_a = grupos[i],
            grupo_b = grupos[j],
            gen_a = genes_a[pos[1]],
            gen_b = genes_b[pos[2]],
            distancia = d,
            genes_intermedios = paste(simbolos, collapse = ",")
        ))
    }
}
close(barra)
caminos <- caminos[order(caminos$distancia), ]
write.table(caminos, file.path(opt$output, "caminos_cortos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nParejas de grupos según la distancia:\n")
print(table(caminos$distancia))

# 7. Heatmap de las distancias (las parejas sin camino quedan en blanco)
d_figura <- dist_grupos
d_figura[is.infinite(d_figura)] <- NA
pdf(file.path(opt$output, "distancias_grupos.pdf"), width = 10, height = 9)
par(mar = c(9, 9, 3, 2))
image(1:n, 1:n, d_figura, col = hcl.colors(10, "YlOrRd", rev = TRUE), axes = FALSE, xlab = "", ylab = "",
      main = paste("Distancia en STRING entre grupos (umbral", opt$umbral, ")"))
axis(1, at = 1:n, labels = grupos, las = 2, cex.axis = 0.6)
axis(2, at = 1:n, labels = grupos, las = 2, cex.axis = 0.6)
# el número de la distancia dentro de cada casilla
text(rep(1:n, times = n), rep(1:n, each = n), d_figura, cex = 0.5)
invisible(dev.off())

cat("\nResultados guardados en:", opt$output, "\n")