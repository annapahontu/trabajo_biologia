#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(igraph)
})

# 1. Argumentos
option_list <- list(
    make_option(c("-c", "--caminos"), type = "character", default = NULL,
                help = "Caminos más cortos del paso 07: results/clusters/umbral_[400|600]/distancias/caminos_cortos.tsv"),
    make_option(c("-g", "--genes"), type = "character", default = NULL,
                help = "Genes semilla del paso 01: results/data/hpo_genes.tsv"),
    make_option(c("-r", "--red"), type = "character", default = NULL,
                help = "Red completa de STRING del paso 02b: results/network/string_completo.rds"),
    make_option(c("-u", "--umbral"), type = "integer", default = 600,
                help = "Umbral de STRING [defecto: 600]"),
    make_option(c("-n", "--n_puentes"), type = "integer", default = 5,
                help = "Número de genes puente que se añaden [defecto: 5]"),
    make_option(c("--max_grado"), type = "integer", default = 300,
                help = "Grado máximo de un gen puente, para descartar hubs [defecto: 300]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Carpeta de salida: results/network/puentes_[400|600]/n_5")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$caminos) || is.null(opt$genes) || is.null(opt$red) || is.null(opt$output)) {
    stop("Uso: Rscript 09_red_puentes.R -c caminos_cortos.tsv -g hpo_genes.tsv -r string_completo.rds -o salida")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# 2. Datos de los pasos anteriores
caminos <- read.delim(opt$caminos, colClasses = "character")
genes <- read.delim(opt$genes, colClasses = "character", quote = "")

g <- readRDS(opt$red)
g <- delete_edges(g, E(g)[combined_score < opt$umbral])
g <- delete_vertices(g, which(degree(g) == 0))

# nos quedamos con las parejas a distancia 2, que tienen un solo gen en medio
caminos <- caminos[caminos$distancia == "2", ]
if (nrow(caminos) == 0) {
    stop("No hay parejas de grupos a distancia 2 en ", opt$caminos)
}

# 3. Genes puente
# miramos cuántos grupos distintos une cada gen intermedio
puentes <- data.frame(gen = unique(caminos$genes_intermedios))
puentes$n_grupos <- sapply(puentes$gen, function(gen) {
    filas <- caminos[caminos$genes_intermedios == gen, ]
    length(unique(c(filas$grupo_a, filas$grupo_b)))
})

# fuera los hubs, que están conectados a casi todo y unen grupos sin que signifique nada
puentes$grado <- degree(g, v = V(g)$name[match(puentes$gen, V(g)$gen)])
puentes <- puentes[puentes$grado <= opt$max_grado, ]

# primero los que unen más grupos y, si empatan, los de menos grado
puentes <- head(puentes[order(-puentes$n_grupos, puentes$grado), ], opt$n_puentes)
cat("Genes puente añadidos:\n")
print(puentes, row.names = FALSE)
write.table(puentes, file.path(opt$output, "genes_puente.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Red de semillas más genes puente
semillas <- V(g)$name[V(g)$gen %in% genes$gene_symbol]
red <- induced_subgraph(g, c(semillas, V(g)$name[match(puentes$gen, V(g)$gen)]))
red <- delete_vertices(red, which(degree(red) == 0))
V(red)$name <- V(red)$gen
cat("Red con puentes:", vcount(red), "nodos y", ecount(red), "aristas\n")

# la guardamos igual que el aristas.csv del paso 02 para poder pasarle el 04 y el 05
write.csv2(as_data_frame(red, what = "edges")[, c("from", "to", "combined_score")],
           file.path(opt$output, "aristas.csv"), row.names = FALSE)

# 5. Figura con las semillas en naranja y los puentes en azul
es_puente <- V(red)$name %in% puentes$gen
set.seed(123)
pdf(file.path(opt$output, "red_puentes.pdf"), width = 8, height = 8)
plot(red, layout = layout_with_fr(red),
     vertex.color = ifelse(es_puente, "steelblue", "orange"),
     vertex.size = 5, vertex.label = ifelse(es_puente, V(red)$name, NA),
     vertex.label.cex = 0.7, edge.width = 0.5,
     main = paste("Semillas y genes puente (umbral", opt$umbral, ")"))
invisible(dev.off())