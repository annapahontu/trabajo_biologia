#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(STRINGdb)
    library(igraph)
})

# 1. Argumentos
option_list <- list(
    make_option(c("-u", "--umbral"), type = "integer", default = 400,
                help = "Umbral mínimo de STRING. Las redes de umbrales mayores se sacan de esta [defecto: 400]"),
    make_option(c("--string_dir"), type = "character", default = NULL,
                help = "Carpeta con los archivos de STRING del paso 02: results/network/string_data"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Archivo de salida: results/network/string_completo.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$string_dir) || is.null(opt$output)) {
    stop("Uso: Rscript 02b_string_completo.R --string_dir string_data -o string_completo.rds")
}
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# 2. Red completa de STRING
# La cargamos una sola vez y la guardamos, para que los pasos que la necesitan
# (distancias y genes puente) no tengan que volver a leer STRING cada vez
cat("Cargando la red completa de STRING, puede tardar unos minutos\n")
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = opt$umbral, input_directory = opt$string_dir)
g <- string_db$get_graph()

# cada interacción viene repetida; al simplificar nos quedamos con la puntuación
# más alta para no perder combined_score, que usamos luego para filtrar por umbral
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
              edge.attr.comb = list(combined_score = "max", "ignore"))

# guardamos el símbolo de cada gen para no depender de STRINGdb en los pasos siguientes
proteinas <- string_db$get_proteins()
V(g)$gen <- proteinas$preferred_name[match(V(g)$name, proteinas$protein_external_id)]

cat("Red de STRING (umbral", opt$umbral, "):", vcount(g), "proteínas y", ecount(g), "interacciones\n")

saveRDS(g, opt$output)
cat("Red guardada en:", opt$output, "\n")
