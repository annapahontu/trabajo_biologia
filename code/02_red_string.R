#!/usr/bin/env Rscript

library(optparse)
library(STRINGdb)
library(igraph)

#1. Argumentos
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Entrada del proyecto: results/data/hpo_genes.tsv (generada por el paso 01)"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Salida del proyecto: results/network (redes, tablas y descargas de STRING)"),
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
    stop(" Uso: Rscript 02_red_string.R -i results/data/hpo_genes.tsv -o results/network")
}
if (!file.exists(opt$input)) {
    stop("No se encuentra el archivo de entrada: ", opt$input, ".El paso 01 debe generar results/data/hpo_genes.tsv antes de ejecutar el paso 02.")
}


