#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
})

# 1. Argumentos 
option_list<- list(
    make_option(c("-d","--data"), type = "character", default = NULL, help = "Carpeta con los datos del paso 01: results/data"),
    make_option(c("-n","--network"), type = "character", default = NULL, help = "Carpeta con las redes del paso 02: results/network"),
    make_option(c("-t", "--topology"), type = "character", default = NULL, help = "Carpeta con las métricas del paso 03: results/topology"),
    make_option(c("-u", "--umbral"), type = "integer", default = 700, help = "Umbral de STRING de la red final [defecto: 700]"),
    make_option(c("-o","--output"), type = "character", default = NULL, help = "Carpeta de salida: results/figuras")
)

opt<- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data) || is.null(opt$network) || is.null(opt$topology) || is.null(opt$output)) {
    stop("Uso: Rscript 08_figuras_resumen.R -d results/data -n results/network -t results/topology -o results/figuras")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

carpeta_umbral<- file.path(opt$network, paste0("umbral_", opt$umbral)) # carpeta de la red con el umbral elegido
cat("Leyendo la red de", carpeta_umbral, "\n")

# 2. Comprobar que existen los archivos de los pasos de antes
archivos <- c(
  hpo_genes= file.path(opt$data, "hpo_genes.tsv"),
  mapeo= file.path(carpeta_umbral, "mapeo_genes.csv"),
  nodos= file.path(carpeta_umbral, "nodos.csv"),
  nodos_exp= file.path(carpeta_umbral, "nodos_expandida.csv"),
  sensibilidad= file.path(opt$network, "sensibilidad_umbrales.csv"),
  metricas= file.path(opt$topology, "metricas_nodos.tsv"),
  metricas_exp= file.path(opt$topology, "expandida", "metricas_nodos.tsv")
)

faltan <- archivos[!file.exists(archivos)]
if (length(faltan) > 0) {
  stop("Faltan archivos de pasos anteriores (ejecuta antes launch.sh):\n  ",
       paste(faltan, collapse = "\n  "))
}

# 3. Leer las tablas 
hpo_genes<- read.delim(archivos["hpo_genes"], stringsAsFactors = FALSE)
mapeo<- read.csv2(archivos["mapeo"], stringsAsFactors = FALSE)
nodos<- read.csv2(archivos["nodos"], stringsAsFactors = FALSE)
nodos_exp<- read.csv2(archivos["nodos_exp"], stringsAsFactors = FALSE)
sensibilidad<- read.csv2(archivos["sensibilidad"], stringsAsFactors = FALSE)
metricas<- read.delim(archivos["metricas"], stringsAsFactors = FALSE)
metricas_exp<- read.delim(archivos["metricas_exp"], stringsAsFactors = FALSE)

cat("Genes HPO:", nrow(hpo_genes), "\n")
cat("Mapeados en STRING:", nrow(mapeo), "\n")
cat("En la red de umbral", opt$umbral, ":", nrow(nodos), "\n")
cat("Semillas en la red expandida:", sum(nodos_exp$tipo == "semilla"), "\n")