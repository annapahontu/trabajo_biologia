#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
})

# 1. Argumentos 
option_list<- list(
    make_option(c("-d","--data"), type = "character", default = NULL, help = "Carpeta con los datos del paso 01: results/data"),
    make_option(c("-n","--network"), type = "character", default = NULL, help = "Carpeta con las redes del paso 02: results/network"),
    make_option(c("-t", "--topology"), type = "character", default = NULL, help = "Carpeta con las métricas del paso 03: results/topology"),
    make_option(c("-u", "--umbrales"), type = "character", default = "400,700,900", help = "Umbrales de STRING separados por comas [defecto: 400,700,900]"),
    make_option(c("-o","--output"), type = "character", default = NULL, help = "Carpeta de salida: results/figuras")
)

opt<- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data) || is.null(opt$network) || is.null(opt$topology) || is.null(opt$output)) {
    stop("Uso: Rscript 08_figuras_resumen.R -d results/data -n results/network -t results/topology -o results/figuras")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# Pasamos "400,700,900" a un vector de números: 400 700 900
umbrales<- as.integer(strsplit(opt$umbrales, ",")[[1]])
cat("Umbrales:", umbrales, "\n")

# 2. Comprobar que existen los archivos de los pasos de antes
# Las métricas por nodo del 03 están en una carpeta por umbral
archivos<- c(
  hpo_genes    = file.path(opt$data, "hpo_genes.tsv"),
  sensibilidad = file.path(opt$network, "sensibilidad_umbrales.csv"),
  file.path(opt$topology, paste0("umbral_", umbrales), "metricas_nodos.tsv")
)

faltan<- archivos[!file.exists(archivos)]
if (length(faltan) > 0) {
  stop("Faltan archivos de pasos anteriores (ejecuta antes launch.sh):\n  ",
       paste(faltan, collapse = "\n  "))
}

# 3. Leer las tablas
hpo_genes<- read.delim(archivos["hpo_genes"], stringsAsFactors = FALSE)
sensibilidad<- read.csv2(archivos["sensibilidad"], stringsAsFactors = FALSE)

# Una tabla de métricas por umbral, guardadas en una lista con el umbral como nombre
metricas<- list()
for (u in umbrales) {
  metricas[[as.character(u)]]<- read.delim(
    file.path(opt$topology, paste0("umbral_", u), "metricas_nodos.tsv"),
    stringsAsFactors = FALSE)
}

cat("Genes de HPO:", nrow(hpo_genes), "\n")
print(sensibilidad[, c("umbral", "genes_mapeados", "nodos_red")])

# 4. Figura 1: embudo de genes
# Cuántos genes del fenotipo quedan en cada paso: los de HPO, los que se mapean en STRING y los que tienen alguna interacción con cada umbral
embudo<- c(nrow(hpo_genes), sensibilidad$genes_mapeados[1], sensibilidad$nodos_red)
names(embudo)<- c("Genes de HPO", "Mapeados en STRING", paste("En la red (umbral", sensibilidad$umbral, ")"))
print(embudo)

write.table(data.frame(paso = names(embudo), genes = embudo), file.path(opt$output, "embudo_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

pdf(file.path(opt$output, "embudo_genes.pdf"), width = 8, height = 4)
par(mar = c(4, 14, 3, 2)) # margen izquierdo amplio para los nombres de los pasos
posiciones<- barplot(rev(embudo), horiz = TRUE, las = 1, col = "lightblue", xlim = c(0, max(embudo) * 1.15), xlab = "Número de genes", main = "Genes del fenotipo en cada paso del análisis")
text(rev(embudo), posiciones, labels = rev(embudo), pos = 4)
invisible(dev.off())