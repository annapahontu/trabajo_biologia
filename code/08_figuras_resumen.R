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
    make_option(c("-o","--output"), type = "character", default = NULL, help = "Carpeta de salida: results/figuras"),
    make_option(c("-k", "--n_etiquetas"), type = "integer", default = 5, help = "Genes con mayor betweenness que se etiquetan en la figura 3 [defecto: 5]")
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

# 5. Figura 2: sensibilidad al umbral de STRING
# Cómo cambia la red al exigir más confianza a las interacciones
pdf(file.path(opt$output, "sensibilidad_umbral.pdf"), width = 10, height = 3.5)
par(mfrow = c(1, 3)) # tres gráficos en una fila

columnas<- c(nodos_red = "Genes en la red", aristas_red = "Interacciones", componentes = "Componentes")
for (col in names(columnas)) {
    plot(sensibilidad$umbral, sensibilidad[[col]], type = "b", pch = 19, col = "darkblue", xaxt = "n", xlab = "Umbral de STRING", ylab = columnas[col], main = columnas[col], xlim = c(300, 1000), ylim = c(0, max(sensibilidad[[col]]) * 1.15))
    axis(1, at = sensibilidad$umbral)
    text(sensibilidad$umbral, sensibilidad[[col]], sensibilidad[[col]], pos = 3)
}
invisible(dev.off())

# 6. Figura 3: grado frente a betweenness en cada umbral
# Los genes puente tienen betweenness alta: muchos caminos más cortos pasan por ellos
pdf(file.path(opt$output, "grado_betweenness.pdf"), width = 5 * length(umbrales), height = 5)
par(mfrow = c(1, length(umbrales)))

puentes<- data.frame()
for (u in umbrales) {
    m<- metricas[[as.character(u)]]
    plot(m$grado, m$betweenness, pch = 19, col = "lightblue", xlab = "Grado", ylab = "Betweenness", main = paste0("Umbral ", u, " (", nrow(m), " genes)"), xlim = c(0, max(m$grado) * 1.15), ylim = c(0, max(m$betweenness) * 1.2)) # xlim e ylim son para dejar margenes

# Etiquetamos los genes con mayor betweenness (solo si es mayor que 0)
    top<- head(m[order(m$betweenness, decreasing = TRUE), ], opt$n_etiquetas)
    top<- top[top$betweenness > 0, ]
    if (nrow(top) > 0) {
        points(top$grado, top$betweenness, pch = 19, col = "darkred")
# Si varios genes caen en el mismo punto, juntamos sus nombres en una sola etiqueta
        etiquetas<- aggregate(gene ~ grado + betweenness, data = top, FUN = function(x) paste(x, collapse = "\n"))
        text(etiquetas$grado, etiquetas$betweenness, etiquetas$gene, pos = rep(c(3, 1), length.out = nrow(etiquetas)), cex = 0.7)
    }
    puentes<- rbind(puentes, data.frame(umbral = u, top[, c("gene", "grado", "betweenness")]))
}
invisible(dev.off())

cat("\nGenes con mayor betweenness en cada umbral:\n")
print(puentes, row.names = FALSE)
write.table(puentes, file.path(opt$output, "top_betweenness.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nFiguras guardadas en:", opt$output, "\n")