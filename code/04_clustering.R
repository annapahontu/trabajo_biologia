#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(igraph)
    library(optparse)
})

# 1. Argumentos 
option_list<- list(
    make_option(c("-i","--input"), type = "character", default = NULL, help = "Aristas de la red: results/network/umbral_700/aristas.csv (paso 02)"),
    make_option(c("-o","--output"), type = "character", default = NULL, help = "Carpeta de salida: results/clusters"),
    make_option(c("-s", "--semilla"), type = "integer", default = 123, help = "Semilla aleatoria [defecto: 123]"),
    make_option(c("-n", "--n_semillas"), type = "integer", default = 100, help = "Número de semillas para la estabilidad de Louvain [defecto: 100]"),
    make_option(c("-m", "--metodo"), type = "character", default = "louvain", help = "Método para la partición final [defecto: louvain]"),
    make_option(c("-r", "--resolucion"), type = "double", default = 1, help = "Resolución de Louvain y Leiden [defecto: 1]")
)

opt<- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
    stop("Uso: Rscript 04_clustering.R -i results/network/umbral_700/aristas.csv -o results/clusters")
}
if (!file.exists(opt$input)) {
    stop("No se encuentra el archivo de entrada: ", opt$input, ". El paso 02 debe ejecutarse antes del paso 04.")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# Comprobamos el método antes de empezar, para no hacer todo el análisis y fallar al final
metodos_validos <- c("edge_betweenness", "walktrap", "fast_greedy", "louvain", "leiden", "infomap", "label_propagation")
if (!(opt$metodo %in% metodos_validos)) {
  stop("Método no válido: ", opt$metodo, ". Opciones: ", paste(metodos_validos, collapse = ", "))
}

# 2. Leer la red 
# el paso 02 guarda las aristas con write.csv2 (separador ";")
aristas<- read.csv2(opt$input, stringsAsFactors = FALSE)
if (!all(c("from", "to") %in% colnames(aristas))) {
  stop("El archivo ", opt$input, " debe tener las columnas 'from' y 'to'")
}

g<- graph_from_data_frame(aristas[, c("from", "to")], directed = FALSE)
g<- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat("Red cargada:", gorder(g), "nodos y", gsize(g), "aristas\n")
cat("Componentes:", components(g)$no, "\n")

# 3. Aplicar los 7 métodos de clustering
comunidades<- list()
# Ponemos la semilla antes de cada método porque varios son aleatorios
set.seed(opt$semilla)
comunidades[["edge_betweenness"]]<- cluster_edge_betweenness(g)
set.seed(opt$semilla)
comunidades[["walktrap"]]<- cluster_walktrap(g)
set.seed(opt$semilla)
comunidades[["fast_greedy"]]<- cluster_fast_greedy(g)
set.seed(opt$semilla)
comunidades[["louvain"]]<- cluster_louvain(g, resolution = opt$resolucion)
set.seed(opt$semilla)
comunidades[["leiden"]]<- cluster_leiden(g, objective_function = "modularity", resolution = opt$resolucion)
set.seed(opt$semilla)
comunidades[["infomap"]]<- cluster_infomap(g)
set.seed(opt$semilla)
comunidades[["label_propagation"]]<- cluster_label_prop(g)

# Tabla con el número de clusters y la modularidad (Q) de cada método.
# Calculamos Q igual para todos, aunque algunos (infomap) no la optimicen
comparacion<- data.frame()
for (metodo in names(comunidades)) {
    comparacion<- rbind(comparacion, data.frame(
        metodo = metodo,
        n_clusters = length(comunidades[[metodo]]),
        tam_cluster_mayor = max(sizes(comunidades[[metodo]])),
        modularidad = modularity(g, membership(comunidades[[metodo]]))
    ))
}

cat("\nComparación de métodos:\n")
print(comparacion)

write.table(comparacion, file.path(opt$output, "comparacion_metodos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 4. NMI entre cada par de métodos
# Como no conocemos las comunidades reales, vemos cuánto coinciden los métodos entre sí (1 = particiones idénticas, 0 = nada en común)
metodos<- names(comunidades)
nmi<- matrix(NA, nrow = length(metodos), ncol = length(metodos), dimnames = list(metodos, metodos))

for(m1 in metodos) {
    for(m2 in metodos) {
        nmi[m1, m2]<- compare(comunidades[[m1]], comunidades[[m2]], method = "nmi")
    }
}

cat("\nNMI entre métodos:\n")
print(round(nmi, 2))

write.table(data.frame(metodo = metodos, round(nmi, 4), check.names = FALSE), file.path(opt$output, "nmi_metodos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap de la NMI con el valor escrito en cada casilla
pdf(file.path(opt$output, "nmi_metodos.pdf"), width = 8, height = 7)
par(mar = c(10, 10, 3, 2)) # pongo márgenes amplios para que quepan los nombres de los métodos
n_met<- length(metodos)
image(1:n_met, 1:n_met, nmi, col = hcl.colors(20, "YlOrRd", rev = TRUE), zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = "NMI entre métodos de clustering")

axis(1, at = 1:n_met, labels = metodos, las = 2)
axis(2, at = 1:n_met, labels = metodos, las = 2)
for(i in 1:n_met) {
    for(j in 1:n_met) {
        text(i, j, round(nmi[i, j], 2), col = "white") # que si no pongo white sale muy oscuro
    }
}
invisible(dev.off())

# 5. Estabilidad de Louvain con varias semillas
# Louvain es aleatorio: lo repetimos y comparamos cada resultado con el de nuestra semilla. Si la NMI es cercana a 1, la partición no depende de la semilla
referencia<- comunidades[["louvain"]]
estabilidad<- replicate(opt$n_semillas, {
    cl<- cluster_louvain(g, resolution = opt$resolucion)
    c(n_clusters =length(cl), modularidad = modularity(cl), nmi_referencia = compare(cl, referencia, method = "nmi"))
})
# la transponemos para tener una fila por repetición
estabilidad<- as.data.frame(t(estabilidad))

cat("\nEstabilidad de Louvain en", opt$n_semillas, "repeticiones:\n")
cat("NMI media con la partición de referencia:", round(mean(estabilidad$nmi_referencia), 3), "\n")
cat("Número de clusters: entre", min(estabilidad$n_clusters), "y", max(estabilidad$n_clusters), "\n")

print(table(estabilidad$n_clusters)) # cuántas repeticiones dan cada número de clusters

write.table(estabilidad, file.path(opt$output, "estabilidad_louvain.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

pdf(file.path(opt$output, "estabilidad_louvain.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2)) # dos gráficos uno al lado del otro
hist(estabilidad$nmi_referencia, main = "NMI con la partición de referencia", xlab = "NMI", col = "lightblue", border = "white")
barplot(table(estabilidad$n_clusters), main = "Número de clusters de Louvain", xlab = "Clusters", ylab = "Repeticiones", col = "lightblue")
invisible(dev.off())

# 6. Partición final
final<- comunidades[[opt$metodo]]
particion_final<- data.frame(
  gene = V(g)$name,
  cluster = as.vector(membership(final))
)
particion_final<- particion_final[order(particion_final$cluster, particion_final$gene), ]

write.table(particion_final, file.path(opt$output, "particion_final.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nPartición final con", opt$metodo, ":", length(final), "clusters\n")
print(split(particion_final$gene, particion_final$cluster))

# 7. Figura de la red coloreada por cluster
set.seed(opt$semilla)
layout<- layout_with_fr(g)

pdf(file.path(opt$output, "red_clusters.pdf"), width = 9, height = 9)
plot(final, g, layout = layout, vertex.size = 3 + 2 * sqrt(degree(g)), vertex.label.cex = 0.6, vertex.label.color = "black", main = paste("Clusters con", opt$metodo, "(resolución", opt$resolucion, ")"))
# el vertex.size es el tamaño del nodo según su grado (sqrt para que los hubs no salgan enormes)
invisible(dev.off())

cat("\nResultados guardados en:", opt$output, "\n")