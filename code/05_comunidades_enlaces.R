#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(igraph)
})

# 1. Argumentos
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, help = "Aristas de la red: results/network/umbral_700/aristas.csv (paso 02)"),
    make_option(c("-o", "--output"), type = "character", default = NULL, help = "Carpeta de salida: results/clusters/enlaces"),
    make_option(c("-s", "--semilla"), type = "integer", default = 123, help = "Semilla aleatoria [defecto: 123]"),
    make_option(c("--particion"), type = "character", default = NULL, help = "Opcional, particion_final.tsv del paso 04 (gene, cluster)"),
    make_option(c("--metricas"), type = "character", default = NULL, help = "Opcional, metricas_nodos.tsv del paso 03 (para la betweenness)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
    stop("Uso: Rscript 05_comunidades_enlaces.R -i results/network/umbral_700/aristas.csv -o results/clusters/enlaces")
}
# los opcionales que no se pasan son NULL y no entran en el vector
for (archivo in c(opt$input, opt$particion, opt$metricas)) {
    if (!file.exists(archivo)) {
        stop("No se encuentra el archivo de entrada: ", archivo)
    }
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# 2. Leer la red y los archivos opcionales
aristas <- read.csv2(opt$input, stringsAsFactors = FALSE)
if (!all(c("from", "to") %in% colnames(aristas))) {
    stop("El archivo ", opt$input, " debe tener las columnas 'from' y 'to'")
}
g <- graph_from_data_frame(aristas[, c("from", "to")], directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
cat("Red cargada:", gorder(g), "nodos y", gsize(g), "aristas\n")

if (!is.null(opt$particion)) {
    particion <- read.delim(opt$particion, stringsAsFactors = FALSE)
}
if (!is.null(opt$metricas)) {
    metricas <- read.delim(opt$metricas, stringsAsFactors = FALSE)
}

# 3. Funciones de comunidades de enlaces (Ahn et al., 2010)

# Similitud entre dos aristas que comparten un nodo (Jaccard de los vecindarios de los otros extremos)
similitud_enlaces <- function(g, e1, e2) {
    n1 <- ends(g, e1)
    n2 <- ends(g, e2)
    k <- intersect(n1, n2)            # el nodo compartido
    if (length(k) != 1) return(0)     # no comparten nodo
    i <- setdiff(n1, k)
    j <- setdiff(n2, k)
    vec_i <- c(i, names(neighbors(g, i)))   # vecindario inclusivo de i
    vec_j <- c(j, names(neighbors(g, j)))   # vecindario inclusivo de j
    length(intersect(vec_i, vec_j)) / length(union(vec_i, vec_j))
}

# Matriz de similitud entre todas las aristas de g
matriz_similitud <- function(g) {
    M <- ecount(g)
    S <- diag(M)
    for (a in 1:(M - 1)) {
        for (b in (a + 1):M) {
            S[a, b] <- S[b, a] <- similitud_enlaces(g, a, b)
        }
    }
    nombres <- apply(as_edgelist(g), 1, paste, collapse = "-")
    dimnames(S) <- list(nombres, nombres)
    S
}

# Comunidades de cada nodo = grupos de todas sus aristas
comunidades_nodos <- function(g, grupos) {
    extremos <- as_edgelist(g)
    tapply(c(grupos, grupos), c(extremos[, 1], extremos[, 2]),
           function(x) sort(unique(x)))
}

# Densidad de partición de un conjunto de grupos de aristas
densidad_particion <- function(g, grupos) {
    D <- 0
    for (c in unique(grupos)) {
        aristas <- which(grupos == c)
        m <- length(aristas)
        n <- length(unique(as.vector(ends(g, aristas))))
        if (n > 2) D <- D + m * (m - (n - 1)) / (n * (n - 1) / 2 - (n - 1))
    }
    D / ecount(g)
}

# 4. Similitud entre aristas y clustering jerárquico
S <- matriz_similitud(g)
# single linkage, como en el método original
hc <- hclust(as.dist(1 - S), method = "single")

# 5. Corte del árbol donde la densidad de partición es máxima
alturas <- sort(unique(hc$height))
D <- sapply(alturas, function(h) densidad_particion(g, cutree(hc, h = h)))
mejor <- alturas[which.max(D)]
grupos <- cutree(hc, h = mejor)

cat("\nAltura de corte elegida:", round(mejor, 3), "\n")
cat("Densidad de partición máxima:", round(max(D), 3), "\n")

# Las comunidades de una sola arista no aportan nada a la densidad, así que las descartamos (quedan como NA).
# Renumeramos las demás de mayor a menor para que la 1 sea la más grande
tam <- sort(table(grupos), decreasing = TRUE)
validas <- names(tam)[tam > 1]
grupos_final <- match(as.character(grupos), validas)
n_com <- length(validas)

cat("Comunidades con más de una arista:", n_com, "\n")
cat("Aristas sueltas descartadas:", sum(is.na(grupos_final)), "\n")
cat("Tamaño de cada comunidad (aristas):", paste(as.vector(tam[tam > 1]), collapse = ", "), "\n")

# Genes de cada comunidad, que usamos para el log y para comparar con Louvain
extremos <- as_edgelist(g)
genes_com <- list()
for (k in 1:n_com) {
    a <- which(grupos_final == k)
    genes_com[[k]] <- sort(unique(c(extremos[a, 1], extremos[a, 2])))
}
cat("\nGenes de cada comunidad\n")
for (k in 1:n_com) {
    cat(k, "-", paste(genes_com[[k]], collapse = ", "), "\n")
}

# 6. Comunidades de cada gen
# sort() quita los NA, así que las aristas descartadas no cuentan
com_nodos <- comunidades_nodos(g, grupos_final)
nodos_com <- data.frame(gene = names(com_nodos), n_comunidades = 0, comunidades = "")
for (i in 1:nrow(nodos_com)) {
    nodos_com$n_comunidades[i] <- length(com_nodos[[i]])
    nodos_com$comunidades[i] <- paste(com_nodos[[i]], collapse = ",")
}
nodos_com <- nodos_com[order(-nodos_com$n_comunidades, nodos_com$gene), ]

write.table(data.frame(from = extremos[, 1], to = extremos[, 2], comunidad = grupos_final),
            file.path(opt$output, "aristas_comunidades.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(nodos_com, file.path(opt$output, "nodos_comunidades.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# genes puente = genes con aristas en más de una comunidad
puentes <- nodos_com[nodos_com$n_comunidades > 1, ]
cat("\nGenes en más de una comunidad:", nrow(puentes), "\n")
print(puentes, row.names = FALSE)

# 7. Figuras
pdf(file.path(opt$output, "densidad_particion.pdf"), width = 7, height = 5)
plot(alturas, D, type = "b", pch = 19, xlab = "Altura de corte", ylab = "Densidad de partición",
     main = paste("Densidad de partición (máximo", round(max(D), 3), "en", round(mejor, 3), ")"))
abline(v = mejor, col = "red", lty = 2)
points(mejor, max(D), col = "red", pch = 19, cex = 1.5)
invisible(dev.off())

# con tantas aristas los nombres de las hojas no se leen, así que no los ponemos
pdf(file.path(opt$output, "dendrograma_enlaces.pdf"), width = 10, height = 6)
plot(hc, labels = FALSE, main = "Dendrograma de aristas", xlab = "", sub = "", ylab = "Altura (1 - S)")
abline(h = mejor, col = "red", lty = 2)
invisible(dev.off())

# aristas coloreadas por comunidad y las descartadas en gris
colores_aristas <- rep("gray85", ecount(g))
colores_aristas[!is.na(grupos_final)] <- rainbow(n_com)[grupos_final[!is.na(grupos_final)]]
colores_nodos <- rep("gray90", vcount(g))

set.seed(opt$semilla)
layout <- layout_with_fr(g)
# el núcleo ciliar es muy denso, con la figura grande y nodos pequeños se leen mejor las etiquetas
pdf(file.path(opt$output, "red_comunidades_enlaces.pdf"), width = 14, height = 14)
plot(g, layout = layout, edge.color = colores_aristas, edge.width = 2, vertex.color = colores_nodos,
     vertex.size = 3, vertex.label.cex = 0.7, vertex.label.color = "black",
     main = paste("Comunidades de enlaces (", n_com, "comunidades)"))
invisible(dev.off())

# 8. Genes puente frente a betweenness (solo con --metricas)
# Queremos ver si los genes que están en varias comunidades son los que hacen de puente en la red
if (!is.null(opt$metricas)) {
    nodos_com$betweenness <- metricas$betweenness[match(nodos_com$gene, metricas$gene)]
    puentes <- nodos_com[nodos_com$n_comunidades > 1, ]
    write.table(puentes[, c("gene", "n_comunidades", "betweenness")], file.path(opt$output, "genes_varias_comunidades.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    top <- nodos_com[order(-nodos_com$betweenness), ]
    top <- top[1:min(5, nrow(top)), ]
    cat("\nGenes con mayor betweenness y su número de comunidades\n")
    print(top[, c("gene", "betweenness", "n_comunidades")], row.names = FALSE)
    cat("De los", nrow(top), "genes con mayor betweenness,", sum(top$n_comunidades > 1), "son genes puente\n")

    pdf(file.path(opt$output, "puentes_betweenness.pdf"), width = 7, height = 6)
    es_puente <- nodos_com$n_comunidades > 1
    plot(nodos_com$n_comunidades, nodos_com$betweenness, pch = 19,
         col = ifelse(es_puente, "red", "gray60"), xlim = c(min(nodos_com$n_comunidades), max(nodos_com$n_comunidades) + 1.5),
         xlab = "Número de comunidades de enlaces", ylab = "Betweenness", main = "Genes puente frente a betweenness")
    # si hay muchos puentes las etiquetas se pisan, así que solo ponemos las de los 15 con más betweenness
    etiquetas <- puentes[order(-puentes$betweenness), ]
    etiquetas <- etiquetas[seq_len(min(15, nrow(etiquetas))), ]
    if (nrow(etiquetas) > 0) {
        text(etiquetas$n_comunidades, etiquetas$betweenness, etiquetas$gene, pos = 4, cex = 0.7)
    }
    invisible(dev.off())
}

# 9. Comparación con la partición de Louvain (solo con --particion)
if (!is.null(opt$particion)) {
    # cuántos genes de cada comunidad de enlaces caen en cada cluster
    filas_com <- c()
    filas_cluster <- c()
    for (k in 1:n_com) {
        filas_com <- c(filas_com, rep(k, length(genes_com[[k]])))
        filas_cluster <- c(filas_cluster, particion$cluster[match(genes_com[[k]], particion$gene)])
    }
    contingencia <- table(comunidad = filas_com, cluster = filas_cluster)
    write.table(data.frame(comunidad = rownames(contingencia), as.data.frame.matrix(contingencia), check.names = FALSE),
                file.path(opt$output, "contingencia_louvain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    mezcladas <- sum(rowSums(contingencia > 0) > 1)
    cat("\nComunidades de enlaces con genes de más de un cluster de Louvain:", mezcladas, "de", n_com, "\n")
}

cat("\nResultados guardados en:", opt$output, "\n")
