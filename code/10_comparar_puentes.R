#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(igraph)
})

# 1. Argumentos
option_list <- list(
    make_option(c("-r", "--results"), type = "character", default = NULL,
                help = "Carpeta principal de resultados"),
    make_option(c("-u", "--umbral"), type = "integer", default = 700,
                help = "Umbral de STRING [defecto: 700]"),
    make_option(c("-n", "--n_puentes"), type = "character", default = "3,5,10",
                help = "N separados por comas [defecto: 3,5,10]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Carpeta de salida")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$results) || is.null(opt$output)) {
    stop("Uso: Rscript 10_comparar_puentes.R -r results -o results/comparacion_puentes")
}

dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

n_puentes <- as.integer(strsplit(opt$n_puentes, "[, ]+")[[1]])

if (any(is.na(n_puentes))) {
    stop("Formato incorrecto en --n_puentes")
}

casos <- c(0, n_puentes)

# 2. Funciones auxiliares
leer_red <- function(archivo) {
    aristas <- read.csv2(archivo, stringsAsFactors = FALSE)
    if (!all(c("from", "to") %in% colnames(aristas))) {
        stop("Faltan columnas from/to en ", archivo)
    }

    g <- graph_from_data_frame(aristas[, c("from", "to")],  directed = FALSE)
    simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
}

leer_particion <- function(archivo) {
    part <- read.delim(archivo,stringsAsFactors = FALSE)

    if (!all(c("gene", "cluster") %in% colnames(part))) {
        stop("Partición inválida: ", archivo)
    }
    part
}

leer_resumen_ora <- function(archivo) {
    if (!file.exists(archivo)) {
        return(data.frame())
    }
    read.csv2(archivo, stringsAsFactors = FALSE)
}

# 3. Datos base
semillas <- read.delim(file.path(opt$results, "data", "hpo_genes.tsv"),stringsAsFactors = FALSE)
semillas <- unique(semillas$gene_symbol)

# red original

particion_original <- leer_particion(
    file.path(opt$results, "clusters", paste0("umbral_", opt$umbral),"particion_final.tsv")
)

g_original <- leer_red(
    file.path(opt$results,"network",paste0("umbral_", opt$umbral),"aristas.csv")
)

# 4. Comparación global
comparacion <- data.frame()
grupos_unidos <- data.frame()

for (N in casos) {
    cat("===Comparando redes con: ", N," puentes===\n")

    if (N == 0) {
        carpeta_red <- file.path(opt$results,"network",paste0("umbral_", opt$umbral))
        carpeta_cluster <- file.path(opt$results, "clusters",paste0("umbral_", opt$umbral))
        carpeta_ora <- file.path(opt$results,"ora",paste0("umbral_", opt$umbral))
    } else {
        carpeta_red <- file.path(opt$results,"network", paste0("puentes_", opt$umbral), paste0("n_", N))
        carpeta_cluster <- file.path(opt$results,"clusters",paste0("puentes_", opt$umbral), paste0("n_", N))
        carpeta_ora <- file.path(opt$results, "ora", paste0("puentes_", opt$umbral),paste0("n_", N))
    }

    red_archivo <- file.path(carpeta_red, "aristas.csv")

    if (!file.exists(red_archivo)) {
        warning("No existe: ", red_archivo)
        next
    }

    g <- leer_red(red_archivo)
    comp <- components(g)
    lcc <- max(comp$csize)
    semillas_conectadas <- sum(semillas %in% V(g)$name)
    semillas_aisladas <- length(semillas) - semillas_conectadas

    # clusters
    archivo_particion <- file.path(carpeta_cluster,"particion_final.tsv")
    part <- leer_particion(archivo_particion)
    # el cluster de cada nodo tiene que ir en el mismo orden que V(g)
    miembros <- part$cluster[match(V(g)$name, part$gene)]
    if (any(is.na(miembros))) {
        stop("Hay nodos de la red sin cluster en ", archivo_particion)
    }
    mod <- modularity(g, as.integer(as.factor(miembros)))

    n_clusters <- length(unique(part$cluster))

    # ORA
    archivo_ora <- file.path(carpeta_ora, "resumen_comunidades.csv")
    clusters_con_ora <- NA
    clusters_analizables <- NA

    if (file.exists(archivo_ora)) {
        ora <- leer_resumen_ora(archivo_ora)
        if (nrow(ora) > 0) {
            clusters_analizables <-  sum(ora$n_genes >= 3)
            clusters_con_ora <- sum(ora$estado_ORA == "con_terminos_significativos")
        }
    }

    # NMI
    genes_comunes <- intersect(particion_original$gene, part$gene)
    m1 <- particion_original$cluster[match(genes_comunes, particion_original$gene)]
    m2 <- part$cluster[match(genes_comunes, part$gene)]
    nmi <- compare(as.numeric(as.factor(m1)), as.numeric(as.factor(m2)), method = "nmi")

    # grupos unidos
    if (N > 0) {
        originales <- split(particion_original$gene,particion_original$cluster)
        ids <- names(originales)
        if (length(ids) > 1) {
            for (i in 1:(length(ids) - 1)) {
                for (j in (i + 1):length(ids)) {
                    c1 <- ids[i]
                    c2 <- ids[j]
                    g1 <- originales[[c1]]
                    g2 <- originales[[c2]]
                    nuevo1 <- unique(part$cluster[part$gene %in% g1])
                    nuevo2 <- unique(part$cluster[part$gene %in% g2])
                    if (length(intersect(nuevo1, nuevo2)) > 0) {
                        grupos_unidos <- rbind(grupos_unidos,data.frame(
                                puentes = N,
                                cluster_original_1 = c1,
                                cluster_original_2 = c2,
                                nuevo_cluster = paste(intersect(nuevo1, nuevo2), collapse = ",")
                            )
                        )
                    }
                }
            }
        }
    }

    comparacion <- rbind(
        comparacion,
        data.frame(
            puentes = N,
            nodos = gorder(g),
            aristas = gsize(g),
            componentes = comp$no,
            componente_gigante = lcc,
            fraccion_gigante = round(lcc / gorder(g), 3),
            semillas_conectadas = semillas_conectadas,
            semillas_aisladas = semillas_aisladas,
            clusters = n_clusters,
            modularidad = round(mod, 4),
            clusters_analizables = clusters_analizables,
            clusters_con_ORA = clusters_con_ora,
            NMI_vs_original = round(nmi, 4)
        )
    )
}

# 5. Guardar tablas
write.table(comparacion, file.path(opt$output, "comparacion_puentes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(grupos_unidos,file.path(opt$output,"grupos_unidos.tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
cat("\nResumen:\n")
print(comparacion)

# 6. Figura principal
pdf(file.path(opt$output,"comparacion_puentes.pdf"), width = 10, height = 8, pointsize = 12)
par(mfrow = c(2, 2))
plot(comparacion$puentes, comparacion$semillas_conectadas, type = "b", pch = 19, col = "#0072B2", 
    xlab = "Genes puente", ylab = "Semillas conectadas", main = "Conectividad")
max_clusters <- max(comparacion$clusters) + 1
plot(comparacion$puentes, comparacion$clusters, type = "b", pch = 19, col = "#0072B2",
    xlab = "Genes puente", ylab = "Clusters", main = "Número de clusters",
    ylim = c(0, max_clusters), yaxt = "n")
axis(2, at = 0:max_clusters)
plot(comparacion$puentes, comparacion$modularidad, type = "b", pch = 19, col = "#0072B2",
    xlab = "Genes puente", ylab = "Modularidad", main = "Modularidad")
# el tope del eje son los clusters que tenían genes suficientes para el ORA
max_ora <- max(comparacion$clusters_analizables, na.rm = TRUE) + 1
plot(comparacion$puentes, comparacion$clusters_con_ORA, type = "b", pch = 19,
    col = "#0072B2", xlab = "Genes puente", ylab = "Clusters con ORA", main = "Enriquecimiento funcional",
    ylim = c(0, max_ora), yaxt = "n")
axis(2, at = 0:max_ora)

invisible(dev.off())

# 7. Red con puentes coloreada por cluster
# Cada color es un cluster y los genes puente van en cuadrado y con su nombre
for(n_red in n_puentes){
    caso <- file.path(paste0("puentes_", opt$umbral), paste0("n_", n_red))
    carpeta_red <- file.path(opt$results, "network", caso)
    carpeta_cluster <- file.path(opt$results, "clusters", caso)

    archivo_aristas <- file.path(carpeta_red, "aristas.csv")
    if (!file.exists(file.path(archivo_aristas))) {
        stop("Falta el archivo de aristas de la red con ", n_red, " puentes: ", archivo_aristas)
    }

    g <- leer_red(archivo_aristas)
    part <- leer_particion(file.path(carpeta_cluster, "particion_final.tsv"))
    puentes <- read.delim(file.path(carpeta_red, "genes_puente.tsv"), stringsAsFactors = FALSE)$gen

    cluster_nodo <- part$cluster[match(V(g)$name, part$gene)]
    # un color por cluster (paletas Okabe-Ito y Paul Tol, sin rojos ni naranjas); si hay más clusters que colores, se interpolan
    colores <- c("#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#332288",
                  "#44AA99", "#999933", "#AA4499", "#88CCEE", "#117733", "#999999")
    n_cl <- length(unique(cluster_nodo))
    if (n_cl > length(colores)) colores <- colorRampPalette(colores)(n_cl)
    colores <- colores[seq_len(n_cl)]
    es_puente <- V(g)$name %in% puentes

    set.seed(123)
    pdf(file.path(opt$output, paste0("red_puentes_clusters_n_", n_red, ".pdf")), width = 9, height = 9, pointsize = 12)
    
    plot(g, layout = layout_with_fr(g),
         vertex.color = colores[as.integer(as.factor(cluster_nodo))],
         vertex.shape = ifelse(es_puente, "square", "circle"),
         vertex.size = ifelse(es_puente, 8, 5),
         vertex.label = V(g)$name,
         # las semillas con el nombre pequeño y en gris, los puentes más grandes y en negrita
         vertex.label.cex = ifelse(es_puente, 0.9, 0.6),
         vertex.label.font = ifelse(es_puente, 2, 1),
         vertex.label.color = ifelse(es_puente, "black", "grey30"),
         vertex.label.dist = 1,
         edge.width = 0.5, edge.color = "grey70")
    legend("bottomleft", legend = c(paste("Cluster", sort(unique(cluster_nodo))), "Gen puente"),
           pch = c(rep(21, length(colores)), 22), pt.bg = c(colores, "white"),
           pt.cex = 1.5, bty = "n", cex = 0.8)
    invisible(dev.off())
}

cat("\nResultados guardados en:", opt$output, "\n")