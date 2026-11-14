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
    make_option(c("--umbral_final"), type = "integer", default = 700,
                help = "Umbral de STRING con el que se construye la red expandida [defecto: 700]"),
    make_option(c("--min_semillas"), type = "integer", default = 2,
                help = "Número mínimo de genes semilla con los que debe interaccionar un vecino para añadirlo [defecto: 2]"),
    make_option(c("--max_anadidos"), type = "integer", default = 50,
                help = "Número máximo de vecinos añadidos a la red expandida [defecto: 50]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
    stop(" Uso: Rscript code/02_red_string.R -i results/data/hpo_genes.tsv -o results/network")
}
if (!file.exists(opt$input)) {
    stop("No se encuentra el archivo de entrada: ", opt$input, ". El paso 01 debe generar results/data/hpo_genes.tsv antes de ejecutar el paso 02.")
}

# 2. Leer los genes del paso 01
genes <- read.delim(opt$input, colClasses = "character", quote = "")
if (!all(c("ncbi_id", "gene_symbol") %in% colnames(genes))) {
    stop("El fichero ", opt$input, " debe contener las columnas 'ncbi_id' y 'gene_symbol'")
}
if (nrow(genes) == 0 || any(is.na(genes$gene_symbol)) ||
    any(genes$gene_symbol == "")) {
    stop("La lista de genes está vacía o contiene símbolos ausentes.")
}
genes <- unique(genes[, c("ncbi_id", "gene_symbol")])

dir.create(opt$output, recursive= TRUE, showWarnings = FALSE)
carpeta_string <- file.path(opt$output, "string_data")
dir.create(carpeta_string, recursive= TRUE, showWarnings = FALSE)

# Configurar STRING para el flujo de trabajo
version_string <- "12.5"
especie <- 9606 # Homo sapiens
semilla <- 123
umbrales <- c(400, 700, 900)
resumen <- data.frame()

# la expansión se hace dentro del bucle, así que el umbral final tiene que ser uno de los tres
if (!(opt$umbral_final %in% umbrales)) {
    stop("El umbral final debe ser uno de: ", paste(umbrales, collapse = ", "))
}

cat("Genes HPO de entrada: ", length(unique(genes$gene_symbol)), "\n")
cat("STRING:", version_string, "Especie:", especie, "\n")

# configuración de STRING con el umbral más bajo, los demás se filtran dentro del bucle
string_db <- STRINGdb$new(
    version= version_string,
    species = especie,
    score_threshold = min(umbrales),
    input_directory = carpeta_string 
)

# Mapeo de genes a ids STRING
string_ids <- string_db$map(
    genes,
    "gene_symbol",
    removeUnmappedRows = TRUE
)
cat("\n") # el aviso de map no termina en salto de línea

if (nrow(string_ids) == 0){
    stop("No se encontraron identificadores STRING para los genes proporcionados. Revisa la entrada.")
}

# Obtención de interacciones para cada umbral
for (umbral in umbrales){
    cat("\n--- Umbral STRING:", umbral, " ---\n")
    carpeta_umbral <- file.path(opt$output, paste0("umbral_", umbral))
    dir.create(carpeta_umbral, recursive = TRUE, showWarnings = FALSE)
    
    # 3. Mapear y obtener interacciones

    # Comprobar genes no mapeados
    no_mapeados <- setdiff(genes$gene_symbol, string_ids$gene_symbol)
    cat("HPO genes: ", nrow(genes), "\n")
    cat("Genes mapeados a STRING:", nrow(string_ids), "\n")

    if (length(no_mapeados) > 0){
        cat("Genes no mapeados a STRING:", length(no_mapeados), "\n")
        cat("Ejemplos de genes no mapeados:", paste(head(no_mapeados, 10), collapse = ", "), "\n")
    } else {
        cat("Todos los genes se mapearon correctamente a identificadores STRING.\n")
    }

    # Obtener interacciones de STRING
    interaction_network <- string_db$get_interactions(string_ids$STRING_id)
    interaction_network <- interaction_network[interaction_network$combined_score >= umbral, ]
    interaction_network <- unique(interaction_network[, 1:2])  # eliminar duplicados

    # sustituir ids STRING por símbolos de genes
    interaction_network$from <- string_ids$gene_symbol[match(interaction_network$from, string_ids$STRING_id)]
    interaction_network$to <- string_ids$gene_symbol[match(interaction_network$to, string_ids$STRING_id)]
    
    if (any(is.na(interaction_network$from)) || 
        any(is.na(interaction_network$to))) {
        stop("Error al mapear los identificadores STRING a símbolos de genes. 
        Revisa el mapeo.")
    }

    # 4. Construir la red de interacciones
    g_umbral <- graph_from_data_frame(interaction_network, directed = FALSE)

    # Genes mapeados a STRING que no estan en la red
    sin_interacciones <- setdiff(string_ids$gene_symbol, V(g_umbral)$name)
    cat("Genes mapeados a STRING sin interacciones:", length(sin_interacciones), "\n")
    if (length(sin_interacciones) > 0){
        cat("Ejemplos de genes sin interacciones:", paste(head(sin_interacciones, 10), collapse = ", "), "\n")
    }

    # 5. Guardar red y resumen descriptivo 
    comp <- components(g_umbral)
    componente_mayor <- 0
    if (gorder(g_umbral) > 0){
        componente_mayor <- max(comp$csize)
    }
    
    resumen <- rbind(resumen, data.frame(
        umbral = umbral,
        genes_mapeados = nrow(string_ids),
        genes_sin_interacciones = length(sin_interacciones),
        interacciones = ecount(g_umbral),
        nodos_red = vcount(g_umbral),
        aristas_red = ecount(g_umbral),
        porcentaje_genes_red = 100 * vcount(g_umbral) / nrow(string_ids),
        componentes = comp$no,
        componente_mayor = componente_mayor
    ))


    # 6. Exportar tablas a TSV como dataframes
    write.csv2(
        string_ids, 
        file.path(carpeta_umbral, "mapeo_genes.csv"), 
        row.names = FALSE
        )
    write.csv2(
        data.frame(gene_symbol = no_mapeados),
        file.path(carpeta_umbral, "genes_no_mapeados.csv"), 
        row.names = FALSE
        )
    write.csv2(
        data.frame(gene_symbol = sin_interacciones),
        file.path(carpeta_umbral, "genes_sin_interacciones.csv"),
        row.names = FALSE
    )
    write.csv2(
        data.frame(gene_symbol = V(g_umbral)$name),
        file.path(carpeta_umbral, "nodos.csv"),
        row.names = FALSE   
    )
    write.csv2(
        as_data_frame(g_umbral, what = "edges"),
        file.path(carpeta_umbral, "aristas.csv"),
        row.names = FALSE   
    )

    # 7. Crear figura con semilla para cada red

    if (gorder(g_umbral) > 0){
        pdf(file.path(carpeta_umbral, paste0("red_umbral_", umbral, ".pdf")),
            width = 8, height = 8)
        set.seed(semilla)
        layout <- layout_with_fr(g_umbral)
       
        plot(g_umbral, layout = layout, vertex.color = "skyblue", edge.color = "deepskyblue",
             vertex.label.color = "black",
            vertex.label.cex = 0.7, vertex.size = 5,
             main = paste("Red STRING (umbral", umbral, ")"))
        dev.off()
    } else {
        cat("No se puede generar la figura de la red STRING para el umbral", umbral, "ya que no hay nodos en la red.\n")
    }

    cat("Genes no mapeados a STRING:", length(no_mapeados), "\n")
    cat("Genes mapeados a STRING sin interacciones:", length(sin_interacciones), "\n")
    cat("Nodos:", vcount(g_umbral), "Aristas:", ecount(g_umbral), "\n")
    cat("Componente mayor: ", componente_mayor, "\n")    

    # 8. Red expandida (solo con el umbral final)
    # Con solo las semillas muchos genes quedan aislados, así que añadimos vecinos de STRING
    # que interaccionan con varias semillas a la vez
    if (umbral == opt$umbral_final) {
        cat("\n--- Red expandida (umbral", umbral, ") ---\n")

        # get_interactions solo da aristas entre las proteínas que le pasamos,
        # para ver los vecinos hace falta el grafo completo de STRING
        g_string <- string_db$get_graph()
        # el objeto se creó con el umbral más bajo, así que quitamos las aristas por debajo
        # del umbral final y las proteínas que se quedan sin ninguna (antes de simplify, que borra combined_score)
        g_string <- delete_edges(g_string, E(g_string)[combined_score < umbral])
        g_string <- delete_vertices(g_string, V(g_string)[degree(g_string) == 0])
        # cada interacción viene repetida, por eso simplificamos
        g_string <- igraph::simplify(g_string, remove.multiple = TRUE, remove.loops = TRUE)
        proteinas <- string_db$get_proteins()

        semillas_ids <- intersect(string_ids$STRING_id, V(g_string)$name)
        N <- vcount(g_string)
        K <- length(semillas_ids)

        # Juntamos los vecinos de todas las semillas, y cuantas más veces sale una proteína,
        # con más semillas interacciona. De paso contamos cuántas semillas toca cada semilla
        vecinos <- c()
        k_semillas <- c()
        for (id_semilla in semillas_ids) {
            vecinos_semilla <- names(neighbors(g_string, id_semilla))
            vecinos <- c(vecinos, vecinos_semilla)
            k_semillas <- c(k_semillas, sum(vecinos_semilla %in% semillas_ids))
        }
        vecinos <- vecinos[!(vecinos %in% semillas_ids)]
        conteo <- table(vecinos)

        candidatos <- data.frame(
            STRING_id = names(conteo),
            k = as.vector(conteo)
        )
        candidatos$grado_string <- as.vector(degree(g_string, candidatos$STRING_id))
        cat("Vecinos de las semillas evaluados:", nrow(candidatos), "\n")

        candidatos <- candidatos[candidatos$k >= opt$min_semillas, ]
        cat("Vecinos que interaccionan con al menos", opt$min_semillas, "semillas:", nrow(candidatos), "\n")

        # Con la hipergeométrica vemos la probabilidad de tocar k o más semillas por azar con su grado.
        # Así no entran los hubs que conectan con todo
        candidatos$p_valor <- phyper(candidatos$k - 1, K, N - K, candidatos$grado_string, lower.tail = FALSE)
        candidatos$p_ajustado <- p.adjust(candidatos$p_valor, method = "BH")
        candidatos <- candidatos[candidatos$p_ajustado < 0.05, ]
        cat("Vecinos con p ajustado < 0.05:", nrow(candidatos), "\n")

        # si siguen siendo muchos nos quedamos con los más significativos (en empate, más semillas)
        candidatos <- candidatos[order(candidatos$p_ajustado, -candidatos$k), ]
        if (nrow(candidatos) > opt$max_anadidos) {
            candidatos <- candidatos[1:opt$max_anadidos, ]
        }
        candidatos$gene <- proteinas$preferred_name[match(candidatos$STRING_id, proteinas$protein_external_id)]
        candidatos$tipo <- rep("anadido", nrow(candidatos))
        cat("Genes añadidos:", nrow(candidatos), "\n")
        cat(paste(candidatos$gene, collapse = ", "), "\n")

        # En la tabla de nodos las semillas van con el símbolo de HPO y los añadidos con el nombre de STRING
        tabla_semillas <- data.frame(
            STRING_id = semillas_ids,
            k = k_semillas,
            grado_string = as.vector(degree(g_string, semillas_ids)),
            p_valor = NA,
            p_ajustado = NA,
            gene = string_ids$gene_symbol[match(semillas_ids, string_ids$STRING_id)],
            tipo = "semilla"
        )
        nodos_exp <- rbind(tabla_semillas, candidatos)

        # Red con las interacciones entre semillas y añadidos
        g_expandida <- induced_subgraph(g_string, nodos_exp$STRING_id)
        V(g_expandida)$name <- nodos_exp$gene[match(V(g_expandida)$name, nodos_exp$STRING_id)]
        # quitamos los nodos sin aristas, igual que en la red de semillas
        g_expandida <- delete_vertices(g_expandida, V(g_expandida)[degree(g_expandida) == 0])
        nodos_exp <- nodos_exp[nodos_exp$gene %in% V(g_expandida)$name, ]

        semillas_red <- nodos_exp$gene[nodos_exp$tipo == "semilla"]
        cat("Semillas que estaban aisladas y ahora quedan conectadas:",
            sum(sin_interacciones %in% semillas_red), "de", length(sin_interacciones), "\n")
        cat("Semillas que siguen fuera de la red:",
            paste(setdiff(string_ids$gene_symbol, semillas_red), collapse = ", "), "\n")
        cat("Nodos:", vcount(g_expandida), "Aristas:", ecount(g_expandida),
            "Componentes:", components(g_expandida)$no, "\n")

        write.csv2(
            as_data_frame(g_expandida, what = "edges"),
            file.path(carpeta_umbral, "aristas_expandida.csv"),
            row.names = FALSE
        )
        write.csv2(
            nodos_exp[, c("gene", "tipo", "k", "grado_string", "p_valor", "p_ajustado")],
            file.path(carpeta_umbral, "nodos_expandida.csv"),
            row.names = FALSE, na = ""
        )
        # lista de semillas para la prueba de módulo de enfermedad del paso 03
        writeLines(semillas_red, file.path(carpeta_umbral, "genes_semilla.txt"))

        # Figura con semillas y añadidos en colores distintos
        colores <- ifelse(V(g_expandida)$name %in% semillas_red, "skyblue", "orange")
        pdf(file.path(carpeta_umbral, paste0("red_expandida_umbral_", umbral, ".pdf")),
            width = 10, height = 10)
        set.seed(semilla)
        layout <- layout_with_fr(g_expandida)
        plot(g_expandida, layout = layout, vertex.color = colores, edge.color = "gray80",
             vertex.label.color = "black", vertex.label.cex = 0.6, vertex.size = 5,
             main = paste("Red expandida (umbral", umbral, ")"))
        legend("bottomleft", legend = c("Semilla (HPO)", "Añadido"),
               fill = c("skyblue", "orange"), bty = "n")
        invisible(dev.off())
    }
}

write.csv2(
    resumen,
    file.path(opt$output, "sensibilidad_umbrales.csv"),
    row.names = FALSE
)
cat("Resumen de los tres umbrales; \n")
print(resumen)
cat("\nResultados guardados en:", file.path(opt$output, "sensibilidad_umbrales.csv"), "\n")

