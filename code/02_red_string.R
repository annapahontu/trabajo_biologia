#!/usr/bin/env Rscript

library(optparse)
library(STRINGdb)
library(igraph)

#1. Argumentos
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Entrada del proyecto: results/data/hpo_genes.tsv (generada por el paso 01)"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Salida del proyecto: results/network (redes, tablas y descargas de STRING)")
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

cat("Genes HPO de entrada: ", length(unique(genes$gene_symbol)), "\n")
cat("STRING:", version_string, "Especie:", especie, "\n")

# Obtención de interacciones para cada umbral
for (umbral in umbrales){
    cat("\n--- Umbral STRING:", umbral, " ---\n")
    carpeta_umbral <- file.path(opt$output, paste0("umbral_", umbral))
    dir.create(carpeta_umbral, recursive = TRUE, showWarnings = FALSE)
    
    # 3. Mapear y obtener interacciones

    # configuración de STRING para el umbral actual
    string_db <- STRINGdb$new(
        version= version_string,
        species = especie,
        score_threshold = umbral,
        input_directory = carpeta_string 
    )

    # Mapeo de genes a ids STRING
    string_ids <- string_db$map(
        genes,
        "gene_symbol",
        removeUnmappedRows = TRUE
    )

    if (nrow(string_ids) == 0){
        stop("No se encontraron identificadores STRING para los genes proporcionados. Revisa la entrada.")
    }

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
             vertex.label.color = "lightcyan",
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
}

write.csv2(
    resumen,
    file.path(opt$output, "sensibilidad_umbrales.csv"),
    row.names = FALSE
)
cat("Resumen de los tres umbrales; \n")
print(resumen)
cat("\nResultados guardados en:", file.path(opt$output, "sensibilidad_umbrales.csv"), "\n")

