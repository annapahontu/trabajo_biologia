#!/usr/bin/env Rscript

paquetes <- c("optparse", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("No se puede cargar el paquete ", p, "
    . Ejecuta bash code/setup.sh y revisa su log.")
  }
}

library(optparse)
library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Argumentos
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Carpeta con los resultados del 01: results/data"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Carpeta de salida para los resultados del ORA: results/ORA"),
    make_option(c("-s", "--semilla"), type = "integer", default = 123,
              help = "Semilla aleatoria para reproducibilidad [defecto: 123]")
)
opt <- parse_args(OptionParser(option_list = option_list))


if (is.null(opt$input)) {
  stop("Error: Debe especificar un archivo de entrada con -i o --input")
}

if (!dir.exists(opt$input)){
  stop("No se encuentra la carpeta de entrada: ", opt$input)
}


if (is.na(opt$semilla)) {
  stop("La semilla debe ser un número entero.")
}

set.seed(opt$semilla) # Establece la semilla para sample_gnm y rewire
# dir.create(opt$output, recursive = TRUE, showWarnings = FALSE) # Directorio de salida

# obtención de tablas de entrada
carpeta_results <- dirname(normalizePath(opt$input))
umbrales <- c(400, 700, 900)

# función para leer las tablas de entrada
leer_tabla <- function(archivo, columnas, sep, permitir_vacia = FALSE){
  if (!file.exists(archivo)) {
    stop("No se encuentra el archivo: ", archivo)
  }
  tabla <- tryCatch(
    read.table(file = archivo, header = TRUE, sep = sep, quote = "\"", comment.char = "", colClasses = "character"),
    error = function(e) stop("Error al leer el archivo ", archivo, ": ", conditionMessage(e))
  )
  if (!all(columnas %in% colnames(tabla))) {
    stop("El archivo ", archivo, " no contiene las columnas requeridas: ", paste(columnas, collapse = ", "))
  }
  if (!permitir_vacia && nrow(tabla) == 0) {
    stop("El archivo ", archivo, " está vacío.")
  }

  for (col in columnas) {
    if (any(is.na(tabla[[col]])) || any(trims(tabla[[col]]) == "")){
      stop("Hay valores ausentes en ", columna, "del archivo", archivo)
    }
  }

  tabla
}

cat("\n=== Paso 06: enriquecimiento funcional Gene Ontology - Biological Process ===\n")
cat("comprobando los resultados de los pasos 01, 02 y 04...\n")
hpo <- leer_tabla(file.path(opt$input), "hpo_genes.tsv"), c("ncbi_id", "gene_symbol"), "\t")

hpo <- unique(hpo[, c("ncbi_id", "gene_symbol")])
if (anyDuplicated(hpo$gene_symbol)|| anyDuplicated(hpo$ncbi_id)){
  stop("Hay identificadores o símbolos repetidos con asociaciones distintas en hpo_genes.tsv.")
}

# Definición del universo de genes a partir de genes_fondo_hpo.txt

archivo_fondo <- file.path(opt$input, "genes_fondo_hpo.txt")
if (!file.exists(archivo_fondo)){
  stop("El paso 01 debe generar: ", archivo_fondo)
}
fondo <- readLines(archivo_fondo, warn= FALSE)

validos <- !is.na(fondo) & trimws(fondo) != "" & trimws(fondo) != "-"
cat("Entradas sin simbolo valido en el fondo HPO:", sum(!validos), "\n")

fondo <- unique(trimws(fondo[validos]))

if (length(fondo) ==0) {
  stop("El fondo HPO no contiene símbolos de genes válidos.")
}
fuera_fondo <- setdiff(hpo$gene_symbol, fondo)
if (length(fuera_fondo) > 0){
  stop("Hay genes de entrada ausentes del fondo HPO: ", paste(fuera_fondo, collapse = ", "),
  ". Revisa la coherencia entre la API y la release del paso 01")
}

redes <- list()

for(umbral in umbrales){
  carpeta_red <- file.path(carpeta_results, "network", paste0("umbra-_", umbral))
  carpeta_cluster <- file.path(carpeta_results, "clusters", paste0("umbral_", umbral))
  nodos <- leer_tabla(file.path(carpeta_red, "nodos.csv"), "gene_symbol", ";", TRUE)$gene_symbol
    no_mapeados <- leer_tabla(file.path(carpeta_red, "genes_no_mapeados.csv"),
                              "gene_symbol", ";", TRUE)$gene_symbol
    sin_interacciones <- leer_tabla(file.path(carpeta_red, "genes_sin_interacciones.csv"),
                                    "gene_symbol", ";", TRUE)$gene_symbol
    particion <- leer_tabla(file.path(carpeta_cluster, "particion_final.tsv"),
                            c("gene", "cluster"), "\t", TRUE)

    # Los tres estados deben ser excluyentes y cubrir los genes del paso 01.
    estados <- c(nodos, no_mapeados, sin_interacciones)
    if (anyDuplicated(estados) || !setequal(estados, hpo$gene_symbol)) {
        stop("Los genes del paso 02 no concuerdan con el paso 01 para el umbral ", umbral,
             ". Comprueba que los archivos proceden de la misma ejecucion y de la red original.")
    }
    if (anyDuplicated(particion$gene) || !setequal(particion$gene, nodos)) {
        stop("La particion del paso 04 no coincide con los nodos del umbral ", umbral)
    }
    if (any(!grepl("^[1-9][0-9]*$", particion$cluster))) {
        stop("Los identificadores de comunidad deben ser enteros positivos: umbral ", umbral)
    }
    redes[[as.character(umbral)]] <- list(
        nodos = nodos, no_mapeados = no_mapeados,
        sin_interacciones = sin_interacciones, particion = particion
    )
}

dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
if(!dir.exists(opt$output)){
  stop("No se puede crear la carpeta de salida: ", opt$output)
}

cat("Genes HPO de entrada: ", nrow(hpo), "\n")
cat("Genes en el fondo HPO proporcionado: ", length(fondo), "\n")
cat("Se analizarán las redes originales de los umbrales: ", umbrales, "\n")


# ORA. Función para el enriquecimiento funcional 
hacer_ora <- function(genes, universo, carpeta, titulo){
  dir.create(carpeta, recursive= TRUE, showWarnings = FALSE)
  cat("\n", titulo, "\n", sep= "")
  cat("Genes de entrada: ", length(genes), "| Fondo proporcionado: ", length(universo), "\n")

  resultados <- tryCatch(
    enrichGO(
      gene          = genes,        # Genes of interest
      universe      = universo,    # Background: all genes measured in the experiment
      OrgDb         = org.Hs.eg.db,       # Gene annotation database
      keyType       = "SYMBOL",         # Type of gene identifier
      ont           = "BP",               # Ontology (BP: biological processes)
      pAdjustMethod = "BH",               # Method for p-value adjustment
      pvalueCutoff  = 0.05,               # Significance threshold
      readable      = TRUE                # Show gene symbols instead of Entrez IDs
    ),
    error = function(e) {
      stop("Error en ", titulo, ": ", conditionMessage(e))
    }
  )

  tabla_resultado_ORA <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(), 
    BgRatio = character(), 
    pvalue = numeric(),
    p.adjust = numeric(),
    geneID = character(),
    Count = integer())

  estado <- "sin_resultado_ORA"
  if (!is.null(resultado)){
    tabla <- as.data.frame(resultado)
    estado z- if (nrow(tabla) > 0) "con_terminos_significativos" else "sin_terminos_significativos"
  }

  write.csv2(tabla, file.path(carpeta, "ORA_BP.csv"), row.names = FALSE)

  archivo_pdf <- file.path(carpeta, "plot_ORA.pdf")
  if (!file.exists(archivo_pdf) && !file.remove(archivo_pdf)){
    stop("No se puede reemplazar la figura: ", archivo_pdf)
  }
  write.csv2(tabla, file.path(carpeta, "ORA_BP.csv"), row.names = FALSE)

  archivo_pdf <- file.path(carpeta, "plot_ORA.pdf")
  if (!file.exists(archivo_pdf) && !file.remove(archivo_pdf)){
    stop("No se puede reemplazar la figura: ", archivo_pdf)
  }

  tabla_simple <- tabla
  if (is.null(resultado)) {
      cat("No se obtuvo un resultado ORA. Revisa los mensajes sobre anotaciones y conjuntos elegibles.\n")
  } else if (nrow(tabla) == 0) {
      cat("No hay terminos que superen los filtros del ORA. No se genera una figura.\n")
  } else {
      cat("Terminos significativos antes de reducir redundancia:", nrow(tabla), "\n")
      go_simple <- simplify(resultado)
      tabla_simple <- as.data.frame(go_simple)
      cat("Terminos despues de simplify():", nrow(tabla_simple), "\n")
      cat("Primeros terminos ordenados por p ajustado, con sus genes:\n")
      tabla_mostrar <- tabla_simple[order(tabla_simple$p.adjust), ]
      print(head(tabla_mostrar[, c("Description", "Count", "GeneRatio", "p.adjust", "geneID")], 5),
            row.names = FALSE)
      if (nrow(tabla_simple) > 0) {
          pdf(archivo_pdf, width = 9, height = 9)
          tryCatch(
              print(dotplot(go_simple, showCategory = 10, title = titulo)),
              finally = dev.off()
          )
      }
  }

  write.csv2(tabla_simple, file.path(carpeta, "ORA_BP_simplificado.csv"), row.names = FALSE)
  writeLines(c(paste("Analisis:", titulo), paste("Estado:", estado),
                paste("Genes de entrada:", length(genes)),
                paste("Genes del fondo proporcionado:", length(universo)),
                "GO BP; keyType SYMBOL; BH; pvalueCutoff 0.05; qvalueCutoff 0.2"),
              file.path(carpeta, "estado_ORA.txt"))
  list(tabla = tabla, tabla_simple = tabla_simple, estado = estado)
}



#grafica_ORA <- dotplot(go_simple, showCategory = 10, title = "GO Enrichment (BP)")

pdf(file.path(opt$output, "plot_ORA.pdf"), width = 9, height = 9)
dotplot(go_simple, showCategory = 10, title = "GO Enrichment (BP)")
invisible(dev.off())
