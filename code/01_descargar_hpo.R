#!/usr/bin/env Rscript
library(optparse)
library(httr)
library(jsonlite)


# Argumentos
option_list <- list(
  make_option(c("--hpo"), type = "character", default = NULL,
              help = "Código del término HPO, por ejemplo HP:0000073"),
  make_option(c("--release"), type = "character", default = NULL,
              help = "Release de HPO, por ejemplo v2026-09-01"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Carpeta de salida"),
  make_option(c("--api"), type = "character",
              default = "https://ontology.jax.org/api/network/annotation/",
              help = "URL base de la API de HPO")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$hpo) || is.null(opt$release) || is.null(opt$output)) {
  stop("Faltan argumentos. Uso: Rscript 01_descargar_hpo.R --hpo HP:... --release v... -o carpeta")
}
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# Descarga una URL. Si el servidor falla, lo vuelve a intentar hasta 3 veces
descargar <- function(url) {
  for (intento in 1:3) {
    response <- tryCatch(GET(url, timeout(120)), error = function(e) NULL)
    if (!is.null(response) && status_code(response) == 200) {
      return(response)
    }
    cat("  Intento", intento, "fallido\n")
    Sys.sleep(5)
  }
  stop("No se ha podido descargar ", url)
}


# 1. Genes y enfermedades del término (API de HPO)
url_api <- paste0(opt$api, opt$hpo)
cat("Descargando genes y enfermedades de", opt$hpo, "\n")
cat(" ", url_api, "\n")

response <- descargar(url_api)
data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))

# Quitamos el prefijo "NCBIGene:" para quedarnos con el número
genes <- data.frame(ncbi_id = sub("NCBIGene:", "", data$genes$id),
                    gene_symbol = data$genes$name)
genes <- unique(genes)
genes <- genes[order(genes$gene_symbol), ]

if (nrow(genes) == 0) {
  stop("La API no ha devuelto ningún gen para ", opt$hpo)
}

# Enfermedades
diseases <- data.frame(disease_id = data$diseases$id,
                       disease_name = data$diseases$name,
                       mondo_id = data$diseases$mondoId)
diseases <- diseases[order(diseases$disease_id), ]

write.table(genes, file.path(opt$output, "hpo_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diseases, file.path(opt$output, "hpo_diseases.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  Genes:", nrow(genes), "\n")
cat("  Enfermedades:", nrow(diseases),
    "(", length(unique(diseases$mondo_id)), "distintas según MONDO )\n")


# 2. genes_to_phenotype.txt de una release fija
url_g2p <- paste0("https://github.com/obophenotype/human-phenotype-ontology/releases/download/",
                  opt$release, "/genes_to_phenotype.txt")
archivo_g2p <- file.path(opt$output, paste0("genes_to_phenotype_", opt$release, ".txt"))

if (file.exists(archivo_g2p)) {
  cat("\ngenes_to_phenotype.txt (", opt$release, ") ya descargado\n")
} else {
  cat("\nDescargando genes_to_phenotype.txt (", opt$release, ")\n")
  cat(" ", url_g2p, "\n")
  response <- descargar(url_g2p)
  writeBin(content(response, as = "raw"), archivo_g2p)
}

g2p <- read.delim(archivo_g2p, colClasses = "character", quote = "")

# Todos los genes con alguna anotación en HPO
fondo <- sort(unique(g2p$gene_symbol))
writeLines(fondo, file.path(opt$output, "genes_fondo_hpo.txt"))
cat("  Genes con alguna anotación en HPO (fondo):", length(fondo), "\n")


# 3. Anotaciones directas del término
# Nos dicen qué gen va con qué enfermedad, y la
# frecuencia con la que aparece el fenotipo en cada enfermedad
anotaciones <- g2p[g2p$hpo_id == opt$hpo,
                   c("ncbi_gene_id", "gene_symbol", "disease_id", "frequency")]
colnames(anotaciones)[1] <- "ncbi_id"
anotaciones <- anotaciones[order(anotaciones$gene_symbol, anotaciones$disease_id), ]

write.table(anotaciones, file.path(opt$output, "hpo_anotaciones_termino.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAnotaciones directas de", opt$hpo, "en genes_to_phenotype:", nrow(anotaciones), "\n")

# Comparamos con la API, ya que la API también incluye los genes de los términos
# hijos, por eso puede devolver algún gen más
genes_directos <- unique(anotaciones$gene_symbol)
solo_api <- setdiff(genes$gene_symbol, genes_directos)
solo_g2p <- setdiff(genes_directos, genes$gene_symbol)
cat("\nGenes de la API:", nrow(genes), "  Genes con anotación directa:", length(genes_directos), "\n")
if (length(solo_api) == 0) solo_api <- "ninguno"
if (length(solo_g2p) == 0) solo_g2p <- "ninguno"
cat("  Solo en la API (vienen de términos hijos):", solo_api, "\n")
cat("  Solo en genes_to_phenotype:", solo_g2p, "\n")

# 3b. Evidencia de cada gen
# HPO indica con qué frecuencia aparece el fenotipo en cada enfermedad, con una
# categoría (HP:0040280 a HP:0040284) o con una fracción de pacientes (por ejemplo 1/13).
# Para cada gen nos quedamos con la frecuencia más alta entre sus enfermedades
categorias <- c("HP:0040280" = "Obligatoria", "HP:0040281" = "Muy frecuente",
                "HP:0040282" = "Frecuente", "HP:0040283" = "Ocasional",
                "HP:0040284" = "Muy rara")
orden <- c("Obligatoria", "Muy frecuente", "Frecuente", "Ocasional", "Muy rara",
           "Pocos casos", "Sin dato", "Sin anotación directa")

# Las fracciones se pasan a la categoría de HPO que les corresponde por porcentaje,
# salvo las de menos de 5 pacientes (1/1, 1/2...), que son casos sueltos
a_categoria <- function(frecuencia) {
  if (frecuencia %in% names(categorias)) return(categorias[[frecuencia]])
  if (grepl("^[0-9]+/[0-9]+$", frecuencia)) {
    partes <- as.numeric(strsplit(frecuencia, "/")[[1]])
    if (partes[2] < 5) return("Pocos casos")
    porcentaje <- 100 * partes[1] / partes[2]
    if (porcentaje == 100) return("Obligatoria")
    if (porcentaje >= 80) return("Muy frecuente")
    if (porcentaje >= 30) return("Frecuente")
    if (porcentaje >= 5) return("Ocasional")
    return("Muy rara")
  }
  "Sin dato"
}
anotaciones$categoria <- sapply(anotaciones$frequency, a_categoria)

evidencia <- data.frame(gene_symbol = genes$gene_symbol)
evidencia$n_enfermedades <- sapply(evidencia$gene_symbol, function(gen) {
  length(unique(anotaciones$disease_id[anotaciones$gene_symbol == gen]))
})
evidencia$frecuencias <- sapply(evidencia$gene_symbol, function(gen) {
  frec <- unique(anotaciones$frequency[anotaciones$gene_symbol == gen])
  if (length(frec) == 0) return("-")
  paste(frec, collapse = ",")
})
# los genes que solo vienen de términos hijos no tienen anotación directa
evidencia$categoria <- sapply(evidencia$gene_symbol, function(gen) {
  cat_gen <- anotaciones$categoria[anotaciones$gene_symbol == gen]
  if (length(cat_gen) == 0) return("Sin anotación directa")
  orden[min(match(cat_gen, orden))]
})

write.table(evidencia, file.path(opt$output, "hpo_evidencia_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nEvidencia de los genes del término (frecuencia más alta de cada gen):\n")
print(table(factor(evidencia$categoria, levels = orden)))

# 4. Resumen de la descarga
resumen <- c(
  paste("Fecha de descarga:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  paste("Término HPO:", opt$hpo),
  paste("URL de la API:", url_api),
  paste("Release de HPO (genes_to_phenotype):", opt$release),
  paste("URL de genes_to_phenotype:", url_g2p),
  paste("Genes (API):", nrow(genes)),
  paste("Enfermedades (API):", nrow(diseases)),
  paste("Genes con anotación directa:", length(genes_directos)),
  paste("Genes solo en la API:", paste(solo_api, collapse = ", ")),
  paste("Genes en el fondo HPO:", length(fondo))
)
writeLines(resumen, file.path(opt$output, "resumen_descarga.txt"))

cat("\nArchivos guardados en", opt$output, "\n")