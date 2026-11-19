#!/usr/bin/env Rscript

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

set.seed(opt$semilla) # Establece la semilla para sample_gnm y rewire
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE) # Directorio de salida

# Obtención de genes para el universo del ORA del paquete DOSE de clusterProfiler
data(geneList, package = "DOSE")
universe_genes <- names(geneList)

# Genes de interés: HPO
hpo_genes = file.path(opt$input, "hpo_genes.tsv")
hpo <- read.delim(hpo_genes, stringsAsFactors = FALSE)


go_enrich <- enrichGO(
  gene          = hpo$ncbi_id,        # Genes of interest
  universe      = names(geneList),    # Background: all genes measured in the experiment
  OrgDb         = org.Hs.eg.db,       # Gene annotation database
  keyType       = "ENTREZID",         # Type of gene identifier
  ont           = "BP",               # Ontology (BP: biological processes)
  pAdjustMethod = "BH",               # Method for p-value adjustment
  pvalueCutoff  = 0.05,               # Significance threshold
  readable      = TRUE                # Show gene symbols instead of Entrez IDs
)

tabla_resultado_ORA <-as.data.frame(go_enrich)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")]

go_simple <- simplify(go_enrich)

#grafica_ORA <- dotplot(go_simple, showCategory = 10, title = "GO Enrichment (BP)")

pdf(file.path(opt$output, "plot_ORA.pdf"), width = 9, height = 9)
dotplot(go_simple, showCategory = 10, title = "GO Enrichment (BP)")
invisible(dev.off())
