#! /usr/bin/env Rscript
library(optparse)
library(igraph)

option_list <- list(
make_option(c("-i", "--input_file"), type="character", default=NULL, help="Input file with graph edges"),
make_option(c("-g", "--graph_file"), type="character", default=NULL, help="Output file with network graph"),
make_option(c("-d", "--dict_file"), type="character", default=NULL, help="File to translate gene symbol to ensembl protein"),
make_option(c("-G", "--gene_subset"), type="character", default=NULL, help="File with gene list associated to HPO")
)
opt <- parse_args(OptionParser(option_list=option_list))

net_table <- read.table(opt$input_file, sep="\t", header=TRUE)
gene_list <- read.table(opt$gene_subset, sep="\t", header=TRUE)
dict_gene <- read.table(opt$dict_file, sep="\t", header=FALSE)

match_vector <- match(gene_list$name, dict_gene$V2)
gene_list <- dict_gene$V1[match_vector]
gene_list <- gene_list[!is.na(gene_list)]

net <- graph_from_data_frame(d=net_table)
nodes_in_graph <- match(gene_list, V(net)$name)
nodes_in_graph <- nodes_in_graph[!is.na(nodes_in_graph)]
gene_list <- V(net)$name[nodes_in_graph]

net <- subgraph(net, gene_list)
print(net) # primera linea da Num nodos y edges
print(transitivity(net))
pdf(file=opt$graph_file)
plot(net, vertex.color = "orange", vertex.size = 5,vertex.label=NA, edge.width = 0.5, main = "Nicely Layout")
dev.off()
