#!/usr/bin/env Rscript

lib_r <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib_r)) {
  stop("Falta la carpeta de instalación. Uso: Rscript instalar_paquetes.R ../software/deps_R")
}
dir.create(lib_r, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_r, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 2)

paquetes_cran <- c("optparse", "httr", "jsonlite", "igraph", "BiocManager")
paquetes_bioc <- c("STRINGdb")

# TRUE si el paquete está instalado en nuestra carpeta
instalado <- function(p) {
  p %in% rownames(installed.packages(lib.loc = lib_r))
}

# Paquetes de CRAN

for (p in paquetes_cran) {
  if (instalado(p)) {
    cat(p, "ya instalado\n")
  } else {
    cat("Instalando", p, "\n")
    install.packages(p, lib = lib_r)
  }
}

# install.packages() solo avisa si falla, así que lo comprobamos aquí
for (p in paquetes_cran) {
  if (!instalado(p)) {
    stop("No se ha podido instalar ", p, ". Revisa la conexión a internet ",
         "y que estén instaladas las librerías del sistema del aviso de setup.sh")
  }
}

# Paquetes de Bioconductor
for (p in paquetes_bioc) {
  if (instalado(p)) {
    cat(p, "ya instalado\n")
  } else {
    cat("Instalando", p, "\n")
    BiocManager::install(p, lib = lib_r, update = FALSE, ask = FALSE)
  }
}

for (p in paquetes_bioc) {
  if (!instalado(p)) {
    stop("No se ha podido instalar ", p)
  }
}

cat("Paquetes de R listos\n")