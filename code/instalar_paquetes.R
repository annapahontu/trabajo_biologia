#!/usr/bin/env Rscript

lib_r <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib_r)) {
  stop("Falta la carpeta de instalación. Uso: Rscript instalar_paquetes.R ../software/deps_R")
}
dir.create(lib_r, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_r, .libPaths()))

# Desactivar bspm en esta sesión
if (requireNamespace("bspm", quietly = TRUE)) {
  bspm::disable()
}

# Repositorios
# En Ubuntu, Debian y sus derivadas, los paquetes de CRAN se descargan ya compilados de
# Posit Package Manager, así no hay que instalar con sudo las librerías para compilarlos
repo_fuente <- "https://cloud.r-project.org"

# campo de /etc/os-release, o "" si no está
campo_os <- function(lineas, campo) {
  linea <- grep(paste0("^", campo, "="), lineas, value = TRUE)
  if (length(linea) == 1) gsub('"', "", sub(paste0("^", campo, "="), "", linea)) else ""
}

codename <- ""
if (file.exists("/etc/os-release")) {
  os <- readLines("/etc/os-release")
  # las distribuciones basadas en Ubuntu (Mint, Pop!_OS...) indican aquí el Ubuntu del que vienen
  codename <- campo_os(os, "UBUNTU_CODENAME")
  if (!nzchar(codename) && campo_os(os, "ID") %in% c("ubuntu", "debian")) {
    codename <- campo_os(os, "VERSION_CODENAME")
  }
}

repo_cran <- repo_fuente
if (nzchar(codename)) {
  repo_posit <- paste0("https://packagemanager.posit.co/cran/__linux__/", codename, "/latest")
  # Posit solo da los paquetes compilados si sabe la versión de R y el sistema
  options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                  paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
  # comprobamos que Posit tiene paquetes para este sistema antes de usarlo
  disponibles <- tryCatch(available.packages(repos = repo_posit), error = function(e) NULL)
  if (!is.null(disponibles) && nrow(disponibles) > 0) repo_cran <- repo_posit
}

if (repo_cran != repo_fuente) {
  cat("Paquetes de CRAN precompilados para", codename, "desde Posit Package Manager\n")
} else {
  cat("No hay paquetes precompilados para este sistema, los de CRAN se compilarán desde el código fuente\n")
}
options(repos = c(CRAN = repo_cran), Ncpus = 2)

paquetes_cran <- c("optparse", "httr", "jsonlite", "igraph", "BiocManager")
paquetes_bioc <- c("STRINGdb")
paquetes_ORA <- c("fgsea", "treeio", "DOSE", "ggtree","clusterProfiler", "org.Hs.eg.db", "enrichplot")

# TRUE si el paquete está instalado en nuestra carpeta
instalado <- function(p) {
  p %in% rownames(installed.packages(lib.loc = lib_r))
}

# TRUE si el paquete carga. Uno precompilado puede instalarse y no cargar si al sistema le falta alguna librería
funciona <- function(p) {
  suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
}

# Si un paquete precompilado no carga, se compila desde CRAN
compilar_si_falla <- function(p) {
  if (instalado(p) && !funciona(p) && repo_cran != repo_fuente) {
    cat(p, "no carga con la versión precompilada, se compila desde el código fuente\n")
    install.packages(p, lib = lib_r, repos = c(CRAN = repo_fuente), type = "source")
  }
}

# Cuando un paquete de Bioconductor no carga, normalmente es por una dependencia
# de CRAN precompilada, así que se compilan las dependencias que no cargan
reparar_dependencias <- function(p) {
  if (instalado(p) && !funciona(p)) {
    instalados <- installed.packages(lib.loc = lib_r)
    deps <- tools::package_dependencies(p, db = instalados, recursive = TRUE)[[1]]
    for (d in intersect(deps, rownames(instalados))) compilar_si_falla(d)
  }
}

# Paquetes de CRAN
for (p in paquetes_cran) {
  if (instalado(p)) {
    cat(p, "ya instalado\n")
  } else {
    cat("Instalando", p, "\n")
    install.packages(p, lib = lib_r)
  }
  compilar_si_falla(p)
}

# install.packages() solo avisa si falla, así que lo comprobamos aquí
for (p in paquetes_cran) {
  if (!instalado(p) || !funciona(p)) {
    stop("No se ha podido instalar ", p, ". Revisa la conexión a internet")
  }
}

# Paquetes de Bioconductor
# BiocManager instala las dependencias de CRAN desde el mismo repositorio, así que también llegan compiladas
for (p in paquetes_bioc) {
  if (instalado(p)) {
    cat(p, "ya instalado\n")
  } else {
    cat("Instalando", p, "\n")
    BiocManager::install(p, lib = lib_r, update = FALSE, ask = FALSE)
  }
  reparar_dependencias(p)
}

for (p in paquetes_bioc) {
  if (!instalado(p) || !funciona(p)) {
    stop("No se ha podido instalar ", p)
  }
}

# Paquetes para el análisis de ORA
for (p in paquetes_ORA) {
  if (instalado(p)) {
    cat(p, "ya instalado\n")
  } else {
    cat("Instalando", p, "\n") 
    BiocManager::install(
      p,
      lib = lib_r,
      update = FALSE,
      ask = FALSE,
      force = TRUE
    )
  } 
  reparar_dependencias(p)
}

for (p in paquetes_ORA) {
  if (!instalado(p) || !funciona(p)) {
    stop("No se ha podido instalar ", p)
  }
}

cat("Paquetes de R listos\n")
