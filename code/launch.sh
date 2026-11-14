#!/usr/bin/env bash

set -Eeuo pipefail   # parar en cuanto un paso falle

# PARÁMETROS DEL ANÁLISIS

# --- Término HPO ---
HPO_ID="HP:0000073"            
HPO_NOMBRE="Ureteral duplication"
# Release fija de HPO para genes_to_phenotype.txt, para que el resultado
# no cambie cuando salga una versión nueva
HPO_RELEASE="v2026-09-01"

# Rutas
CODE_DIR="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$CODE_DIR")"
RES="$RAIZ/results"
cd "$CODE_DIR"

if [ ! -d "$RAIZ/software/deps_R" ]; then
    echo "ERROR: no están instaladas las dependencias. Ejecuta primero ./setup.sh"
    exit 1
fi

# Paquetes de R instalados por setup.sh
export R_LIBS_USER="$RAIZ/software/deps_R"
# Para que R lea bien las tildes de los scripts
export LC_ALL=C.UTF-8

mkdir -p "$RES"/{data,logs}
LOGS="$RES/logs"

# Funciones auxiliares
PASO_ACTUAL="inicio"

# Si algo falla, decimos en qué paso y dónde está su log
trap 'echo ""; echo "ERROR en el paso: $PASO_ACTUAL"; echo "Mira el log en results/logs/"; exit 1' ERR

# paso NOMBRE_LOG comando argumentos...
# Ejecuta el comando, lo muestra en pantalla y lo guarda en el log
paso() {
    local nombre_log="$1"
    shift
    PASO_ACTUAL="$nombre_log"
    "$@" 2>&1 | tee "$LOGS/$nombre_log.log"
}

echo "=============================================="
echo " Análisis de $HPO_ID ($HPO_NOMBRE)"
echo "=============================================="

# Genes y enfermedades del término HPO
echo ""
echo "Descargando genes de HPO"
paso 01_descargar_hpo Rscript 01_descargar_hpo.R \
    --hpo "$HPO_ID" \
    --release "$HPO_RELEASE" \
    -o "$RES/data"
    
# Construir las redes de STRING y comparar umbrales
paso 02_red_string Rscript 02_red_string.R \
    -i "$RES/data/hpo_genes.tsv" \
    -o "$RES/network"

# Analizar topología de la red (de umbral 700)
paso 03_topologia Rscript 03_topologia.R \
    -i "$RES/network/umbral_700/aristas.csv" \
    -m completo \
    -n 1000 \
    -s 123 \
    -o "$RES/topology"

# Detectar comunidades en la red (de umbral 700)
paso 04_clustering Rscript 04_clustering.R \
    -i "$RES/network/umbral_700/aristas.csv" \
    -m louvain \
    -r 1 \
    -n 100 \
    -s 123 \
    -o "$RES/clusters"

echo ""
echo "=============================================="
echo " Análisis terminado en $((SECONDS / 60)) min $((SECONDS % 60)) s"
echo " Resultados en: $RES"
echo "=============================================="