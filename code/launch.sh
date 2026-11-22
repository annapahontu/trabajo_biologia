#!/usr/bin/env bash

set -Eeuo pipefail   # parar en cuanto un paso falle

# PARÁMETROS DEL ANÁLISIS

# --- Término HPO ---
HPO_ID="HP:0000073"            
HPO_NOMBRE="Ureteral duplication"
# Release fija de HPO para genes_to_phenotype.txt, para que el resultado
# no cambie cuando salga una versión nueva
HPO_RELEASE="v2026-09-01"

# --- Red expandida del paso 02 ---
# true para construirla, false para saltarla (es la parte que más tarda del paso 02)
EXPANDIR=false

# --- Genes puente del paso 09 ---
# umbral de la red sobre la que se añaden los puentes
UMBRAL_PUENTES=700
# cuántos genes puente se prueban a añadir
N_PUENTES="3 5 10"
# los genes con más interacciones que esto se descartan por ser hubs
MAX_GRADO=150

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
    local inicio=$SECONDS
    shift
    PASO_ACTUAL="$nombre_log"
    "$@" 2>&1 | tee "$LOGS/$nombre_log.log"
    echo "Tiempo de $nombre_log: $((SECONDS - inicio)) s"
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
    --expandir "$EXPANDIR" \
    -o "$RES/network"

# Red completa de STRING, cargada una sola vez para las distancias entre clusters
paso 02b_string_completo Rscript 02b_string_completo.R \
    -u 400 \
    --string_dir "$RES/network/string_data" \
    -o "$RES/network/string_completo.rds"

# Topología, clustering, comunidades de enlaces y distancias entre clusters de la red de cada umbral,
# para ver en cuál se distinguen mejor los clusters
for umbral in 400 700 900; do
    paso "03_topologia_umbral_$umbral" Rscript 03_topologia.R \
        -i "$RES/network/umbral_$umbral/aristas.csv" \
        -g "$RES/network/umbral_$umbral/mapeo_genes.csv" \
        -r "$RES/network/string_completo.rds" \
        -u "$umbral" \
        -m completo \
        -n 1000 \
        -s 123 \
        -o "$RES/topology/umbral_$umbral"

    paso "04_clustering_umbral_$umbral" Rscript 04_clustering.R \
        -i "$RES/network/umbral_$umbral/aristas.csv" \
        -m louvain \
        -r 1 \
        -n 100 \
        -s 123 \
        -o "$RES/clusters/umbral_$umbral"

    paso "05_comunidades_enlaces_umbral_$umbral" Rscript 05_comunidades_enlaces.R \
        -i "$RES/network/umbral_$umbral/aristas.csv" \
        -s 123 \
        --particion "$RES/clusters/umbral_$umbral/particion_final.tsv" \
        --metricas "$RES/topology/umbral_$umbral/metricas_nodos.tsv" \
        -o "$RES/clusters/umbral_$umbral/enlaces"
done

paso 06_ora Rscript 06_ora.R \
    -i "$RES/data" \
    -o "$RES/ora" \
    -u 400,700,900 \
    -m 3 \
    -s 123

for umbral in 400 700 900; do
    paso "07_distancias_clusters_umbral_$umbral" Rscript 07_distancias_clusters.R \
        -p "$RES/clusters/umbral_$umbral/particion_final.tsv" \
        -g "$RES/data/hpo_genes.tsv" \
        -u "$umbral" \
        -r "$RES/network/string_completo.rds" \
        -d 3 \
        -o "$RES/clusters/umbral_$umbral/distancias"
done

# Añadimos unos pocos genes puente a la red, sacamos de nuevo los clusters
# y repetimos el ORA solo con las semillas para ver cómo cambian
for n in $N_PUENTES; do
    paso "09_red_puentes_n_$n" Rscript 09_red_puentes.R \
        -c "$RES/clusters/umbral_$UMBRAL_PUENTES/distancias/caminos_cortos.tsv" \
        -g "$RES/data/hpo_genes.tsv" \
        -r "$RES/network/string_completo.rds" \
        -u "$UMBRAL_PUENTES" \
        -n "$n" \
        --max_grado "$MAX_GRADO" \
        -o "$RES/network/puentes_$UMBRAL_PUENTES/n_$n"

    paso "04_clustering_puentes_n_$n" Rscript 04_clustering.R \
        -i "$RES/network/puentes_$UMBRAL_PUENTES/n_$n/aristas.csv" \
        -m louvain \
        -r 1 \
        -n 100 \
        -s 123 \
        -o "$RES/clusters/puentes_$UMBRAL_PUENTES/n_$n"

    paso "06_ora_puentes_n_$n" Rscript 06_ora.R \
        -i "$RES/data" \
        -o "$RES/ora/puentes_$UMBRAL_PUENTES/n_$n" \
        -m 3 \
        -s 123 \
        --particion "$RES/clusters/puentes_$UMBRAL_PUENTES/n_$n/particion_final.tsv" \
        --puentes "$RES/network/puentes_$UMBRAL_PUENTES/n_$n/genes_puente.tsv"
done

echo ""
echo "=============================================="
echo " Análisis terminado en $((SECONDS / 60)) min $((SECONDS % 60)) s"
echo " Resultados en: $RES"
echo "=============================================="
