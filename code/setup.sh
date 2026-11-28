#!/usr/bin/env bash

set -euo pipefail

# Rutas absolutas, para que funcione, se lance desde donde se lance
CODE_DIR="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$CODE_DIR")"
LIB_R="$RAIZ/software/deps_R"
mkdir -p "$LIB_R"
# Para que R muestre bien las tildes de los mensajes
export LC_ALL=C.UTF-8

echo "=============================================="
echo " Instalando paquetes de R en $LIB_R"
echo "=============================================="

# 1. Comprobar que R está instalado
if ! command -v Rscript > /dev/null; then
    echo "ERROR: no se encuentra R. Instálalo con: sudo apt install r-base r-base-dev"
    exit 1
fi
echo "$(R --version | head -n 1)"

# 2. Avisar si el R es antiguo
# Con R anterior a 4.5, BiocManager usa una versión antigua de Bioconductor que no es compatible con los paquetes actuales de CRAN y falla el ORA
version_r=$(Rscript -e 'cat(as.character(getRversion()))')
if [ "$(printf '%s\n' "4.5.0" "$version_r" | sort -V | head -n 1)" != "4.5.0" ]; then
    echo ""
    echo "AVISO: tienes R $version_r. Los paquetes del ORA (paso 06) necesitan R 4.5 o posterior."
    echo "Instrucciones para actualizar: https://cloud.r-project.org/bin/linux/ubuntu/"
    echo ""
fi

# 3. Avisar si faltan librerías del sistema para compilar paquetes de R
# Solo hacen falta si algún paquete precompilado no carga y hay que compilarlo
librerias="build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev libuv1-dev libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev libpng-dev libjpeg-dev libtiff5-dev"
faltan=""
for lib in $librerias; do
    if ! dpkg -s "$lib" > /dev/null 2>&1; then
        faltan="$faltan $lib"
    fi
done
if [ -n "$faltan" ]; then
    echo ""
    echo "AVISO: faltan estas librerías del sistema:$faltan"
    echo "Solo hacen falta si algún paquete precompilado no carga y hay que compilarlo."
    echo "En ese caso, el administrador del equipo puede instalarlas con este comando"
    echo "  sudo apt install -y$faltan"
    echo ""
fi

# 4. Paquetes de R
export R_LIBS_USER="$LIB_R"
Rscript "$CODE_DIR/instalar_paquetes.R" "$LIB_R"


echo ""
echo "=============================================="
echo " Instalación terminada. Ejecuta ./launch.sh"
echo "=============================================="