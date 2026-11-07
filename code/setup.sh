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

# 2. Avisar si faltan herramientas del sistema que necesitan algunos
# paquetes de R para compilar 
faltan=""
for programa in gcc g++ gfortran make; do
    if ! command -v "$programa" > /dev/null; then
        faltan="$faltan $programa"
    fi
done
if [ ! -f /usr/include/curl/curl.h ] && [ ! -f /usr/include/x86_64-linux-gnu/curl/curl.h ]; then
    faltan="$faltan libcurl(dev)"
fi
if [ ! -f /usr/include/png.h ] && [ ! -f /usr/include/libpng16/png.h ]; then
    faltan="$faltan libpng(dev)"
fi
if [ -n "$faltan" ]; then
    echo ""
    echo "AVISO: puede que falten estas herramientas del sistema:$faltan"
    echo "Si la instalación falla, ejecuta (una sola vez):"
    echo "  sudo apt install r-base-dev libcurl4-openssl-dev libssl-dev libpng-dev"
    echo ""
fi

# 3. Paquetes de R
Rscript "$CODE_DIR/instalar_paquetes.R" "$LIB_R"

echo ""
echo "=============================================="
echo " Instalación terminada. Ejecuta ./launch.sh"
echo "=============================================="