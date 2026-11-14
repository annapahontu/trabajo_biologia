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

# 2. Paquetes de R
Rscript "$CODE_DIR/instalar_paquetes.R" "$LIB_R"

echo ""
echo "=============================================="
echo " Instalación terminada. Ejecuta ./launch.sh"
echo "=============================================="