# Análisis de fenotipo basado en HPO y redes STRING

Pipeline para el análisis de fenotipos utilizando Human Phenotype Ontology (HPO), redes de interacción proteína-proteína de STRING, análisis topológico, detección de comunidades y enriquecimiento funcional (ORA).

## Requisitos

- Ubuntu 22.04 o superior
- R ≥ 4.5
- Conexión a Internet

### Dependencias del sistema

```bash
sudo apt update

sudo apt install -y \
build-essential \
gfortran \
libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
libglpk-dev \
libuv1-dev \
libfontconfig1-dev \
libfreetype-dev \
libharfbuzz-dev \
libfribidi-dev \
libpng-dev \
libjpeg-dev \
libtiff5-dev
```

Comprobar la versión de R:

```bash
R --version
```

## Instalación

Desde la carpeta raíz del proyecto:

```bash
chmod +x code/setup.sh
./code/setup.sh
```

Las dependencias de R se instalarán automáticamente en:

```text
software/deps_R/
```

## Ejecución

Ejecutar la pipeline completa:

```bash
chmod +x launch.sh
./launch.sh
```

Todos los parámetros del análisis (término HPO, release de HPO, umbrales de STRING, etc.) se configuran en `launch.sh`.

## Resultados

Todos los resultados se guardan dentro de:

```text
results/
```

Estructura general:

```text
results/
├── data/
├── network/
├── topology/
├── clusters/
├── ora/
├── figuras/
├── comparacion_puentes/
└── logs/
```

### Directorios principales

| Carpeta | Contenido |
|----------|-----------|
| `results/data` | Genes y enfermedades asociadas al término HPO |
| `results/network` | Redes STRING y análisis de sensibilidad |
| `results/topology` | Métricas topológicas y modelos nulos |
| `results/clusters` | Particiones y comunidades |
| `results/ora` | Resultados del enriquecimiento funcional |
| `results/figuras` | Figuras resumen |
| `results/comparacion_puentes` | Comparación entre redes con distintos números de genes puente |
| `results/logs` | Logs de ejecución |

## Reproducibilidad

Todos los scripts utilizan una semilla aleatoria por defecto (`123`) para facilitar la reproducibilidad de los resultados.
