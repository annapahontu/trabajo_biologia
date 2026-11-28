# Análisis de redes para el fenotipo HP:0000073 (duplicación ureteral)

Este repositorio contiene el análisis en R del término HPO HP:0000073 (Ureteral duplication). A partir de los genes asociados al término en la Human Phenotype Ontology se construye su red de interacciones en STRING y se analiza en varias fases.

- Topología de la red y comparación con modelos nulos.
- Búsqueda de módulos con Louvain y con comunidades de enlaces.
- Enriquecimiento funcional (ORA) de Gene Ontology de cada módulo.
- Adición de unos pocos genes intermedios ("genes puente") para ver si unen grupos de genes que en la red original aparecen separados.

Todo se ejecuta con dos scripts, `setup.sh` y `launch.sh`, sin tocar nada a mano.

## Requisitos

- Linux. Se ha probado en Ubuntu, y también funciona en Windows con WSL.
- R 4.5 o posterior, con sus herramientas de compilación (en Ubuntu o Debian, `sudo apt install r-base r-base-dev`). Con versiones anteriores, los paquetes del ORA fallan.
- Conexión a internet, porque los paquetes y los datos se descargan durante la instalación y la ejecución.

Con R instalado no hace falta nada más ni permisos de administrador. En Ubuntu, Debian y las distribuciones basadas en ellas (Linux Mint, Pop!_OS...), los paquetes de CRAN se descargan ya compilados de Posit Package Manager, así que no hay que instalar librerías de desarrollo del sistema. En otras distribuciones se compilan desde CRAN, y entonces sí hacen falta esas librerías.

## Instalación y ejecución paso a paso

Todos los comandos se escriben en una terminal de Linux (en Windows, en la terminal de Ubuntu de WSL).

### 1. Comprobar la versión de R

Tiene que ser la 4.5 o posterior.

```bash
R --version
```

### 2. Clonar el repositorio

```bash
git clone https://github.com/annapahontu/trabajo_biologia.git
cd trabajo_biologia
```

### 3. Instalar las dependencias

Solo hace falta la primera vez.

```bash
bash code/setup.sh
```

El script instala todos los paquetes de R en `software/deps_R/`, dentro del propio proyecto, y comprueba que cada uno carga. Si un paquete precompilado no carga, lo vuelve a instalar compilándolo desde CRAN. Termina con el mensaje `Paquetes de R listos`. Si se corta a medias, se puede volver a lanzar, porque salta los paquetes ya instalados.

### 4. Ejecutar el análisis

```bash
bash code/launch.sh
```

`launch.sh` ejecuta todos los pasos en orden, muestra por pantalla lo que va haciendo y guarda la salida de cada paso en `results/logs/`. Si un paso falla, se detiene y dice cuál ha sido.

### 5. Revisar los resultados

Al terminar aparece el mensaje `Análisis terminado` con el tiempo total, y todo queda en `results/`. 

## Parámetros

Los parámetros principales están en `launch.sh`.

| Parámetro | Valor | Qué controla |
|---|---|---|
| `HPO_ID` | `HP:0000073` | Término HPO analizado |
| `HPO_RELEASE` | `v2026-09-01` | Release fija de `genes_to_phenotype.txt` de HPO |
| `UMBRALES` | `400 500 600 700 800` | Umbrales de confianza de STRING que se exploran |
| `UMBRAL_PUENTES` | `600` | Umbral de la red a la que se añaden los genes puente |
| `N_PUENTES` | `3 5 10` | Número de genes puente que se prueban |
| `MAX_GRADO` | `150` | Grado a partir del cual un gen se descarta como puente por ser hub |

Todos los pasos aleatorios usan la semilla 123.

## Pasos del análisis

| Script | Qué hace |
|---|---|
| `01_descargar_hpo.R` | Descarga los genes y enfermedades del término y las anotaciones de HPO |
| `02_red_string.R` | Construye la red de STRING para cada umbral |
| `02b_string_completo.R` | Guarda la red completa de STRING, que se usa como referencia |
| `03_topologia.R` | Métricas de la red, modelos nulos, test de módulo de enfermedad y robustez |
| `04_clustering.R` | Clustering con varios métodos y estabilidad de Louvain |
| `05_comunidades_enlaces.R` | Comunidades de enlaces |
| `06_ora.R` | Enriquecimiento funcional (GO) de todos los genes y de cada cluster |
| `07_distancias_clusters.R` | Distancias y caminos más cortos entre clusters en STRING |
| `08_figuras_resumen.R` | Figuras de resumen (embudo de genes, sensibilidad al umbral, grado frente a betweenness) |
| `09_red_puentes.R` | Añade a la red los genes puente |
| `10_comparar_puentes.R` | Compara la red original con las redes con puentes |

El lanzador primero explora todos los umbrales y después hace el análisis completo con 600 como umbral principal y 400 como referencia. A 400, Louvain agrupa en un solo cluster los 17 genes del bloque denso, mientras que a 600 esos genes se reparten en dos grupos de 11 y 6, una partición que se repite en las 100 ejecuciones de Louvain y coincide con la de Leiden.

Los genes puente se eligen entre los genes que, en STRING, quedan justo en medio de dos grupos de la red (a distancia 2). Se descartan los hubs y se añaden los que unen más grupos. Con cada red con puentes se repiten el clustering y el ORA, quitando antes los puentes para que el resultado dependa solo de los genes del fenotipo.

## Carpeta de resultados

Todo lo que genera el análisis se guarda en `results/`, que no se sube al repositorio.

```text
results/
├── data/                  genes, enfermedades y anotaciones de HPO
├── network/               redes de STRING (una por umbral y una por número de puentes)
├── topology/              métricas, modelos nulos y robustez
├── clusters/              clustering, comunidades de enlaces y distancias
├── ora/                   enriquecimiento funcional
├── comparacion_puentes/   comparación entre la red original y las redes con puentes
├── figuras/               figuras de resumen
├── logs/                  salida de cada paso
└── session_info.txt       versiones de R y de los paquetes
```

Dentro de `network/`, `topology/`, `clusters/` y `ora/` hay una subcarpeta `umbral_<u>/` por umbral y `puentes_600/n_<N>/` para las redes con puentes.

## Ejecutar un solo paso

Cada script acepta sus parámetros por línea de comandos (`--help` muestra cuáles son). Para lanzar uno suelto hay que estar en `code/` e indicar dónde están los paquetes. Los comandos exactos de cada paso están en `launch.sh`.

```bash
cd code
export R_LIBS_USER=../software/deps_R
export LC_ALL=C.UTF-8
RES=../results

Rscript 04_clustering.R -i $RES/network/umbral_600/aristas.csv -m louvain -r 1 -n 100 -s 123 -o $RES/clusters/umbral_600
```

Cada paso lee lo que generan los anteriores, así que solo se puede lanzar suelto si esos resultados ya existen.

## Problemas frecuentes

Si un paso falla, el error aparece al final de su log, en `results/logs/<paso>.log`.

| Problema | Solución |
|---|---|
| `ERROR: no están instaladas las dependencias` | Ejecutar `bash code/setup.sh` antes que el lanzador |
| `/usr/bin/env: 'bash\r'` o `$'\r': command not found` | Los scripts tienen saltos de línea de Windows, algo que pasa si el repositorio se clonó desde Windows. Se arregla con `sed -i 's/\r$//' code/*.sh` |
| `Permission denied` al lanzar `./code/launch.sh` | Lanzarlo con `bash code/launch.sh` |
| `setup.sh` avisa de librerías del sistema que faltan | No es un problema si la instalación termina con `Paquetes de R listos`. La lista de librerías y el comando son de Ubuntu y Debian |
| Un paquete falla al compilar después del mensaje `no carga con la versión precompilada` | Al sistema le falta alguna librería. Hay que instalarla con el comando que imprime `setup.sh` (requiere permisos de administrador) y volver a lanzarlo |
| El paso 06 falla al cargar clusterProfiler, enrichplot o DOSE | Versión de R anterior a 4.5. Hay que actualizar R, borrar `software/deps_R/` y volver a ejecutar `setup.sh` |
| `No se ha podido descargar ...` o errores de conexión | El servidor de HPO, GitHub o STRING no responde. Suele ser temporal, basta con volver a lanzar |
| `Hay genes de entrada ausentes del fondo HPO` en el paso 06 | La API de HPO ha añadido un gen que no está en la release fija. Hay que cambiar `HPO_RELEASE` por una más reciente |
| Avisos de `Arguments in ... must be used` con `by = "Count"` | No es un error, se puede ignorar |
| El paso 03 tarda mucho | Es normal, porque genera 1000 redes aleatorias. Para hacer pruebas se puede bajar con `-n 100` |
| Un paso no encuentra un archivo de otro paso | Se ha lanzado un paso suelto antes de tiempo o una ejecución anterior se cortó. Hay que volver a lanzar `launch.sh` entero |
