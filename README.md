# adcp-pipeline

Procesamiento de relevamientos batimétricos con ADCP **SonTek RiverSurveyor
M9/S5** para la Delegación Alto Valle (DPA / AIC), Río Negro, Patagonia
Argentina.

Toma los `.mat` exportados por RiverSurveyor Live y entrega secciones
transversales con **cotas de fondo absolutas** en el datum vertical IGN SRVN16,
CRS POSGAR07 faja 2 (**EPSG:5344**), referidas al sistema de progresivas
oficiales de cada río.

## Estado

| versión | archivo | estado |
|---|---|---|
| v6 | `process_adcp_bathimetric_v6.py` | **producción** |
| v5 | `process_adcp_bathimetric_v5.py` | anterior |
| v4 | `process_adcp_bathimetric_v4.py` | anterior |

## Qué hace

- Extrae del `.mat` el haz vertical, los 4 haces laterales, la trayectoria GNSS,
  el rumbo y las velocidades medias.
- Reconstruye las huellas de fondo de cada haz en `(E, N, profundidad)`.
- Proyecta la nube sobre el eje de la sección y construye un perfil de fondo en
  grilla regular, con rampas de margen desde las distancias declaradas en campo.
- Ubica cada transecta sobre una red multi-río (`RiverNetwork`), le asigna
  progresiva oficial y rotula los brazos alrededor de islas.
- Interpola el pelo de agua desde GNSS y/o estaciones limnimétricas y convierte
  a cotas absolutas.
- **Detecta las repeticiones de un aforo y las promedia en un único perfil**,
  fusionando las nubes de haces sobre un eje común.
- Exporta CSV, PNG y shapefiles por transecta, por grupo y consolidados.

## Uso rápido

```bash
python process_adcp_bathimetric_v6.py ./crudo/adcp \
    --centerline ../ejes/R_LI-NE-NE.shp \
    --water-surface-csv ./crudo/gnss/PA_NXYZD.csv \
    --outdir ./procesado --survey-name relevamiento-20260415
```

En producción conviene un archivo de campaña:

```bash
python process_adcp_bathimetric_v6.py --config campanha.ini
```

Referencia completa de flags, entradas y salidas: **[`docs/USO_v6.md`](docs/USO_v6.md)**.

## Instalación

Requiere Python 3.10 o superior.

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
```

`geopandas` y `shapely` son opcionales en sentido estricto: sin ellos el script
corre igual pero no genera shapefiles ni progresivas.

## Documentación

- [`docs/USO_v6.1.md`](docs/USO_v6.1.md) — uso, flags, entradas y salidas
- [`docs/CAMBIOS_v6_vs_v5.md`](docs/CAMBIOS_v6_vs_v5.md) — cambios de v6
- [`docs/CAMBIOS_v4_vs_v3.md`](docs/CAMBIOS_v4_vs_v3.md) — cambios de v4
- [`docs/USO_v4.md`](docs/USO_v4.md) — referencia histórica

## Referencias

El tratamiento de profundidades sigue el modelo de **USGS QRev**
(Mueller, D.S., 2020, *QRev*, U.S. Geological Survey software release,
<https://doi.org/10.5066/P9OZ8QDL>): referencia primaria seleccionable,
profundidades compuestas, promedio IDW de los 4 haces y filtrado de picos por
haz. QRev calcula caudal; acá el objetivo es batimetría, así que se agrega la
opción de conservar las huellas espaciales de cada haz en vez de colapsarlas.

## Datos de campaña

Los datos crudos **no** se versionan en este repositorio. Se guardan aparte
(`D:\Datos\ADCP\campanhas\`) y se referencian desde el `.ini` de cada campaña.
