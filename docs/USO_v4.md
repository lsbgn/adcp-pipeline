# Procesador batimétrico ADCP v4 — Uso, entradas y salidas

`process_adcp_bathimetric_v4.py` — extracción de secciones batimétricas 1D
perpendiculares al flujo desde archivos `.mat` de SonTek M9/S5 (RiverSurveyor
Live). No exporta velocidades: sólo el lecho y los productos de control.

---

## Requisitos

Python 3 con: `numpy`, `scipy`, `matplotlib`, `pyproj`, y para shapefiles /
progresivas: `geopandas` + `shapely`. Sin `geopandas` corre igual pero omite
shapefiles y el cálculo de progresivas.

---

## Uso

### Una transecta

```bash
python process_adcp_bathimetric_v4.py transecta.mat \
    --centerline r_nqn_inf.shp \
    --water-surface-csv PA_NXYZD.csv
```

### Una salida completa (todas las transectas)

```bash
python process_adcp_bathimetric_v4.py carpeta_salida/ \
    --centerline r_nqn_inf.shp \
    --water-surface-csv PA_NXYZD.csv \
    --survey-name salida_2026-04-15
```

`matfiles` acepta un archivo, **varios archivos**, una **carpeta** (toma todos los
`*.mat`) o un **patrón glob**. Con 2 o más transectas se activa el modo
relevamiento. `--centerline` y `--water-surface-csv` son opcionales pero se
recomiendan (dan progresiva/brazo y cotas absolutas IGN SRVN16).

### Con archivo de configuración (recomendado por campaña)

En vez de tipear todas las opciones, se puede usar un `campania.ini`:

```bash
python process_adcp_bathimetric_v4.py --config campania.ini
```

El archivo tiene dos secciones: `[procesamiento]` (entradas y parámetros) y
`[campania]` (metadatos de trazabilidad). Las **rutas relativas se resuelven
respecto de la ubicación del `.ini`**, así que conviene dejarlo dentro de la
carpeta de la campaña. **La línea de comandos tiene prioridad** sobre el archivo
(por ejemplo `--config campania.ini --dx 1.0` usa dx=1.0). Ver la plantilla
`campania_ejemplo.ini`. Cada corrida deja un registro `procesamiento.txt` con la
versión del script, los parámetros y los metadatos usados.

---

## Entradas

| Entrada | Formato | Obligatoria | Descripción |
|---|---|---|---|
| Transecta(s) | `.mat` SonTek M9/S5 | Sí | Uno o varios / carpeta / glob. |
| Centerline | `.shp` (u otro vector) | No | Eje del río; puede tener brazos (isla). Da progresiva + brazo. |
| Pelo de agua | `.csv` (X, Y, cota) | No | Puntos de superficie del agua. Da cotas absolutas en el datum. |

**CSV de pelo de agua:** columnas por defecto `East`, `North`, `H_correg`
(configurables). CRS por defecto **EPSG:5344** (POSGAR07 faja 2). Los puntos se
proyectan a progresiva y se interpola la cota; **cubre sólo el tramo relevado**,
fuera de él se extrapola con aviso.

### Opciones principales

| Opción | Def. | Qué hace |
|---|---|---|
| `--config` | — | Archivo INI (campania.ini) con entradas, parámetros y metadatos. |
| `--outdir` | — | Carpeta de salida. |
| `--survey-name` | `salida_<perfil>` | Nombre de la carpeta del relevamiento. |
| `--dx` | `0.5` | Paso de grilla del perfil [m]. |
| `--bin-half-width` | `2·dx` | Semiventana de muestreo en s [m]. |
| `--centerline` | — | Vector del eje del río (progresivas + brazos). |
| `--chainage-offset` | `0.0` | Progresiva base sumada (origen oficial del eje). |
| `--chainage-reverse` | off | Invierte el sentido de progresivas (se ignora si hay pelo de agua). |
| `--water-surface-csv` | — | CSV de pelo de agua → cotas absolutas. |
| `--ws-east-col` / `--ws-north-col` / `--ws-elev-col` | `East` / `North` / `H_correg` | Nombres de columnas del CSV. |
| `--ws-crs` | `EPSG:5344` | CRS de las coordenadas del CSV. |
| `--datum-name` | `IGN SRVN16` | Etiqueta del datum vertical. |
| `--water-surface-elev` | `0.0` | Cota constante de pelo de agua si NO se da CSV. |
| `--blend-beams` | off | Vuelve al lecho mezclado VB+laterales (comportamiento v3). |
| `--vb-min` | `1` | Mínimo de puntos de VB en una celda para usar solo VB. |
| `--no-edge-extrapolation` | off | Desactiva la rampa de margen a profundidad 0. |
| `--offset-scale` | `3.0` | Escala [m] de la penalización por offset perpendicular. |
| `--no-offset-weighting` | off | Desactiva esa penalización. |

---

## Salidas por transecta

En `./out_<perfil>/` (una transecta) o `<salida>/<perfil>/` (relevamiento):

| Archivo | Contenido |
|---|---|
| `bathymetric_profile.csv` | Perfil del lecho grillado (columnas abajo). |
| `raw_bed_points.csv` | Cada huella de beam proyectada (sin pesos). |
| `plan_view_map.png` | Vista en planta: trayectoria, eje de la sección, huellas por frecuencia. |
| `cross_section.png` | Sección transversal: nube + lecho suavizado (en cota absoluta si hay pelo de agua). |
| `axis.shp` | Eje de la sección (línea). |
| `profile_points.shp` | Puntos grillados del perfil. |
| `raw_beam_points.shp` | Huellas de beams proyectadas. |
| `process_log.txt` | Log de la corrida de esa transecta. |

## Salidas del relevamiento (`<salida>/_resumen/`)

| Archivo | Contenido |
|---|---|
| `survey_index.csv` | Una fila por transecta: progresiva, brazo, cotas, ancho, avisos. |
| `survey_profiles_all.csv` | Todos los puntos grillados de todas las secciones. |
| `survey_plan_view.png` | Traza del río + ejes de secciones, coloreados por brazo. |
| `survey_long_profile.png` | Perfil longitudinal: pelo de agua + thalweg vs. progresiva. |
| `survey_axes.shp` | Ejes de todas las secciones (EPSG:5344). |
| `survey_profile_points.shp` | Puntos de lecho de todas las secciones (EPSG:5344). |

---

## Esquema de columnas (CSV)

**`bathymetric_profile.csv`**
```
perfil, progresiva_m, brazo, s_m, depth_m, ws_elev_m, bed_elev_m,
latitude_deg, longitude_deg, x_utm, y_utm, x_posgar07, y_posgar07
```
- `s_m`: distancia transversal (0 = margen izquierda, MI).
- `depth_m`: profundidad bajo el pelo de agua.
- `ws_elev_m`: cota del pelo de agua interpolada (vacío si no hay CSV).
- `bed_elev_m`: **cota absoluta del lecho** (= `ws_elev_m` − `depth_m`) si hay pelo
  de agua; si no, relativa (0 en la superficie).

**`survey_index.csv`**
```
perfil, progresiva_m, progresiva, brazo, ws_elev_m, thalweg_elev_m,
max_depth_m, width_m, n_samples, dist_axis_m, theta_flow_deg, ws_note
```
- `progresiva`: en notación argentina `k+mmm,mm`.
- `thalweg_elev_m`: cota mínima del lecho de la sección.
- `dist_axis_m`: distancia del centroide de la transecta al eje del río.
- `ws_note`: aviso si el pelo de agua se extrapoló fuera del rango relevado.

**`survey_profiles_all.csv`**
```
perfil, progresiva_m, brazo, s_m, depth_m, ws_elev_m, bed_elev_m,
x_posgar07, y_posgar07
```

---

## Sistemas de coordenadas

- **Geometría de la transecta** (perfil, gráficos, shapefiles por transecta):
  UTM automático según longitud (acá **WGS84 UTM 19S, EPSG:32719**) y también
  **POSGAR07 faja** (acá **EPSG:5344**).
- **Progresivas y pelo de agua**: se calculan en el **CRS del centerline**
  (EPSG:5344).
- **Shapefiles consolidados del relevamiento**: **EPSG:5344**.
- **Cotas absolutas**: datum vertical **IGN SRVN16** (etiqueta configurable).

---

## Convención de la sección

- `s = 0` → **margen izquierda (MI)**; `s = s_max` → **margen derecha (MD)**
  (hidrográfico, mirando aguas abajo).
- El **lecho lo define el beam vertical (VB)** donde existe; los laterales sólo
  rellenan huecos interiores. Las **márgenes** son rectas desde el último VB de
  cada lado hasta la orilla (a la distancia de margen del *edge*).
