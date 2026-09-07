# Uso de `process_adcp_bathimetric_v6.py`

Procesa relevamientos batimétricos con SonTek RiverSurveyor M9/S5 desde los
`.mat` exportados por RiverSurveyor Live, y entrega cotas de fondo absolutas en
el datum vertical IGN SRVN16, CRS POSGAR07 faja 2 (EPSG:5344).

---

## Invocación

```bash
# Una transecta
python process_adcp_bathimetric_v6.py 20260415110139.mat

# Campaña completa (carpeta, lista o glob) -> modo survey
python process_adcp_bathimetric_v6.py ./crudo/adcp --centerline ../ejes/R_LI-NE-NE.shp

# Con archivo de configuración (recomendado en producción)
python process_adcp_bathimetric_v6.py --config campanha.ini
```

Dos o más `.mat` disparan el **modo survey**: carpeta por transecta, más
`_resumen/` con los productos consolidados.

Precedencia: **línea de comandos > archivo `--config` > valores por defecto**.

---

## Entradas

| entrada | formato | obligatoria |
|---|---|---|
| transectas | `.mat` de RiverSurveyor Live | sí |
| eje del río | shapefile / gpkg / geojson con `river`, `role`, `seg_id` (opcional `label`, `flow_dir`) | no, pero sin él no hay progresivas |
| pelo de agua GNSS | CSV con `East`, `North`, `H_correg` | no |
| registro de estaciones | CSV | no |
| lecturas de nivel | CSV puntual o serie temporal (autodetectado) | no |

El eje admite varios ríos en un mismo archivo. Los `role = anabranch`
(brazos alrededor de islas) se resuelven solos y se rotulan MI / M2… / MD de
izquierda a derecha mirando aguas abajo.

---

## Flags

### Malla y lecho

| flag | default | qué hace |
|---|---|---|
| `--dx` | `0.5` | paso de la grilla del perfil [m] |
| `--bin-half-width` | `2*dx` | semiancho de la ventana en `s` para muestrear la nube |
| `--no-edge-extrapolation` | — | desactiva las rampas de margen desde `Setup.Edges_*` |

### Referencia de profundidad (modelo QRev)

| flag | valores | default |
|---|---|---|
| `--depth-ref` | `auto`, `vb`, `bt` | `auto` |
| `--composite` | `on`, `off` | `on` |
| `--bt-avg` | `idw`, `simple` | `idw` |
| `--bt-geometry` | `footprints`, `ensemble` | `footprints` |
| `--edge-anchor` | `ref`, `cloud` | `ref` |
| `--vb-min` | entero | `1` |
| `--depth-filter` | `smooth`, `off` | `smooth` |

`auto` respeta `Setup.depthReference`, o sea lo que eligió el operador en campo.
`footprints` mantiene los 4 haces laterales en su huella real (mejor cobertura
cerca de las márgenes); `ensemble` los colapsa a un punto promediado en la
posición del bote, que es lo que hace QRev para caudal. `--vb-min` es el mínimo
de puntos de la referencia primaria en la ventana de un nodo para usarla sola.

### Ponderación y acimut de la sección

| flag | valores | default |
|---|---|---|
| `--density-weighting` | `on`, `off` | `on` |
| `--flow-weighting` | `discharge`, `density`, `none` | `discharge` |
| `--use-edge-ensembles` | — | desactivado |
| `--section-orientation` | `flow`, `centerline` | `flow` |
| `--offset-scale` | metros | `3.0` |
| `--no-offset-weighting` | — | — |

`--flow-weighting none` reproduce el cálculo de v5.

### Agrupación de aforos

| flag | default | qué hace |
|---|---|---|
| `--no-group-repeats` | — | desactiva la agrupación (comportamiento v5) |
| `--group-tol` | `15.0` | máxima diferencia de progresiva entre repeticiones [m] |
| `--group-angle-tol` | `25.0` | máximo ángulo entre ejes de sección [°], y tolerancia QA de dispersión |
| `--group-max-span` | `0.0` | tope opcional de duración de una ocupación [min]; 0 = desactivado |

### Progresivas y pelo de agua

| flag | default |
|---|---|
| `--centerline` | — |
| `--chainage-offset` | `0.0` |
| `--chainage-reverse` | — |
| `--water-surface-csv` | — |
| `--water-surface-elev` | `0.0` |
| `--ws-east-col` / `--ws-north-col` / `--ws-elev-col` | `East` / `North` / `H_correg` |
| `--ws-crs` | `EPSG:5344` |
| `--stations` / `--readings` | — |
| `--ws-qc-tol` | `0.10` |
| `--datum-name` | `IGN SRVN16` |

---

## Salidas

### Por transecta — `<survey>/<perfil>/`

| archivo | contenido |
|---|---|
| `bathymetric_profile.csv` | perfil en grilla: `s`, profundidad, cota, `E/N`, lat/lon, POSGAR07 |
| `raw_bed_points.csv` | nube cruda de haces |
| `cross_section.png` | sección con la nube y el pelo de agua |
| `plan_view_map.png` | planta con trayectoria, huellas y eje |
| `axis.shp`, `profile_points.shp`, `raw_beam_points.shp` | vectoriales EPSG:5344 |
| `process_log.txt` | log completo con los QA |

### Por grupo de aforo — `<survey>/<primera>_xN/`

Lo mismo, más:

| archivo | contenido |
|---|---|
| `grupo_repetibilidad.png` | repeticiones superpuestas, perfil promediado y banda ±1σ |

### Consolidado — `<survey>/_resumen/`

| archivo | contenido |
|---|---|
| `survey_index.csv` | una fila por perfil **o grupo**, con caudal medido (`q_m3s`) |
| `survey_profiles_all.csv` | todos los perfiles apilados |
| `grupos_qc.csv` | σ del lecho y total, RMS por repetición, dispersiones y **caudal del aforo** (media, desvío, CV, rango, valores) |
| `survey_plan_view.png` | planta de la campaña |
| `survey_long_profile.png` | perfil longitudinal: pelo de agua y thalweg |
| `survey_axes.shp`, `survey_profile_points.shp`, `survey_raw_beam_points.shp` | vectoriales consolidados |
| `procesamiento.txt` | registro de trazabilidad de la corrida |

**Los agregados de campaña ven sólo el perfil de grupo.** Un aforo de 4
repeticiones aparece una vez. Las carpetas individuales se conservan igual.

---

## Archivo de configuración

Ver `campanha_ejemplo_v6.ini`. Secciones reconocidas:

- `[campanha]` (o `[campania]`, `[campaña]`) — metadatos de trazabilidad
- `[procesamiento]` — cualquier flag, sin los guiones iniciales
- `[progresivas]` — offset de km oficial por río
- `[grupos]` — agrupación manual de aforos

Las rutas relativas se resuelven respecto de la ubicación del `.ini`. Se admiten
comentarios en línea con `;`.

---

## Lectura de los QA

```
[QA] dispersión de θ entre repeticiones: máx 4.2° (tol 25°) PASS
[QA] repetibilidad del lecho medido (s=3.5–44.5 m, donde las 4 repeticiones se solapan): σ media = 0.107 m, σ máx = 1.108 m
[QA] CV del caudal entre repeticiones: 13.3% (referencia USGS ≤5%) WARN
[QA] orientación: recorrido vs eje = +12.5°  (ancho proyectado pierde 2.4%)
[QA] orientación: sección ⟂flujo vs ⟂traza = -8.3° PASS
[QA] axis perpendicular to flow within ±10°? PASS
[QA] left bank at s = 0?  PASS
```

- **σ del lecho medido** es la repetibilidad real del fondo, calculada sólo
  donde todas las repeticiones se solapan. La σ que incluye las rampas de margen
  es mayor y está dominada por cuánto se acercó cada pasada al banco, no por
  desacuerdo sobre el fondo. La σ pica en los taludes, donde un error chico de
  posición horizontal se traduce en una diferencia vertical grande.
- **CV del caudal** compara las repeticiones entre sí. La referencia USGS es
  ≤5 %; por encima conviene revisar lecho móvil, referencia de track o largo de
  las transectas. Es independiente de la repetibilidad batimétrica.
- **recorrido vs eje** mide cuánto ancho de sección se pierde por cruzar
  oblicuamente (`1 − cos θ`). Con 10° se pierde 1.5 %; con 30°, 13 %.
- **⟂flujo vs ⟂traza** responde "¿mis secciones son transversales al río?". Un
  desvío sistemático y de un solo signo apunta al cálculo o a la traza; uno que
  dispersa y sigue curvas, islas y bifurcaciones indica oblicuidad real del
  flujo.
- **WS GNSS vs estación** compara las dos fuentes de pelo de agua donde ambas
  existen, con tolerancia `--ws-qc-tol`.

---

## Notas

- Las profundidades del `.mat` **ya incluyen** el calado del transductor
  (`Setup.sensorDepth`). No se le suma nada.
- `Edges_0__DistanceToBank` es la margen **izquierda** y `Edges_1` la derecha,
  siempre, sin importar `startEdge`.
- `System.Time` son segundos desde 2000-01-01.
- Sin `--centerline` no hay progresivas: la agrupación cae al desplazamiento
  perpendicular entre centroides, y el perfil longitudinal no se genera.
