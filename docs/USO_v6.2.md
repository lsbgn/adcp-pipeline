# process_adcp_bathimetric.py v6.2 — uso

Procesa relevamientos batimétricos con ADCP SonTek RiverSurveyor M9/S5.
Toma los `.mat` exportados por RiverSurveyor Live y entrega secciones
transversales con cotas de fondo absolutas en el datum vertical IGN SRVN16,
referidas al sistema de progresivas oficiales de cada río.

---

## 1. Invocación

```bash
# transecta única
python process_adcp_bathimetric.py ./crudo/adcp/20260825145258.mat \
    --centerline ../ejes/R_LI-NE-NE.shp \
    --water-surface-csv ./crudo/gnss/PA_20260825.csv \
    --outdir ./procesado

# relevamiento (2 o más .mat) -> modo survey
python process_adcp_bathimetric.py ./crudo/adcp \
    --centerline ../ejes/R_LI-NE-NE.shp \
    --water-surface-csv ./crudo/gnss/PA_20260825.csv \
    --outdir ./procesado --survey-name relevamiento-20260825

# producción: archivo de campaña
python process_adcp_bathimetric.py --config campanha.ini
```

Las rutas **relativas del `.ini`** se resuelven respecto de la ubicación del
`.ini`. Las de la línea de comandos, respecto del directorio actual. La línea de
comandos siempre tiene prioridad sobre el `.ini`.

Toda opción de la línea de comandos es una clave válida del `.ini`
(`--group-tol` se escribe `group-tol`). Una clave desconocida se avisa; un valor
fuera de las opciones corta la corrida.

---

## 2. Sistemas de coordenadas (6.2)

Esta es la parte que cambió más respecto de v6.1, y la que hay que entender si
el relevamiento cruza de faja.

### Entradas — cada una declara la suya

| Clave | Qué declara | Default |
|---|---|---|
| `centerline-crs` | CRS del eje, si le falta el `.prj` | el `.prj` del archivo |
| `ws-crs` | CRS del CSV de pelo de agua | `EPSG:5344` |
| `stations-crs` | CRS del registro de estaciones **en CSV** | el valor de `ws-crs` |

Un registro de estaciones **vectorial** (`.shp`, `.gpkg`) usa su propio CRS y no
necesita `stations-crs`. Esto importa si tenés estaciones en faja 3 y el pelo de
agua en faja 2: en v6.1 el registro estaba obligado a usar `ws-crs`.

Antes de procesar, el script verifica que el CRS declarado sea creíble: mide la
mediana de la distancia de los puntos al eje más cercano y, si supera 1 km,
prueba las fajas POSGAR07 y las zonas UTM y nombra la que pone los datos sobre
el cauce.

```
[error] CSV de pelo de agua: sus puntos quedan a 712 km de la red con el CRS
        declarado (EPSG:5345). Con EPSG:5344 la mediana baja a 0 m — revisá --ws-crs.
```

### Salida — una sola, para todo

```ini
out-crs = EPSG:5344
```

Se resuelve una vez, al inicio: `--out-crs` → CRS del `--centerline` →
`EPSG:5344` con aviso. Debe ser **proyectado y métrico**: un eje geográfico
cortaría la corrida, porque las progresivas y las tolerancias quedarían en
grados.

**Todas** las capas salen en ese CRS: las de campaña, las de cada transecta, las
de cada grupo y la nube de recorridos. Las columnas de los CSV se llaman `x_out`
/ `y_out` (en v6.1 eran `x_posgar07` / `y_posgar07`, que no decían la faja), y el
EPSG efectivo queda en `procesamiento.txt` y en el rótulo de los ejes de las
figuras.

**Relevamiento a caballo de dos fajas.** El UTM de cada transecta se sigue usando
como marco interno de cálculo, pero no llega a ninguna salida. Cuando la faja
POSGAR natural de una transecta no coincide con `out-crs`, su log lo dice:

```
[info] la faja POSGAR07 natural de esta transecta sería EPSG:5345;
       se escribe en 5344 como todo el relevamiento
```

El límite faja 2/3 cae en lon −67,5: de Chichinales para abajo, el Negro es
faja 3.

---

## 3. Entradas

### `.mat` de transectas

`inputs` acepta una carpeta, varios `.mat` (uno por línea o separados por coma) o
un glob. Dos o más disparan el modo survey.

**Ensambles sin fix GNSS.** RiverSurveyor escribe `0/0` cuando pierde el fix. Se
descartan, se reportan, y por encima de `gps-max-gap` (0,30) la transecta no se
procesa salvo `--allow-gps-gaps`:

```
[warn] GNSS: 3/90 ensambles sin fix (3.3%) — descartados (0/0 de RiverSurveyor)
```

### Red de ejes (`centerline`)

Un feature por tramo río + rol, con atributos `river` / `role`
[`/ seg_id / label / flow_dir`]. La topología entre ríos y los recorridos
cabecera → desembocadura se detectan solos.

### Pelo de agua GNSS (`water-surface-csv`)

Columnas E/N/H (configurables con `ws-east-col`, `ws-north-col`,
`ws-elev-col`). Una columna opcional `rio` fuerza el río de un punto, útil en la
confluencia.

### Escalas (`stations` + `readings`)

Registro con `station_id`, `name`, `river`, `gauge_zero` y X/Y (o geometría), más
las lecturas de la campaña: spot (`station_id,nivel`) o serie temporal
(`station_id,datetime,nivel`), autodetectado.

### Recorridos (`tracks`, 6.2)

Los `.mat` que **no** son transectas: trayectos longitudinales entre secciones,
aproximaciones, reposicionamientos, mediciones estáticas.

```ini
tracks = ./crudo/adcp/recorridos
```

No hay que listarlos: la carpeta es la declaración y el clasificador etiqueta
cada archivo por geometría. Ver §6.

---

## 4. Pelo de agua

GNSS como fuente primaria; escalas como puente, relleno y control. La cota de
cada confluencia se estima desde cada tributario con datos y se promedia con
peso 1/brecha.

Fuera del rango de datos el modelo **continúa con la pendiente local**. En 6.2
esa pendiente tiene guardas:

```ini
ws-slope-min-span = 200      ; [m] separación mínima para definirla
ws-slope-max      = 2.0      ; [m/km] banda física admisible
ws-default-slope  = auto     ; auto | <m/km>
```

Orden: ventana del extremo (si abarca ≥ `ws-slope-min-span`) → ajuste sobre el
recorrido completo → `ws-default-slope`. En los tres pasos se rechaza una
pendiente fuera de `[-ws-slope-max, 0]`. El motivo queda en la nota y en
`water_surface_qc.csv`.

Esto existe porque v6.1 aceptaba cualquier par separado más de 1 m: dos escalas a
1,7 m con ceros 229,15 y 229,56 daban −241 m/km, o sea 482 m de error a 2 km.

`ws_source` por transecta:

| Valor | Qué significa |
|---|---|
| `survey` | interpolación entre puntos GNSS |
| `bridge` | interpolación que mezcla GNSS, escala o confluencia |
| `station` | interpolación entre escalas |
| `extrap` | continuación con pendiente fuera del rango de datos |
| `none` | sin pelo de agua (cota relativa) |

Si hay fuentes configuradas y ninguna transecta puede usarlas, la corrida se
detiene antes de procesar. `allow-relative = true` permite seguir igual.

---

## 5. Opciones

### Malla y lecho

| Clave | Default | Qué hace |
|---|---|---|
| `dx` | `0.5` | espaciamiento de la grilla del perfil [m] |
| `bin-half-width` | `2*dx` | semiancho de muestreo por nodo [m] |
| `no-edge-extrapolation` | `false` | desactiva las rampas de margen |

### Referencia de profundidad (modelo QRev)

| Clave | Default | Opciones |
|---|---|---|
| `depth-ref` | `auto` | `auto` \| `vb` \| `bt` |
| `composite` | `on` | `on` \| `off` |
| `bt-avg` | `idw` | `idw` \| `simple` |
| `bt-geometry` | `footprints` | `footprints` \| `ensemble` |
| `edge-anchor` | `ref` | `ref` \| `cloud` |
| `vb-min` | `1` | mín. puntos de la primaria por nodo |
| `depth-filter` | `smooth` | `smooth` \| `off` |

### Ponderación y acimut

| Clave | Default | Notas |
|---|---|---|
| `density-weighting` | `on` | |
| `flow-weighting` | `discharge` | `discharge` \| `density` \| `none` |
| `section-orientation` | `flow` | `centerline` no está implementado; avisa |
| `offset-scale` | `3.0` | escala gaussiana de la penalización perpendicular [m] |
| `no-offset-weighting` | `false` | |

### Ubicación y agrupación

| Clave | Default | Notas |
|---|---|---|
| `transect-locate` | `crossing` | `crossing` \| `centroid` (regla v6.0) |
| `chainage-offset` | `0` | para ríos sin `[progresivas]` |
| `group-tol` | `15.0` | máx. Δ progresiva entre repeticiones [m] |
| `group-angle-tol` | `25.0` | máx. ángulo entre ejes de sección [°] |
| `group-max-span` | `0` | tope de duración de una ocupación [min], 0 = off |
| `no-group-repeats` | `false` | `true` = comportamiento v5 |

### GNSS y tiempo (6.2)

| Clave | Default | Notas |
|---|---|---|
| `gps-max-gap` | `0.30` | fracción sin fix que aborta la transecta |
| `allow-gps-gaps` | `false` | seguir igual |
| `time-epoch` | `auto` | `auto` \| `sontek` \| `unix` \| `datenum` |

`auto` construye los candidatos plausibles de `System.Time` y elige el que caiga
a menos de 7 días del sello del nombre de archivo (`YYYYMMDDhhmmss`); si no hay
sello o discrepan, usa SonTek. Los rangos de las codificaciones se solapan de
verdad, así que ningún orden de ramas los resuelve solo.

---

## 6. Recorridos (6.2)

### Clasificación

El discriminante es el tramo de progresiva que cubre la traza (`dkm`) contra su
extensión transversal (`dofs`): una sección no avanza progresiva y cruza el
cauce; un longitudinal hace lo contrario.

| clase | criterio |
|---|---|
| `estatico` | traza de punta a punta < 30 m |
| `longitudinal` | `dkm > max(50 m, 2 × dofs)` |
| `aproximacion` | `dkm > 30 m` o `dofs > 30 m` |
| `otro` | el resto |

La sección `[recorridos]` del INI es **opcional**, sólo para corregir un caso
puntual o forzar el río de un archivo medido en una confluencia:

```ini
[recorridos]
; <id de archivo> = <clase> [| <río>]
20260825151230 = longitudinal | Negro
20260825153044 = aproximacion
```

Una transecta de `inputs` que cubra mucha más progresiva que traza deja un
`[warn]` en el pre-escaneo. Nunca se reclasifica sola.

### Opciones

| Clave | Default | Qué hace |
|---|---|---|
| `tracks` | — | carpeta, `.mat` sueltos o glob |
| `tracks-beams` | `all` | `all` \| `vb` \| `slant` |
| `tracks-decimate` | `1` | conservar 1 de cada N ensambles |
| `tracks-thin` | `0` | raleo por grilla [m], 0 = off |
| `tracks-csv` | `on` | también escribir el CSV |
| `tracks-per-file` | `false` | además, un shp/CSV por archivo |
| `tracks-locate-tol` | `250` | [m] más allá: `OFF_AXIS`, sin progresiva |
| `tracks-feed-sections` | `false` | no implementado a propósito; avisa |

### Atributos de la nube

Geometría `PointZ` con **Z = `bed_elev`**: entra directo a un TIN o a un IDW sin
mapear campos.

| Campo | Contenido |
|---|---|
| `pt_id` | correlativo global de la campaña |
| `file` | stem del `.mat` |
| `class` | `longitudinal` / `aproximacion` / `estatico` / `otro` |
| `ens` | índice de ensamble dentro del archivo |
| `beam_id` | `VB`, `B1`…`B4` |
| `beam_type` | `vert` / `slant` |
| `src` | `VB` / `BT` |
| `freq_khz` | 500 (VB) / 1000 / 3000 |
| `depth` | profundidad desde el transductor (incluye draft) |
| `ws_elev` | pelo de agua en ese punto |
| `ws_src` | `survey` / `bridge` / `station` / `extrap` / `none` |
| `bed_elev` | cota de fondo SRVN16 |
| `river`, `brazo` | río y brazo asignados |
| `km_ofic` | progresiva oficial [m] |
| `off_axis` | offset perpendicular firmado [m], + = margen izquierda |
| `dist_axis` | distancia al eje [m] |
| `t_utc` | ISO-8601 del ensamble |
| `spike` | 0/1, rechazado por el filtro de picos (marcado, no descartado) |
| `flag` | `REL`, `AMBIGUOUS`, `OFF_AXIS`, concatenados |

La progresiva y el offset se calculan **por huella**, no por ensamble: un haz
lateral aterriza a `profundidad × tan 25°` del bote (≈2,1 m con 4,5 m de agua), y
a la escala a la que se grilla eso importa.

El pelo de agua también es **por punto**: un longitudinal de 2 km atraviesa unos
2 m de carga.

### Filtrar al haz vertical

```
"beam_id" = 'VB'
```

en QGIS, o directamente `tracks-beams = vb`.

---

## 7. Salidas

### Por transecta (`<id>/`)

```
bathymetric_profile.csv     perfil en grilla regular; columnas x_out / y_out
raw_bed_points.csv          nube cruda; columnas beam_index, beam_id
plan_view.png               traza, nube y eje de la sección
cross_section.png           sección transversal
axis.shp                    eje de la sección (LineString)
profile_points.shp          nodos del perfil (Point)
raw_beam_points.shp         huellas de haz (Point)
process_log.txt             log de la transecta
```

### Por grupo de aforo (`<id>_xN/`)

Lo mismo, más `group_qc.png` con la dispersión entre repeticiones.

### De campaña (`_resumen/`)

```
survey_index.csv            un renglón por transecta / grupo
survey_profiles_all.csv     todos los perfiles concatenados
survey_plan_view.png        vista en planta
survey_long_profile.png     perfil longitudinal, un panel por recorrido
survey_axes.shp             ejes de todas las secciones
survey_profile_points.shp   nodos de todos los perfiles
survey_raw_beam_points.shp  nube fusionada de todas las transectas (PointZ)
grupos_qc.csv               estadísticas de los aforos
water_surface_qc.csv        cada punto GNSS, escala y confluencia, con su rol
tracks_beam_points.shp      nube de los recorridos (PointZ)      [6.2]
tracks_beam_points.csv      la misma, en CSV                     [6.2]
tracks_index.csv            un renglón por archivo de recorrido  [6.2]
tracks_plan_view.png        nube de recorridos sobre los ejes    [6.2]
procesamiento.txt           parámetros efectivos, red, avisos
procesamiento.log           transcripción completa de la consola
```

Las capas muy grandes caen a GeoPackage (`.gpkg`) pasando ~1,5 M de features.

---

## 8. Chequeos al terminar una corrida

1. `CRS de salida` en consola: que sea el que esperabas.
2. `pelo de agua por perfil`: que no haya `none`.
3. Los `[warn] GNSS:` por transecta.
4. `water_surface_qc.csv`: puntos `AMBIGUOUS`, escalas de control fuera de
   tolerancia, y las notas de extrapolación que digan `pendiente por defecto`.
5. `grupos_qc.csv`: que la agrupación automática coincida con lo de campo.
6. `tracks_index.csv`: la clase asignada a cada recorrido y sus avisos.

---

## 9. Instalación

Python 3.10 o superior.

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
```

`geopandas` y `shapely` son opcionales en sentido estricto: sin ellos el script
corre pero no genera shapefiles ni progresivas. `scipy` es necesario para la
línea de lecho PCHIP del perfil longitudinal (sin él cae a una polilínea recta).

En Windows con el `.venv` activo, usar siempre `python`, no `py`.
