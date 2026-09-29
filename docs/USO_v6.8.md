# process_adcp_bathimetric.py v6.8 — guía de uso

Extrae secciones transversales batimétricas de archivos `.mat` de SonTek RiverSurveyor (M9/S5). Genera el perfil del lecho perpendicular al flujo, georreferenciado, en progresiva oficial y en cota IGN SRVN16. No exporta velocidades.

## 1. Requisitos

- Python 3.10 o superior, con el `.venv` del repo activo. En Windows, usar `python` y no `py`.
- `numpy`, `scipy`, `matplotlib`, `pyproj`, `geopandas`, `shapely` (ver `requirements.txt`). Sin geopandas corre igual, pero no escribe shapefiles.

## 2. Formas de correrlo

```powershell
# Campaña completa desde su INI (recomendado)
python process_adcp_bathimetric.py --config D:\Datos\ADCP\campanhas\2026-04-15\campanha.ini

# Una transecta suelta, sin red ni pelo de agua (cotas relativas)
python process_adcp_bathimetric.py 20260415122618.mat --outdir prueba_122618

# Carpeta de .mat con red y pelo de agua, sin INI
python process_adcp_bathimetric.py transectas\ --centerline R_LI-NE-NE.shp `
    --water-surface-csv pelo_de_agua.csv --survey-name neuquen_0415

# Pisar una clave del INI desde la línea de comandos
python process_adcp_bathimetric.py --config campanha.ini --position-fill off --outdir salida_v66
```

Prioridad de los valores: **CLI > INI > default**. Las rutas del INI se resuelven respecto de la carpeta del INI; las de la CLI, respecto de la carpeta desde donde se corre.

## 3. El INI de campaña

Hay un ejemplo completo y comentado en `campanha_ejemplo_v6.ini`. Secciones:

| Sección | Para qué |
|---|---|
| `[campanha]` | Metadatos libres (nombre, fecha, instrumental). Solo van al registro. |
| `[procesamiento]` | Entradas y parámetros; cualquier opción de la CLI sin los `--`. |
| `[progresivas]` | Offset oficial por río: `Neuquen = 0`. |
| `[grupos]` | Forzar un aforo: `NQ31 = <id1> <id2> ...` |
| `[rios]` | Forzar el río de una transecta: `<id> = Negro` |
| `[ubicacion]` | Forzar río, progresiva y brazo: `<id> = Neuquen \| <progresiva> \| MI` |
| `[orientacion]` | Forzar la sección al recorrido o a un azimut: `<id> = recorrido` o `<id> = 135` |
| `[saltos]` | Discontinuidades del pelo de agua: `<nombre> = Neuquen \| <progresiva>` o `<nombre> = <x> <y>` |
| `[recorridos]` | Corregir la clase de un recorrido: `<id> = longitudinal \| Neuquen` |
| `[descartes]` | Sacar del ajuste una fuente de profundidad en todo o parte de una sección (6.8, ver apartado 7): `<id> = vb s 131.5-134.0` |

Las claves desconocidas y los valores fuera de las opciones válidas se avisan o cortan la corrida; no se ignoran en silencio.

## 4. Claves de `[procesamiento]`

### Entradas y salida

| Clave | Default | |
|---|---|---|
| `inputs` | — | Carpetas, archivos o globs de `.mat` de secciones |
| `tracks` | — | `.mat` que no son secciones: se exportan como nube de puntos |
| `outdir` | `./<survey-name>/` | |
| `survey-name` | `salida_<primer perfil>` | |

### Red, ubicación y CRS

| Clave | Default | |
|---|---|---|
| `centerline` | — | Red de ejes (`river`, `role`, `seg_id`, `label`, `flow_dir`) |
| `centerline-crs` | del `.prj` | Solo si falta el `.prj` |
| `out-crs` | CRS del centerline, si no EPSG:5344 | Un único CRS para todas las salidas |
| `transect-locate` | `crossing` | `crossing`: donde la sección cruza el eje. `centroid`: regla v6.0 |
| `chainage-offset` | 0 | Offset global para los ríos que no están en `[progresivas]` |
| `chainage-reverse` | false | |

### Pelo de agua

| Clave | Default | |
|---|---|---|
| `water-surface-csv` | — | Puntos GNSS del pelo de agua (fuente primaria) |
| `ws-east-col`, `ws-north-col`, `ws-elev-col` | `East`, `North`, `H_correg` | |
| `ws-crs` | EPSG:5344 | |
| `stations`, `readings` | — | Registro de escalas y lecturas de la campaña (respaldo y QC) |
| `stations-crs` | = `ws-crs` | |
| `ws-qc-tol` | 0,10 m | Escala contra superficie GNSS; dispersión en la confluencia |
| `ws-axis-tol` | 300 m | Punto o escala más lejos del eje: no se usa |
| `ws-extrap-tol` | 200 m | Extrapolación más larga: `[warn]` |
| `ws-slope-min-span` | 200 m | |
| `ws-slope-max` | 2,0 m/km | |
| `ws-default-slope` | `auto` | |
| `water-surface-elev` | 0 | Cota constante, solo sin CSV |
| `datum-name` | `IGN SRVN16` | |
| `allow-relative` | false | Seguir con cotas relativas si nadie obtiene pelo de agua |

### Posiciones GNSS del M9 (nuevo en 6.7)

| Clave | Default | |
|---|---|---|
| `position-fill` | `bt` | `bt`: descarta por HDOP y por salto GGA-BT y reconstruye con bottom track. `interp`: solo HDOP, interpolación en el tiempo. `off`: sin filtro (v6.6) |
| `gga-hdop-max` | 4 | HDOP mayor: inválido. 0 = sin filtro |
| `gga-hdop-change` | 2 | Desvío del HDOP respecto de la media de los válidos. 0 = sin filtro |
| `gga-bt-tol` | 1,5 m | Distancia de un fijo a la tendencia local GGA − BT. 0 = sin control de saltos |
| `gga-bt-window` | 30 s | Semiventana de esa tendencia |
| `s-fold-tol` | 2 m | QA: retroceso de s en el cruce |
| `edge-drift-tol` | 5 m | QA: deriva del bote en la medición de borde |
| `gps-max-gap` | 0,30 | Fracción de ensambles sin fix que descarta la transecta |
| `allow-gps-gaps` | false | |

### Montaje y actitud

| Clave | Default | |
|---|---|---|
| `beam-layout` | `auto` | Numeración de haces leída de la matriz del archivo |
| `pitch-roll` | `off` | `on` / `inv` inclinan las huellas |
| `pitch-sign`, `roll-sign` | `bow-up`, `stbd-down` | |
| `antenna-offset` | — | `proa,estribor` [m], solo si no se cargó en RSL |
| `antenna-height` | 0 | |
| `gnss-lag` | 0 s | |
| `time-epoch` | `auto` | |

### Profundidad

| Clave | Default | |
|---|---|---|
| `depth-ref` | `auto` | `vb` / `bt`; `auto` = lo elegido en RSL |
| `composite` | `on` | Relleno QRev con la otra fuente |
| `bt-avg` | `idw` | |
| `bt-geometry` | `footprints` | |
| `depth-filter` | `smooth` | Filtro de picos por haz |
| `vb-min` | 1 | |

### Lecho y márgenes

| Clave | Default | |
|---|---|---|
| `dx` | 0,5 m | Paso de la grilla del perfil |
| `bed-fit` | `loess` | `median` reproduce v6.3 |
| `bin-half-width` | 0,75 m (loess) | |
| `edge-shape` | `auto` | Forma cargada en RSL (triangular / rectangular) |
| `edge-extent` | `mean` | |
| `edge-anchor` | `ref` | |
| `no-edge-extrapolation` | false | |
| `offset-scale` | 3,0 m | Escala gaussiana de la penalización perpendicular |
| `no-offset-weighting` | false | |

### Sección y aforos

| Clave | Default | |
|---|---|---|
| `flow-weighting` | `discharge` | Azimut de la sección = dirección del caudal unitario total |
| `density-weighting` | `on` | |
| `use-edge-ensembles` | false | |
| `section-orientation` | `flow` | |
| `group-tol` | 15 m | |
| `group-angle-tol` | 25° | |
| `group-max-span` | 0 (off) | |
| `no-group-repeats` | false | |

### Recorridos y figuras

| Clave | Default | |
|---|---|---|
| `tracks-beams` | `all` | `vb` / `slant` |
| `tracks-decimate` | 1 | |
| `tracks-thin` | 0 m | |
| `tracks-csv` | `on` | |
| `tracks-per-file` | false | |
| `tracks-locate-tol` | 250 m | |
| `planview-labels` | `auto` | |

## 5. Salidas

**Por transecta o grupo** (`<salida>/<perfil>/`):

- `bathymetric_profile.csv`: s, profundidad, pelo de agua, cota del lecho, coordenadas.
- `raw_bed_points.csv`: cada huella de haz proyectada. Desde 6.8 termina en `ens` (número de ensamble de RSL) y `descarte` (la regla que lo sacó; vacío si se usó).
- `cross_section.png` y `plan_view_map.png`.
- `axis.shp`, `profile_points.shp`, `raw_beam_points.shp`.
- `process_log.txt`.

**Campaña** (`<salida>/_resumen/`):

- `survey_index.csv`, `survey_profiles_all.csv`, `grupos_qc.csv`, `water_surface_qc.csv`.
- `survey_plan_view.png`, `survey_long_profile.png`.
- `survey_axes.shp`, `survey_profile_points.shp`, `survey_raw_beam_points.shp`.
- `tracks_index.csv`, `tracks_beam_points.shp/.csv`, `tracks_plan_view.png`.
- `procesamiento.txt` (parámetros efectivos y todas las advertencias) y `procesamiento.log` (la consola completa).

## 6. Cómo revisar las posiciones (6.7)

En `survey_index.csv`, conviene ordenar por `pos_shift_max_m` y mirar primero:

- **`pos_rebuilt > 0` con `pos_shift_max_m > 1`.** El filtro movió posiciones de forma apreciable. Abrí `plan_view_map.png`: el GGA descartado aparece en gris, unido a su posición nueva en violeta. La traza naranja tiene que ser continua y suave.
- **`pos_rebuilt` mayor que la mitad de `n_samples`.** La georreferencia se apoya en pocos fijos. El log lo marca con `[warn]`.
- **`s_fold_m > 2`.** La traza vuelve sobre sí misma más de 2 m. Revisá si el corte muestra dos profundidades en el mismo s.
- **`edge_drift_m > 5`.** El bote derivó durante la medición de borde. Los retornos lejos de la línea ya no pesan en el ajuste, pero la distancia a la orilla cargada en RSL se midió desde otro lugar, así que conviene revisarla.

En `process_log.txt` de cada perfil:

```
[info] GNSS: calidad GGA 1:63 (1 autónomo, 2 diferencial, 4 RTK); HDOP 1.2/2.0/22.7 (mín/med/máx), satélites 4-7; relleno=bt
[warn] posiciones: 20/63 ensambles reconstruidos (bt) — 8 por HDOP, 12 por salto GGA-BT > 1.5 m, 0 sin fix; corrimiento máx 10.14 m
[QA] retroceso de s en el cruce: 1.25 m (tol 2.0 m) PASS
[QA] deriva en el borde inicial (10 ens.): 0.9 m transversal a la sección, 2.3 m a lo largo (tol 5.0 m) PASS
```

**Comparar contra v6.6.** Corré con `--position-fill off --outdir salida_v66` y compará los `survey_index.csv`. Las transectas limpias son idénticas; solo cambian las que tenían posiciones inválidas.

**Si el filtro corrige de más** (poco probable con el autónomo del M9): subí `gga-bt-tol` a 2–3 m, o usá `interp` para confiar solo en el HDOP. Si hubo fondo móvil fuerte, la tendencia GGA − BT deriva más rápido; la ventana de 30 s lo absorbe en cruces normales. No la bajes de 20 s.

## 7. Descartes manuales — `[descartes]` (6.8)

El filtro de picos trabaja sobre la serie de cada haz y no rechaza una tanda de varios ensambles en una barranca, porque ahí el rango intercuartil es grande. Con `[descartes]` se saca una fuente donde se sabe que está mal. El composite rellena esos nodos con la fuente secundaria, igual que QRev cuando falta la primaria.

```ini
[descartes]
; <id> = <regla>, <regla>, ...        (o una regla por renglón, con sangría)
; <regla> = <fuente>[+<fuente>] [<condición> <valor>] ...
20260825153036 = vb s 131.5-134.0
```

**`<id>`**: el nombre del `.mat` sin extensión, o solo sus 14 dígitos (`20260825153036` vale para `20260825153036r`). También puede ser el id de un grupo (`<id>_x4`) o un nombre de `[grupos]`; entonces las condiciones se leen sobre el eje del grupo.

**Fuentes**

| Fuente | Qué saca |
|---|---|
| `vb` | Haz vertical |
| `b1` … `b4` | Un haz inclinado (el número de RSL) |
| `bt` | Los cuatro inclinados |
| `todos` | Todos los haces |

Se combinan con `+`: `vb+b1`.

**Condiciones** (todas deben cumplirse; sin condición se saca la fuente en toda la sección)

| Condición | Sobre qué | Ejemplo |
|---|---|---|
| `s` | Distancia transversal [m], como en `cross_section.png` | `s 131.5-134.0` |
| `ens` | Número de ensamble de RSL (desde 1; columna `ens` de `raw_bed_points.csv`) | `ens 27-30`, `ens 28` |
| `prof` | Profundidad del haz [m] | `prof <3` |
| `cota` | Cota del punto [m]; necesita pelo de agua | `cota >249.4` |
| `dif` | Diferencia con la otra fuente en el mismo ensamble [m]: para el VB, contra la mediana de los inclinados; para un inclinado, contra el VB | `dif >1.5` |

Valores: `A-B`, `<X`, `>X` (bordes incluidos), o un número suelto para `ens`. **El separador decimal es el punto**; la coma separa reglas.

**Ejemplos**

```ini
[descartes]
; el raigón del VB en el Negro 2+018
20260825153036 = vb s 131.5-134.0
; lo mismo, por ensambles
20260825153036 = vb ens 27-30
; lo mismo, sin decir dónde: el VB que difiere más de 1.5 m de los inclinados
20260825153036 = vb dif >1.5
; solo lo que sube por encima de 249.4 m en ese tramo
20260825153036 = vb cota >249.4 s 125-136
; varias reglas en la misma transecta
20260825153036 = vb s 131.5-134.0, b1 ens 28
; un haz que anduvo mal en toda la sección
20260415125906 = b4
; en un aforo, sobre el eje del grupo
NQ31 = vb s 40-42 prof <1.5
```

**Cómo leer los rangos.** `s` se lee en la sección **sin descartes**, o sea el `cross_section.png` que miraste para escribir la regla. Después la extensión medida se recalcula con los puntos que quedan. Si eso corre el origen de s (solo pasa si descartás los puntos de un extremo), el log lo dice.

**Qué queda registrado**

- `process_log.txt`: una línea por regla con cuántos puntos sacó, de qué haz, y el rango de s y de ensambles. Si una regla no coincide con ningún punto, sale un `[warn]`.
- `cross_section.png`: los puntos descartados como cruces rojas.
- `raw_bed_points.csv`: los descartados siguen en el archivo, con la regla en `descarte`.
- `survey_index.csv`: columna `descartes` al final, con la cantidad de puntos (en un grupo, los de los miembros más los del grupo).
- La nube de la campaña (`survey_raw_beam_points.shp`) no incluye los descartados.

Una regla mal escrita corta la corrida al empezar, con el renglón y el motivo.

**Cuidado con `dif`.** En una barranca los haces inclinados ven fondo a otra profundidad que el VB, y la diferencia es real. En el Negro 2+018, `vb dif >1.0` saca, además del raigón (ens. 27–30), dos ensambles del pie de la barranca (21–22). Conviene acotarlo con `s` o `ens`, o subir el umbral.

**En aforos.** Una regla con el id de una transecta se aplica a esa repetición, y el grupo ya no usa esos puntos. Una regla con el id del grupo se aplica a la nube fusionada, sobre el eje del grupo. `ens` no vale en una regla de grupo, porque cada repetición numera sus ensambles.

## 8. Problemas frecuentes

| Síntoma | Qué mirar |
|---|---|
| Lecho en serrucho o escalones que RSL no muestra | Columnas `pos_*` y `plan_view_map.png`. Si `pos_rebuilt = 0`, el problema no es el GNSS: mirar `edge_drift_m` y la dispersión de la nube. |
| "puntos flotando" sobre o bajo el lecho cerca de una margen | Deriva en el borde (`[warn] borde ...`) |
| Perfil con cota relativa | `ws_source = none`: falta pelo de agua en esa progresiva |
| La transecta cae en otro río | `[rios]` o `[ubicacion]` |
| Sección cruzada a una obra | `[orientacion]` con `recorrido` o un azimut |
| Un pico del VB (raigón, ramas, peces) que el filtro no saca | `[descartes]` con `vb` y el tramo (apartado 7) |

## 9. Comparar dos corridas — `tools/comparar_corridas.py`

```powershell
python tools\comparar_corridas.py salida_v67 salida_v68            # A = referencia, B = nueva
python tools\comparar_corridas.py salida_v67 salida_v68 --plots all
python tools\comparar_corridas.py salida_v66 salida_off --check    # regresión: exit 1 si algo cambia
```

Escribe `<B>/_comparacion/` con:

- `comparacion.csv`: una fila por perfil, ordenada del cambio más grande al más chico. Tiene el `estado`, RMS / sesgo / máx |Δz| y el s del máximo, Δancho, Δprof. máx., Δthalweg, Δpelo de agua, Δprogresiva, cuánto se movió el eje, si cambió la nube cruda (`nube_igual`, sobre las 12 columnas de v6.7), cuántos puntos descartó `[descartes]` en cada corrida (`descartes_a`, `descartes_b`) y las columnas `pos_*` de la corrida B.
- `parametros_diff.txt`: las claves con distinto valor entre los dos `procesamiento.txt`, más la versión del script de cada corrida.
- `<perfil>.png`: lechos A (gris, a trazos) y B (negro) superpuestos, con Δz(s) abajo. Por defecto, solo para los que cambian.

| `estado` | Significado |
|---|---|
| `igual` | `bathymetric_profile.csv` idéntico byte a byte |
| `igual_num` | Archivo distinto pero máx \|Δz\| ≤ `--tol-z` (5 mm) y mismo ancho: ruido numérico |
| `cambia` | Cambio real |
| `solo_A` / `solo_B` | El perfil está en una sola corrida (p. ej., cambió la agrupación de aforos) |
| `sin_solape` | Las dos secciones no se superponen en el terreno |

La comparación se hace **sobre el terreno**, no sobre s: los nodos de B se proyectan con sus coordenadas sobre la línea de sección de A. Así, un s = 0 corrido (sección más ancha o más angosta) o un eje apenas rotado no aparecen como un cambio del lecho. Se compara la cota cuando las dos corridas tienen pelo de agua; si no, la profundidad.

Para cambiar parámetros, una corrida por variante y la misma referencia:

```powershell
python process_adcp_bathimetric.py --config campanha.ini --outdir salida_ref
python process_adcp_bathimetric.py --config campanha.ini --outdir salida_tol2 --gga-bt-tol 2
python tools\comparar_corridas.py salida_ref salida_tol2
```

Para ver qué cambia un descarte, la corrida de referencia es la misma campaña sin la regla (comentá el renglón con `;`):

```powershell
python process_adcp_bathimetric.py --config campanha.ini --outdir salida_con_descartes
python tools\comparar_corridas.py salida_ref salida_con_descartes
```
