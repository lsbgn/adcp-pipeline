# Uso de `process_adcp_bathimetric.py` — v6.1

Procesa relevamientos batimétricos con SonTek RiverSurveyor M9/S5 desde los
`.mat` exportados por RiverSurveyor Live. Entrega cotas de fondo absolutas en
el datum vertical IGN SRVN16, CRS POSGAR07 faja 2 (EPSG:5344).

El nombre del archivo es estable. La versión está en `__version__` y queda
registrada en `procesamiento.txt` junto con el commit de git usado en la
corrida. Para que ese registro sea fiel, **correr siempre con los cambios
commiteados**.

Los cambios respecto de v6.0 están en `docs/CAMBIOS_v6.1.md`.

---

## Invocación

```bash
# Con archivo de configuración (recomendado en producción), desde la carpeta de la campaña
python D:/SIG/adcp-pipeline/process_adcp_bathimetric.py --config campanha.ini

# Una transecta
python process_adcp_bathimetric.py 20260415110139.mat --centerline ../ejes/R_LI-NE-NE.shp

# Campaña completa (carpeta, lista o glob) -> modo survey
python process_adcp_bathimetric.py ./crudo/adcp --centerline ../ejes/R_LI-NE-NE.shp
```

Dos o más `.mat` disparan el **modo survey**: una carpeta por transecta, más
`_resumen/` con los productos consolidados.

Precedencia: **línea de comandos > archivo `--config` > valores por defecto**.

Con el `.venv` activo, usar `python` y no `py`, que puede resolver al
intérprete global.

---

## Entradas

| entrada | formato | obligatoria |
|---|---|---|
| transectas | `.mat` de RiverSurveyor Live | sí |
| red de ejes | shapefile / gpkg / geojson con `river`, `role`, `seg_id` (opcional `label`, `flow_dir`) | no, pero sin ella no hay progresivas ni pelo de agua |
| pelo de agua GNSS | CSV con `East`, `North`, `H_correg` (opcionales `Punto`, `rio`) | no |
| registro de estaciones | CSV con X/Y, o capa de puntos | no |
| lecturas de nivel | CSV puntual o serie temporal (autodetectado) | sólo si se usa el registro |

### Red de ejes

- **Un feature por tramo río + rol**, con estos atributos:

  | atributo | contenido | obligatorio |
  |---|---|---|
  | `river` | nombre del río (Neuquen, Limay, Negro…) | sí |
  | `role` | `main` o `anabranch` | sí |
  | `seg_id` | identificador estable | sí |
  | `label` | fuerza el rótulo de un brazo | no |
  | `flow_dir` | `US2DS` o `DS2US` | no |

- **Brazos alrededor de islas.** Los `role = anabranch` se resuelven solos y se rotulan MI / M2… / MD, de izquierda a derecha mirando aguas abajo. La progresiva de un brazo es su proyección sobre el cauce principal.
- **Sentido de cada río** (km 0 aguas arriba), en este orden de prioridad:
  1. `flow_dir`.
  2. La pendiente del pelo de agua GNSS de ese río (≥3 puntos a menos de 50 m).
  3. `--chainage-reverse`.
  4. El orden de digitalización.
- **Topología entre ríos.** Se detecta sola por los extremos de las líneas. El extremo aguas abajo de cada tributario tiene que caer sobre la línea del río en el que desemboca, a menos de 2 m (enganchado con snap). Si dos ríos terminan en un nodo y un tercero empieza ahí (Neuquén + Limay → Negro), eso es una confluencia.

La consola informa la topología y los recorridos:

```
[info] topology: 'Neuquen' -> 'Negro' at its km 0 m
[info] topology: 'Limay' -> 'Negro' at its km 0 m
[info] topology: 'Negro' -> outlet
[info] path Neuquen>Negro: 77.03 km
[info] path Limay>Negro: 69.73 km
```

Si un tributario aparece como `outlet`, falta el enganche o el sentido está
invertido. En QGIS hay dos formas de corregirlo:

- Activar el snapping (`View > Toolbars > Snapping Toolbar`) y mover el extremo sobre la otra línea.
- Usar `Processing Toolbox > Snap geometries to layer`.

### Pelo de agua GNSS

- **Columnas.**
  - Coordenadas y cota: por defecto `East`, `North`, `H_correg`. Se cambian con `--ws-east-col` / `--ws-north-col` / `--ws-elev-col`.
  - Coordenadas en `--ws-crs` (default EPSG:5344).
  - Identificador: columna `Punto` (o `id` / `name`); si no existe, el número de fila.
  - Se acepta el BOM UTF-8 de Excel.
- **Río de cada punto.**
  - Por defecto, el eje más cercano.
  - Una columna opcional `rio` fuerza el río de ese punto; si queda vacía, rige la asignación automática.
- **Puntos marcados o descartados.**
  - Un punto con dos ríos candidatos a menos de 25 m de diferencia se marca `AMBIGUOUS` en `water_surface_qc.csv`. Es típico de las confluencias.
  - Se descartan, con aviso, los puntos a más de 250 m de su eje y los que caen fuera del tramo digitalizado.
- **Buenas prácticas.**
  - Sacar del CSV los puntos que no representan el pelo de agua del río. Un punto sobre un brazo que comunica dos ríos en la confluencia toma el nivel del más bajo y distorsiona la cota de confluencia.
  - En la confluencia misma, completar `rio` cuando la asignación no sea evidente.

### Registro de estaciones

- **Formato.**
  - CSV con `station_id`, `name`, `river`, `X`, `Y` (en `--ws-crs`) y `gauge_zero` (cero de escala, SRVN16).
  - Opcional: `chainage`, progresiva interna en m que reemplaza a la proyección sobre el eje.
  - Otras columnas (`operador`, …) se ignoran.
  - También sirve una capa de puntos con los mismos atributos, en su propio CRS.
- **Proyección.** Cada escala se proyecta sobre su **río declarado**. El nombre se compara sin tildes ni mayúsculas.
- **Exclusiones.** Quedan fuera, con el motivo anotado en consola y en `water_surface_qc.csv` (`excluded`):
  - escalas sin `gauge_zero` (no niveladas),
  - escalas a más de 250 m de su eje,
  - escalas fuera del tramo digitalizado,
  - escalas en un río que no está en la red.

### Lecturas de nivel

- **Formatos.**

  | modo | columnas |
  |---|---|
  | puntual | `station_id,nivel` |
  | serie temporal | `station_id,datetime,nivel` |

  La serie se detecta sola cuando hay fecha y más de una lectura por estación.
- **Cota.** `ws_elev = gauge_zero + nivel`.
- **Fechas.** Se aceptan `AAAA-MM-DD HH:MM[:SS]`, `AAAA-MM-DDTHH:MM[:SS]`, `DD/MM/AAAA HH:MM[:SS]` y `AAAA-MM-DD`.
- **Serie temporal.** Cada transecta usa la lectura interpolada a su hora (mediana de `System.Time`). Las horas tienen que estar en la misma referencia que el reloj del ADCP: no se convierte el huso.

---

## Cómo se ubica cada transecta

El río, el brazo y la progresiva salen del punto donde la **sección cruza un
eje**, como en HEC‑RAS.

1. **La sección es la sección real**: el eje ⟂ al flujo medio sobre el que se arma el perfil, abarcando la extensión de la traza proyectada sobre él, prolongado un 15 % (mínimo 15 m) a cada lado. **No es la traza del bote.** En una confluencia el bote puede subir por un tributario y bajar por el otro mientras la sección que representa cruza sólo el cauce de aguas abajo.
2. Si la sección cruza varios ejes, gana el cruce más cercano al centroide de la traza.
3. Si no cruza ninguno, se prueba el eje principal de la traza (`loc_method = crossing-track`, con nota en el log).
4. Si tampoco, se usa el eje más cercano al centroide (`nearest`, con `[warn]`).

`loc_method`, en `survey_index.csv` y en el log, dice cuál de los cuatro casos
se aplicó.

El criterio de v6.0, el eje más cercano al centroide, falla en las confluencias:
una sección sobre el Negro a pocos metros del nodo, con la traza cargada a una
margen, tiene el centroide más cerca del final del eje del Limay. El log lo
avisa:

```
[info] 20260825145258: Negro km 0+005.00 (crossing) — nearest axis to the track centroid is Limay (64 m) — the crossing wins
```

Hay dos formas de intervenir:

- `--transect-locate centroid` vuelve a la regla de v6.0, para reproducir progresivas viejas.
- La sección `[rios]` del INI fuerza el río de una transecta puntual. Si la sección no cruza el eje de ese río, la progresiva sale de proyectar el centroide sobre él, y el log lo aclara.

En transectas simétricas y perpendiculares al eje, los tres criterios coinciden.
Difieren en trazas curvas, asimétricas u oblicuas.

El **pre-escaneo por GPS** que corre antes de procesar no conoce todavía la
dirección del flujo, así que usa el eje de la traza. Su conteo por río es
provisional y el log lo dice; cada transecta se reubica sobre su sección real
al procesarse.

---

## Pelo de agua: red y recorridos

No hay tronco. Cada **recorrido** va de una cabecera a la desembocadura
(`Neuquen>Negro`, `Limay>Negro`) con progresiva continua. El tramo bajo una
confluencia pertenece a todos los recorridos que pasan por él: es el río común
donde se encuentran los pelos de agua de los tributarios. Cada río tiene su
propio pelo de agua.

1. **Anclas.** Son los puntos GNSS y las escalas que **no** quedan entre puntos GNSS.
   - Una escala entre puntos GNSS es sólo **control**: se compara contra el GNSS con tolerancia `--ws-qc-tol`.
   - Una escala fuera del rango GNSS es ancla. Por ejemplo, Isla Jordan y María Elvira cuando el GNSS cubre sólo el Neuquén.
2. **Cota de confluencia.** Se estima desde **cada tributario con datos**.
   - Para cada uno se interpola, a lo largo de su recorrido, entre su último valor aguas arriba del nodo y el primero aguas abajo. Es la pendiente de cada río continuada en el río común.
   - Los estimados se promedian con peso 1/brecha y dan una sola cota.
   - La diferencia entre ramas se controla con `--ws-qc-tol`.
3. **Dentro del rango de datos**, el pelo de agua se interpola linealmente en progresiva entre valores consecutivos.
4. **Fuera del rango**, se continúa con la **pendiente local**: mínimos cuadrados sobre los valores más cercanos a ese extremo, hasta abarcar al menos 1 km.
   - Un tributario sin datos propios toma la cota de la confluencia y sigue aguas arriba con la pendiente del tramo común.
   - Una extrapolación más larga que `--ws-extrap-tol` (200 m) deja la transecta con WARN.

Cada transecta registra `ws_source` y `ws_note`:

| `ws_source` | significa |
|---|---|
| `survey` | interpolación entre puntos GNSS |
| `bridge` | interpolación que mezcla GNSS, escala o cota de confluencia |
| `station` | interpolación entre escalas |
| `extrap` | continuación con pendiente fuera del rango (`EXTRAP` en mayúsculas = más allá de `ws-extrap-tol`) |
| `constant` | cota constante `--water-surface-elev` |
| `none` | sin pelo de agua: cotas **relativas**, siempre con `[warn]` |

Si hay fuentes de pelo de agua configuradas y **ninguna** transecta puede
usarlas, la corrida se detiene antes de procesar. Con `--allow-relative` sigue
igual.

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

- `auto` respeta `Setup.depthReference`, o sea lo que eligió el operador en campo.
- `footprints` mantiene los 4 haces laterales en su huella real, con mejor cobertura cerca de las márgenes. `ensemble` los colapsa a un punto promediado en la posición del bote, que es lo que hace QRev para caudal.
- `--vb-min` es el mínimo de puntos de la referencia primaria en la ventana de un nodo para usarla sola.
- `--blend-beams` (v5) se sigue aceptando pero es obsoleto. Equivale a `composite = on` + `bt-geometry = footprints`.

### Ponderación y acimut de la sección

| flag | valores | default |
|---|---|---|
| `--density-weighting` | `on`, `off` | `on` |
| `--flow-weighting` | `discharge`, `density`, `none` | `discharge` |
| `--use-edge-ensembles` | — | desactivado |
| `--section-orientation` | `flow` | `flow` |
| `--offset-scale` | metros | `3.0` |
| `--no-offset-weighting` | — | — |

- `--flow-weighting none` reproduce el cálculo de v5.
- `--section-orientation centerline` **no está implementada**. Ya era inerte en v6.0, donde además no avisaba; ahora emite un `[warn]`. El acimut de la sección siempre sale ⟂ al flujo medio.

### Agrupación de aforos

| flag | default | qué hace |
|---|---|---|
| `--no-group-repeats` | — | desactiva la agrupación (comportamiento v5) |
| `--group-tol` | `15.0` | máxima diferencia de progresiva entre repeticiones [m] |
| `--group-angle-tol` | `25.0` | máximo ángulo entre ejes de sección [°], y tolerancia QA de dispersión |
| `--group-max-span` | `0.0` | tope opcional de duración de una ocupación [min]; 0 = desactivado |

### Red y progresivas

| flag | default | qué hace |
|---|---|---|
| `--centerline` | — | red de ejes |
| `--transect-locate` | `crossing` | `crossing` (cruce de la sección) o `centroid` (regla v6.0) |
| `--chainage-offset` | `0.0` | offset de km oficial para los ríos que no figuran en `[progresivas]` |
| `--chainage-reverse` | — | medir desde el último vértice (sólo ríos sin `flow_dir` ni pelo de agua) |

### Pelo de agua

| flag | default | qué hace |
|---|---|---|
| `--water-surface-csv` | — | puntos GNSS de pelo de agua |
| `--ws-east-col` / `--ws-north-col` / `--ws-elev-col` | `East` / `North` / `H_correg` | columnas del CSV |
| `--ws-crs` | `EPSG:5344` | CRS del CSV de pelo de agua y del registro de estaciones en CSV |
| `--stations` / `--readings` | — | registro de estaciones y lecturas de la campaña |
| `--ws-qc-tol` | `0.10` | [m] control escala vs GNSS y dispersión de la cota de confluencia |
| `--ws-extrap-tol` | `200` | [m] extrapolación con pendiente más larga que esto → WARN |
| `--allow-relative` | — | seguir aunque ninguna transecta pueda tener pelo de agua |
| `--water-surface-elev` | `0.0` | cota constante para transectas sin otra fuente; 0 = desactivado |
| `--datum-name` | `IGN SRVN16` | rótulo del datum vertical |

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
| `process_log.txt` | log completo: ubicación (método y avisos), fuente del pelo de agua, QA |

### Por grupo de aforo — `<survey>/<primera>_xN/`

Lo mismo, más:

| archivo | contenido |
|---|---|
| `grupo_repetibilidad.png` | repeticiones superpuestas, perfil promediado y banda ±1σ |

### Consolidado — `<survey>/_resumen/`

| archivo | contenido |
|---|---|
| `survey_index.csv` | una fila por perfil **o grupo**: río, km oficial, brazo, cotas, ancho, `ws_source`, `ws_note`, caudal medido (`q_m3s`), `recorrido`, `loc_method` |
| `survey_profiles_all.csv` | todos los perfiles apilados |
| `grupos_qc.csv` | σ del lecho y total, RMS por repetición, dispersiones y caudal del aforo (media, desvío, CV, rango, valores) |
| `water_surface_qc.csv` | cada punto GNSS, escala y confluencia (ver abajo) |
| `survey_plan_view.png` | planta de la campaña |
| `survey_long_profile.png` | perfil longitudinal, **un panel por recorrido** |
| `survey_axes.shp`, `survey_profile_points.shp`, `survey_raw_beam_points.shp` | vectoriales consolidados |
| `procesamiento.txt` | registro de trazabilidad: versión, commit, parámetros efectivos, red, pelo de agua, advertencias |
| `procesamiento.log` | la consola completa de la corrida |

- **Los agregados de campaña ven sólo el perfil de grupo.** Un aforo de 4 repeticiones aparece una vez. Las carpetas individuales se conservan igual.
- **`survey_long_profile.png`.** Cada panel muestra:
  - el pelo de agua del modelo,
  - los puntos GNSS,
  - las escalas (ancla o control),
  - la cota estimada de cada confluencia.

  El tramo común bajo una confluencia se repite, con los mismos valores, en los paneles que lo necesitan. El eje x es la progresiva continua del recorrido, y cada confluencia está rotulada con el km oficial del río receptor.
- **`survey_chainage_m`** (en `survey_index.csv`) es la progresiva continua sobre el recorrido en el que se grafica la transecta. La progresiva oficial de cada río está en `km_oficial_m`.

### `water_surface_qc.csv`

Columnas:

- `kind`: `gnss`, `gauge` o `junction`.
- `id` y `river`.
- `km_internal_m` y `km_oficial_m`.
- `dist_axis_m`.
- `alt_river` y `alt_dist_m`: el segundo río más cercano y su distancia.
- `elev_m`, `role`, `note`.
- `x`, `y`.

Valores de `role`:

| `role` | significa |
|---|---|
| `anchor` | valor usado para construir el pelo de agua |
| `qc` | escala entre puntos GNSS: sólo control, con la diferencia en `note` |
| `junction` | cota estimada de una confluencia, con el estimado de cada rama en `note` |
| `dropped` | punto GNSS no usado (lejos del eje, fuera del tramo, fila ilegible) |
| `excluded` | escala fuera del cálculo (sin cero, lejos del eje, fuera del tramo) |
| `no reading` | escala válida sin lectura en esta campaña |

---

## Archivo de configuración

Ver `campanha_ejemplo_v6.ini`. Secciones reconocidas:

- `[campanha]` (o `[campania]`, `[campaña]`): metadatos de trazabilidad.
- `[procesamiento]`: cualquier flag, sin los guiones iniciales (`group-tol = 15`).
- `[progresivas]`: offset de km oficial por río. El nombre se compara sin tildes ni mayúsculas contra el eje; uno que no coincide se avisa.
- `[grupos]`: agrupación manual de aforos, `<nombre> = <id> <id> …`. Lo declarado acá sale del criterio automático, así que **un grupo de un solo miembro fija esa transecta como suelta** y sirve para deshacer una agrupación espuria:

  ```ini
  [grupos]
  suelta_113623 = 20260825113623
  suelta_113717 = 20260825113717
  ```

  El nombre del grupo es libre; un grupo de un miembro conserva el nombre original de la transecta.
- `[rios]`: río forzado por transecta, `<id> = <río>`. Se avisa si el id no está entre las entradas o si el río no está en la red.

Reglas:

- Las rutas relativas se resuelven respecto de la ubicación del `.ini`.
- `;` inicia un comentario en **todas** las secciones, incluida `[campanha]`: todo lo que sigue se descarta. Para separar ítems dentro de un valor, usar ` | `.
- Una clave o sección desconocida se avisa en consola y en `procesamiento.txt`.
- Un valor fuera de las opciones válidas (por ejemplo `depth-ref = vv`) corta la corrida con un error.

---

## Lectura de los QA y avisos

### Del aforo y la sección

```
[QA] dispersión de θ entre repeticiones: máx 4.2° (tol 25°) PASS
[QA] repetibilidad del lecho medido (s=3.5–44.5 m, donde las 4 repeticiones se solapan): σ media = 0.107 m, σ máx = 1.108 m
[QA] CV del caudal entre repeticiones: 13.3% (referencia USGS ≤5%) WARN
[QA] orientación: recorrido vs eje = +12.5°  (ancho proyectado pierde 2.4%)
[QA] orientación: sección ⟂flujo vs ⟂traza = -8.3° PASS
[QA] axis perpendicular to flow within ±10°? PASS
[QA] left bank at s = 0?  PASS
```

- **σ del lecho medido.** Es la repetibilidad real del fondo, calculada sólo donde todas las repeticiones se solapan.
  - La σ que incluye las rampas de margen es mayor y está dominada por cuánto se acercó cada pasada al banco, no por desacuerdo sobre el fondo.
  - La σ pica en los taludes, donde un error chico de posición horizontal se traduce en una diferencia vertical grande.
- **CV del caudal.** Compara las repeticiones entre sí. La referencia USGS es ≤5 %; por encima conviene revisar lecho móvil, referencia de track o largo de las transectas. Es independiente de la repetibilidad batimétrica. El caudal oficial del aforo se calcula en QRev.
- **Recorrido vs eje.** Mide cuánto ancho de sección se pierde por cruzar oblicuamente (`1 − cos θ`). Con 10° se pierde 1.5 %; con 30°, 13 %.
- **⟂flujo vs ⟂traza.** Responde "¿mis secciones son transversales al río?". Se evalúa contra el eje del río de la transecta.
  - Un desvío sistemático y de un solo signo apunta al cálculo o a la traza.
  - Uno que dispersa y sigue curvas, islas y bifurcaciones indica oblicuidad real del flujo.

### De la ubicación y el pelo de agua

```
[info] 20260825145258: Negro km 0+005.00 (crossing) — nearest axis to the track centroid is Limay (64 m) — the crossing wins
[info] confluence into 'Negro' at km 0 (Limay + Neuquen): stage 254.565 m — Neuquen: 254.604 (interp Pt149->NEG-01 over 1767 m); Limay: 254.526 (interp Pt150->NEG-01 over 1742 m); adopted 254.565 (interp)
[QA] confluence into 'Negro': the tributaries' estimates of the junction stage differ by 0.078 m (tol 0.10) PASS
[QA] gauge NEU-04 (Toma Cipolleti) between GNSS points: 263.740 m vs GNSS 263.788 m, Δ=-0.048 m (tol 0.10) PASS
[QA] water-surface interpolation: WARN — EXTRAP 500 m upstream of Pt1 at -1.04 m/km (9 values over 1000 m)
[info] pelo de agua por perfil: bridge 3, extrap 1, station 3, survey 8
```

(Valores de la prueba sintética, a modo de ejemplo.)

- **Ubicación.** Indica con qué criterio se ubicó la transecta. Cuando el cruce contradice al centroide, lo dice.
- **Confluencia.** Lista el estimado de cada tributario y la cota adoptada.
  - Una diferencia entre ramas mayor que `ws-qc-tol` suele venir de un punto mal asignado, un punto sobre un brazo de conexión o un cero de escala dudoso.
  - `stage EXTRAPOLATED` avisa que falta pelo de agua de un lado del nodo.
- **Escala entre puntos GNSS.** Controla escala contra GNSS. Un desvío grande apunta al cero, a la lectura o a la hora de la lectura.
- **EXTRAP.** La transecta está más de `ws-extrap-tol` fuera del rango de datos. La cota usa la pendiente local, pero conviene cubrir ese tramo con puntos GNSS o una escala.
- **Pelo de agua por perfil.** Es lo primero que hay que mirar: no debería aparecer `none`.
- **`water surface RISES downstream`.** Entre dos valores consecutivos el pelo de agua sube más de 5 cm hacia aguas abajo. Apunta a un punto mal ubicado, al cero de una escala o a una lectura.

---

## Qué revisar después de cada corrida

1. Consola, bloque inicial:
   - la topología y los recorridos son los esperados,
   - la cantidad de transectas por río es la esperada,
   - no hay avisos de ubicación inesperados.
2. La línea `pelo de agua por perfil`: sin `none`, y con `extrap` sólo donde no hay datos.
3. Los QA de confluencia y de escalas de control están en PASS.
4. En `water_surface_qc.csv` no aparecen puntos `AMBIGUOUS` o `dropped`, ni escalas `excluded`, que no se esperaban.
5. En `procesamiento.txt`, revisar las secciones `[advertencias de la corrida]` y `[advertencias QA por transecta]`.
6. En `survey_long_profile.png`, el pelo de agua es continuo en las confluencias y no sube hacia aguas abajo.

---

## Problemas frecuentes

| síntoma | causa probable | qué hacer |
|---|---|---|
| un tributario figura como `outlet` | extremo no enganchado al río receptor, o sentido invertido | snap en QGIS; revisar `flow_dir` |
| `the section line crosses no centerline` | transecta parcial o eje desplazado respecto del cauce | revisar el eje; si hace falta, forzar con `[rios]` |
| transecta asignada al río equivocado | geometría del eje en la confluencia | `[rios]` para esa transecta |
| punto GNSS `AMBIGUOUS` o en el río equivocado | punto en la confluencia o sobre un brazo de conexión | columna `rio`, o sacarlo del CSV |
| escala `excluded` | sin `gauge_zero`, lejos del eje o fuera del tramo digitalizado | completar el registro o extender el eje |
| `N transect(s) will get NO water surface` | río sin valores en ninguno de sus recorridos | revisar `water_surface_qc.csv` y la topología |
| dos transectas distintas agrupadas como aforo | su diferencia de progresiva quedó bajo `group-tol`; sobre un brazo, la progresiva es la proyección en el cauce principal y comprime la separación real | revisar `grupo_repetibilidad.png`; si son distintas, declarar cada una como grupo de un miembro en `[grupos]` |
| `unknown key` / `unknown section` en el INI | error de tipeo | corregir el nombre (`--help` lista las claves) |
| `[progresivas] 'X' matches no river` | el nombre no coincide con el atributo `river` del eje | corregir el nombre |

---

## Notas

- Las profundidades del `.mat` **ya incluyen** el calado del transductor (`Setup.sensorDepth`). No se le suma nada.
- `Edges_0__DistanceToBank` es la margen **izquierda** y `Edges_1` la derecha, siempre, sin importar `startEdge`.
- `System.Time` son segundos desde 2000-01-01.
- Sin `--centerline` no hay progresivas ni pelo de agua absoluto. La agrupación cae al desplazamiento perpendicular entre centroides, y el perfil longitudinal usa el eje simple.
- La profundidad de cada perfil no depende de la ubicación ni del pelo de agua: entre v6.0 y v6.1, `depth_m` es idéntica. Cambian la cota absoluta (`bed_elev_m`) y, en trazas asimétricas u oblicuas, la progresiva.
