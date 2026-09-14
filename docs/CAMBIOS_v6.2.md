# process_adcp_bathimetric.py — cambios en v6.2

Dos bloques: **correcciones** sobre nueve bugs auditados en v6.1, y la **nube de
haces de los recorridos**, los `.mat` de la campaña que no son transectas.

El núcleo numérico no cambia: `build_beam_cloud`, `build_profile`, el acimut por
caudal y los filtros de profundidad son los de v6.0/v6.1. Lo que sí cambia de
resultado está listado en *Qué resultados cambian*.

---

## A. Correcciones

### A.1 Un solo CRS para toda la corrida (bugs 1, 4, 6, 7)

Los cuatro eran la misma falla: no había un CRS de salida definido en un lugar.

`process_one` elegía la faja POSGAR de cada transecta con
`auto_posgar07_crs(data["lon"])`, mientras `export_survey_shapefiles`,
`export_survey_beam_cloud` y `plot_survey_planview` estampaban
`CRS.from_epsg(5344)` fijo. El límite faja 2/3 de esa función cae en lon −67,5:
**todo el Negro aguas abajo de Chichinales es faja 3**, así que esas transectas
se escribían 739 km al este, sin un solo `[warn]`. Los shapefiles por transecta,
además, se rotulaban con el UTM interno, de modo que una campaña entregaba capas
en dos o tres CRS distintos.

Ahora hay una clave nueva y una regla única:

```ini
out-crs = EPSG:5344      ; default: el CRS del centerline; si no tiene, 5344 + [warn]
```

- Se resuelve **una vez**, al inicio: `--out-crs` → CRS del `--centerline` →
  `EPSG:5344` con aviso.
- Debe ser **proyectado y métrico**. Un CRS geográfico corta la corrida
  (bug 7: no había un solo `is_geographic` en el script, y un eje en lat/lon
  daba progresivas en grados con las tolerancias de 250 m, 15 m y
  `TRACK_LOC_MARGIN_MIN` interpretadas también como grados).
- `--centerline-crs` declara el CRS del eje cuando falta el `.prj`.
- **Toda** capa de salida va en ese CRS: `axis.shp`, `profile_points.shp`,
  `raw_beam_points.shp` (por transecta y por grupo), `survey_axes.shp`,
  `survey_profile_points.shp`, `survey_raw_beam_points.shp/gpkg` y
  `tracks_beam_points.shp/gpkg`. No queda ningún `5344` literal en el código.
- El UTM por transecta sigue siendo el marco **interno** donde se arman las
  huellas de haz y se proyecta la sección; se transforma una sola vez, al
  exportar. Así la corrección no toca el núcleo validado.
- `auto_posgar07_crs` queda como diagnóstico: si la faja natural de una
  transecta no coincide con `out-crs`, el log de esa transecta lo dice.
- Las columnas `x_posgar07` / `y_posgar07` de los CSV pasan a `x_out` / `y_out`,
  y `procesamiento.txt` registra el EPSG **efectivo**. El nombre viejo no decía
  en qué faja estaba.

**Aforos (bug 6).** `merge_group_clouds` reproyecta la nube de cada repetición
al UTM del miembro de referencia antes de concatenar. v6.1 concatenaba `E`/`N`
crudos: dos repeticiones a distinto lado de lon −66 daban 517 km de corrimiento.
Notar la asimetría que tenía v6.1: las *trazas* sí se reproyectaban por miembro,
sólo la nube no.

### A.2 Pendiente de continuación del pelo de agua (bug 2)

`_end_slope` ajustaba con los valores más cercanos al extremo y el único piso
era `np.ptp(ps) < 1.0` — un metro. Dos escalas a 1,7 m con ceros 229,15 y
229,56 daban **−241 m/km**, o sea 482 m de corrimiento a 2 km. El guard
existente (`slope > 0`) sólo atrapaba el caso que sube; el que baja pasaba
entero. El caso muerde cuando las anclas de un recorrido abarcan menos que
`WS_SLOPE_WINDOW` (1000 m): con una tercera ancla a 4 km el ajuste se recupera.

Claves nuevas:

```ini
ws-slope-min-span = 200      ; [m] separación mínima para definir una pendiente
ws-slope-max      = 2.0      ; [m/km] banda física admisible
ws-default-slope  = auto     ; auto | <m/km>
```

Orden de resolución, con el motivo anotado en la nota y en
`water_surface_qc.csv`:

1. Ventana del extremo, si abarca ≥ `ws-slope-min-span`.
2. Ajuste sobre el recorrido completo.
3. `ws-default-slope`. Con `auto`, la pendiente del recorrido cuyos valores
   abarcan más; si no hay ninguno, 0 con aviso.

En los tres pasos se rechaza una pendiente fuera de `[-ws-slope-max, 0]`. Toda
extrapolación que no haya podido usar la ventana del extremo queda marcada
`EXTRAP` aunque no supere `ws-extrap-tol`: la duda es la pendiente, no la
distancia.

`2.0 m/km` es el default porque el Neuquén, el Limay y el Negro corren entre
0,55 y 1,0 m/km.

### A.3 Dropouts GNSS 0/0 (bug 3)

RiverSurveyor escribe `0/0` —no `NaN`— cuando pierde el fix, y `0/0` es finito,
así que pasaba `np.isfinite(...).all()` y se usaba como posición del bote. El
alcance era mayor que eso: los `0/0` de `Latitude`/`Longitude` entraban también
en el promedio que elige la zona UTM y la faja POSGAR (un solo dropout entre dos
muestras buenas movía 19S → 23S y faja 2 → faja 7).

`_mask_gnss_dropouts` invalida el ensamble en `lat`, `lon` y `GPS.UTM` a la vez,
en `extract_data`, antes de que nada lo lea. El mismo criterio se aplica en el
pre-escaneo, para que la transecta se ubique con las mismas posiciones con las
que se procesa. Se reporta `[warn] GNSS: n/m ensambles sin fix (…%)` y, por
encima de `gps-max-gap` (0,30), la transecta no se procesa salvo
`--allow-gps-gaps`.

Donde `GPS.UTM` falta pero `Lat`/`Lon` son válidos, el ensamble se reconstruye
en vez de descartar la matriz entera como hacía v6.1.

### A.4 CRS declarado sin verificar (bug 5)

Dos cosas distintas.

**a)** El registro de estaciones estaba atado a `ws-crs`
(`crs_override=args.ws_crs`). Un registro *vectorial* usa su `.prj`, pero uno en
CSV no podía declarar el suyo — un problema real con estaciones en faja 3 y
pelo de agua en faja 2. Clave nueva `stations-crs` (default: el valor de
`ws-crs`, por compatibilidad).

**b)** Nadie verificaba que el CRS declarado fuera el real. `check_crs_plausibility`
mide la mediana de la distancia de los puntos al eje más cercano; si supera
1 km prueba las fajas POSGAR07 y las zonas UTM 19S/20S y nombra la que pone los
datos sobre el cauce, antes de cortar:

```
[error] CSV de pelo de agua: sus puntos quedan a 712 km de la red con el CRS
        declarado (EPSG:5345). Con EPSG:5344 la mediana baja a 0 m — revisá --ws-crs.
```

En v6.1 esto ya no era del todo silencioso —la regla de los 250 m descartaba los
puntos y §5 cortaba la corrida— pero el mensaje apuntaba al síntoma
("ninguna transecta tiene pelo de agua"), tres pasos después de la causa.

### A.5 `representative_time` (bug 8)

Un tiempo Unix de 2026 (1,79e9) caía en la rama SonTek y se leía como 2056: la
rama Unix era inalcanzable por debajo de 2e9.

Pero el **orden de ramas no estaba mal**: los rangos se solapan de verdad
(2e8–2e9 es SonTek 2006-2063 y Unix 1976-2033), así que ningún orden lo resuelve.
`_decode_epoch` construye todos los candidatos plausibles y, con
`time-epoch = auto`, elige el que caiga a menos de 7 días del sello del nombre
de archivo (`YYYYMMDDhhmmss`); si no hay sello o discrepan, usa SonTek, que es lo
que escribe RiverSurveyor Live y lo que `sontek_time` verificó contra los nombres
de archivo. `--time-epoch sontek|unix|datenum` fuerza la interpretación.
`sontek_time` usa el mismo decodificador.

### A.6 Colores de brazo (bug 9)

`hash()` de un `str` está salado por proceso: `MI` salía `tab:green`, `tab:blue`
y `tab:brown` en tres corridas consecutivas, y en una de ellas `MD` y `M3`
colisionaron en el mismo color. Ahora el color sale de la **posición en el
conjunto ordenado de brazos presentes en la corrida**, registrado con
`set_brazo_order()` antes de graficar: estable entre corridas y sin colisiones
dentro de una.

---

## B. Nube de haces de los recorridos

Los `.mat` que no son transectas —trayectos longitudinales entre secciones,
aproximaciones, reposicionamientos, mediciones estáticas— tienen fondo y GNSS
válidos y se descartaban. Ahora sus huellas de haz salen como **una nube de
puntos categorizada**, pensada para interpolar entre secciones y para el MDT del
cauce.

```ini
[procesamiento]
tracks = ./crudo/adcp/recorridos     ; carpeta, .mat sueltos o glob
```

No entran en `survey_index.csv`, ni en la agrupación de aforos, ni en el perfil
longitudinal, ni en el ranking de recorridos de la red. Es una salida paralela.

### B.1 Clasificación

No hay que listar nada: la carpeta es la declaración y el clasificador etiqueta
cada archivo solo. El discriminante es el tramo de **progresiva** que cubre la
traza (`dkm`) contra su extensión **transversal** (`dofs`, el rango de los
offsets al eje): una sección no avanza progresiva y cruza el cauce; un
longitudinal hace lo contrario.

| clase | criterio |
|---|---|
| `estatico` | traza de punta a punta < 30 m |
| `longitudinal` | `dkm > max(50 m, 2 × dofs)` |
| `aproximacion` | `dkm > 30 m` o `dofs > 30 m` |
| `otro` | el resto |

La sección `[recorridos]` del INI es sólo para el caso que salga mal, o para
forzar el río de un archivo medido en una confluencia:

```ini
[recorridos]
; <id de archivo> = <clase> [| <río>]
20260825151230 = longitudinal | Negro
```

Además, una transecta de `inputs` que cubra mucha más progresiva que traza deja
un `[warn]` en el pre-escaneo (*"parece un recorrido longitudinal… ¿va en
'tracks'?"*). Nunca se reclasifica sola.

### B.2 Procesamiento

1. `extract_data` con la máscara GNSS de A.3.
2. Filtro de picos con el mismo `--depth-filter`, pero **marcando en vez de
   descartar**: los puntos rechazados salen con `spike = 1`. Para MDT importa
   ver qué se filtró.
3. `build_beam_cloud` con `bt-geometry = footprints` forzado: colapsar al
   ensamble tira 4 puntos de 5, que es justo la cobertura lateral útil.
4. Cada ensamble se ubica en la red: río, brazo, progresiva oficial, distancia
   al eje y **offset perpendicular firmado** (+ hacia margen izquierda mirando
   aguas abajo).
5. **Progresiva y offset por huella, no por ensamble.** Un haz lateral aterriza
   a `profundidad × tan 25°` del bote — unos 2,1 m con 4,5 m de agua. Llevar el
   km del ensamble a sus cuatro huellas embarraría la nube justo a la escala a
   la que se la va a grillar. Se transfiere con el marco local del eje.
6. **Pelo de agua punto por punto.** Un longitudinal de 2 km atraviesa ~2 m de
   carga. `WaterSurfaceModel.resolve_series` evalúa el modelo sobre una grilla
   de progresivas de 1 m por río, cacheada, y agrega las progresivas de las
   anclas: el resultado es exacto, no una aproximación, y cuesta O(km) en vez de
   O(puntos).
7. `bed_elev = ws(punto) − depth(punto)`. Sin pelo de agua queda relativo y el
   punto se marca `REL`.

### B.3 Atributos

Geometría `PointZ` con **Z = `bed_elev`**, para que entre directo a un TIN o a un
IDW sin mapear campos. Nombres ≤ 10 caracteres por el límite del DBF.

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
| `depth` | profundidad desde el transductor (ya incluye draft) |
| `ws_elev` | pelo de agua en ese punto |
| `ws_src` | `survey` / `bridge` / `station` / `extrap` / `none` |
| `bed_elev` | cota de fondo SRVN16 |
| `river`, `brazo` | río y brazo asignados |
| `km_ofic` | progresiva oficial [m] |
| `off_axis` | offset perpendicular firmado al eje [m] |
| `dist_axis` | distancia al eje [m] |
| `t_utc` | ISO-8601 del ensamble |
| `spike` | 0/1 — rechazado por el filtro de picos |
| `flag` | `REL`, `AMBIGUOUS`, `OFF_AXIS`, concatenados |

`km_ofic` está a propósito: en `survey_raw_beam_points` se excluyó porque ahí
todos los puntos de una transecta comparten progresiva; acá varía punto a punto
y es la mitad de la información útil.

Para quedarse sólo con el haz vertical: `"beam_id" = 'VB'` en QGIS, o
`tracks-beams = vb`.

### B.4 Ubicación en confluencias

Un recorrido no tiene sección perpendicular, así que la regla del cruce de 6.1
no aplica: cada ensamble se ubica por eje más cercano. Si dos ríos quedan a menos
de 25 m de diferencia, el punto se marca `AMBIGUOUS`; el río se puede forzar por
archivo desde `[recorridos]`. Un punto a más de `tracks-locate-tol` (250 m) de
todo eje se **conserva**, con `km_ofic` vacío y `flag = OFF_AXIS`: para MDT un
punto en un brazo no digitalizado sigue sirviendo.

### B.5 Lo que queda explícitamente afuera

`--tracks-feed-sections` **no está implementado a propósito**. Meter puntos de
recorrido en las nubes de las secciones transversales es el paso siguiente, pero
cambia perfiles ya validados. La clave se acepta para que no rompa un INI y
emite un `[warn]`.

---

## C. Figuras

### `cross_section.png`

Los recuadros MI/MD estaban dentro de los ejes, en `y_floor + 0.98*span`, con
bbox blanco y `fontsize=11` en negrita, y la leyenda iba en `loc="lower right",
ncol=3` sobre el relleno del lecho. Regla nueva: **adentro del área de datos no
va nada que no sea dato**. Los márgenes son rótulos grises sobre el eje x, la
leyenda baja afuera en una fila sin marco, los rellenos y la nube bajan de
opacidad, el lecho pasa de `lw=2.0` a `1.6` y la grilla queda sólo horizontal.
Los números que había que buscar en el log (ancho, prof. máx, thalweg, pelo de
agua) están en la línea gris del subtítulo.

### `survey_plan_view.png`

Tres causas distintas:

- **Colores inestables** — es el bug 9, resuelto en A.6.
- **Leyenda con lo que no se ve.** El loop recorría *todos* los `network.axes` y
  *todos* los anabranches sin mirar la ventana de zoom, y el rótulo del brazo no
  tenía deduplicación (`seenr` sólo se chequeaba para el cauce principal). Por
  eso aparecían el Limay y el Negro con el relevamiento entero en el Neuquén, y
  `Neuquen MI (eje)` se repetía una vez por isla. Ahora la ventana se calcula
  primero, la geometría se clipea contra ella y a la leyenda sólo entra lo
  dibujado.
- **`62k` repetida.** `format_progresiva(kmo).split("+")[0] + "k"` tiraba todo lo
  que seguía al `+`: 62+150, 62+480 y 62+905 imprimían los tres `62k`. Ahora va
  la progresiva completa, con colocación por colisión real: se mide la caja de
  texto con el renderer y se prueba contra la geometría muestreada cada 3 px
  —cauce, ejes de brazo, **todas** las secciones, la escala gráfica—, contra las
  cajas ya puestas y contra el borde del área de dibujo, en doce posiciones
  candidatas a los dos extremos de la sección. Una sección con otra a menos de
  28 px arranca por los corrimientos largos y lleva siempre una guía fina: en un
  racimo, la cercanía sola no dice qué etiqueta es de qué sección.

Se agregan escala gráfica y norte, el rótulo de los ejes sale de `out-crs` (v6.1
decía `POSGAR07 f2` hardcodeado) y el tamaño de figura sale del aspecto de los
datos, con la ventana expandida para calzar la caja: con `figsize=(11, 12)` fijo
y `set_aspect("equal")`, un tramo de 50 km con secciones de 140 m dejaba media
figura en blanco.

### `survey_long_profile.png`

Línea continua de lecho por **interpolación monótona PCHIP**, no spline: no
sobrepasa entre puntos, o sea nunca dibuja un lecho más bajo que el thalweg
medido. `lw=1.3`, una línea por (recorrido, brazo) —MI y MD no se unen—, cortada
donde el hueco entre secciones supera `max(3 × espaciamiento mediano, 2 km)`.
Los marcadores bajan a `ms=3.8` y las verticales thalweg→pelo de agua a
`lw=0.6`. Las etiquetas de brazo, que se apilaban todas en `thal_min(thal)`,
ahora rotulan sólo la primera sección de cada corrida de brazo.

### `tracks_plan_view.png` (nueva)

Nube de los recorridos sobre los ejes, coloreada por clase, con las secciones
transversales de fondo.

---

## Claves nuevas

| Clave (INI / CLI) | Default | Qué hace |
|---|---|---|
| `out-crs` | CRS del centerline | CRS único de todas las salidas |
| `centerline-crs` | — | Declara el CRS del eje si falta el `.prj` |
| `stations-crs` | valor de `ws-crs` | CRS del registro de estaciones en CSV |
| `ws-slope-min-span` | `200` | [m] separación mínima para definir pendiente |
| `ws-slope-max` | `2.0` | [m/km] banda física admisible |
| `ws-default-slope` | `auto` | Pendiente de reserva |
| `gps-max-gap` | `0.30` | Fracción sin fix que aborta la transecta |
| `allow-gps-gaps` | `false` | Seguir igual |
| `time-epoch` | `auto` | `auto` \| `sontek` \| `unix` \| `datenum` |
| `tracks` | — | `.mat` que no son transectas |
| `tracks-beams` | `all` | `all` \| `vb` \| `slant` |
| `tracks-decimate` | `1` | Conservar 1 de cada N ensambles |
| `tracks-thin` | `0` | [m] raleo por grilla, 0 = off |
| `tracks-csv` | `on` | `on` \| `off` |
| `tracks-per-file` | `false` | Además, un shp/CSV por archivo |
| `tracks-locate-tol` | `250` | [m] más allá: `OFF_AXIS`, sin progresiva |
| `tracks-feed-sections` | `false` | No implementado a propósito; avisa |
| `[recorridos]` | — | `<id> = <clase> [\| <río>]` |

## Salidas que cambian

- **Nuevo:** `_resumen/tracks_beam_points.shp` (o `.gpkg`),
  `tracks_beam_points.csv`, `tracks_index.csv`, `tracks_plan_view.png`.
- **Renombrado:** en `bathymetric_profile.csv`, `raw_bed_points.csv`,
  `profile_points.shp` y `raw_beam_points.shp`, `x_posgar07` / `y_posgar07`
  pasan a `x_out` / `y_out`. Se agrega `beam_id` (`VB`, `B1`…`B4`) junto a
  `beam_index`.
- **CRS:** todas las capas en `out-crs`. Las por transecta dejan de estar en
  el UTM interno.
- **`procesamiento.txt`:** registra el EPSG efectivo de salida y las claves
  nuevas.

## Qué resultados cambian respecto de v6.1

- **Transectas con dropouts GNSS:** cambian nube, centroide, acimut y
  progresiva. Es la corrección, no una regresión.
- **Cualquier campaña al este de lon −67,5:** las capas cambian de rótulo y de
  coordenadas. El tramo del Negro aguas abajo de Chichinales cae íntegro en
  faja 3.
- **Recorridos con anclas de pelo de agua a menos de 200 m:** cambia la
  extrapolación.
- **La campaña 2026-08-25** está en faja 2 de punta a punta (confluencia
  Neuquén–Limay ≈ −68,06), así que el bug 1 no la tocó. Si tiene dropouts, sí
  cambia por A.3.

## Validación

Campaña sintética armada para pegarle a los bugs: eje de 52 km entre lon −67,8 y
−67,2 (cruza el límite faja 2/3), 9 transectas, dos con dropouts `0/0`
inyectados, dos escalas separadas 1,7 m con ceros 229,15 y 229,56, dos
recorridos longitudinales de 2,5 km.

| Chequeo | Resultado |
|---|---|
| CRS de las 7 capas de salida | EPSG:5344 en todas; la nube dentro del bbox del eje |
| Transecta en faja 3 natural | `[info]` en su log; se escribe en 5344 como el resto |
| Dropouts GNSS | `3/90 ensambles sin fix (3,3%) — descartados`; zona y faja sin alterar |
| `_end_slope`, par a 1,7 m | −0,90 m/km (pendiente por defecto), antes −241,2 m/km |
| `_end_slope`, par a 50 m | rechazado; a 500 m aceptado (−0,90 m/km) |
| `_end_slope`, 12 anclas / 50 km | −1,00 m/km exacto |
| Eje en EPSG:4326 | la corrida corta con el mensaje de CRS geográfico |
| `ws-crs` declarado 5345 siendo 5344 | corta nombrando 5344 y la mediana resultante |
| Aforo forzado por `[grupos]` | nubes reproyectadas al UTM de referencia; σ lecho 0,026 m |
| Pelo de agua en los recorridos | error máximo **0,3 mm** contra la pendiente analítica de 1,0 m/km |
| Huellas de un ensamble a 4,6 m | ±2,14 m del VB en proa/popa/babor/estribor = `4,6·tan 25°` |
| `--tracks-beams vb --tracks-decimate 3` | 200 puntos, sólo `VB` |
| Transecta única (sin `--config`) | capas en 5344 |

## Al reprocesar una campaña de v6.1

1. Correr `test_synthetic.py`.
2. Agregar `out-crs` al `.ini` si la campaña cruza fajas, o dejar que lo tome
   del eje. Confirmar en consola la línea `CRS de salida`.
3. Revisar los `[warn] GNSS:` por transecta. Un porcentaje alto de dropouts en
   una transecta que antes "salía bien" explica diferencias de progresiva.
4. En `water_surface_qc.csv`, revisar las notas de extrapolación: las que digan
   `pendiente por defecto` son las que v6.1 resolvía con un ajuste de dos
   valores pegados.
5. Si tenés estaciones en otra faja que el pelo de agua, declarar
   `stations-crs`.
6. Poner los `.mat` que no son transectas en una carpeta y apuntar `tracks` a
   ella. No hay que listarlos.
