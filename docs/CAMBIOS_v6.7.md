# Cambios en v6.7 respecto de v6.6

**Tema:** filtro de las posiciones GNSS del bote y controles de calidad de la traza.

## Por qué

El M9 se posiciona con su GPS interno, casi siempre en modo autónomo (GGA calidad 1). En la campaña 2026-04-15 aparecieron secciones con el lecho en serrucho o con escalones que no existen, aunque RiverSurveyor Live las mostraba bien.

La causa no es el modo autónomo en sí. Su sesgo de ~2 m es estable mientras el HDOP y la constelación no cambian: corre la traza entera y no deforma la sección. El problema aparece cuando entra o sale un satélite. Ahí la posición salta varios metros **a lo largo de la línea de sección**, y s se estira o se pliega. La penalización por distancia perpendicular (v6.4) no ve ese error, porque es paralelo a la sección.

| Archivo | Qué pasaba en v6.6 |
|---|---|
| 20260415122618 | HDOP 22,6 en el borde inicial y 16,5 a mitad del cruce. Los 10 ensambles de borde (prof. ~2,0 m) quedaban 10 m adentro del cauce, mezclados con fondo de 2,6 m. Lecho en serrucho y prof. máx. 2,56 m contra 2,74 m medidos con VB. |
| 20260415135202 | HDOP 4,3–7,1 en los primeros 25 ensambles: la rampa de la margen de inicio quedaba comprimida unos 3 m. |
| 20260415131650 | El GNSS está bien. El bote derivó 21 m transversal a la sección durante la medición de borde. El ajuste ya lo descartaba, pero no quedaba ningún aviso. |

RiverSurveyor Live no muestra el problema porque su eje "Track" es la distancia recorrida acumulada, que siempre crece. No proyecta sobre una línea de sección.

## Qué cambia

### 1. Filtro de posiciones GGA (`position-fill = bt`, por defecto)

Se aplica a cada ensamble, antes de cualquier uso de la posición: azimut de la sección, nube de haces, ubicación en la red.

1. **HDOP, como QRev.** Un fijo con HDOP > `gga-hdop-max` (4) es inválido. Después, iterativamente, también lo es el que se aparta de la media de los válidos más de `gga-hdop-change` (2).
2. **Consistencia GGA–bottom track.** Se integra la trayectoria del bottom track (velocidad del bote = −`BottomTrack.BT_Vel`, en ENU). La diferencia GGA − BT tiene que variar despacio: sesgo del receptor, desalineación de la brújula, fondo móvil. Se ajusta una tendencia lineal local de ±`gga-bt-window` (30 s). Un fijo a más de `gga-bt-tol` (1,5 m) de esa tendencia es un salto aunque su HDOP sea normal (en 122618 el HDOP era 2,0 y la traza retrocedía mientras el bote avanzaba).
3. **Reconstrucción.** Cada tramo inválido se rearma con los incrementos del bottom track entre los fijos válidos de ambos lados, repartiendo el cierre linealmente en el tiempo (compensación de poligonal). Un tramo al principio o al final de la traza cuelga del fijo válido más cercano.
4. **El GGA válido no se toca.** Una transecta limpia da el mismo perfil que en v6.6, byte a byte.

La calidad GGA (1/2/4) **no se filtra**: el default de QRev (≥ 2) dejaría sin posición a casi toda la campaña. Se informa en el log y en `survey_index.csv`.

Modos:

| `position-fill` | Qué hace |
|---|---|
| `bt` (por defecto) | Filtro de HDOP y de saltos; lo inválido se reconstruye con bottom track |
| `interp` | Solo filtro de HDOP; lo inválido se interpola linealmente en el tiempo |
| `off` | Sin filtro; reproduce v6.6 exacto |

Si `Setup.coordinateSystem` no es ENU (2), o el bottom track no sirve, se pasa a `interp` con un `[warn]`.

Los ensambles sin fix (0/0) ya no se descartan: se reconstruyen con bottom track como cualquier otro fijo inválido. `gps-max-gap` sigue rechazando la transecta que tiene demasiados.

El filtro corre también sobre los **recorridos** (`tracks`). En el **pre-escaneo** de ubicación solo corre el filtro de HDOP, porque ahí no hay bottom track a mano.

### 2. QA: retroceso de s durante el cruce

Es el mayor retroceso de s en el sentido de avance, sobre los ensambles Step 3 y con las posiciones finales. Se marca `[warn]` por encima de `s-fold-tol` (2 m). Si difiere del GGA sin filtrar, el log muestra los dos valores.

Un bote que retrocede sobre el mismo fondo no hace daño. Un retroceso grande significa que puede haber dos profundidades en el mismo s.

> En el mensaje del 28/9 había propuesto 1 m. Con las posiciones corregidas, 122618 todavía retrocede 1,25 m al comenzar el cruce, y es real: el bottom track también lo ve, y el fondo es plano ahí. Por eso el umbral queda en 2 m.

### 3. QA: deriva en los bordes

Para cada margen, mide la dispersión de las posiciones de los ensambles de borde (Step 2/4), transversal y a lo largo de la línea de sección. Mide cuánto se movió el bote mientras debía estar quieto; la oblicuidad del cruce no entra. Se marca `[warn]` por encima de `edge-drift-tol` (5 m). Es solo un aviso: la penalización perpendicular ya saca del ajuste los retornos alejados.

### 4. Salidas

- `plan_view_map.png`: si se reconstruyó algo, muestra el GGA registrado (gris a trazos), los fijos descartados (×), la línea que une cada fijo descartado con su posición nueva y las posiciones reconstruidas (violeta).
- `process_log.txt`: una línea `[info] GNSS` (calidad, HDOP mín/med/máx, satélites), otra `[info]` o `[warn] posiciones` (cuántos ensambles se reconstruyeron y por qué, y el corrimiento máximo) y las líneas `[QA]` de retroceso y deriva.
- `survey_index.csv` suma estas columnas al final (las anteriores no cambian de posición):

| Columna | Contenido |
|---|---|
| `gga_q` | Cantidad de ensambles por calidad GGA (`1:93 2:22`) |
| `pos_hdop` | Fijos invalidados por HDOP |
| `pos_jump` | Fijos invalidados por salto GGA-BT |
| `pos_rebuilt` | Posiciones reconstruidas en total (incluye las que no tenían fix) |
| `pos_shift_max_m` | Mayor distancia entre el GGA registrado y la posición reconstruida |
| `s_fold_m` | Retroceso de s en el cruce |
| `edge_drift_m` | Mayor deriva en los bordes (transversal o a lo largo) |

  En las filas de grupo (aforos), los conteos se suman y los máximos se toman entre repeticiones.
- `tracks_index.csv` suma `pos_reconstruidas` y `pos_corrimiento_max_m`.

### Claves nuevas (CLI e INI)

| Clave | Por defecto | |
|---|---|---|
| `position-fill` | `bt` | `bt` · `interp` · `off` |
| `gga-hdop-max` | 4 | 0 = sin filtro |
| `gga-hdop-change` | 2 | 0 = sin filtro |
| `gga-bt-tol` | 1,5 m | 0 = sin control de saltos |
| `gga-bt-window` | 30 s | |
| `s-fold-tol` | 2 m | |
| `edge-drift-tol` | 5 m | |

Un INI v6.6 funciona sin cambios y toma estos valores por defecto.

## Validación

**Sintética** (respuesta conocida): traza recta con sesgo autónomo de 2,5 m, deriva lenta y ruido de 0,15 m. Se inyecta un salto de 5 m durante 6 s (los primeros 4 con HDOP 12) y, a continuación, otro de 3,2 m durante 4 s con HDOP normal. Se detectan los 10 ensambles; el error de las posiciones reconstruidas es de 0,16 m como máximo (sin filtro: 5,1 m). Una traza limpia no se modifica: 0 posiciones reconstruidas y los fijos válidos idénticos.

**Cinco `.mat` de la campaña 2026-04-15, sin red ni pelo de agua:**

| Perfil | Reconstruidas | Corrim. máx. | Ancho v6.6 → v6.7 | Prof. máx. v6.6 → v6.7 | Perfil |
|---|---|---|---|---|---|
| 20260415121329 | 4/68 (HDOP) | 0,13 m | 57,14 → 57,15 m | 2,056 → 2,052 m | cambio mínimo |
| 20260415122618 | 20/63 (8 HDOP + 12 saltos) | 10,14 m | 31,02 → 31,95 m | 2,561 → 2,695 m | **corregido**: V suave como en RSL; VB máx. 2,74 m |
| 20260415125906 | 0/90 | — | 48,18 → 48,18 m | 1,600 → 1,600 m | idéntico |
| 20260415131650 | 0/115 | — | 63,87 → 63,87 m | 1,965 → 1,965 m | idéntico; `[warn]` de deriva en el borde inicial (21,2 m) |
| 20260415135202 | 25/98 (24 HDOP + 1 salto) | 2,26 m | 49,11 → 48,50 m | 3,153 → 3,158 m | margen de inicio sin el falso escalón |

Con `position-fill = off`, `bathymetric_profile.csv` y `raw_bed_points.csv` de los cinco salen byte a byte iguales a v6.6.

Los umbrales se probaron en un rango: HDOP-change de 1 a 3, `gga-bt-tol` de 1 a 2 m y ventana de 15 a 30 s. Con los valores por defecto, los archivos limpios no se tocan. Con HDOP-change = 1, 125906 perdía 11 fijos sanos, y con una ventana de 15 s 122618 empeoraba.

## Pendiente

- Correr la campaña completa (64 archivos) y revisar las filas con `pos_rebuilt > 0` o `edge_drift_m > 5`.
- 135202 cambió la forma de la margen de inicio. Conviene mirarla contra el punto de borde medido en campo, si lo hay.
