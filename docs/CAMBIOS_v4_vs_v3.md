# Procesador batimétrico ADCP — Cambios de la v4 respecto de la v3

Resumen de las correcciones y mejoras incorporadas en `process_adcp_bathimetric_v4.py`
(internamente documentado como v7.0) frente a la v3.

---

## 1. Progresivas sobre centerline real, con brazos

- Nueva clase `RiverRoute`: arma la traza del río a partir del shapefile de eje
  (uno o varios tramos) y **detecta la topología sola**.
  - Une los tramos que comparten extremos (`linemerge`) → **camino principal**.
  - Un tramo cuyos **dos extremos caen sobre el camino principal** (tolerancia
    `BRANCH_SNAP_TOL = 15 m`) se reconoce como **brazo** (isla). Calcula dónde
    **abre** y dónde **cierra** ese brazo en progresiva.
- **Progresivas "del conjunto"**: los dos brazos comparten el mismo origen y el
  **mismo valor de progresiva en el punto de bifurcación**; aguas abajo cada uno
  acumula su propia progresiva. Cada sección se rotula además con su brazo.
- **Rótulo Norte / Sur automático** (por latitud del brazo vs. el cauce principal
  en el tramo de isla).
- **Orientación automática** de la traza (qué extremo es aguas abajo) usando la
  **pendiente del pelo de agua**: la progresiva crece hacia aguas abajo sin
  necesidad de intervención. Si no se pasa CSV de pelo de agua, se usa el sentido
  de digitalización o la bandera `--chainage-reverse`.
- Detectado para `r_nqn_inf.shp`: camino principal **75.159 m**; **Brazo Norte =
  `id=2`** (canal corto al norte de la isla), abre en **66+565,72** y cierra en
  **69+284,11**; el cauce principal en ese tramo queda como **Brazo Sur**.

> Reemplaza a la función `compute_chainage` de la v3, que sólo tomaba la línea más
> larga y no distinguía brazos.

## 2. Cotas absolutas sobre nivel del mar (IGN SRVN16)

- Nueva clase `WaterSurfaceProfile`: toma el CSV de puntos de pelo de agua
  (X, Y, cota), **proyecta cada punto a su progresiva** sobre la traza y arma el
  **perfil longitudinal del pelo de agua** (cota vs. progresiva), interpolado
  linealmente. Es la forma físicamente correcta de transportar la superficie del
  agua a lo largo del río (respeta la pendiente y los meandros).
- `bed_elev_m` pasa a ser **cota absoluta** = pelo de agua − profundidad. Se
  agregan las columnas **`ws_elev_m`** (cota de pelo de agua interpolada) y
  **`brazo`** a todas las salidas.
- **Fuera del rango relevado** de pelo de agua, extrapola con el valor del extremo
  más cercano (clamp constante, no recta) y lo **marca con un aviso** (`ws_note` /
  `[WARN EXTRAP …]`), por seguridad.
- La **sección transversal se dibuja en cota absoluta** (eje "Cota IGN SRVN16 [m]",
  superficie del agua en su cota real) cuando se provee el pelo de agua; si no,
  mantiene el modo relativo de la v3.

## 3. Modo relevamiento (multi-archivo)

- Se le pasa una **carpeta**, **varios `.mat`** o un **patrón (glob)** y procesa
  **todas las transectas de la salida**.
- Cada transecta va a su subcarpeta, más una carpeta **`_resumen/`** con los
  productos consolidados (ver el MD de uso): índice, nube total, vista en planta
  del relevamiento, **perfil longitudinal (pelo de agua + thalweg)** y shapefiles
  combinados.

## 4. Corrección del lecho: prioridad del beam vertical (VB)

- **Problema (v3):** en cada celda del perfil entraba 1 punto de VB (peso 1,0)
  contra ~4 laterales (0,5 c/u = 2,0). Por mayoría, **los laterales le ganaban al
  VB y aplanaban la sección** (pendiente de margen medida: VB −0,090 m/m vs.
  dibujada −0,065 m/m).
- **Corrección (v4, por defecto):** el **beam vertical (nadir) define el lecho
  donde existe**; los laterales sólo **rellenan huecos interiores** donde no hay
  VB. La línea del lecho pasa a seguir el VB en todo el ancho.

## 5. Corrección de las márgenes: rampa desde el último VB

- **Problema:** las márgenes se anclaban al punto **más externo de toda la nube**
  (casi siempre un lateral), y los laterales que sobresalen hacia la costa metían
  un pozo/quiebre en la margen.
- **Corrección (v4, por defecto):** las márgenes se anclan a la **extensión del
  VB**. Desde el **último VB disponible** de cada lado sale una **recta hasta la
  orilla**, a la distancia de margen del *edge* (`Setup.Edges_*`). Los laterales
  ya **no extrapolan** la costa. Es coherente con el instrumento: el *edge
  distance* del SonTek se mide desde el nadir (VB) del último ensemble a la costa.

## 6. Muestreo del perfil

- `build_profile` ahora usa **`bin_half_width = 2·dx`** por defecto
  (`bin_half_width=None` → `2*dx`), con la firma acordada.
- **`dx` por defecto pasa a 0,50 m** (antes 0,25 m).

## 7. Nuevas opciones de línea de comandos

- Posicional **`matfiles`** admite 1 o varios archivos / carpeta / glob.
- Nuevas: `--survey-name`, `--water-surface-csv`, `--ws-east-col`,
  `--ws-north-col`, `--ws-elev-col`, `--ws-crs`, `--datum-name`, `--blend-beams`
  (vuelve al lecho mezclado v3), `--vb-min` (mínimo de puntos de VB para usar solo
  VB en una celda).
- Cambian los valores por defecto de `--dx` (0,5) y `--bin-half-width` (2·dx).
- Se mantienen: `--centerline`, `--chainage-offset`, `--chainage-reverse`,
  `--offset-scale`, `--no-offset-weighting`, `--no-edge-extrapolation`,
  `--water-surface-elev`, `--outdir`.

## 8. Otros

- **Archivo de configuración INI** (`--config campania.ini`): un solo archivo por
  campaña con `[procesamiento]` (entradas, eje, pelo de agua, parámetros) y
  `[campania]` (metadatos: operador, equipo, nivel, observaciones). Evita tipear
  las opciones a mano. Las **rutas relativas se resuelven respecto del `.ini`**
  (pensado para dejarlo dentro de la carpeta de la campaña). **La línea de comandos
  tiene prioridad** sobre el archivo (CLI > config > default).
- **Registro de trazabilidad** `procesamiento.txt` (en `_resumen/` o en la carpeta
  de salida): versión del script, commit de git si está disponible, fecha,
  metadatos de campaña, entradas, parámetros usados y resumen por transecta.
- Constante **`__version__`** (v4.0), informada en el log y el registro.
- **Depurado**: se quitaron imports sin uso (`os`, `unary_union`, `nearest_points`);
  el núcleo numérico validado (perfil, ruta, pelo de agua, gráficos) se dejó
  intacto para no introducir regresiones.
- El **log** informa el modo de lecho usado (`bed=VB-priority (vb_min=1)` o
  `blended`) y el ancho de banda efectivo; se guarda un **`process_log.txt`** por
  transecta.
- Nuevos QA por transecta: **distancia del centroide al eje** del río y **aviso de
  extrapolación** de pelo de agua.
- Compatibilidad hacia atrás: sin `--centerline` no calcula progresiva; sin
  `--water-surface-csv` las cotas quedan relativas al pelo de agua (como la v3).

---

### Confirmado

1. **Brazo Norte = `id=2`** (canal corto al norte de la isla) — coincide con la
   detección automática N/S.
2. **Origen de progresivas** en el extremo NW del eje (`--chainage-offset 0`, el
   valor por defecto).
