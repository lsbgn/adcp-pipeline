# process_adcp_bathimetric.py v6.4 — cambios respecto de v6.3

Origen: en varias secciones de 2026-04-15 el lecho dibujado se apartaba de los
puntos del haz vertical hasta 0,3 m. Las dos transectas revisadas
(20260415124913 y 20260415124412) tenían `depth reference: VB` y 39 nodos VB,
así que no era un problema de referencia: era el estimador.

## 1. Por qué se apartaba

Tres causas, todas dentro de `build_profile()` de v6.3:

1. **La ventana y el suavizado colgaban de `dx`.** `bin-half-width` valía `2*dx`
   y la ventana del Savitzky-Golay era de 5 nodos. Con `dx = 1.0` (el valor del
   INI de campaña) eso son 4 m de mediana más 4 m de ajuste cuadrático: se
   recortan el thalweg y los resaltos de menos de ~5 m.
2. **La mediana ponderada ignora dónde cae cada punto dentro de su ventana.**
   Sólo mira profundidades y pesos, así que cuando la velocidad del bote varía
   (los pesos por densidad varían con ella) el nodo toma la profundidad del lado
   "pesado" y el lecho se corre horizontalmente.
3. **El Savitzky-Golay corría sobre toda la grilla**, incluidos los ceros
   forzados de las orillas. En el empalme entre lecho medido y rampa el ajuste
   cuadrático sobrepasa: ése es el pozo de 20260415124412 en s = 38 m.

Y una cuarta, menor: el ancho salía cuantizado a `dx`, porque el último nodo
caía en un múltiplo de `dx` (de ahí los anchos de 48,0 y 40,0 m exactos).

## 2. Qué cambia

### `bed-fit = loess` (nuevo default; `median` reproduce v6.3 exacto)

En cada nodo, regresión lineal local ponderada sobre los puntos de la
referencia primaria:

- núcleo tricúbico de semiancho `bin-half-width` = **0,75 m fijos**,
  independiente de `dx`; sólo se ensancha hasta el doble para juntar 3 puntos;
- recta con 2 puntos dentro de un semiancho de los extremos del lecho medido,
  donde la ventana es de un solo lado (si no, el promedio local tira hacia
  adentro: probado, 0,29 → 0,10 m en un caso);
- dos iteraciones robustas LOWESS (bicuadrada), pero **sólo donde la ventana
  tiene 5 puntos o más**: con dos o tres no hay forma de distinguir un pico de
  un thalweg angosto. La variante "dejando uno afuera" se probó y se descartó:
  cortaba thalwegs reales (máximo contra el lecho real 0,06 → 0,35 m);
- relleno compuesto con la otra fuente donde la primaria no llega, e
  interpolación lineal en lo que quede;
- **sin Savitzky-Golay global**: las rampas se pegan después del ajuste, así que
  el empalme con la margen ya no puede deformar el lecho medido.

### Márgenes exactas y forma declarada en RSL

- Nodos exactos en los extremos del lecho medido y en las dos orillas; el ancho
  deja de ser múltiplo de `dx`.
- `edge-shape = auto` lee `Setup.Edges_0__Method` / `Edges_1__Method`
  (2 triangular, 1 rectangular, 0 Q de usuario; mapeo verificado contra el
  lector SonTek de QRevPy). Una margen rectangular mantiene la profundidad del
  extremo hasta la orilla y cierra con una pared vertical de 1 cm de ancho en
  planta, para que `s` siga siendo estrictamente creciente.
- Los nodos de margen llevan código de fuente 5 (`SRC_EDGE`).

### Figura de la sección

Lecho medido en trazo lleno, márgenes extrapoladas en trazo discontinuo, tramos
interpolados punteados, y la leyenda dice cuál es la referencia primaria.

### Pitch y roll en las huellas (`pitch-roll = off | on | inv`)

Cada haz, el vertical incluido, se inclina con la actitud del ensamble
(R = Rz(rumbo)·Ry(pitch)·Rx(roll)) antes de ubicar su huella. Con `off`
reproduce la geometría nivelada de v6.3 hasta 1e-14 m.

Queda en **off** por default: SonTek no documenta el signo de pitch/roll del M9,
y con el signo cambiado el error de la huella se duplica en vez de anularse.

### Escalas hidrométricas en `survey_plan_view.png`

Cuadrado rojo relleno si la escala tuvo lectura en la campaña, hueco si no.
Se dibujan clipeadas a la ventana (una escala a 30 km no la agranda) y sus
rótulos se ubican antes que las progresivas, que los esquivan. El registro se
construye aunque no haya archivo de lecturas.

### Otros

- `survey_index.csv` suma al final `forma_mi`, `forma_md`, `bed_fit`.
- Vista en planta: los decimales de los ejes salen del paso de ticks exacto. Con
  paso de 0,25 km, v6.3 imprimía "5687.2" dos veces.

## 3. Validación

### Regresión contra v6.3

Campaña sintética georreferenciada (10 transectas, un aforo de 4 pasadas, un
recorrido, pelo de agua GNSS, registro de escalas): v6.4 con `bed-fit = median`
da **206 archivos idénticos byte a byte** contra v6.3. El único que cambia es
`survey_index.csv`, por las tres columnas nuevas al final.

### Contra el lecho real (misma campaña sintética, `dx = 1.0`)

| | v6.3 | v6.4 |
|---|---|---|
| RMS del lecho contra el real | 0,067 m | **0,021 m** |
| Máximo | 0,250 m | **0,064 m** |
| VB a más de 10 cm de la línea | 8,9 % | **0,3 %** |
| RMS de los VB contra la línea | 0,062 m | **0,023 m** |
| Error de posición de la orilla derecha | 0,25 m | **0,04 m** |

### Escenarios de transecta suelta (6 semillas cada uno)

| Escenario | VB > 10 cm: v6.3 → v6.4 | RMS real | Máximo real |
|---|---|---|---|
| bote a 0,5 m/s | 13,4 % → 0,3 % | 0,083 → 0,015 m | 0,29 → 0,04 m |
| bote a 1,0 m/s | 13,2 % → 0,0 % | 0,113 → 0,021 m | 0,30 → 0,06 m |
| bote a 1,5 m/s | 16,7 % → 0,5 % | 0,082 → 0,031 m | 0,25 → 0,09 m |
| ruido de 4 cm | 17,5 % → 0,6 % | 0,115 → 0,037 m | 0,30 → 0,10 m |
| aforo de 4 pasadas | 12,8 % → 0,1 % | 0,099 → 0,016 m | 0,28 → 0,05 m |
| velocidad alternada 0,2/1,2 m | 20,5 % → 0,8 % | 0,155 → 0,063 m | 0,44 → 0,28 m |

### Semiancho del núcleo

Sobre la campaña sintética: 0,75 m da RMS 0,021 m; 1,5 m da 0,026 m. Por eso el
default es 0,75 m.

### Pitch/roll

Campaña sintética generada con la convención aeronáutica, inclinación p95 3,9°,
profundidad mediana 1,3 m:

| | RMS del lecho | Máximo |
|---|---|---|
| `off` | 0,021 m | 0,064 m |
| `on` (signo correcto) | 0,018 m | 0,047 m |
| `inv` (signo cambiado) | 0,026 m | 0,085 m |

### Límite conocido

El filtro de picos estilo QRev (`depth-filter = smooth`, multiplicador 15)
descarta puntos VB reales del fondo de una fosa angosta cuando el bote pasa
rápido: en el escenario "bote rápido sobre la fosa" descarta entre 2 y 4 puntos VB del fondo en la mitad de las corridas
y el máximo contra el lecho real queda en 0,295 m (v6.3: 0,329 m). Sin filtro,
v6.4 da 0,104 m y v6.3 0,239 m. No se tocó el filtro en esta versión.

## 4. Cómo reprocesar 2026-04-15

1. En el `.ini` de la campaña agregar `bed-fit = loess` y `edge-shape = auto`
   (o copiar los bloques nuevos de `campanha_ejemplo.ini`).
2. Correr primero con `bed-fit = median`: tiene que dar exactamente lo mismo que
   la corrida v6.3 (control de que nada más se movió).
3. Correr con `bed-fit = loess` y comparar con `compara_corridas.py`. Las
   profundidades cambian en toda la sección, no sólo en los extremos: es el
   cambio de estimador.
4. `python tools\diagnostico_mat.py <carpeta de .mat> > diag_mat.txt` para
   confirmar la forma de margen que se cargó en RSL, la distribución de
   pitch/roll y la numeración de haces.
