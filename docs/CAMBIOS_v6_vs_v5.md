# Cambios de v6 respecto de v5

`process_adcp_bathimetric_v6.py` — abril 2026

El núcleo numérico validado en v4/v5 (extracción del `.mat`, geometría de
huellas de haz, mediana ponderada por bin, rampas de margen, exportaciones por
transecta, red multi-río `RiverNetwork`, pelo de agua multi-estación) se
mantiene. v6 agrega cuatro capacidades y corrige dos sesgos, todo validado
contra una campaña real del Río Neuquén inferior del 15/04/2026 y contra el
código de referencia de USGS QRev (`QRevPy`, `Classes/TransectData.py`,
`Classes/DepthData.py`, `Classes/DepthStructure.py`).

---

## 1. Agrupación y promediado de aforos

**Problema.** En un aforo se repite la transecta un número par de veces (mínimo
4, alternando margen de partida). v5 entregaba cuatro perfiles separados en la
misma progresiva.

**Solución.** Las repeticiones se detectan y se promedian en un único perfil.

### Criterio de agrupación

Cuatro condiciones simultáneas:

1. mismo río y mismo brazo;
2. `|Δ progresiva| ≤ --group-tol` (default 15 m). Sin `--centerline` se usa la
   componente **perpendicular al eje de sección** del desplazamiento entre
   centroides, que es la distancia río arriba/abajo real. No se usa la
   distancia recta: una repetición que cubrió un tramo algo distinto de la
   *misma* recta de sección tiene un offset a lo largo de la sección que no
   debe romper el grupo;
3. ejes de sección paralelos dentro de `--group-angle-tol` (default 25°, módulo
   180°);
4. **contigüidad temporal**: el grupo debe ser una corrida ininterrumpida en la
   secuencia ordenada por tiempo.

La condición 4 reemplaza la ventana temporal absoluta que se evaluó primero.
Una ventana no discrimina: secciones consecutivas separadas 50 m suelen estar a
menos de 5 min una de otra, así que cualquier ventana lo bastante amplia para
contener un aforo real también las contiene. La contigüidad no necesita
parámetro y cubre el caso que la progresiva sola no ve: volver a pasar por el
mismo lugar horas después con otro nivel, que **no** debe fusionarse porque el
pelo de agua cambió.

`--group-max-span` (minutos, default 0 = desactivado) queda como red de
seguridad opcional.

La sección `[grupos]` del INI permite forzar la agrupación a mano; lo que se
declara ahí tiene prioridad y el resto pasa por el criterio automático.

### Método de promediado

Se **fusionan las nubes de haces** de las repeticiones sobre un eje común y se
corre `build_profile` **una sola vez**. No se promedian los perfiles ya
construidos.

Motivo: cada punto se proyecta desde su `(E, N)` real, así que repeticiones con
distinto ancho, distinta margen de partida y distinta trayectoria se combinan
sin ninguna decisión de remuestreo ni de alineación. En un aforo sintético de 4
repeticiones sobre un lecho conocido, la fusión bajó el RMSE de 0.046 m (media
de los perfiles individuales) a 0.028 m, un 40 %.

### Eje del grupo

Dirección del vector caudal unitario **agrupado** de todos los ensembles
in-transect de todas las repeticiones (§4). No es un promedio de ángulos por
transecta: una vez excluidos los ensembles de margen, una repetición con más
ensembles aportó más travesía real sobre el mismo campo de flujo, y esa
información se conserva. El caso patológico —una repetición que divagó a otro
campo de flujo— se detecta y se reporta, no se diluye.

### Salidas

- Las carpetas por transecta individual se siguen escribiendo intactas.
- Se agrega `<grupo>/` con perfil, sección, shapefiles y
  `grupo_repetibilidad.png` (las repeticiones superpuestas, el perfil promedio y
  la banda ±1σ).
- **Todos los agregados de campaña** (`survey_index.csv`, perfil longitudinal,
  planta, shapefiles, nube consolidada) ven **sólo el perfil de grupo**.
- `_resumen/grupos_qc.csv`: σ nodo a nodo, RMS de cada repetición contra el
  promedio, dispersión de progresiva y de acimut, y el **caudal del aforo**
  (media, desvío, CV, rango y valor de cada repetición).

### Alineación de las repeticiones para el QC

Las repeticiones se comparan proyectando los nodos de cada perfil sobre el
**eje de la sección fusionada**, usando sus `(E, N)` reales, igual que se
proyectó la nube. Un nodo que está en un punto del terreno cae en el mismo `s`
en ambos, así que la comparación es punto a punto.

Una versión anterior normalizaba por **fracción de ancho**, o sea estiraba cada
repetición hasta el ancho del grupo. Con anchos de 44.0, 43.5, 42.5 y 45.5 m
contra 46.0 m del promedio, eso corría los rasgos hasta 1.15 m en horizontal, y
sobre un talud ese corrimiento se convierte en varios decímetros de cota
aparente: bastaba para que el perfil promediado se dibujara **fuera** de la
envolvente de las cuatro repeticiones de las que salió. El perfil estaba bien;
la comparación no. Con la alineación geométrica, sobre el aforo real los nodos
fuera de la envolvente pasaron de 63/86 a 9/83 y el excedente máximo de 0.26 m
a 0.046 m.

Un pequeño excedente residual es legítimo y no un error: el perfil del grupo es
la mediana ponderada de la nube agrupada y después se suaviza, no es la media
aritmética de cuatro perfiles ya construidos, cada uno suavizado por su cuenta.

### σ del lecho vs σ de las rampas

Se reportan dos repetibilidades. La **σ de la zona común** (donde las N
repeticiones se solapan) es la del lecho efectivamente medido. La **σ sobre toda
la grilla** incluye las rampas de margen, y ahí domina cuánto se acercó cada
pasada al banco, no el desacuerdo sobre el fondo: la nube agrupada llega más
lejos que cualquier pasada individual, así que su rampa arranca donde algunas
repeticiones no tienen con qué compararse.

---

## 2. Referencia de profundidad al modelo QRev

v5 tenía un booleano opaco (`--blend-beams`) más `--vb-min`. v6 adopta el
vocabulario de QRev:

| flag | valores | default | qué hace |
|---|---|---|---|
| `--depth-ref` | `auto`, `vb`, `bt` | `auto` | referencia primaria. `auto` respeta `Setup.depthReference`, o sea lo que eligió el operador en campo |
| `--composite` | `on`, `off` | `on` | rellena los nodos que la primaria no alcanzó con la otra fuente |
| `--bt-avg` | `idw`, `simple` | `idw` | cómo colapsan los 4 haces a una profundidad |
| `--bt-geometry` | `footprints`, `ensemble` | `footprints` | huellas reales de los 4 haces, o un punto promediado en la posición del bote (lo que hace QRev) |
| `--edge-anchor` | `ref`, `cloud` | `ref` | a qué se anclan las rampas de margen |

`--blend-beams` y `--vb-min` se siguen aceptando como alias.

**IDW, no media simple.** QRev pondera los 4 haces por el inverso de su alcance
desde el transductor, así los retornos más cortos pesan más; en un lecho en
pendiente el haz largo mira ladera abajo y sesga la media hacia profundidades
mayores. Portado de `DepthData.average_depth`.

**Regla de ≥2 haces.** Con menos de 2 haces válidos el promedio BT se descarta y
pasa a la alternativa (`DepthStructure.composite_depths`).

**Código de fuente por nodo.** `build_profile` devuelve además `src_grid` con el
código QRev de la fuente efectivamente usada (1 = BT, 2 = VB, 3 = DS,
4 = interpolado), y el log reporta el reparto. Permite auditar qué parte del
perfil vino del nadir, cuál de los laterales y cuál es relleno.

**`build_profile` generalizado.** La regla `vb_priority` de v5 pasó a operar
sobre *la referencia seleccionada* en vez de sobre el VB fijo. Con
`--depth-ref vb` (el default heredado) el resultado es equivalente a v5.

---

## 3. Filtro de picos de profundidad

v5 **no tenía ningún** rechazo de outliers: una detección de fondo espuria
entraba directo a la nube y corría el bin.

v6 agrega `--depth-filter {smooth, off}`, default `smooth`. Inspirado en
`DepthData.filter_smooth` de QRev: residuos contra un suavizado robusto,
contrastados con un IQR corrido, y umbral de rechazo igual al **máximo** entre
el criterio IQR, el 5 % de la profundidad medida y 0.10 m — así nunca rechaza
por ruido menor al que el instrumento puede resolver. Cada haz se filtra por
separado, como en QRev: un retorno malo en un haz no debe descartar los otros
tres.

Sobre la transecta `20260415104958.mat` rechaza 1 muestra VB de 91 y ninguna
lateral: conservador.

---

## 4. Dirección de la sección: ponderación y ensembles de margen

`compute_flow_direction` ya era una **suma vectorial**, así que las velocidades
mayores siempre pesaron más. v6 agrega dos correcciones y un filtro.

**Densidad.** Los ensembles no están equiespaciados; una media plana
sobre-representa donde el bote se frenó.

**Profundidad.** Un ensemble sobre el thalweg mueve mucho más agua que uno sobre
una planicie a la misma velocidad.

Las dos juntas dan `peso = intervalo × profundidad`, y el resultado es la
dirección del **vector caudal unitario total de la sección**, que es exactamente
la cantidad a la que una sección transversal debería ser perpendicular.
Controlado por `--flow-weighting {discharge, density, none}`, default
`discharge`; `none` reproduce v5.

**Ensembles de margen excluidos.** RiverSurveyor marca cada ensemble con
`System.Step`: 3 = travesía, 2 y 4 = margen (bote detenido midiendo el borde).
QRev usa `Step == 3` como conjunto in-transect. v5 no lo miraba. Medido sobre
`20260415104958.mat`:

| tramo | ens. | extensión | densidad | \|V_bote\| |
|---|---|---|---|---|
| Step 2 (margen, inicio) | 9 | 1.46 m | 6.18 ens/m | 0.185 m/s |
| Step 3 (travesía) | 77 | 41.91 m | 1.84 ens/m | 0.581 m/s |
| Step 4 (margen, fin) | 5 | 0.31 m | 16.02 ens/m | 0.118 m/s |

En ese archivo el sesgo que introducen es de sólo +0.04°, pero replicando los
ensembles de margen para simular una demora en campa (motor trabado, salida
demorada) el error satura en +0.47°: acotado por la diferencia angular entre el
flujo de margen y el de centro de cauce, que en un río con recirculación fuerte
puede ser de decenas de grados. `--use-edge-ensembles` los reincorpora.

---

## 5. Peso por densidad en la nube

`--density-weighting {on, off}`, default `on`. Cada ensemble recibe un peso
proporcional al intervalo que representa a lo largo del recorrido (media
distancia a sus vecinos, normalizado a media 1), multiplicado con el peso de haz
y la penalización de offset ya existentes.

Los ensembles de margen **no se descartan** de la nube: traen información de
fondo cerca del banco, que es donde escasea. Lo que se les saca es el peso
espurio por amontonamiento. Un ensemble a 0.1 m del anterior pesa ~1/6 de uno a
0.59 m. Resuelve además el bote que se frena en medio del cauce, y preserva el
beneficio batimétrico de la divagación (más puntos), que la penalización
perpendicular ya degrada suavemente según el offset.

---

## 6. QC de orientación de la sección

Nuevo `section_orientation_qc`, y `RiverNetwork.tangent_azimuth` para
soportarlo. Reporta por transecta:

- **recorrido vs eje**: qué tan oblicuamente cruzó el bote respecto del plano de
  sección, y el porcentaje de ancho proyectado que se pierde (`= 1 − cos θ`).
- **sección ⟂flujo vs ⟂traza**: el ángulo entre la perpendicular al flujo y la
  perpendicular al eje del río en esa progresiva.

Es el diagnóstico para "mis secciones no se ven transversales al río". Si el
desvío es sistemático y de un solo signo en toda la campaña, sospechar del
cálculo o de la digitalización de la traza. Si dispersa y se correlaciona con
curvas, islas y bifurcaciones, el flujo realmente es oblicuo ahí y las secciones
están bien.

`--section-orientation {flow, centerline}` queda en `flow` (comportamiento v5).
Cambiar el default movería todos los resultados previos, así que se deja al
usuario decidir con los números del QC sobre sus propias campañas.

---

## 7. Dos supuestos verificados, sin cambio de código

Ambos se habían señalado como posibles bugs y resultaron **correctos en v5**.

**Mapeo de márgenes.** QRev asigna `Edges_0__DistanceToBank` a la margen
**izquierda** y `Edges_1` a la **derecha**, siempre, sin mirar `startEdge`; ese
campo sólo decide qué código de `System.Step` cuenta como ensembles de cada
borde. Comprobado sobre `20260415104958.mat` (`startEdge = 1`, `Edges_0 = 2.0`,
`Edges_1 = 5.0`): la geometría confirma que el bote arrancó por la margen
derecha, pero eso no implica que `Edges_0` sea la margen de partida. El mapeo de
v5 es correcto.

**Calado del transductor.** `Setup.sensorDepth = 0.11 m`. El comentario de
`DepthData.depth_orig_m` en QRev dice explícitamente que las profundidades del
archivo *ya incluyen* el draft; `draft_in` se guarda sólo para poder cambiarlo
después por diferencia. v5 hace bien en no sumarlo. v6 lo lee igual porque el
promedio IDW necesita el alcance desde el transductor.

De paso queda confirmado que `depthReference = 0` es VB: `Summary.Depth` es
idéntico a `VB_Depth` en los 91 ensembles, y QRev hace
`if depthReference < 0.5: selected = 'vb_depths'`.

---

## 8. Caudal del aforo

`Summary.Total_Q` del `.mat` es acumulado a lo largo de la travesía, así que su
último valor es el caudal medido de la transecta; coincide con la suma de las
cinco componentes (Top / Middle / Bottom / Left / Right), que se usa como
control de consistencia. Se conservan las estimaciones de margen que hizo
RiverSurveyor: es el caudal medido, no un recálculo.

- `survey_index.csv` gana la columna `q_m3s` por transecta o grupo.
- `grupos_qc.csv` gana media, desvío, CV, mínimo, máximo y el valor de cada
  repetición.
- El título de `grupo_repetibilidad.png` muestra Q y su CV.

El **coeficiente de variación** es el número a mirar. La práctica USGS es que
cada transecta caiga dentro de ~5 % de la media; por encima de eso la medición
suele necesitar más repeticiones o tiene un problema (lecho móvil, referencia de
track inestable, transectas demasiado cortas). Se reporta como QA con esa
referencia.

Ojo: es **independiente** de la repetibilidad batimétrica. Un lecho estable con
un caudal ruidoso es perfectamente posible, y los dos números se leen por
separado.

---

## Validación

Campaña real del 15/04/2026, Río Neuquén inferior, 7 transectas y el eje
multi-río `R_LI-NE-NE.shp` (Limay 65.7 km, Neuquén 75.2 km con 3 brazos,
Negro 3.3 km), EPSG:5344.

Agrupación detectada: **1 grupo de 4 repeticiones** (`20260415110139`,
`110321`, `110522`, `110702`) más 3 transectas sueltas.

| | resultado |
|---|---|
| dispersión de progresiva dentro del grupo | 1.60 m |
| dispersión de acimut entre repeticiones | 4.2° (tol 25°) PASS |
| nube fusionada | 2057 puntos, 413 de la referencia primaria |
| ancho del perfil promediado | 46.00 m |
| σ media del lecho medido (zona común, s = 3.5–44.5 m) | 0.107 m |
| σ media incluyendo rampas de margen | 0.132 m |
| RMS de cada repetición contra el promedio | 0.223 / 0.356 / 0.472 / 0.093 m |
| **caudal del aforo** | **59.132 m³/s** |
| desvío / CV del caudal | 7.885 m³/s / **13.3 %** |
| caudales por repetición | 55.006 / 64.018 / 50.229 / 67.276 m³/s |

`20260415111031` quedó **fuera** del grupo, correctamente: está 42.9 m río abajo
(km 59+174.17 contra 59+131.29) y sus distancias a margen son 1.0/3.0 m contra
3.0/4.0 m de las cuatro repeticiones.

Comparación de referencias sobre ese mismo aforo:

| `--depth-ref` / `--bt-geometry` | ancho | σ media | σ máx |
|---|---|---|---|
| `auto` (= `vb`) / `footprints` | 46.00 m | 0.0704 m | 0.1862 m |
| `bt` / `footprints` | 47.00 m | 0.0620 m | 0.1709 m |
| `bt` / `ensemble` | 46.00 m | 0.0767 m | 0.2809 m |
| `vb` / `footprints` | 46.00 m | 0.0704 m | 0.1862 m |

`auto` reproduce exactamente `vb`, confirmando la lectura de
`Setup.depthReference`. Ojo con la interpretación: **la repetibilidad no es
exactitud**. `bt/footprints` repite mejor porque promedia 4× más puntos, lo que
da una respuesta más suave, no necesariamente más correcta; el nadir sigue
siendo la medición físicamente más confiable del fondo. `bt/ensemble` empeora la
σ máxima porque al colapsar los haces se pierde la cobertura espacial cerca de
las márgenes.

Todos los tests de aceptación siguen pasando (perpendicularidad al flujo,
margen izquierda en s = 0, perfil sin NaN).
