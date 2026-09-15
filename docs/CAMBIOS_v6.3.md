# process_adcp_bathimetric.py — cambios en v6.3

Versión corta: un cambio de fondo en cómo se define el ancho de una sección, y
dos de legibilidad en la vista en planta.

**Cambia números en toda la campaña**, no sólo en los aforos. Ver *Qué cambia*
al final antes de reprocesar.

---

## 1. El ancho de la sección (`edge-extent`)

### El problema

Antes de trazar las rampas de margen hay que decidir hasta dónde llega el
**lecho medido**. v6.2 tomaba el punto más externo de toda la nube, o sea un
máximo sobre un conjunto. Un máximo está sesgado hacia afuera y **crece con la
cantidad de muestras**.

En un aforo las muestras son las repeticiones: con 4 pasadas hay 4
oportunidades de que alguna se haya arrimado más a la margen, y el ancho del
grupo terminaba siendo el de la pasada más audaz. Con 6 habría sido peor. Eso
hacía que el perfil promediado fuera sistemáticamente **más ancho que
cualquiera de las secciones que promedia**, lo cual no es promediar.

Dentro de una transecta sola pasa lo mismo sobre los 5 haces. Cada haz deja de
ver fondo en un ensamble distinto, y además la huella de un haz lateral cae a
`r = profundidad · tan(25°)` del bote: el haz de proa gana siempre el máximo y
corre el extremo `+r` hacia afuera. En la margen, con 0,5 a 1,5 m de agua, son
0,2 a 0,7 m por lado.

### La regla nueva

```ini
edge-extent = mean      ; mean (6.3) | max (comportamiento 6.2)
```

Con `mean`, para cada combinación **(repetición × haz)** se toma el punto válido
más externo de cada lado, y se promedian esas estimaciones. Una sola definición
que sirve igual para una transecta suelta (una "repetición") y para un aforo.

Dos propiedades que importan:

- **No depende de cuántas pasadas se hicieron.** El promedio de los alcances es
  el mismo con 2 repeticiones que con 6.
- **Los corrimientos de huella se cancelan solos.** El bote cruza a lo largo de
  la sección, así que los haces de proa y popa se corren `+r` y `−r` en `s`, y
  los de babor y estribor se corren perpendicular, casi sin tocar `s`. Al
  promediar, el sesgo desaparece sin tener que elegir entre la posición del bote
  y la de la huella.

Un haz necesita al menos `EXTENT_MIN_PTS` (3) puntos válidos para aportar una
estimación. Si ninguna combinación califica, se cae al criterio viejo.

Los puntos que queden fuera del rango promediado **no se pierden**: siguen en
`raw_beam_points`, simplemente no extienden la sección.

### Lo que no cambia

- `edge-anchor = ref | cloud` sigue decidiendo si el extremo se mide sobre la
  referencia primaria o sobre toda la nube. Es ortogonal a `edge-extent`.
- `Edges_0/1__DistanceToBank` se sigue usando igual, y sigue siendo la mediana
  entre repeticiones en un grupo.
- La **forma** de la margen (triangular / rectangular, `Setup.Edges_*`) sigue sin
  leerse: la rampa es siempre recta a profundidad 0, que equivale a asumir
  margen triangular. Queda para una versión posterior.

### Trazabilidad

- El log de cada transecta y de cada grupo dice la extensión adoptada y la
  envolvente, para que la diferencia sea visible:

  ```
  [info] extensión medida (promedio de 20 combinaciones repetición×haz):
         s = -71.88 a 71.88 m (envolvente: -75.40 a 75.40 m)
  ```

- En un grupo se informa además la dispersión del extremo entre repeticiones, y
  se emite `[warn]` si una pasada cubrió más de 5 m que otra.
- `grupos_qc.csv` suma seis columnas: `s_min_medio_m`, `s_max_medio_m`,
  `disp_s_min_m`, `disp_s_max_m`, `s_min_por_repeticion_m`,
  `s_max_por_repeticion_m`. Si una pasada se aparta, ahora es una columna y no
  algo que había que deducir de la forma de la banda gris.
- `merge_group_clouds` agrega un array `rep` con la repetición de origen de cada
  punto, que es lo que hace posible el cálculo por pasada.

---

## 2. Etiquetas de la vista en planta (`planview-labels`)

Con 56 secciones, v6.2 intentaba colocarlas todas y sólo descartaba una cuando
no encontraba lugar libre. Resultado: los tramos abiertos no ganaban nada y los
codos quedaban ilegibles, con una maraña de guías cruzadas.

```ini
planview-labels = auto      ; auto | all | every:N | none
```

`auto` ralea **por densidad en píxeles**, no "una de cada N". Se recorren las
secciones en orden de progresiva y se etiqueta una cada ~30 px; las salteadas
quedan con un tick y siguen en `survey_index.csv`. Una sección que abre un
racimo después de un hueco se etiqueta siempre, para que ningún grupo de
secciones quede sin ninguna referencia.

La ventaja sobre "1 cada 3" es que se adapta a la escala dibujada: los tramos
sueltos conservan todas sus etiquetas, los racimos se ralean, y un recorte de un
tramo muestra más etiquetas sin tocar la configuración. Efecto secundario: al
haber menos etiquetas apretadas, muchas menos necesitan guía.

`all` reproduce v6.2.

---

## 3. Ejes en miles

Las coordenadas POSGAR se dibujaban completas: `2.582.000`, siete dígitos para
distinguir tramos de 2 km. Ahora se dividen por mil y se marca `×10³ m` en el
vértice de cada eje.

- Los decimales salen del **paso de los ticks**, no son fijos: con 2 km entre
  ticks quedan enteros; con 50 m harían falta dos decimales.
- Por debajo de 2 km de extensión se vuelve a metros lisos, donde el `×10³`
  estorba más de lo que ayuda.
- El valor escalado se imprime **sin separador de miles**: `2610`, no `2.610`,
  que al lado de un `×10³` se lee como 2,61.

Aplica a `survey_plan_view.png` y a `tracks_plan_view.png`.

---

## Claves nuevas

| Clave (INI / CLI) | Default | Qué hace |
|---|---|---|
| `edge-extent` | `mean` | `mean` (6.3) \| `max` (regla v6.2) |
| `planview-labels` | `auto` | `auto` \| `all` \| `every:N` \| `none` |

## Salidas que cambian

- `grupos_qc.csv`: seis columnas nuevas al final del bloque de σ.
- `survey_plan_view.png` y `tracks_plan_view.png`: ejes en miles con `×10³ m`;
  etiquetas raleadas.
- `survey_index.csv`: `width_m` cambia de valor. La lista de perfiles no.

## Qué cambia respecto de v6.2

**Esto sí mueve números en toda la campaña.** Los anchos y las rampas de margen
cambian en todas las transectas, no sólo en los grupos:

- **Aforos:** la corrección grande. El ancho baja del de la pasada más externa
  al promedio de las pasadas.
- **Transectas simples:** corrección chica, de centímetros a un par de
  decímetros por margen, según la profundidad del último retorno válido. Con
  `depth-ref` resolviendo al haz vertical y `edge-anchor = ref`, el extremo ya
  salía de puntos bajo el bote, así que ahí el cambio es sólo de máximo a
  promedio sobre los ensambles.
- **Profundidades en la zona medida:** no cambian. Lo que cambia es dónde
  terminan los datos y dónde arrancan las rampas.

Para separar efectos al reprocesar, `edge-extent = max` reproduce v6.2 exacto.
Corriendo las dos y comparando con `compara_corridas.py` ves qué movió cada
cosa.

---

## Validación

`measured_extent` sobre una nube construida con 4 repeticiones que alcanzan
distinta distancia de la margen (70,0 / 72,0 / 74,5 / 71,0 m por lado) y 5 haces
con corrimiento de huella de ±0,9 m:

| Caso | `max` (v6.2) | `mean` (v6.3) | Real |
|---|---|---|---|
| Grupo de 4 repeticiones | 150,80 m | **143,75 m** | 143,75 m |
| Transecta sola, 5 haces | 141,80 m | **140,00 m** | 140,00 m |

El `mean` del grupo reproduce al centímetro el promedio de los anchos de las
cuatro pasadas. El `max` los supera en 7 m, y ese exceso crecería con una quinta
repetición.

Raleo de etiquetas sobre un trazado de 43 secciones con tramos abiertos y dos
racimos: `auto` deja 20 etiquetadas y 23 con tick; `all` las 43; `every:3`
deja 15.

El resto del pipeline sin cambios: la campaña sintética de v6.2 corre completa,
con las mismas progresivas, el mismo pelo de agua y las mismas capas en
EPSG:5344.
