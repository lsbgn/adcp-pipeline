# process_adcp_bathimetric.py v6.5 — cambios respecto de v6.4

Origen: el diagnóstico de los 64 `.mat` de la campaña 2026-04-15
(`diag_mat.txt`).

## 1. Numeración de haces: el haz 2 está a babor

Desde v3 el script suponía que, visto desde arriba, los haces se numeran en
sentido horario: haz 1 a proa, haz 2 a **estribor**. En el M9 es al revés.
Hay tres evidencias independientes:

| Evidencia | Resultado |
|---|---|
| Manual de RSL, "XYZ Coordinates" | XYZ dextrógiro, +X hacia el haz 1, +Z arriba: +Y es babor |
| `Transformation_Matrices` de los 64 archivos, 3 MHz y 1 MHz | el haz 2 a +90° del haz 1 en ese marco, o sea a babor, en todos |
| Profundidad de cada haz inclinado contra el perfil VB en su huella | antihorario gana en 63 de 64 archivos; RMS mediano 0,084 m contra 0,150 m del horario |

Qué estaba mal:

- En 3 MHz, los haces 1 y 3 (proa y popa) estaban bien. Los haces 2 y 4
  estaban cambiados de lado.
- En 1 MHz, los cuatro haces estaban espejados respecto de la línea
  proa-popa. En la campaña del 15/04, 1 MHz es la frecuencia dominante en 42
  de los 55 archivos con datos suficientes para separarla.
- Cada huella espejada queda hasta `2 · profundidad · tan 25°` fuera de lugar:
  0,9 m con 1 m de agua y 2,8 m con 3 m.

Dónde pega el error:

- **Perfil de sección con referencia VB:** casi no cambia. Sólo lo afectan los
  nodos de relleno compuesto cerca de las márgenes.
- **Nubes de haces:** cambian. Eso incluye `raw_bed_points`,
  `survey_raw_beam_points` y la nube de **recorridos**, que es la que alimenta
  el MDT del cauce.
- En un recorrido longitudinal, los haces laterales caían del lado equivocado
  del eje del río.

La corrección es la clave nueva `beam-layout`:

- `auto` (default) lee el sentido de numeración de la matriz de cada archivo.
  Usa el signo del producto vectorial entre haz 1 y haz 2, que no depende del
  signo global de la matriz.
- Sin matriz, usa `ccw`.
- `cw` reproduce v6.4 byte a byte.

Cada transecta registra en su log qué numeración usó y de dónde la tomó.

## 2. Actitud

### Las profundidades no se tocan

Según el manual de RSL, `VB_Depth` y `BT_Beam_Depth` ya salen verticales,
"including ... compensation for tilt". QRev también se apoya en eso.

La sección [3] del diagnóstico v1 decía "SIN CORREGIR" en varios archivos.
Estaba contaminada por la numeración espejada y por el corrimiento de huella
que no modelaba. **No usar esas líneas.** Esa sección se eliminó del
diagnóstico.

### Lo que falta es dónde cae cada huella

El 15/04, 24 de las 64 transectas tienen `|pitch|` p95 por encima de 8° (hasta
12°), y el pitch mediano es positivo en 62 de 64. En esas 24, a la profundidad
mediana, el haz vertical pega entre 0,2 y 0,55 m adelante del transductor; en
lo hondo, más.

### Claves

- `pitch-roll` sigue en `off` por default. Los signos pasan a ser dos claves
  separadas:
  - `pitch-sign = bow-up | bow-down`
  - `roll-sign = stbd-down | port-down | none` (`none`: sólo pitch)

  En v6.4 sólo se podían invertir los dos a la vez. `inv` sigue aceptándose.
- `antenna-height`: la antena GNSS montada sobre un mástil se inclina con el
  bote, y el transductor queda desplazado `altura · sen(inclinación)` respecto
  de ella.
- `gnss-lag`: retardo de la posición GNSS. Corre cada posición
  `retardo · velocidad` hacia adelante.
- Con `pitch-roll = off`, una transecta cuya inclinación corre la huella del VB
  más de 0,20 m deja un `[warn] actitud sin corregir`.

### Corrección de lectura

`Compass.Pitch` y `Compass.Roll` figuran en el manual de RSL como NS × 3. v6.4
los aplanaba a 3·NS, y en ese caso `pitch-roll` se desactivaba sin avisar.
Ahora se toma la primera columna.

## 3. Diagnóstico v2 (`tools/diagnostico_mat.py`)

- **[3] Numeración de haces.** Compara la matriz y los datos, y marca cuál es
  la regla del M9 y cuál la de v3–v6.4.
- **[4] Direcciones de haz.** Los haces se normalizan apuntando hacia abajo, así
  que desaparece la ambigüedad de 180°. Se omite la matriz del haz vertical.
- **[5] Nueva: signo de pitch/roll, altura de antena y retardo GNSS.**
  - Agrupa las pasadas repetidas sobre una misma sección.
  - Para cada combinación de signos, con altura de antena (0 a 1,5 m) y
    retardo (−2 a +2 s) ajustados, mide cuánto coinciden los perfiles VB de
    pasadas en sentidos opuestos. La inclinación corre las huellas hacia la
    proa, y al invertirse el sentido ese corrimiento se invierte sobre la
    sección.
  - Estima antes el poder de la prueba: el desacuerdo que produciría una
    inclinación real (pendiente del lecho × corrimiento de huella) contra el
    ruido entre pasadas del mismo sentido. Por debajo de 0,5 dice que no
    puede decidir.
  - Decide el signo del pitch y el del roll por separado. Si sólo el pitch
    queda determinado, recomienda `roll-sign = none`.
  - Toma el retardo GNSS sólo si mejora de verdad el modelo físico, y avisa
    cuando un retardo sin actitud explicaría casi lo mismo.
  - Termina con la línea a copiar en el `.ini`, o con la recomendación de
    dejar `off`.

## 4. Validación

Campañas sintéticas escritas como `.mat` de RSL: 4 transectas sueltas y 2
aforos de 4 pasadas alternadas, con estas características:

- geometría de haces conocida;
- 1 MHz en 2 de cada 3 ensambles;
- bote cangrejeando 30°;
- cabeceo proa-arriba que crece con la velocidad (p95 ≈ 9,5°);
- roll con desvío de 2°;
- profundidades hasta 4,4 m.

### Diagnóstico

| Verdad | Numeración | Signo detectado | Antena | Retardo |
|---|---|---|---|---|
| ccw, proa-arriba / estribor-abajo, antena 0,6 m, sin retardo | ccw 12/12 | igual a la verdad | 0,50 m | 0,0 s |
| ccw, proa-abajo / babor-abajo, sin antena, retardo 0,5 s | ccw 12/12 | igual a la verdad | 0,00 m | +0,5 s |
| cw, proa-arriba / estribor-abajo, antena 1,0 m | cw 12/12 | igual a la verdad | 0,75 m | 0,0 s |

- El signo del pitch se separa con claridad: el signo equivocado da RMS 0,093
  contra 0,050 m.
- El del roll se separa por poco (0,008 m), porque el roll es chico.
- La altura de antena sale con la resolución de la grilla. Conviene medirla.

### Pipeline contra el lecho real (primera campaña de la tabla)

| Corrida | RMS del perfil | Nube de haces inclinados: RMS | p95 |
|---|---|---|---|
| v6.4 | 0,119 m | 0,255 m | 0,621 m |
| v6.5, `beam-layout = cw` | 0,119 m | 0,255 m | 0,621 m |
| v6.5 default (ccw, actitud off) | 0,119 m | 0,106 m | 0,269 m |
| v6.5 ccw + actitud con signo correcto, antena 0,6 m | **0,040 m** | **0,030 m** | **0,058 m** |
| v6.5 ccw + actitud con la antena que estimó el diagnóstico (0,5 m) | 0,040 m | 0,030 m | 0,058 m |
| v6.5 ccw + actitud con signo invertido | 0,221 m | 0,195 m | 0,512 m |

Regresión: v6.5 con `beam-layout = cw` y `pitch-roll = off` da **245 archivos
idénticos** byte a byte contra v6.4.

## 5. Resultado sobre las campañas reales

| | 2026-04-15 | 2026-08-25 |
|---|---|---|
| Numeración: ccw mejor que cw | 64 de 64 | 62 de 62 |
| Matrices | ccw en todas | ccw en todas |
| Transectas con `\|pitch\|` p95 > 8° | 24 de 64 | 35 de 61 |
| Ocupaciones con pasadas opuestas | 4 (dos de sólo 2 pasadas) | 2 |
| RMS entre opuestas, nivelado / mismo sentido | 0,185 / 0,155 m | 0,127 / 0,131 m |
| Mejor hipótesis con signo (sin retardo) | 0,182 m | 0,127 m |

Ninguna hipótesis de actitud mejora la coincidencia entre pasadas, en ninguna
de las dos campañas. Eso **no** prueba que el pitch sea irrelevante:

- En 2026-04-15, las dos ocupaciones de 4 pasadas tienen `|pitch|` p95 de 1 a 2°
  y no pueden decir nada. Sólo el par 125238/125430 tiene cabeceo fuerte.
- En 2026-08-25 el aforo de la tarde sí tiene cabeceo (pitch mediano de 3,7 a
  8,6° según el sentido), pero las opuestas no discrepan más que las del mismo
  sentido. Lo más probable es que la sección de aforo sea casi plana, y ahí un
  corrimiento horizontal no cambia la profundidad.

La primera versión del diagnóstico no estimaba el poder de la prueba. La v2.1
lo hace y da la tabla por ocupación; hay que volver a correrla para confirmar.

Mientras tanto: **`pitch-roll = off`**. El signo se puede fijar en banco en un
minuto (ver `USO_v6.5.md` §9), y con el signo fijo ya no hace falta que los
aforos lo resuelvan.

## 6. Cómo reprocesar 2026-04-15

1. Correr el diagnóstico v2.1 sobre la carpeta de transectas y mirar la
   sección [5], incluida la tabla por ocupación y el poder de la prueba:
   ```powershell
   python .\tools\diagnostico_mat.py <carpeta de transectas> > diag_mat_v2.txt
   ```
2. Fijar el signo de pitch y roll en banco, y medir la altura de la antena
   GNSS sobre la cara del transductor en el catamarán.
3. En el `.ini`, dejar `beam-layout = auto`. Encender `pitch-roll` con el signo
   de banco (o el de la sección [5], si tuvo poder) y la `antenna-height`
   medida.
4. Reprocesar y comparar con `compara_corridas.py` contra la corrida v6.4:
   - Los perfiles cambian en toda la sección si se enciende la actitud.
   - Si no se enciende, los perfiles sólo cambian donde hubo relleno compuesto.
   - Las nubes de haces cambian siempre.
