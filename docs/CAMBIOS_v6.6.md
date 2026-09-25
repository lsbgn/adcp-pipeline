# process_adcp_bathimetric.py v6.6 — cambios respecto de v6.5

Origen: campaña 2026-09-04, Dique Ballester. Una sección recorrida al pie de las
compuertas, paralela a ellas, se dibujó atravesando la estructura, y la escala
NEU-03 dejó de anclar el pelo de agua cuando se sumó el punto Pt1 del embalse.

## 1. `[orientacion]`: sección a lo largo de la trayectoria

La sección es la normal al flujo medido. Al pie de las compuertas el flujo es un
chorro con recirculación y su dirección no es la del río: la normal giró y cortó
la obra.

```ini
[orientacion]
20260904140552 = recorrido     ; eje de la trayectoria de la embarcación
20260904141510 = 17.5          ; o un azimut fijo de la línea
```

- `recorrido`: eje principal (componentes principales) de las posiciones GNSS
  en transecta, sin los ensambles estacionarios de margen. Aviso si la traza no
  es recta (linealidad < 0,80).
- De los dos sentidos de la línea se elige el más cercano al del flujo, así que
  `s = 0` sigue siendo la margen izquierda. Si la línea queda casi paralela al
  flujo, aviso: el sentido es ambiguo.
- Todo lo que sigue usa la línea forzada: proyección de la nube, penalización
  por distancia a la sección, ubicación sobre el eje (progresiva), traza en
  planta, shapefiles.
- El control `axis perpendicular to flow` pasa a `SKIPPED` con el motivo.
- En un grupo alcanza con forzar un miembro o el nombre del grupo; se usa la
  traza de todas las repeticiones.

Validación sintética (flujo girado 35° en una transecta de 47,8 m de ancho):
sin override el ancho sale 41,2 m (la sección gira 35°); con `recorrido`,
47,9 m y azimut de la línea 16,8° contra 16,9° de la verdad.

## 2. `[saltos]`: discontinuidades del pelo de agua

Qué pasaba con NEU-03: una escala es de control (no ancla) cuando hay puntos
GNSS aguas arriba y aguas abajo en el mismo recorrido. Pt1, en el embalse, la
dejó entre puntos: pasó a control y el pelo de agua entre Pt1 y el siguiente
punto aguas abajo se interpoló a través de un salto de más de 2 m.

```ini
[saltos]
Dique Ballester = Neuquen | 29+490
```

Cada recorrido de la red se corta en tramos en los saltos declarados. Usan sólo
valores del mismo tramo:

- interpolación y extrapolación del pelo de agua;
- ajuste de la pendiente de continuación (ventana del extremo, tramo completo,
  pendiente por defecto `auto`);
- cota de confluencias;
- la decisión escala-ancla / escala-control.

Un tramo sin valores queda sin pelo de agua y la nota dice por qué; no se
extrapola desde el otro lado. `survey_long_profile.png` dibuja el salto
(línea vertical punteada violeta) y el escalón del pelo de agua.

También por coordenadas: `Compuertas = <X> <Y>` en el CRS del eje.

Validación sintética (salto de 2,00 m en km 60+000, escala G1 100 m aguas abajo,
puntos GNSS a ambos lados):

| | sin `[saltos]` | con `[saltos]` | verdad |
|---|---|---|---|
| G1 | control, Δ = −0,855 m | ancla | — |
| pelo en km 59+800 | 267,914 | 268,201 | 268,200 |
| pelo en km 60+333 | 265,855 | 265,663 | 265,667 |

## 3. Aviso de salto sin declarar

Una escala de control a más de 0,5 m del pelo interpolado entre puntos GNSS
avisa en consola: salto sin declarar o cero de escala mal cargado.

## 4. Compatibilidad

Sin `[orientacion]` ni `[saltos]` la salida es idéntica a v6.5
(`survey_profiles_all.csv` byte a byte en el caso sintético de referencia).
