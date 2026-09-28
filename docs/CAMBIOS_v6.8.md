# Cambios en v6.8 respecto de v6.7

**Tema:** descartes manuales de profundidad por sección, con la sección nueva `[descartes]` del INI.

## Por qué

20260825153036 es la sección del río Negro en la progresiva 2+018. En esa sección el haz vertical ve algo angosto 2.5–3.7 m por encima del fondo, cerca de la MD. RSL también lo muestra.

| s [m] | Ensambles de RSL | VB | Haces inclinados cerca |
|---|---|---|---|
| 135.6–134.6 | 25–26 | sin fondo (0) | 4.3–5.0 |
| 132.9–131.7 | 27–30 | 1.50 · 2.48 · 2.55 · 2.73 | B3 4.3–4.9; B2 3.9–4.1; B4 3.3; B1 4.0 / 2.84 / 5.24 |
| a ambos lados | | 5.2–5.3 | 5.0–5.3 |

B1 pasa por el mismo s dos veces con un segundo de diferencia y lee 2.84 m la primera y 5.24 m la segunda: lo que ve es un objeto (un raigón o un tronco), no un escalón del lecho.

El filtro de picos de v6 no lo saca, por dos razones:

- Trabaja sobre la serie de cada haz, y aquí son 4 ensambles seguidos.
- Están en plena rampa de la barranca, donde el rango intercuartil móvil es grande.

La única salida en v6.7 era `depth-ref = bt` para toda la transecta. Eso cambia 252 de los 292 nodos para arreglar 7.

## Qué cambia

### 1. Sección `[descartes]` del INI

```ini
[descartes]
; <id> = <regla>, <regla>, ...
; <regla> = <fuente>[+<fuente>] [<condición> <valor>] ...
20260825153036 = vb s 131.5-134.0
```

- **Fuentes:** `vb`, `b1` a `b4`, `bt` (los cuatro inclinados) y `todos`. Se combinan con `+`.
- **Condiciones.** Se cumplen todas a la vez. Sin condiciones, se descarta la fuente en toda la sección.
  - `s`: distancia transversal [m], leída en la sección sin descartes.
  - `ens`: número de ensamble de RSL, `System.Sample`, desde 1.
  - `prof`: profundidad del haz [m].
  - `cota`: cota del punto [m]. Necesita pelo de agua.
  - `dif`: para el VB, la diferencia con la mediana de los inclinados del mismo ensamble; para un inclinado, la diferencia con el VB.
- **Valores:** `A-B`, `<X` o `>X`, con los bordes incluidos. Para `ens` vale además un número suelto.
- **`<id>`:**
  - El nombre del `.mat` o sus 14 dígitos.
  - El id de un grupo o un nombre de `[grupos]`. En ese caso la regla se aplica a la nube fusionada, sobre el eje del grupo, y no admite `ens`.
  - Una regla con el id de una transecta se aplica a esa repetición y el grupo hereda el descarte. Nunca se vuelve a aplicar sobre el eje del grupo, aunque `20260825153036` coincida con `20260825153036r_x4`.

Los puntos descartados salen del ajuste del lecho. Donde falta la fuente primaria, el composite completa con la secundaria, igual que QRev.

**Marco de s.** Las condiciones `s` se evalúan en la sección sin descartes, que es el `cross_section.png` que se miró para escribir la regla. Después, la extensión medida se recalcula con los puntos que quedan. Si eso corre el origen de s, el log lo informa.

### 2. Salidas

| Salida | Cambio |
|---|---|
| `process_log.txt` | Una línea por regla: puntos descartados por haz, rango de s y de ensambles. `[warn]` si la regla no coincide con ningún punto |
| `cross_section.png` | Los puntos descartados como cruces rojas, con su entrada en la leyenda. El subtítulo pasa a punto decimal |
| `raw_bed_points.csv` | Dos columnas al final: `ens` (ensamble de RSL) y `descarte` (la regla, vacía si el punto se usó). Los descartados siguen en el archivo |
| `survey_index.csv` | Columna `descartes` al final (en un grupo, suma la de los miembros y la del grupo) |
| Nube de la campaña (`survey_raw_beam_points.shp`) | Sin los puntos descartados |
| `procesamiento.txt` | Lista de reglas leídas |

Una regla mal escrita corta la corrida al arrancar, con el renglón y el motivo. Por ejemplo, con coma decimal (`vb s 131,5-134`) la corrida corta y avisa que el separador es el punto.

### 3. `tools/comparar_corridas.py`

- `nube_igual` compara solo las 12 columnas de v6.7 de `raw_bed_points.csv`. Así, una corrida v6.7 y una v6.8 de los mismos datos dan `si`.
- Columnas nuevas `descartes_a` y `descartes_b`, con la cantidad de puntos descartados en cada corrida.

## Validación

**Sin reglas, v6.8 reproduce v6.7.**

- Los 5 `.mat` de la campaña 2026-04-15 dan `bathymetric_profile.csv` idénticos byte a byte, y `nube_igual = si`.
- 20260825153036 también sale idéntico.

**Con `20260825153036 = vb s 131.5-134.0`:**

- Se descartan los 4 puntos del VB de los ensambles 27–30 y nada más.
- Cambian 7 nodos, entre s = 131.0 y 134.0 m. El perfil coincide con el parche manual hecho antes de la implementación (diferencia máxima 0.000 m).

| s [m] | 131.5 | 132.0 | 132.5 | 133.0 | 133.5 | 134.0 |
|---|---|---|---|---|---|---|
| v6.7 | 3.64 | 2.61 | 2.29 | 1.29 | 1.79 | 3.22 |
| v6.8 | 5.31 | 5.30 | 5.30 | 4.22 | 4.07 | 5.17 |

La loma de ~1 m que queda en s = 133 es la de los haces inclinados (B2, B3, B4 leen 3.3–4.5 m ahí).

**Las otras formas de escribir la regla dan el mismo resultado.** Estas tres descartan exactamente esos 4 puntos:

- `vb ens 27-30`
- `vb dif >1.5`
- `vb cota >249.4 s 125-136`

`vb cota >249.4 s 125-136` se probó con pelo de agua constante y `bed-fit = median`.

**`dif` necesita acotarse.** Umbrales de `dif` en esa sección:

| Umbral | Ensambles que saca |
|---|---|
| 0.5 m | 3, 4, 21, 22, 27–30 |
| 1.0 m | 21, 22, 27–30 |
| 1.5 m y 2.0 m | 27–30 |

Los ensambles 21–22 están al pie de la barranca: ahí la diferencia es real. Por eso conviene acotar `dif` con `s` o `ens`.

**Grupos.** Se armó un grupo de prueba con dos copias de la transecta. La regla de un miembro quedó solo en ese miembro. La regla del grupo, `b1 s 133-134 prof <3.5`, sacó un punto por repetición. `survey_index.csv` suma los dos: 6 puntos.

**Reglas inválidas.** Cortan con un mensaje claro. Se probaron:

- coma decimal;
- fuente `b5`;
- condición desconocida;
- `s 5` sin rango.

**Reglas que no aplican.** No cortan la corrida, avisan con `[warn]`:

- una regla que no coincide con ningún punto;
- `cota` sin pelo de agua.
