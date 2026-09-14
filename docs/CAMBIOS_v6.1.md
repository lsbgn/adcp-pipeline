# process_adcp_bathimetric.py — cambios en v6.1

Versión de corrección. El núcleo numérico (armado de la nube de haces, acimut por caudal, filtros de profundidad y `build_profile`) es idéntico a v6.0: 32 de 33 funciones coinciden a nivel AST, y la restante (`group_repeat_stats`) sólo suprime un `RuntimeWarning` inofensivo.

## Por qué

En la campaña 2026-08-25 (Neuquén inferior + confluencia con el Negro) v6.0 terminó sin errores, pero con **las 56 transectas en cotas relativas**. Una cadena de fallas lo explica:

1. Había datos en el Neuquén y en el Limay, o sea dos cabeceras. `build_survey_chain` exigía una sola, así que la cadena continua quedaba indefinida.
2. Sin cadena, los 40 puntos GNSS de pelo de agua se descartaban.
3. `estaciones.csv` se leía con `gpd.read_file`, que para un CSV devuelve un `DataFrame` sin `.crs`. Las escalas tampoco entraban.
4. Sin pelo de agua, cada transecta caía a cota relativa y el aviso no quedaba en ningún archivo.

Además, la transecta de la confluencia (20260825145258), cuya sección cruza el eje del Negro, se asignaba al Limay: su centroide quedaba más cerca del final del eje del Limay.

## Cambios

### 1. Localización de transectas por cruce de sección

El río y la progresiva de cada transecta salen del punto donde la **línea de sección cruza un eje**, como en HEC-RAS.

La línea de sección es la **sección real: el eje ⟂ al flujo sobre el que se arma el perfil**, con la extensión de la traza proyectada sobre él, prolongado un 15 % (mínimo 15 m) a cada lado. No es la traza del bote: en una confluencia el bote puede subir por un tributario y bajar por el otro mientras la sección que representa cruza sólo el cauce de aguas abajo.

Orden de resolución:

1. La sección ⟂ al flujo. Si cruza varios ejes, gana el cruce más cercano al centroide de la traza.
2. Si no cruza ninguno, el eje principal de la traza (`loc_method = crossing-track`, con nota).
3. Si tampoco, el eje más cercano al centroide (`nearest`, con `[warn]`).

El pre-escaneo por GPS no conoce todavía la dirección del flujo, así que usa el eje de la traza y su conteo por río es **provisional**: cada transecta se reubica sobre su sección real al procesarse.

- `transect-locate = centroid` reproduce la regla de v6.0.
- La sección `[rios]` del INI fuerza el río de una transecta puntual.
- El log de cada transecta dice cómo se ubicó. Ejemplo: `nearest axis to the track centroid is Limay (64 m) — the crossing wins`.

### 2. Red de ríos y recorridos, sin tronco

La red se lee como tal. El extremo aguas abajo de cada río se enlaza con el río en el que desemboca, en su km 0 o a mitad de tramo. Cada **recorrido cabecera → desembocadura** (`Neuquen>Negro`, `Limay>Negro`) tiene una progresiva continua.

El tramo bajo una confluencia pertenece a todos los recorridos que pasan por él, así que no hay que elegir un tronco. El orden de los recorridos, por cantidad de transectas, se usa sólo para los gráficos y el índice.

### 3. Pelo de agua por río, unido en la confluencia (`WaterSurfaceModel`)

- **Anclas.** Son los puntos GNSS más las escalas que no quedan entre puntos GNSS. Una escala entre puntos GNSS es sólo de control: se compara contra `ws-qc-tol`.
- **Cota de confluencia.** Se estima desde cada tributario con datos. Para cada uno se interpola, a lo largo de su recorrido, entre su último valor aguas arriba y el primero aguas abajo del nodo; es decir, la pendiente de cada río continúa en el río común. Los estimados se promedian con peso 1/brecha y se obtiene una sola cota. La dispersión entre ramas se controla con `ws-qc-tol`.
- **Fuera del rango de datos.** Se **continúa con la pendiente local**: mínimos cuadrados sobre los valores más cercanos, abarcando al menos 1 km. v5 y v6.0 extendían un valor plano, que a unos 0.9 m/km erra 0.5 m cada 500 m.
  - Un tributario sin datos propios toma la cota de la confluencia y la pendiente del tramo común.
  - Una extrapolación más larga que `ws-extrap-tol` (200 m) es WARN.
- **`ws_source`.** Toma uno de estos valores:

  | Valor | Qué significa |
  |---|---|
  | `survey` | Interpolación entre puntos GNSS |
  | `bridge` | Interpolación que mezcla GNSS, escala o confluencia |
  | `station` | Interpolación entre escalas |
  | `extrap` | Continuación con pendiente fuera del rango de datos |
  | `none` | Sin pelo de agua (cota relativa) |

- **`_resumen/water_surface_qc.csv`.** Lista cada punto GNSS, escala y confluencia con río, km, distancia al eje, río alternativo, rol (anchor / qc / dropped / excluded / junction) y motivo.
- **Puntos ambiguos.** Un punto GNSS con dos ríos candidatos a menos de 25 m de diferencia se marca `AMBIGUOUS`. Una columna opcional `rio` en el CSV de pelo de agua fuerza el río de ese punto.
- **Puntos lejanos.** Un punto GNSS a más de 250 m de su eje, o fuera del tramo digitalizado, se descarta con aviso.

### 4. Estaciones

- El registro se lee desde CSV (columnas X/Y en `ws-crs`) o desde un vectorial.
- Cada escala se proyecta sobre su **río declarado**, no sobre el eje más cercano de cualquier río.
- Se excluyen, con el motivo anotado, las escalas:
  - sin `gauge_zero`,
  - a más de 250 m del eje,
  - fuera del tramo digitalizado.
- Las lecturas en serie temporal se interpolan a la hora de cada transecta. `representative_time` ahora conoce la época SonTek (2000-01-01).

### 5. Nunca más cotas relativas en silencio

- Una transecta sin pelo de agua deja un `[warn]` y queda con `ws_source = none`.
- Si hay fuentes de pelo de agua configuradas y ninguna transecta puede usarlas, la corrida se detiene antes de procesar. `allow-relative = true` permite seguir igual.

### 6. Configuración

- El esquema de claves del INI sale del parser de la línea de comandos. En v6.0 las claves nuevas (`depth-ref`, `composite`, `bt-avg`, `group-tol`, …) **se ignoraban sin aviso**.
- Se validan las opciones y se avisan las claves y secciones desconocidas.
- `[progresivas]` compara los nombres sin tildes ni mayúsculas (`Neuquén` = `Neuquen`) y avisa los que no coinciden con ningún río.
- `;` es comentario en todas las secciones, también en `[campanha]`.

### 7. Registro de la corrida

- `procesamiento.log` guarda la consola completa, incluidos warnings de Python y traceback si la corrida se cae.
- `procesamiento.txt` registra todos los parámetros efectivos, la red, los recorridos, el pelo de agua, cada `[warn]` y las advertencias QA de cada transecta.

### 8. `--section-orientation centerline` deja de ser silencioso

La opción nunca estuvo implementada: `flow_direction_for` no lee `args.section_orientation`, y el acimut de la sección siempre sale ⟂ al flujo medio. En v6.0 elegirla no hacía nada y no avisaba. Ahora emite un `[warn]` y lo dice la ayuda. El comportamiento numérico no cambió.

### 9. Correcciones menores

- `format_progresiva` imprimía `1+1000.00` para 1999.999 m.
- `survey_index.csv` se escribe con el módulo `csv`, porque las notas pueden contener comas.
- El QA de orientación usa la tangente del eje del río de la transecta, no la del eje más cercano.
- La ayuda de `--blend-beams` lo declara obsoleto.
- El docstring de uso nombraba `process_adcp_bathimetric_v4.py`.

## Claves nuevas

| Clave (INI / CLI) | Default | Qué hace |
|---|---|---|
| `transect-locate` | `crossing` | `crossing` o `centroid` (regla v6.0) |
| `ws-extrap-tol` | `200` | [m] extrapolación con pendiente más larga → WARN |
| `allow-relative` | `false` | Seguir aunque ninguna transecta tenga pelo de agua |
| `[rios]` | — | `<id de transecta> = <río>` |
| columna `rio` en el CSV de PA | — | Río forzado para ese punto GNSS |

`--mainstem` y `--junction-tol`, del borrador intermedio, no existen en esta versión.

## Salidas que cambian

- **Nuevo:** `_resumen/water_surface_qc.csv` y `_resumen/procesamiento.log`.
- **`survey_index.csv`:**
  - Se agregan al final las columnas `recorrido` y `loc_method`; las posiciones previas no cambian.
  - `survey_chainage_m` es la progresiva continua sobre el recorrido en el que se grafica la transecta.
- **`survey_long_profile.png`:**
  - Un panel por recorrido con transectas propias.
  - Muestra el pelo de agua del modelo, los puntos GNSS, las escalas (ancla o control) y la cota estimada de cada confluencia.
- **`ws_source`:** desaparece `clamp` y aparece `extrap`.

## Compatibilidad

- `RiverNetwork.build_survey_chain` / `survey_chainage`, `WaterSurfaceProfile` y `StationWS.elev_at` se conservan tal cual (API legacy, para `test_synthetic.py`). El pipeline ya no los usa.
- La progresiva por cruce coincide con la de v6.0 en transectas simétricas. Puede diferir unos metros en trazas asimétricas u oblicuas; `transect-locate = centroid` reproduce v6.0.
- En el INI de 2026-08-25, las claves v6 ignoradas por v6.0 tenían su valor por defecto. Ese bug no alteró los perfiles de esa campaña.

## Validación sobre red sintética

Red sintética armada alrededor de los 40 puntos reales de `PA_20260825.csv`, con las escalas reales y 21 transectas `.mat` sintéticas: 12 en el Neuquén (con un aforo ×4) y 9 en el Negro (con un aforo ×4).

| Chequeo | Resultado |
|---|---|
| Río y km de las 15 salidas (13 simples + 2 aforos) | Río correcto en todas; \|Δkm\| ≤ 0.38 m |
| Sección sobre el Negro 5 m bajo la confluencia, traza cargada a una margen | Centroide: Limay km 65+590.8 (el error real). Cruce: **Negro km 0+005** |
| Traza que recorre ambos tributarios con sección ⟂ al flujo sobre el otro cauce | La sección manda sobre el eje de la traza; el respaldo se rotula `crossing-track` |
| Regresión de la sección ⟂ al flujo sobre A y B | Profundidades, ríos, progresivas y pelo de agua idénticos a la corrida validada |
| Cota de confluencia | Neuquén 254.604, Limay 254.526, adoptada 254.565 (Δ 0.078 m, PASS) |
| Pelo de agua en Negro km 5 / 500 / 1200 / 2000 / 3800 | Igual al cálculo manual, al mm |
| NEU-04 (control entre puntos GNSS) | −0.048 m, PASS |
| Transecta 500 m aguas arriba de Pt1 | 265.795 con pendiente −1.04 m/km (antes: 265.273 plano) |
| Limay sin datos propios (Pt150 forzado al Neuquén) | Toma la cota de la confluencia y continúa con la pendiente del Negro |
| Lecturas en serie temporal | Pelo de agua interpolado a la hora de cada transecta; igual al cálculo manual |
| API legacy vs v6.0 | Mismos resultados |

## Al reprocesar 2026-08-25

1. Correr `test_synthetic.py`.
2. Reprocesar con tu `campanha.ini`. Confirmar en consola:
   - la topología (`Neuquen -> Negro`, `Limay -> Negro`),
   - la línea `confluence into 'Negro'`,
   - `pelo de agua por perfil` sin `none`.
3. Revisar en la consola y en `process_log.txt` cómo se ubicó 20260825145258: debe decir Negro por `crossing`.
4. En `water_surface_qc.csv`:
   - Ver qué puntos GNSS quedaron en el Limay y si alguno está marcado `AMBIGUOUS`. Si un punto de la desembocadura del Neuquén cayó en el Limay, completar la columna `rio`.
   - Revisar la dispersión entre ramas de la cota de confluencia.
5. Pt150 queda fuera del CSV: se tomó sobre un brazo que comunica el Neuquén con el Limay en la confluencia, así que no representa a ningún río. Sin él, la cota de confluencia sale sólo de la rama Neuquén (Pt149 → Isla Jordan). Dejar constancia en `observaciones` del INI.
6. Isla Jordan sigue a verificar. Desde Pt149, la caída hasta la escala da unos 0.9–1.0 m/km, contra 0.55 m/km entre Isla Jordan y María Elvira. Puede ser un quiebre real de pendiente o un cero de escala bajo. Un punto GNSS de pelo de agua junto a la escala en la próxima campaña lo resuelve.
7. El CV de los aforos queda fuera de alcance: los caudales se calculan con QRev.
