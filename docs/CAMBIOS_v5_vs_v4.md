# Procesador batimétrico ADCP — Cambios de la v5 respecto de la v4

`process_adcp_bathimetric_v5.py`. Resumen de lo nuevo frente a la v4. **El
núcleo numérico validado de la v4 queda intacto** (lectura del `.mat`,
`build_beam_cloud`, `build_profile` con prioridad de VB y márgenes desde el
último VB, gráficos y exports por transecta, QA). Los cambios son en las capas de
ruteo, pelo de agua, agregación del relevamiento y config.

---

## 1. Centerline como RED multi-río (`RiverNetwork`, reemplaza `RiverRoute`)

El shapefile del eje se lee ahora como una **red de canales**: una entidad por
tramo continuo de `río` + `rol`. Atributos:

| Campo | Tipo | Obligatorio | Valores |
|---|---|---|---|
| `river` | texto | sí (para modo multi-río) | `Neuquen`, `Limay`, `Negro` (ASCII, sin tilde) |
| `role` | texto | sí | `main` / `anabranch` |
| `seg_id` | entero | recomendado | id estable de la entidad |
| `label` | texto | opcional | fuerza el rótulo del brazo (ej. `MI`) |
| `flow_dir` | texto | opcional | `US2DS` / `DS2US` (fuerza el sentido) |

- Un relevamiento puede **abarcar más de un río** (Neuquén inferior → confluencia
  → Río Negro). Cada transecta queda etiquetada con su `river` y se reporta en el
  **km oficial de ese río**.
- **Compatibilidad hacia atrás:** si el shapefile **no** tiene columnas
  `river`/`role`, se comporta como la v4 (un solo río `río`, la línea más larga es
  el tronco y las partes-isla que pegan al tronco en ambos extremos son brazos).
  El `r_nqn_inf.shp` de la v4 sigue funcionando sin tocar.

### Preparación en QGIS (para el nuevo esquema)
- Cortar el tronco en la **confluencia** (ahí termina `Neuquen` y arranca `Negro`).
- **Snapear exactamente** todos los extremos que se tocan (confluencia, aperturas
  y cierres de brazos). La red los reconoce como el mismo nodo con tolerancia de
  `NODE_SNAP_TOL = 2 m`; el trazado tiene que venir snapeado, esos 2 m solo
  absorben redondeo sub-métrico.
- Incluir el **Limay** (Arroyito → confluencia) como `river=Limay`, `role=main`:
  queda dibujado pero **inerte** (no lleva transectas en esta campaña).

## 2. Progresiva POR RÍO + progresiva continua del relevamiento

- **Interna**: corre sobre el tronco de cada río desde su extremo de aguas arriba
  (= km 0). km0: Neuquén en El Chañar, Limay en Arroyito, Negro en la confluencia.
- **Oficial** (lo que se reporta): `km_oficial = offset_río + interna`, con los
  offsets por río en el INI, sección `[progresivas]`. Si cada línea se dibuja
  desde su mojón km0, los offsets van en 0.
- **Progresiva en brazos**: proyección al **eje principal** (reporta el km del
  cauce principal a esa altura); estar en un brazo solo cambia el rótulo, no el km.
- **Eje continuo del relevamiento**: para el perfil longitudinal se arma una
  progresiva continua **por el camino relevado** a través de la confluencia. Como
  solo el Neuquén y el Negro llevan transectas ahí, la continuidad Neuquén→Negro
  se resuelve **sola** (el Limay se ignora). Si algún día se relevaran a la vez
  dos ríos que entran a la confluencia, el caso es ambiguo: el perfil queda
  **por río** con aviso, hasta declarar un `mainstem` (ver punto 6).

## 3. Brazos para N canales (MI / M2… / MD, reemplaza Norte/Sur)

- Las islas se detectan por topología (≥2 canales que comparten apertura y
  cierre) y se rotulan **de izquierda a derecha mirando aguas abajo**: `MI`, `M2`,
  …, `MD`. Con 2 canales queda `MI` (brazo) y `MD` (el tronco en ese tramo), o al
  revés según la geometría. `label` explícito en la entidad tiene prioridad.
- **Cambio respecto de v4**: se reemplazan los rótulos `Brazo Norte` / `Brazo Sur`
  por `MI`/`MD`, que son determinísticos y generalizan a cualquier N. Las viejas
  confirmaciones de la v4 (Brazo Norte = `id=2`, origen de progresivas por
  `--chainage-offset`) quedan **sin efecto** bajo la v5.

## 4. Pelo de agua MULTI-ESTACIÓN (GNSS primario, estaciones secundarias)

- **GNSS (PA) sigue siendo la fuente primaria**: el CSV de pelo de agua se lleva a
  la progresiva continua y se interpola igual que en la v4.
- **Registro de estaciones** (shapefile de puntos): `station_id`, `name`, `river`,
  `gauge_zero` (cero de escala, cota SRVN16) y `chainage` opcional. Cada estación
  se ubica en la red.
- **Lecturas por campaña** (CSV): `spot` (`station_id,nivel`) o serie temporal
  (`station_id,datetime,nivel`), autodetectado. `ws_elev = gauge_zero + nivel`.
- **Rol de las estaciones**: relleno cuando el GNSS no cubre (en vez del clamp de
  la v4), única fuente si no hay GNSS, y **QC** (si hay ambas, se compara
  `|WS_gnss − WS_estación|` contra `--ws-qc-tol`, default 0,10 m, y se marca).
- Cada transecta guarda de dónde salió el pelo de agua: `ws_source` =
  `survey` / `station` / `clamp` / `constant`.
- Interpolación entre estaciones: **lineal en progresiva** entre las dos que
  encierran el punto **del mismo río**; fuera del rango, clamp con aviso. Se
  chequea (y avisa) que el nivel decrezca aguas abajo.

## 5. Nueva salida del relevamiento: nube de puntos unificada

`_resumen/survey_raw_beam_points.shp` (o `.gpkg` si supera los límites del
shapefile): **todos los retornos de fondo de todas las transectas** en un solo
`PointZ`, EPSG:5344, con **Z = cota de lecho SRVN16**. Atributos (el km se excluye
a propósito): `pt_id`, `transect`, `ens`, `beam_id`, `beam_type` (`slant`/`vert`),
`bed_elev`, `depth`, `river`, `ws_elev`, `flag` (fuente del pelo de agua; `REL` si
la cota es relativa por falta de PA).

## 6. `mainstem` (tronco canónico) — NO cableado

El gancho para declarar el tronco canónico (Limay→Negro, cuando se releve el
Limay) queda **deliberadamente sin conectar** en la v5. Hoy no hace falta: la
cadena Neuquén→Negro se resuelve por el camino relevado. Se cablea cuando toque.

## 7. INI y línea de comandos

- Nueva sección de INI **`[progresivas]`**: `Neuquen = 0`, `Limay = 0`,
  `Negro = 0` (offsets oficiales por río).
- Nuevas keys/opciones: `stations` / `--stations`, `readings` / `--readings`,
  `ws_qc_tol` / `--ws-qc-tol`.
- `--chainage-offset` pasa a ser el **offset global de reserva** para ríos que no
  figuren en `[progresivas]` (los de `[progresivas]` mandan). `--centerline`
  ahora es la red multi-río (misma opción).
- `survey_index.csv`, `survey_profiles_all.csv` y los shapefiles consolidados
  suman columnas `river`, `km_oficial` y `ws_source`; el índice además trae
  `survey_chainage_m` (progresiva continua). El perfil longitudinal usa esa
  progresiva continua como eje X.

---

## Estado / validación

- **Núcleo numérico**: sin cambios respecto de la v4 validada.
- **Capas nuevas** (red multi-río, progresiva por río + continua, brazos MI/MD,
  pelo de agua multi-estación, nube unificada): **compilan y pasan un test de
  geometría sintético** (`test_v5_synthetic.py`), pero **faltan validar contra el
  shapefile multi-río real y una campaña real** antes de producción.
- **Serie temporal de estaciones**: la conversión de `System.Time` del `.mat` a
  fecha/hora (`representative_time`) es **heurística y no probada** contra tiempo
  real; solo importa para el matcheo por hora en modo serie (esta campaña no tiene
  lecturas). Confirmar codificación y huso del instrumento antes de usarla.
