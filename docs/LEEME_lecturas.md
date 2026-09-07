# Lecturas de nivel de estaciones — formato de entrada (v5)

> Los valores de `X`, `Y` y `gauge_zero` en `estaciones_ejemplo.csv` son
> PLACEHOLDERS: reemplazalos por los relevados. Los `station_id` sí deben
> coincidir entre el registro de estaciones y las lecturas.

En un relevamiento real usás UN archivo de lecturas por campaña. El modo
(spot / serie) se detecta por la presencia de la columna `datetime` y por la
cantidad de registros por estación. Acá van los dos formatos por separado.

## Registro de estaciones (capa de puntos, EPSG:5344)

Se arma una vez y se reutiliza entre campañas. `estaciones_ejemplo.csv` es la
tabla de atributos; podés generar la capa en QGIS con "Agregar capa de texto
delimitado" a partir de X,Y, o usar tu propio `.shp`.

| campo      | descripción                                          |
|------------|------------------------------------------------------|
| station_id | id corto y estable (lo referencian las lecturas)     |
| name       | nombre para el informe                               |
| river      | Neuquen \| Limay \| Negro (ASCII, sin tildes)        |
| X, Y       | POSGAR07 faja 2 (EPSG:5344)                          |
| gauge_zero | cero de escala, cota SRVN16 [m]                      |
| chainage   | opcional; si va vacío se proyecta X,Y sobre el eje   |

## Lecturas por campaña

Cota de pelo de agua en la estación:

    ws_elev = gauge_zero + nivel        # nivel = altura de escala [m]

### Modo spot — una lectura por estación

Altura al momento del relevamiento.

    station_id, nivel

### Modo serie — varias lecturas por estación (recomendado para ADCP)

    station_id, datetime, nivel

- Más de un registro por estación ⇒ se trata como serie.
- Cada transecta toma el nivel por timestamp (interpolado entre registros).
- `datetime` en horario local (Argentina, UTC-3, sin DST).

> IMPORTANTE: el `datetime` debe estar en la MISMA convención horaria que los
> ensembles del `.mat`. Según la config del equipo, SonTek puede loguear en
> UTC o en hora local — confirmalo por equipo antes de matchear series.
