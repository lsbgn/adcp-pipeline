"""End-to-end synthetic check of v6.2: builds a two-river network that STRADDLES
the POSGAR faja 2/3 boundary (lon -67.5), fake .mat transects with GNSS
dropouts, a longitudinal track, gauges and a WS CSV, then runs the pipeline."""
import os, shutil, subprocess, sys, math, datetime
import numpy as np
import scipy.io as sio
from pyproj import CRS, Transformer
import geopandas as gpd
from shapely.geometry import LineString

from pathlib import Path
HERE = Path(__file__).resolve().parent

ROOT = str(HERE / "prueba")
shutil.rmtree(ROOT, ignore_errors=True)
for d in ("adcp", "tracks", "gnss", "ejes", "niveles"):
    os.makedirs(f"{ROOT}/{d}", exist_ok=True)

OUT = CRS.from_epsg(5344)                     # deliberate: data crosses into f3
to_out = Transformer.from_crs(4326, OUT, always_xy=True)
to_ll = Transformer.from_crs(OUT, 4326, always_xy=True)
SONTEK_EPOCH = datetime.datetime(2000, 1, 1)

# --- axis: Neuquen from lon -67.8 to -67.2, i.e. across the faja boundary ----
lonA = np.linspace(-67.80, -67.20, 200)
latA = -39.00 + 0.004 * np.sin(np.linspace(0, 2.4, 200))
xA, yA = to_out.transform(lonA, latA)
gpd.GeoDataFrame({"river": ["Neuquen"], "role": ["main"]},
                 geometry=[LineString(np.column_stack([xA, yA]))],
                 crs=OUT).to_file(f"{ROOT}/ejes/eje.shp")
axis = LineString(np.column_stack([xA, yA]))
L = axis.length
print(f"[test] eje: {L / 1000:.1f} km, lon {lonA[0]} a {lonA[-1]} "
      f"(cruza el límite faja 2/3 en -67,5)")

# --- water surface: 1.0 m/km, plus gauges -----------------------------------
WS0, SLOPE = 300.0, -1.0e-3


def ws_at(s):
    return WS0 + SLOPE * s


with open(f"{ROOT}/gnss/PA.csv", "w") as f:
    f.write("Pto,East,North,H_correg\n")
    for k, s in enumerate(np.linspace(1500, L - 1500, 12)):
        p = axis.interpolate(s)
        f.write(f"Pt{k+1},{p.x:.3f},{p.y:.3f},{ws_at(s):.3f}\n")

# two gauges 1.7 m apart, downstream of every GNSS point: the v6.1 slope trap
with open(f"{ROOT}/ejes/estaciones.csv", "w") as f:
    f.write("station_id,name,river,gauge_zero,X,Y\n")
    for sid, s, zero in (("NEG-04", L - 300.0, 229.15), ("NEG-05", L - 298.3, 229.56)):
        p = axis.interpolate(s)
        f.write(f"{sid},Escala {sid},Neuquen,{zero:.3f},{p.x:.3f},{p.y:.3f}\n")
with open(f"{ROOT}/niveles/niveles.csv", "w") as f:
    f.write("station_id,nivel\n")
    for sid, s, zero in (("NEG-04", L - 300.0, 229.15), ("NEG-05", L - 298.3, 229.56)):
        f.write(f"{sid},{ws_at(s) - zero:.3f}\n")


def bed_depth(s, off):
    """Analytic channel: 4.5 m thalweg, parabolic banks, half-width 70 m."""
    return max(0.0, 4.5 * (1.0 - (off / 70.0) ** 2) - 0.25 * math.sin(s / 800.0))


def write_mat(path, s_center, n=90, half=70.0, longitudinal=False,
              dropouts=(), t0=None, span=1200.0):
    """One fake RiverSurveyor .mat. A transect crosses the axis; a longitudinal
    run travels ALONG it."""
    if longitudinal:
        ss = np.linspace(s_center - span / 2, s_center + span / 2, n)
        offs = 12.0 * np.sin(np.linspace(0, 3, n))
    else:
        ss = np.full(n, s_center)
        offs = np.linspace(-half, half, n)
    E, N, dep = [], [], []
    for s, o in zip(ss, offs):
        p = axis.interpolate(float(np.clip(s, 0, L)))
        q = axis.interpolate(float(np.clip(s + 5, 0, L)))
        tx, ty = q.x - p.x, q.y - p.y
        m = math.hypot(tx, ty) or 1.0
        nx, ny = -ty / m, tx / m
        E.append(p.x + o * nx); N.append(p.y + o * ny)
        dep.append(bed_depth(s, o))
    E, N, dep = np.array(E), np.array(N), np.array(dep)
    lon, lat = to_ll.transform(E, N)
    utmc = auto_utm(lat, lon)
    ex, ny_ = Transformer.from_crs(4326, utmc, always_xy=True).transform(lon, lat)
    utm = np.column_stack([ex, ny_])
    for i in dropouts:                                   # RiverSurveyor 0/0
        utm[i] = 0.0
        lat[i] = 0.0
        lon[i] = 0.0
    hdg = np.degrees(np.arctan2(np.gradient(ex), np.gradient(ny_))) % 360.0
    t0 = t0 or datetime.datetime(2026, 8, 25, 14, 52, 58)
    t = (t0 - SONTEK_EPOCH).total_seconds() + np.arange(n) * 2.0
    mdic = dict(
        BottomTrack=dict(VB_Depth=dep,
                         BT_Beam_Depth=np.column_stack([dep + 0.05 * np.sin(np.arange(n) + k)
                                                        for k in range(4)]),
                         BT_Frequency=np.full(n, 3000.0), BT_Depth=dep),
        GPS=dict(UTM=utm, Latitude=lat, Longitude=lon),
        System=dict(True_North_ADP_Heading=hdg, Time=t,
                    Step=np.full(n, 3.0)),
        Compass=dict(Pitch=np.zeros(n), Roll=np.zeros(n)),
        Setup=dict(Edges_0__DistanceToBank=3.0, Edges_1__DistanceToBank=3.0,
                   startEdge=1, sensorDepth=0.1, depthReference=0.0),
        Summary=dict(Depth=dep, Mean_Vel=np.column_stack([np.full(n, 0.8),
                                                          np.full(n, 0.1)]),
                     Total_Q=np.full(n, 420.0)),
    )
    sio.savemat(path, mdic, oned_as="column")


def auto_utm(lat, lon):
    z = int(math.floor((float(np.nanmean(lon)) + 180) / 6) + 1)
    return CRS.from_epsg(32700 + z)


# 9 transects spread along the axis (so some land in faja 3), one pair 6 m apart
# to exercise the aforo grouping; two files get GNSS dropouts.
ts = list(np.linspace(1800, L - 1800, 8)) + [float(np.linspace(1800, L - 1800, 8)[4]) + 6.0]
for k, s in enumerate(ts):
    write_mat(f"{ROOT}/adcp/2026082514{k:02d}00.mat", s,
              dropouts=(10, 11, 12) if k in (2, 5) else (),
              t0=datetime.datetime(2026, 8, 25, 14, 0, 0) + datetime.timedelta(minutes=6 * k))
# two longitudinal runs between sections
for k, s in enumerate(np.linspace(4000, L - 4000, 2)):
    write_mat(f"{ROOT}/tracks/2026082516{k:02d}00.mat", s, n=300,
              longitudinal=True, span=2500.0,
              t0=datetime.datetime(2026, 8, 25, 16, 0, 0) + datetime.timedelta(minutes=30 * k))

with open(f"{ROOT}/campanha.ini", "w", encoding="utf-8") as f:
    f.write("""[campanha]
campanha = 2026-08-25 prueba sintetica v6.2
[procesamiento]
inputs            = ./adcp
tracks            = ./tracks
outdir            = ./salida
survey-name       = prueba62
centerline        = ./ejes/eje.shp
water-surface-csv = ./gnss/PA.csv
ws-crs            = EPSG:5344
stations          = ./ejes/estaciones.csv
readings          = ./niveles/niveles.csv
dx                = 1.0
group-tol         = 15.0
[progresivas]
Neuquen = 0
""")
print("[test] entradas listas")
r = subprocess.run([sys.executable, str(HERE / "process_adcp_bathimetric.py"),
                    "--config", f"{ROOT}/campanha.ini"],
                   cwd=str(HERE), capture_output=True, text=True)
print(r.stdout[-9000:])
if r.returncode != 0:
    print("STDERR:", r.stderr[-4000:])
sys.exit(r.returncode)
