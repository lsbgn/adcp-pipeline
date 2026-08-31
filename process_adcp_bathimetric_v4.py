#!/usr/bin/env python3
"""
process_adcp_bathimetric.py
================================================================================
ADCP CROSS-SECTION EXTRACTION — Bathymetry-only pipeline (v4.0)

Builds a hydraulically consistent 1D bathymetric cross-section perpendicular
to the mean flow, from a SonTek M9/S5 .mat file (RiverSurveyor Live export).

Velocities are NOT exported. Only the bed profile and supporting QA artifacts.

Changes vs v3:
    * River chainage (progresiva) now understands a BRANCHED centerline: a river
      that splits into two brazos around an island keeps a single, continuous
      set of progresivas (measured along the whole river), and each transect is
      additionally tagged "Brazo Norte" / "Brazo Sur". Branch topology (split /
      rejoin nodes, N-S identity) is detected automatically and reported.
    * Absolute bed elevation in a real datum (default IGN SRVN16). Pass a CSV of
      water-surface points (X, Y, elevation); the script converts them to the
      river's longitudinal water-surface profile (elevation vs progresiva) and
      interpolates the water-surface elevation at each transect. bed_elev_m then
      becomes an ABSOLUTE elevation (m a.s.l.), not a value relative to the
      water surface. The route is auto-oriented so progresiva increases
      downstream, using the water-surface slope itself.
    * Multi-file SURVEY mode: point the script at several .mat files (or a folder)
      and it processes every transect, then writes a consolidated survey index,
      a merged point cloud, a whole-survey plan view, a longitudinal thalweg
      profile, and combined shapefiles — all keyed by progresiva + brazo.
    * build_profile() sampling half-window defaults to 2*dx (bin_half_width=None),
      and the default grid spacing is dx = 0.50 m.

Kept from v3: perfil identifier on every output; unified slant-beam weights;
perpendicular-offset penalty; frequency-coloured plan view; Spanish plot text;
startEdge=0 bug fix; HEC-RAS CSV removed; per-beam weights kept internal.

Per-transect outputs (in ./out_<basename>/ or <survey>/<perfil>/):
    bathymetric_profile.csv   perfil, progresiva, brazo, s, depth, ws_elev, bed_elev, ...
    raw_bed_points.csv        every projected beam footprint (no weights)
    plan_view_map.png         track + axis + beam cloud (by frequency)
    cross_section.png         raw points + smoothed bed (absolute datum if WS given)
    axis.shp / profile_points.shp / raw_beam_points.shp

Extra survey-level outputs (multi-file mode), in <survey>/_resumen/:
    survey_index.csv          one row per transect (progresiva, brazo, cotas, ancho…)
    survey_profiles_all.csv   every gridded profile point of every transect
    survey_plan_view.png      centerline + all section axes, coloured by brazo
    survey_long_profile.png   water surface + thalweg vs progresiva
    survey_axes.shp / survey_profile_points.shp

Usage:
    # single transect
    python process_adcp_bathimetric_v4.py <input.mat> [--outdir <dir>]
        [--centerline <river_axis.shp>] [--water-surface-csv <ws.csv>]
        [--chainage-offset <m>] [--offset-scale <m>] [--no-offset-weighting]

    # whole survey (folder of .mat, or several files)
    python process_adcp_bathimetric_v4.py <survey_dir>/  --centerline <axis.shp>
        --water-surface-csv <ws.csv> [--survey-name <name>] [--outdir <dir>]

Tested with: SonTek RiverSurveyor M9, .mat v5 export
================================================================================
"""

from __future__ import annotations

import argparse
import configparser
import csv
import datetime
import glob
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyproj import CRS, Transformer

# Optional shapefile output + chainage (geopandas / shapely)
try:
    import geopandas as gpd
    from shapely.geometry import LineString, Point
    from shapely.ops import linemerge
    HAS_GPD = True
except Exception as e:
    HAS_GPD = False
    print(f"[warn] geopandas not available — shapefiles/chainage will be skipped ({e})")

# Default vertical datum label for absolute bed elevations.
DATUM_NAME_DEFAULT = "IGN SRVN16"

# Script version (reported in logs and the traceability record).
__version__ = "4.0"

# Tolerance [m] for deciding that a line's endpoint lies ON the main path,
# i.e. that the line is a side branch (island split) rather than a disjoint reach.
BRANCH_SNAP_TOL = 15.0


# ============================================================================ #
#  STEP 0 — file inspection helpers
# ============================================================================ #

def load_mat(path: str | Path) -> dict:
    """Load a SonTek .mat file with struct unwrapping."""
    return sio.loadmat(path, struct_as_record=False, squeeze_me=True)


def get_field(struct, name, default=None):
    """Safely retrieve a field from a scipy mat_struct."""
    if struct is None:
        return default
    if hasattr(struct, "_fieldnames") and name in struct._fieldnames:
        return getattr(struct, name)
    return default


# ============================================================================ #
#  STEP 1 — data extraction
# ============================================================================ #

def extract_data(mat: dict) -> dict:
    """Pull canonical fields out of the .mat structure."""
    bt   = mat.get("BottomTrack")
    gps  = mat.get("GPS")
    sysd = mat.get("System")
    comp = mat.get("Compass")
    setp = mat.get("Setup")
    summ = mat.get("Summary")

    vb_depth   = np.asarray(get_field(bt, "VB_Depth"), dtype=float)
    bt_beams   = np.asarray(get_field(bt, "BT_Beam_Depth"), dtype=float)  # (n, 4)
    bt_freq    = np.asarray(get_field(bt, "BT_Frequency"), dtype=float)   # (n,)

    utm        = get_field(gps, "UTM")
    lat        = np.asarray(get_field(gps, "Latitude"), dtype=float)
    lon        = np.asarray(get_field(gps, "Longitude"), dtype=float)

    heading    = np.asarray(
        get_field(sysd, "True_North_ADP_Heading",
                  default=get_field(sysd, "Heading")),
        dtype=float,
    )
    pitch      = np.asarray(get_field(comp, "Pitch", default=np.zeros_like(heading)), dtype=float)
    roll       = np.asarray(get_field(comp, "Roll",  default=np.zeros_like(heading)), dtype=float)

    time       = np.asarray(get_field(sysd, "Time", default=np.arange(len(vb_depth))), dtype=float)

    mean_vel   = get_field(summ, "Mean_Vel")     # (n, 2) east, north
    if mean_vel is not None:
        mean_vel = np.asarray(mean_vel, dtype=float)

    # Edge distances from Setup (used for bank extrapolation).
    # NOTE: `or 0.0` here is safe for the edge distances because a 0.0 m edge
    # is meaningless anyway; but for startEdge a legitimate 0 (=Left) must NOT
    # be turned into 1, so it is handled explicitly below (was a bug in v5).
    edge_left  = float(get_field(setp, "Edges_0__DistanceToBank", default=0.0) or 0.0)
    edge_right = float(get_field(setp, "Edges_1__DistanceToBank", default=0.0) or 0.0)
    _se        = get_field(setp, "startEdge", default=1)
    start_edge = int(_se) if _se is not None else 1   # 0 = Left, 1 = Right

    return dict(
        vb_depth=vb_depth,
        bt_beams=bt_beams,
        bt_freq=bt_freq,
        utm=np.asarray(utm, dtype=float) if utm is not None else None,
        lat=lat, lon=lon,
        heading=heading, pitch=pitch, roll=roll,
        time=time,
        mean_vel=mean_vel,
        edge_left=edge_left, edge_right=edge_right, start_edge=start_edge,
    )


def auto_utm_crs(lat: np.ndarray, lon: np.ndarray) -> CRS:
    """Auto-detect UTM zone from mean lat/lon."""
    lat_m = float(np.nanmean(lat))
    lon_m = float(np.nanmean(lon))
    zone = int(math.floor((lon_m + 180.0) / 6.0) + 1)
    south = lat_m < 0
    epsg = (32700 if south else 32600) + zone
    return CRS.from_epsg(epsg)


def auto_posgar07_crs(lon: np.ndarray) -> CRS:
    """
    Pick the appropriate POSGAR 2007 / Argentina zone (EPSG 5343..5349)
    from longitude. Central meridians:
        Zona 1: -72°  -> EPSG 5343
        Zona 2: -69°  -> EPSG 5344
        Zona 3: -66°  -> EPSG 5345
        Zona 4: -63°  -> EPSG 5346
        Zona 5: -60°  -> EPSG 5347
        Zona 6: -57°  -> EPSG 5348
        Zona 7: -54°  -> EPSG 5349
    Each zone is 3° wide (fajas Gauss-Krüger).
    """
    lon_m = float(np.nanmean(lon))
    cms   = [-72, -69, -66, -63, -60, -57, -54]
    epsgs = [5343, 5344, 5345, 5346, 5347, 5348, 5349]
    idx = int(np.argmin([abs(lon_m - cm) for cm in cms]))
    return CRS.from_epsg(epsgs[idx])


# ============================================================================ #
#  STEP 2 — flow direction
# ============================================================================ #

def compute_flow_direction(mean_vel: np.ndarray) -> tuple[float, float]:
    """
    Return (theta_flow, theta_section) in radians, math convention
    (CCW from +East, +X axis).

    theta_section = theta_flow + 90 deg  (perpendicular axis).
    """
    if mean_vel is None or mean_vel.size == 0:
        raise ValueError("Cannot compute flow direction: no Mean_Vel available.")
    ve = mean_vel[:, 0]
    vn = mean_vel[:, 1]
    mask = np.isfinite(ve) & np.isfinite(vn) & ((ve != 0) | (vn != 0))
    if not np.any(mask):
        raise ValueError("Cannot compute flow direction: all Mean_Vel are zero/NaN.")
    theta_flow    = math.atan2(np.nanmean(vn[mask]), np.nanmean(ve[mask]))
    theta_section = theta_flow + math.pi / 2.0
    return theta_flow, theta_section


# ============================================================================ #
#  STEP 3 — beam footprints
# ============================================================================ #

# SonTek M9 azimuth sets (instrument frame, deg from FORWARD, CW)
AZ_3000 = np.array([0.0, 90.0, 180.0, 270.0])     # F, S, A, P
AZ_1000 = np.array([45.0, 135.0, 225.0, 315.0])   # FS, AS, AP, FP
TILT_DEG = 25.0

# Nominal frequencies [kHz] — used for the plan-view frequency colouring
FREQ_VB   = 500.0     # vertical beam (nadir)
FREQ_HIGH = 3000.0    # high-frequency slant set (0/90/180/270)
FREQ_LOW  = 1000.0    # low-frequency  slant set (45/135/225/315)
FREQ_SPLIT = 2500.0   # threshold to classify a sample's active slant set

# --- Bed-detection weights ------------------------------------------------- #
# For BATHYMETRY the vertical beam (VB, nadir) is the reference measurement.
# All four slant beams share the same 25° tilt, so their nadir-referenced
# reliability is the same regardless of instrument azimuth or frequency. The
# previous azimuth/frequency-specific weights (0.8 / 0.4 / 0.6) were not
# physically justified for bed detection: instrument azimuth does not map
# cleanly onto the section plane once boat heading / ferry angle are taken into
# account. That geometry is now handled separately by the perpendicular-offset
# penalty in main() (OFFSET_SCALE_DEFAULT).
#
# >>> Tune these if you want a different VB-vs-slant balance. <<<
W_VB    = 1.0     # vertical beam (500 kHz, nadir)
W_SLANT = 0.5     # every slant beam (25° tilt), any frequency / azimuth

# Default Gaussian scale [m] for the perpendicular-offset penalty (see main()).
OFFSET_SCALE_DEFAULT = 3.0


def beam_weight(az_deg: float, freq_khz: float) -> float:
    """
    Base bed-detection weight for a slant beam.

    All four slant beams share the same 25° tilt, so they receive the same base
    weight for bathymetry (see the note next to W_SLANT). The arguments are kept
    so per-beam / per-frequency tuning stays a one-line change if ever needed.
    """
    return W_SLANT


def build_beam_cloud(data: dict, utm: np.ndarray) -> dict:
    """
    For every sample (ensemble), generate the VB footprint and the 4 active
    slant-beam footprints in (E, N, depth) using vessel heading.

    Note on terminology:
        * SAMPLE (ensemble) = one row in the .mat arrays. It is the average
          of several individual acoustic pings. This is what we iterate over
          here — it is the unit of bottom-detection output.
        * PING = one individual acoustic shot. System.Pings stores how many
          pings were averaged into each sample (typically 7–35 in this data).

    Returns columns of equal length, one entry per (sample, beam) pair.
    """
    n = utm.shape[0]                           # number of SAMPLES
    heading  = data["heading"]                 # deg, true north, +CW
    bt_beams = data["bt_beams"]                # (n_samples, 4)
    bt_freq  = data["bt_freq"]                 # (n_samples,)
    vb_depth = data["vb_depth"]

    pts_e, pts_n, pts_z, pts_w, pts_b, pts_f = [], [], [], [], [], []

    tilt = math.radians(TILT_DEG)

    for i in range(n):                         # loop over SAMPLES
        E0, N0 = utm[i]
        if not (np.isfinite(E0) and np.isfinite(N0)):
            continue
        hdg = heading[i] if np.isfinite(heading[i]) else 0.0
        hdg_rad = math.radians(hdg)            # ADP +CW from true north

        # --- vertical beam ---
        if np.isfinite(vb_depth[i]) and vb_depth[i] > 0:
            pts_e.append(E0); pts_n.append(N0); pts_z.append(float(vb_depth[i]))
            pts_w.append(W_VB); pts_b.append(0); pts_f.append(FREQ_VB)

        # --- 4 slant beams (active set chosen by per-sample frequency) ---
        f = float(bt_freq[i]) if np.isfinite(bt_freq[i]) else FREQ_HIGH
        az_set = AZ_3000 if f >= FREQ_SPLIT else AZ_1000

        for k in range(4):
            d = bt_beams[i, k]
            if not np.isfinite(d) or d <= 0:
                continue
            az_inst = math.radians(az_set[k])         # instrument frame, +CW from forward
            # World azimuth (true north, +CW) = heading + instrument azimuth
            az_world = hdg_rad + az_inst
            r = d * math.tan(tilt)                    # horizontal radius of footprint
            # Compass azimuth (CW from N): dE = r*sin(az), dN = r*cos(az)
            dE = r * math.sin(az_world)
            dN = r * math.cos(az_world)
            pts_e.append(E0 + dE)
            pts_n.append(N0 + dN)
            pts_z.append(float(d))
            pts_w.append(beam_weight(az_set[k], f))
            pts_b.append(k + 1)
            pts_f.append(f)

    return dict(
        E=np.asarray(pts_e),
        N=np.asarray(pts_n),
        depth=np.asarray(pts_z),
        weight=np.asarray(pts_w),
        beam_index=np.asarray(pts_b, dtype=int),
        beam_freq_khz=np.asarray(pts_f),
    )


# ============================================================================ #
#  STEP 4 — project to cross-section axis
# ============================================================================ #

def project_to_axis(
    E: np.ndarray, N: np.ndarray,
    centroid: tuple[float, float],
    theta_section: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Project (E,N) onto the cross-section axis.
    The axis passes through the centroid and has direction (cos θ, sin θ)
    in (E,N) using math convention (θ from +East, CCW).

    Returns:
        s        — along-section distance (signed)
        d_offset — perpendicular offset (signed). Its MAGNITUDE is the distance
                   of the point from the section plane and is invariant to any
                   later axis flip; main() uses |d_offset| for the offset penalty.
    """
    Ec, Nc = centroid
    ux, uy = math.cos(theta_section), math.sin(theta_section)   # axis unit vector
    nx, ny = -uy, ux                                            # left-perpendicular
    dE = E - Ec
    dN = N - Nc
    s        = dE * ux + dN * uy
    d_offset = dE * nx + dN * ny
    return s, d_offset


# ============================================================================ #
#  STEP 4b — river chainage (progresiva)
# ============================================================================ #

def format_progresiva(p_m: float) -> str:
    """Argentine chainage notation: metres -> 'k+mmm.mm' (e.g. 12340.5 -> '12+340.50')."""
    sign = "-" if p_m < 0 else ""
    p = abs(float(p_m))
    km = int(p // 1000)
    m  = p - km * 1000.0
    return f"{sign}{km}+{m:06.2f}"


def _pick_col(cols, preferred, aliases):
    """Case-insensitive column resolver for the water-surface CSV."""
    low = {c.lower(): c for c in cols}
    if preferred and preferred.lower() in low:
        return low[preferred.lower()]
    for a in aliases:
        if a in low:
            return low[a]
    return None


class RiverRoute:
    """
    A river centerline that may split into branches (brazos) around islands.

    The route keeps a single, continuous chainage (progresiva) measured along a
    'main path' from its upstream end. Where the river splits, ONE channel is the
    main path and the OTHER is a branch; both share the progresiva at the split
    node, and every transect is additionally tagged with a brazo label
    ('Brazo Norte' / 'Brazo Sur' / '').  This realises the requirement of keeping
    the progresivas of the whole river ('del conjunto') while telling the two
    brazos apart.

    Chainage is computed in `crs` (the centerline's own projected CRS unless
    overridden).  Callers pass points already reprojected to `crs`.
    """

    def __init__(self, centerline_path, crs_override=None,
                 reverse: bool = False, offset: float = 0.0,
                 ws_xy_elev=None, log=None):
        if not HAS_GPD:
            raise RuntimeError("RiverRoute needs geopandas/shapely")
        self.offset = float(offset)
        self._reverse = bool(reverse)
        self._log = log if log is not None else []

        gdf = gpd.read_file(centerline_path)
        if gdf.crs is None:
            if crs_override is not None:
                gdf = gdf.set_crs(crs_override)
                self._log.append("[warn] centerline has no .prj; assuming "
                                 f"{CRS.from_user_input(crs_override).to_string()}")
            else:
                raise ValueError("centerline has no CRS (.prj) and no override given")
        elif crs_override is not None:
            gdf = gdf.to_crs(crs_override)
        self.crs = gdf.crs

        geoms = [g for g in gdf.geometry if g is not None and not g.is_empty]
        if not geoms:
            raise ValueError("centerline file has no usable geometry")

        # linemerge on the RAW lines: it stitches lines that share endpoints but
        # (unlike unary_union) does NOT node them where a branch merely touches
        # another line's interior — so an island brazo stays a separate part
        # instead of chopping the main channel into pieces.
        merged = linemerge(geoms)
        parts = [merged] if merged.geom_type == "LineString" else list(merged.geoms)
        parts.sort(key=lambda g: g.length, reverse=True)
        main = self._orient(parts[0], ws_xy_elev)
        self.main = main
        self.main_len = main.length

        # Classify the remaining parts. A part whose BOTH endpoints snap onto the
        # main path (within BRANCH_SNAP_TOL) is a branch (island split); anything
        # else is a disjoint reach we cannot chain, and we warn.
        self.branches = []
        for g in parts[1:]:
            d0 = main.distance(Point(g.coords[0]))
            d1 = main.distance(Point(g.coords[-1]))
            if d0 <= BRANCH_SNAP_TOL and d1 <= BRANCH_SNAP_TOL:
                self._add_branch(g)
            else:
                self._log.append(
                    f"[warn] centerline part of length {g.length:.0f} m is not "
                    f"connected to the main path (endpoint gaps {d0:.0f}/{d1:.0f} m) "
                    "— ignored for chainage")

        self._log.append(
            f"[info] centerline: main path {self.main_len:.0f} m, "
            f"{len(self.branches)} branch(es)")
        for b in self.branches:
            self._log.append(
                f"[info]   {b['label']}: {b['geom'].length:.0f} m, splits at "
                f"{format_progresiva(self.offset + b['split_s'])}, rejoins at "
                f"{format_progresiva(self.offset + b['rejoin_s'])} "
                f"(main path there = {b['main_label']})")

    # ------------------------------------------------------------------ #
    def _orient(self, line, ws_xy_elev):
        """Return `line` oriented so chainage increases DOWNSTREAM."""
        if ws_xy_elev:
            xy = np.asarray([(p[0], p[1]) for p in ws_xy_elev], float)
            zs = np.asarray([p[2] for p in ws_xy_elev], float)
            s = np.asarray([line.project(Point(x, y)) for x, y in xy])
            ok = np.isfinite(s) & np.isfinite(zs)
            if ok.sum() >= 3 and float(np.ptp(s[ok])) > 1.0:
                slope = float(np.polyfit(s[ok], zs[ok], 1)[0])
                if slope > 0:      # elevation rises with raw chainage -> raw dir is upstream
                    self._log.append(
                        "[info] route oriented from water surface: raw direction runs "
                        f"UPSTREAM, reversed so progresiva grows downstream "
                        f"(WS slope was +{slope*1000:.2f} m/km)")
                    return LineString(list(line.coords)[::-1])
                self._log.append(
                    "[info] route orientation confirmed from water surface "
                    f"(WS slope {slope*1000:.2f} m/km, decreasing downstream)")
                return line
        if self._reverse:
            self._log.append("[info] route reversed by --chainage-reverse")
            return LineString(list(line.coords)[::-1])
        return line

    def _add_branch(self, g):
        main = self.main
        s0 = main.project(Point(g.coords[0]))
        s1 = main.project(Point(g.coords[-1]))
        if s1 < s0:                                   # branch digitised upstream
            g = LineString(list(g.coords)[::-1])
            s0, s1 = s1, s0
        bmid = g.interpolate(0.5 * g.length)
        mmid = main.interpolate(0.5 * (s0 + s1))
        if bmid.y >= mmid.y:
            blabel, mlabel = "Brazo Norte", "Brazo Sur"
        else:
            blabel, mlabel = "Brazo Sur", "Brazo Norte"
        self.branches.append(dict(geom=g, split_s=s0, rejoin_s=s1,
                                  label=blabel, main_label=mlabel))

    # ------------------------------------------------------------------ #
    def chainage(self, x: float, y: float):
        """(progresiva_m, brazo_label, dist_to_axis_m) for a point in `crs`."""
        p = Point(x, y)
        best_seg, best_d, best_i = "main", self.main.distance(p), None
        for i, b in enumerate(self.branches):
            d = b["geom"].distance(p)
            if d < best_d:
                best_seg, best_d, best_i = f"branch{i}", d, i
        if best_seg == "main":
            s = self.main.project(p)
            brazo = ""
            for b in self.branches:
                if b["split_s"] <= s <= b["rejoin_s"]:
                    brazo = b["main_label"]
                    break
            return self.offset + s, brazo, best_d
        b = self.branches[best_i]
        s = b["split_s"] + b["geom"].project(p)
        return self.offset + s, b["label"], best_d

    def snapped_xy(self, x: float, y: float):
        """Point on the route nearest to (x, y), in `crs` — for QA plots."""
        p = Point(x, y)
        line = self.main
        best_d = self.main.distance(p)
        for b in self.branches:
            d = b["geom"].distance(p)
            if d < best_d:
                best_d, line = d, b["geom"]
        q = line.interpolate(line.project(p))
        return q.x, q.y


def read_ws_csv(csv_path, east_col="East", north_col="North", elev_col="H_correg"):
    """Read a water-surface CSV -> list of (x, y, elevation) in the file's CRS."""
    rows = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rd = csv.DictReader(f)
        cols = rd.fieldnames or []
        ec = _pick_col(cols, east_col, ["east", "x", "este", "x_utm", "x_posgar07"])
        nc = _pick_col(cols, north_col, ["north", "y", "norte", "y_utm", "y_posgar07"])
        hc = _pick_col(cols, elev_col, ["h_correg", "z", "elev", "cota", "h", "altura"])
        if not (ec and nc and hc):
            raise ValueError(f"WS CSV: could not resolve E/N/elev columns among {cols}")
        for r in rd:
            try:
                rows.append((float(r[ec]), float(r[nc]), float(r[hc])))
            except (TypeError, ValueError):
                continue
    if len(rows) < 2:
        raise ValueError("WS CSV: need at least 2 valid points")
    return rows


class WaterSurfaceProfile:
    """
    Longitudinal water-surface profile: elevation as a function of progresiva.

    Built from a CSV of surveyed water-surface points (easting, northing,
    elevation).  Each point is projected onto the river route to get its
    progresiva; the profile is the piecewise-linear elevation(progresiva).  This
    is the physically correct way to carry a water surface along a river: the
    surface follows the longitudinal slope, so interpolating in progresiva (not
    in 2-D space) honours meanders and avoids mixing far-apart banks.
    """

    def __init__(self, route: "RiverRoute", points, log=None):
        """`points` = list of (x, y, elevation) ALREADY in route.crs."""
        self._log = log if log is not None else []
        self.route = route
        rows = list(points)
        if len(rows) < 2:
            raise ValueError("WS: need at least 2 valid points")

        prog, elev, dist = [], [], []
        for x, y, h in rows:
            s, _b, d = route.chainage(x, y)
            prog.append(s); elev.append(h); dist.append(d)
        order = np.argsort(prog)
        prog = np.asarray(prog)[order]
        elev = np.asarray(elev)[order]
        dist = np.asarray(dist)[order]

        # de-duplicate near-identical progresivas (average their elevation)
        up, ue = [], []
        i = 0
        while i < len(prog):
            j = i
            while j + 1 < len(prog) and (prog[j + 1] - prog[i]) < 0.01:
                j += 1
            up.append(float(prog[i:j + 1].mean()))
            ue.append(float(elev[i:j + 1].mean()))
            i = j + 1
        self.prog = np.asarray(up)
        self.elev = np.asarray(ue)
        self.raw_prog = prog
        self.raw_elev = elev
        self.max_dist = float(np.nanmax(dist)) if len(dist) else float("nan")
        self.slope_m_per_km = float(np.polyfit(self.prog, self.elev, 1)[0] * 1000.0)

        self._log.append(
            f"[info] water surface: {len(self.prog)} pts, progresiva "
            f"{self.prog.min():.0f}–{self.prog.max():.0f} m, elev "
            f"{self.elev.min():.3f}–{self.elev.max():.3f} m, slope "
            f"{self.slope_m_per_km:.2f} m/km, max offset from axis "
            f"{self.max_dist:.0f} m")

    def elev_at(self, progresiva):
        """(ws_elev, note). `note` flags extrapolation beyond the surveyed range."""
        lo, hi = float(self.prog[0]), float(self.prog[-1])
        if progresiva < lo:
            return float(self.elev[0]), f"EXTRAP {progresiva - lo:.0f} m (below WS range)"
        if progresiva > hi:
            return float(self.elev[-1]), f"EXTRAP +{progresiva - hi:.0f} m (above WS range)"
        return float(np.interp(progresiva, self.prog, self.elev)), ""


# ============================================================================ #
#  STEP 5 — build final bed profile
# ============================================================================ #

def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Weighted median of `values` with positive `weights`."""
    if values.size == 0:
        return np.nan
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    if cw[-1] <= 0:
        return float(np.nanmedian(v))
    cutoff = 0.5 * cw[-1]
    return float(v[np.searchsorted(cw, cutoff)])


def build_profile(
    s_pts: np.ndarray, depth_pts: np.ndarray, w_pts: np.ndarray,
    s_data_min: float, s_data_max: float, dx: float = 0.5,
    bin_half_width: float | None = None,
    edge_left_dist: float = 0.0,
    edge_right_dist: float = 0.0,
    beam_index: np.ndarray | None = None,
    vb_priority: bool = True,
    vb_min: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build the single-valued bed profile depth = f(s) on a regular grid.

    `bin_half_width` is the half-width [m] of the s-window used to sample the
    beam cloud at each grid node. When None it defaults to 2*dx, i.e. every node
    is informed by points within ±2 grid steps — wide enough to always capture
    neighbours at dx = 0.5 m, narrow enough not to over-smooth real bedforms.

    Beam selection at each node (when `beam_index` is supplied and `vb_priority`):
        * The vertical beam (VB, nadir, beam_index == 0) is the reliable bed
          measurement. Where a node's window holds at least `vb_min` VB points,
          the bed is the weighted median of the VB points ALONE.
        * Only where the window has too few VB points (the near-bank / edge zones
          the nadir beam never reached) do the slant beams take over, extending
          coverage toward the banks.
    This stops the more numerous slant beams (4 per sample vs 1 VB) from out-
    voting the nadir beam and flattening the bank slope. Set vb_priority=False to
    restore the older behaviour (weighted median of all beams together).
    """
    if bin_half_width is None:
        bin_half_width = 2 * dx

    # Reference the bank ramps to the VERTICAL-BEAM (nadir) extent, not the full
    # cloud. The slant beams reach further toward the banks but are unreliable on
    # the bank slope, so they must NOT define the bank: the bed is measured only
    # across the VB span, and from the outermost VB on each side a straight ramp
    # descends to depth 0 at the bank (edge_left/right metres beyond the last VB).
    ref_min, ref_max = s_data_min, s_data_max
    if vb_priority and beam_index is not None:
        _bi = np.asarray(beam_index)
        _ok = (np.isfinite(s_pts) & np.isfinite(depth_pts)
               & (depth_pts > 0) & (_bi == 0))
        if np.any(_ok):
            ref_min = float(np.min(s_pts[_ok]))
            ref_max = float(np.max(s_pts[_ok]))

    # Shift so the outermost LEFT reference point sits at s = edge_left_dist
    # (i.e., left bank at s = 0).
    shift = edge_left_dist - ref_min
    sp_shift = s_pts + shift
    s_data_min_new = edge_left_dist                       # left edge of measured bed
    s_data_max_new = ref_max + shift                      # right edge of measured bed

    # Total cross-section length: left bank (0) → right bank (data_max_new + edge_right)
    grid_min = 0.0
    grid_max = s_data_max_new + edge_right_dist
    s_grid = np.arange(grid_min, grid_max + dx / 2.0, dx)
    depth_grid = np.full_like(s_grid, np.nan)

    finite = np.isfinite(sp_shift) & np.isfinite(depth_pts) & (depth_pts > 0)
    sp = sp_shift[finite]; dp = depth_pts[finite]; wp = w_pts[finite]
    if beam_index is not None:
        bi = np.asarray(beam_index)[finite]
        is_vb = bi == 0
    else:
        is_vb = None

    # Bin-and-weighted-median in the data zone
    for i, sc in enumerate(s_grid):
        if sc < s_data_min_new - dx or sc > s_data_max_new + dx:
            continue
        mask = (sp >= sc - bin_half_width) & (sp <= sc + bin_half_width)
        if not np.any(mask):
            continue
        if vb_priority and is_vb is not None:
            vb_mask = mask & is_vb
            if int(np.count_nonzero(vb_mask)) >= vb_min:
                depth_grid[i] = weighted_median(dp[vb_mask], wp[vb_mask])
                continue
        depth_grid[i] = weighted_median(dp[mask], wp[mask])

    # --- force depth = 0 at the banks ---------------------------------------
    depth_grid[0] = 0.0       # left bank (s = 0)
    depth_grid[-1] = 0.0      # right bank (s = s_grid[-1])

    # --- linear gap-fill across the whole grid ------------------------------
    valid = np.isfinite(depth_grid)
    if valid.sum() >= 2:
        depth_grid = np.interp(s_grid, s_grid[valid], depth_grid[valid])

    # --- light smoothing (Savitzky-Golay) -----------------------------------
    win = max(5, int(round(2.0 / dx)))
    if win % 2 == 0:
        win += 1
    if len(depth_grid) >= win:
        try:
            depth_grid = savgol_filter(depth_grid, window_length=win, polyorder=2)
            depth_grid = np.clip(depth_grid, 0.0, None)
        except Exception:
            pass

    # Re-pin banks to exactly 0 after smoothing (smoothing can drift them slightly)
    depth_grid[0]  = 0.0
    depth_grid[-1] = 0.0

    return s_grid, depth_grid, shift


# ============================================================================ #
#  STEP 6 — exports
# ============================================================================ #

def export_csvs(outdir: Path, perfil_id: str, progresiva_m,
                s_grid: np.ndarray, depth_grid: np.ndarray,
                axis_origin_utm: tuple[float, float],
                theta_section: float,
                utm_crs: CRS, posgar_crs: CRS,
                water_surface_elev: float = 0.0,
                brazo: str = "", ws_elev=None):
    """
    Write bathymetric_profile.csv (HEC-RAS CSV intentionally not produced).

    If `ws_elev` is given (interpolated from the water-surface profile) the bed
    elevation is ABSOLUTE (m in the vertical datum): bed_elev = ws_elev - depth.
    Otherwise it falls back to the constant `water_surface_elev` (default 0, i.e.
    relative to the water surface) and the ws_elev_m column is left blank.
    """
    Ec, Nc = axis_origin_utm
    ux, uy = math.cos(theta_section), math.sin(theta_section)

    # Reconstruct grid locations in UTM. Caller provides axis_origin_utm so that
    # s = 0 corresponds to the left bank on the axis.
    E_grid = Ec + s_grid * ux
    N_grid = Nc + s_grid * uy

    # Reproject UTM → lat/lon (WGS84) and UTM → POSGAR07
    tr_to_ll      = Transformer.from_crs(utm_crs, CRS.from_epsg(4326),  always_xy=True)
    tr_to_posgar  = Transformer.from_crs(utm_crs, posgar_crs,           always_xy=True)
    lon_g, lat_g  = tr_to_ll.transform(E_grid, N_grid)
    xp,    yp     = tr_to_posgar.transform(E_grid, N_grid)

    wse = float(ws_elev) if ws_elev is not None else float(water_surface_elev)
    bed_elev = wse - depth_grid
    prog_str = "" if progresiva_m is None else f"{progresiva_m:.3f}"
    ws_str   = "" if ws_elev is None else f"{wse:.3f}"

    path1 = outdir / "bathymetric_profile.csv"
    with open(path1, "w", encoding="utf-8") as f:
        f.write("perfil,progresiva_m,brazo,s_m,depth_m,ws_elev_m,bed_elev_m,"
                "latitude_deg,longitude_deg,x_utm,y_utm,x_posgar07,y_posgar07\n")
        for i in range(len(s_grid)):
            f.write(f"{perfil_id},{prog_str},{brazo},{s_grid[i]:.4f},{depth_grid[i]:.4f},"
                    f"{ws_str},{bed_elev[i]:.4f},"
                    f"{lat_g[i]:.8f},{lon_g[i]:.8f},"
                    f"{E_grid[i]:.4f},{N_grid[i]:.4f},"
                    f"{xp[i]:.4f},{yp[i]:.4f}\n")

    return path1, (E_grid, N_grid, lon_g, lat_g, xp, yp, bed_elev)


def export_raw_points(outdir: Path, perfil_id: str, cloud: dict, s_pts: np.ndarray,
                       utm_crs: CRS, posgar_crs: CRS):
    """Write raw_bed_points.csv (per-beam weights intentionally omitted)."""
    tr_to_ll     = Transformer.from_crs(utm_crs, CRS.from_epsg(4326), always_xy=True)
    tr_to_posgar = Transformer.from_crs(utm_crs, posgar_crs,          always_xy=True)
    lon, lat = tr_to_ll.transform(cloud["E"], cloud["N"])
    xp, yp   = tr_to_posgar.transform(cloud["E"], cloud["N"])

    path = outdir / "raw_bed_points.csv"
    with open(path, "w", encoding="utf-8") as f:
        f.write("perfil,s_m,depth_m,beam_index,beam_freq_khz,"
                "x_utm,y_utm,latitude_deg,longitude_deg,x_posgar07,y_posgar07\n")
        for i in range(len(s_pts)):
            f.write(f"{perfil_id},{s_pts[i]:.4f},{cloud['depth'][i]:.4f},"
                    f"{int(cloud['beam_index'][i])},{cloud['beam_freq_khz'][i]:.0f},"
                    f"{cloud['E'][i]:.4f},{cloud['N'][i]:.4f},"
                    f"{lat[i]:.8f},{lon[i]:.8f},"
                    f"{xp[i]:.4f},{yp[i]:.4f}\n")
    return path


# ============================================================================ #
#  STEP 7 — QA plots
# ============================================================================ #

def _subtitle(perfil_id: str, progresiva_m, brazo: str = "") -> str:
    """Second title line: transect identifier, chainage and (if any) brazo."""
    parts = [f"Perfil: {perfil_id}"]
    if progresiva_m is not None:
        parts.append(f"Progresiva {format_progresiva(progresiva_m)}")
    if brazo:
        parts.append(brazo)
    return "   —   ".join(parts)


def plot_plan_view(outdir: Path, perfil_id: str, progresiva_m,
                    utm: np.ndarray, cloud: dict,
                    axis_origin: tuple[float, float], theta_section: float,
                    s_grid_extent: tuple[float, float],
                    theta_flow: float, brazo: str = ""):
    """
    Plan view with:
        - boat track (orange)
        - beam footprints coloured by frequency band (VB 500 / 3000 / 1000 kHz)
          so the 1000-vs-3000 kHz beam-set detection can be checked visually
        - cross-section axis (black) with MI (s=0) and MD (s=s_max) markers
        - flow arrow emanating FROM the transect, pointing downstream
    (All on-plot text in Spanish.)
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # --- beam footprints coloured by frequency band -----------------------
    freq = cloud["beam_freq_khz"]
    m_vb   = np.isclose(freq, FREQ_VB)
    m_high = freq >= FREQ_SPLIT
    m_low  = (~m_vb) & (~m_high) & (freq > 0)
    ax.scatter(cloud["E"][m_low],  cloud["N"][m_low],  s=6, c="tab:blue",
               alpha=0.45, label="Beams 1000 kHz")
    ax.scatter(cloud["E"][m_high], cloud["N"][m_high], s=6, c="tab:red",
               alpha=0.55, label="Beams 3000 kHz")
    ax.scatter(cloud["E"][m_vb],   cloud["N"][m_vb],   s=8, c="black",
               alpha=0.65, label="VB (500 kHz)")

    # --- boat track --------------------------------------------------------
    ax.plot(utm[:, 0], utm[:, 1], "-", color="tab:orange", lw=1.5,
            label="Trayectoria del bote")
    ax.plot(utm[:, 0], utm[:, 1], "o", ms=3, color="tab:orange", alpha=0.7)
    ax.plot(utm[0, 0], utm[0, 1], "^", ms=10, color="tab:orange",
            markeredgecolor="black", label="Inicio de trayectoria")

    # --- cross-section axis: s=0 (MI) → s=s_max (MD) -----------------------
    Ec, Nc = axis_origin
    ux, uy = math.cos(theta_section), math.sin(theta_section)
    s0, s1 = s_grid_extent
    Ea, Na = Ec + s0 * ux, Nc + s0 * uy        # LEFT bank  (s = 0)
    Eb, Nb = Ec + s1 * ux, Nc + s1 * uy        # RIGHT bank (s = s_max)
    ax.plot([Ea, Eb], [Na, Nb], "k-", lw=2, label="Eje de la sección")
    ax.plot(Ea, Na, "ks", ms=10)
    ax.plot(Eb, Nb, "ko", ms=10)
    ax.annotate("MI  (s = 0)", xy=(Ea, Na),
                xytext=(8, 8), textcoords="offset points",
                fontsize=11, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85))
    ax.annotate(f"MD  (s = {s1:.1f} m)", xy=(Eb, Nb),
                xytext=(8, -16), textcoords="offset points",
                fontsize=11, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85))

    # --- flow arrow: tail on the transect (axis midpoint), head downstream --
    mE = 0.5 * (Ea + Eb); mN = 0.5 * (Na + Nb)
    track_dx = float(np.nanmax(utm[:, 0]) - np.nanmin(utm[:, 0]))
    track_dy = float(np.nanmax(utm[:, 1]) - np.nanmin(utm[:, 1]))
    track_extent = max(track_dx, track_dy)
    L = max(2.0, 0.18 * track_extent)
    fx = math.cos(theta_flow); fy = math.sin(theta_flow)   # +flow = downstream
    hE, hN = mE + L * fx, mN + L * fy
    ax.annotate(
        "", xy=(hE, hN), xytext=(mE, mN),
        arrowprops=dict(arrowstyle="-|>", color="navy", lw=3, mutation_scale=25),
    )
    ax.text(mE + 1.15 * L * fx, mN + 1.15 * L * fy, "FLUJO",
            color="navy", fontsize=11, fontweight="bold",
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="navy", alpha=0.9))

    # Lock axis limits to the data extent (with a small pad).
    pad = 0.08 * track_extent + 1.0
    all_E = np.concatenate([utm[:, 0], cloud["E"], [Ea, Eb, hE]])
    all_N = np.concatenate([utm[:, 1], cloud["N"], [Na, Nb, hN]])
    ax.set_xlim(np.nanmin(all_E) - pad, np.nanmax(all_E) + pad)
    ax.set_ylim(np.nanmin(all_N) - pad, np.nanmax(all_N) + pad)

    ax.set_xlabel("Este [m]")
    ax.set_ylabel("Norte [m]")
    ax.set_aspect("equal")
    ax.set_title("Vista en planta — trayectoria, eje de la sección y huellas de beams\n"
                 + _subtitle(perfil_id, progresiva_m, brazo))
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = outdir / "plan_view_map.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_cross_section(outdir: Path, perfil_id: str, progresiva_m,
                        s_pts: np.ndarray, depth_pts: np.ndarray,
                        beam_idx: np.ndarray,
                        s_grid: np.ndarray, depth_grid: np.ndarray,
                        brazo: str = "", ws_elev=None,
                        datum_name: str = DATUM_NAME_DEFAULT):
    """Cross-section: slant beams in one translucent colour; VB distinct.

    If `ws_elev` is given, the vertical axis is ABSOLUTE elevation in the datum
    (bed and beams plotted as ws_elev - depth, water surface drawn at ws_elev).
    Otherwise the axis is relative to the water surface (0 = surface).
    (All on-plot text in Spanish.)"""
    fig, ax = plt.subplots(figsize=(11, 5.5))

    absolute = ws_elev is not None
    e0 = float(ws_elev) if absolute else 0.0        # water-surface level on the plot
    y_bed  = e0 - depth_grid
    y_surf = e0
    y_floor = float(np.nanmin(y_bed)) - 0.5

    # Bed fill + water column
    ax.fill_between(s_grid, y_bed, y_floor,
                    color="#7a5230", alpha=0.5, zorder=1, label="Lecho (relleno)")
    ax.fill_between(s_grid, y_surf, y_bed,
                    color="#cfe6f5", alpha=0.6, zorder=0, label="Columna de agua")

    # Raw beam footprints: slant beams unified (single colour, more transparent),
    # vertical beam kept distinct.
    m_lat = beam_idx > 0
    m_vb  = beam_idx == 0
    if np.any(m_lat):
        ax.scatter(s_pts[m_lat], e0 - depth_pts[m_lat], s=8, c="#5b7fa6",
                   alpha=0.20, label="Beams laterales", zorder=2)
    if np.any(m_vb):
        ax.scatter(s_pts[m_vb], e0 - depth_pts[m_vb], s=16, c="black",
                   alpha=0.75, label="Beam vertical (VB)", zorder=3)

    # Smoothed bed
    ax.plot(s_grid, y_bed, "k-", lw=2.0, label="Lecho suavizado", zorder=4)

    # Bank annotations (MI = margen izquierda, MD = margen derecha)
    s0, s1 = float(s_grid[0]), float(s_grid[-1])
    span = y_surf - y_floor
    ax.axvline(s0, color="dimgray", lw=1, ls=":", alpha=0.7)
    ax.axvline(s1, color="dimgray", lw=1, ls=":", alpha=0.7)
    ax.text(s0, y_floor + 0.98 * span, "MI\n(margen izq.)",
            ha="left", va="top", fontsize=11, fontweight="bold", color="black",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85))
    ax.text(s1, y_floor + 0.98 * span, "MD\n(margen der.)",
            ha="right", va="top", fontsize=11, fontweight="bold", color="black",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85))

    ax.set_xlim(s0 - 0.5, s1 + 0.5)
    ax.set_ylim(y_floor, y_surf + max(0.15, 0.05 * span))
    ax.set_xlabel("Distancia transversal, s [m]")
    if absolute:
        ax.set_ylabel(f"Cota {datum_name} [m]")
        surf_lbl = f"Superficie del agua ({y_surf:.2f} m {datum_name})"
    else:
        ax.set_ylabel("Cota relativa a la superficie del agua [m]")
        surf_lbl = "Superficie del agua"
    ax.set_title("Sección transversal — nube de puntos y lecho suavizado\n"
                 + _subtitle(perfil_id, progresiva_m, brazo))
    ax.axhline(y_surf, color="navy", lw=1, ls="--", alpha=0.7, label=surf_lbl)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8, ncol=3)
    fig.tight_layout()
    path = outdir / "cross_section.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ============================================================================ #
#  STEP 8 — shapefiles
# ============================================================================ #

def export_axis_shp(outdir: Path, perfil_id: str, progresiva_m,
                     axis_origin: tuple[float, float], theta_section: float,
                     s_grid_extent: tuple[float, float],
                     utm_crs: CRS, brazo: str = ""):
    """Write axis.shp (single LineString covering the full s-range)."""
    if not HAS_GPD:
        return None
    Ec, Nc = axis_origin
    ux, uy = math.cos(theta_section), math.sin(theta_section)
    s0, s1 = s_grid_extent
    line = LineString([(Ec + s0 * ux, Nc + s0 * uy),
                       (Ec + s1 * ux, Nc + s1 * uy)])
    gdf = gpd.GeoDataFrame(
        {"perfil": [perfil_id],
         "progr_m": [np.nan if progresiva_m is None else float(progresiva_m)],
         "brazo": [brazo],
         "theta_deg": [math.degrees(theta_section)],
         "length_m": [s1 - s0]},
        geometry=[line], crs=utm_crs,
    )
    path = outdir / "axis.shp"
    gdf.to_file(path)
    return path


def export_profile_points_shp(outdir: Path, perfil_id: str, progresiva_m,
                              geo: dict, utm_crs: CRS, posgar_crs: CRS,
                              brazo: str = "", ws_elev=None):
    if not HAS_GPD:
        return None
    s_m       = geo["s_m"]
    depth_m   = geo["depth_m"]
    bed_elev  = geo["bed_elev"]
    E         = geo["E"]; N = geo["N"]
    lon       = geo["lon"]; lat = geo["lat"]
    xp        = geo["xp"];  yp  = geo["yp"]
    prog_val  = np.nan if progresiva_m is None else float(progresiva_m)
    ws_val    = np.nan if ws_elev is None else float(ws_elev)
    geom = [Point(E[i], N[i]) for i in range(len(s_m))]
    gdf = gpd.GeoDataFrame({
        "perfil": perfil_id,
        "progr_m": prog_val,
        "brazo": brazo,
        "s_m": s_m, "depth_m": depth_m,
        "ws_elev_m": ws_val, "bed_elev_m": bed_elev,
        "lat": lat, "lon": lon,
        "x_utm": E, "y_utm": N,
        "x_posgar07": xp, "y_posgar07": yp,
    }, geometry=geom, crs=utm_crs)
    path = outdir / "profile_points.shp"
    gdf.to_file(path)
    return path


def export_raw_beam_points_shp(outdir: Path, perfil_id: str, cloud: dict, s_pts: np.ndarray,
                                utm_crs: CRS, posgar_crs: CRS, brazo: str = ""):
    """raw_beam_points.shp (per-beam weights intentionally omitted)."""
    if not HAS_GPD:
        return None
    tr_to_ll     = Transformer.from_crs(utm_crs, CRS.from_epsg(4326), always_xy=True)
    tr_to_posgar = Transformer.from_crs(utm_crs, posgar_crs,          always_xy=True)
    lon, lat = tr_to_ll.transform(cloud["E"], cloud["N"])
    xp, yp   = tr_to_posgar.transform(cloud["E"], cloud["N"])

    geom = [Point(cloud["E"][i], cloud["N"][i]) for i in range(len(s_pts))]
    gdf = gpd.GeoDataFrame({
        "perfil": perfil_id,
        "brazo": brazo,
        "s_m": s_pts,
        "depth_m": cloud["depth"],
        "beam_idx": cloud["beam_index"],
        "freq_khz": cloud["beam_freq_khz"],
        "lat": lat, "lon": lon,
        "x_utm": cloud["E"], "y_utm": cloud["N"],
        "x_posgar07": xp, "y_posgar07": yp,
    }, geometry=geom, crs=utm_crs)
    path = outdir / "raw_beam_points.shp"
    gdf.to_file(path)
    return path


# ============================================================================ #
#  STEP 9 — acceptance tests
# ============================================================================ #

def acceptance_tests(theta_flow: float, theta_section: float,
                      s_grid: np.ndarray, depth_grid: np.ndarray,
                      edge_left: float,
                      log: list[str]):
    # Perpendicularity test (works modulo 180° because the section is invariant
    # under axis flip). Reduce the angle to [-90, 90], then the deviation from
    # perpendicular is |90 − |a||.
    raw_angle = math.degrees(theta_section - theta_flow)
    a = (raw_angle + 180) % 360 - 180     # to [-180, 180]
    if a > 90:  a -= 180
    if a < -90: a += 180
    perp_err = abs(90.0 - abs(a))
    ok_perp = perp_err <= 10.0
    log.append(f"[QA] axis perpendicular to flow within ±10°? "
               f"{'PASS' if ok_perp else 'FAIL'} "
               f"(deviation = {perp_err:.1f}°)")

    ok_left = abs(s_grid[0]) < 1e-3
    log.append(f"[QA] left bank at s = 0?  "
               f"{'PASS' if ok_left else 'FAIL'} "
               f"(s_grid[0] = {s_grid[0]:.4f})")

    ok_cont = np.all(np.isfinite(depth_grid))
    log.append(f"[QA] profile continuous (no NaN)? "
               f"{'PASS' if ok_cont else 'FAIL'}")

    return ok_perp and ok_cont and ok_left


# ============================================================================ #
#  PER-TRANSECT DRIVER
# ============================================================================ #

def process_one(matpath: Path, args, route=None, wsp=None,
                datum_name: str = DATUM_NAME_DEFAULT,
                survey_outdir: Path | None = None, quiet: bool = False):
    """
    Full pipeline for ONE .mat transect. Returns a result dict (or None on
    failure). When `survey_outdir` is given, outputs go to <survey_outdir>/<perfil>/.
    """
    perfil_id = matpath.stem
    if survey_outdir is not None:
        outdir = survey_outdir / perfil_id
    elif args.outdir:
        outdir = Path(args.outdir)
    else:
        outdir = Path.cwd() / f"out_{perfil_id}"
    outdir.mkdir(parents=True, exist_ok=True)

    log: list[str] = []
    log.append(f"[info] input  : {matpath}")
    log.append(f"[info] perfil : {perfil_id}")
    log.append(f"[info] outdir : {outdir}")

    # ------------------------------------------------------------ STEP 0/1
    try:
        mat = load_mat(matpath)
        data = extract_data(mat)
    except Exception as e:
        log.append(f"[error] could not read/extract {matpath.name}: {e}")
        if not quiet:
            print("\n".join(log))
        return None

    log.append(f"[info] samples: {len(data['vb_depth'])}  "
               f"(each sample averages several acoustic pings)")
    log.append(f"[info] BT freq mix: "
               f"{dict(zip(*np.unique(data['bt_freq'], return_counts=True)))}")
    log.append(f"[info] depths VB : min={np.nanmin(data['vb_depth']):.2f}m, "
               f"max={np.nanmax(data['vb_depth']):.2f}m")
    log.append(f"[info] edges   : L={data['edge_left']:.2f}m, R={data['edge_right']:.2f}m "
               f"(startEdge={data['start_edge']}; 0=Left, 1=Right)")

    utm_crs    = auto_utm_crs(data["lat"], data["lon"])
    posgar_crs = auto_posgar07_crs(data["lon"])
    log.append(f"[info] UTM CRS    : {utm_crs.name} ({utm_crs.to_epsg()})")
    log.append(f"[info] POSGAR07   : {posgar_crs.name} ({posgar_crs.to_epsg()})")

    if data["utm"] is not None and np.isfinite(data["utm"]).all():
        utm = data["utm"]
        log.append("[info] using GPS.UTM directly")
    else:
        tr = Transformer.from_crs(CRS.from_epsg(4326), utm_crs, always_xy=True)
        Ex, Ny = tr.transform(data["lon"], data["lat"])
        utm = np.column_stack([Ex, Ny])
        log.append("[info] derived UTM from Lat/Lon")

    # ------------------------------------------------------------ STEP 2
    theta_flow, theta_section = compute_flow_direction(data["mean_vel"])
    log.append(f"[info] theta_flow    = {math.degrees(theta_flow):+.1f}°  (math, CCW from East)")
    log.append(f"[info] theta_section = {math.degrees(theta_section):+.1f}°  (perpendicular)")

    # ------------------------------------------------------------ STEP 3
    cloud = build_beam_cloud(data, utm)
    log.append(f"[info] beam cloud points: {len(cloud['E'])}")

    # ------------------------------------------------------------ STEP 4
    centroid = (float(np.nanmean(utm[:, 0])), float(np.nanmean(utm[:, 1])))
    s_pts_TL, d_off_TL = project_to_axis(cloud["E"], cloud["N"], centroid, theta_section)
    s_track_TL, _      = project_to_axis(utm[:, 0], utm[:, 1], centroid, theta_section)
    s_pts   = -s_pts_TL
    s_track = -s_track_TL
    ux_LR = -math.cos(theta_section)
    uy_LR = -math.sin(theta_section)
    s_data_min = float(np.nanmin(s_pts))
    s_data_max = float(np.nanmax(s_pts))

    edge_l = data["edge_left"]  if not args.no_edge_extrapolation else 0.0
    edge_r = data["edge_right"] if not args.no_edge_extrapolation else 0.0
    log.append(f"[info] edges (final mapping): LEFT (s=0)={edge_l:.2f}m, "
               f"RIGHT (s=s_max)={edge_r:.2f}m")

    Ec, Nc = centroid
    theta_section_final = math.atan2(uy_LR, ux_LR)
    # axis_origin (geographic location of profile s = 0) is set AFTER build_profile,
    # from the shift it returns, so it stays consistent with the VB-referenced banks.

    # --- perpendicular-offset penalty -------------------------------------
    if args.no_offset_weighting or args.offset_scale <= 0:
        offset_penalty = np.ones_like(cloud["weight"])
        log.append("[info] perpendicular-offset penalty: DISABLED")
    else:
        d_perp = np.abs(d_off_TL)
        offset_penalty = np.exp(-0.5 * (d_perp / args.offset_scale) ** 2)
        log.append(f"[info] perpendicular-offset penalty: scale={args.offset_scale:.2f} m, "
                   f"mean factor={float(np.nanmean(offset_penalty)):.2f}, "
                   f"max |offset|={float(np.nanmax(d_perp)):.2f} m")
    eff_weight = cloud["weight"] * offset_penalty

    # ------------------------------------------------------------ STEP 4b
    # River chainage (progresiva) + brazo, from the route (if given).
    progresiva_m = None
    brazo = ""
    dist_axis = float("nan")
    if route is not None:
        tr_c = Transformer.from_crs(utm_crs, route.crs, always_xy=True)
        cx_r, cy_r = tr_c.transform(centroid[0], centroid[1])
        progresiva_m, brazo, dist_axis = route.chainage(cx_r, cy_r)
        log.append(f"[info] progresiva: {progresiva_m:.2f} m "
                   f"({format_progresiva(progresiva_m)})"
                   + (f"  [{brazo}]" if brazo else "")
                   + f"   dist to axis = {dist_axis:.1f} m")
    else:
        log.append("[info] progresiva: not computed (no --centerline)")

    # ------------------------------------------------------------ STEP 4c
    # Water-surface elevation at this progresiva -> absolute datum.
    ws_elev = None
    ws_note = ""
    if wsp is not None and progresiva_m is not None:
        ws_elev, ws_note = wsp.elev_at(progresiva_m)
        msg = (f"[info] water surface @ progresiva = {ws_elev:.3f} m {datum_name}"
               f"  -> bed elevations are ABSOLUTE ({datum_name})")
        if ws_note:
            msg += f"   [WARN {ws_note}]"
        log.append(msg)
    elif args.water_surface_elev:
        ws_elev = float(args.water_surface_elev)
        log.append(f"[info] water surface (constant) = {ws_elev:.3f} m")

    # ------------------------------------------------------------ STEP 5
    s_grid, depth_grid, s_shift = build_profile(
        s_pts, cloud["depth"], eff_weight,
        s_data_min, s_data_max, dx=args.dx,
        bin_half_width=args.bin_half_width,
        edge_left_dist=edge_l, edge_right_dist=edge_r,
        beam_index=cloud["beam_index"],
        vb_priority=not args.blend_beams, vb_min=args.vb_min,
    )
    bhw = args.bin_half_width if args.bin_half_width is not None else 2 * args.dx
    s_pts_final = s_pts + s_shift
    # Geographic origin of the profile (s = 0 = left bank), consistent with the
    # shift build_profile actually applied (VB-referenced when vb_priority).
    axis_origin = (Ec - s_shift * ux_LR, Nc - s_shift * uy_LR)
    mode = "blended" if args.blend_beams else f"VB-priority (vb_min={args.vb_min})"
    log.append(f"[info] grid: {len(s_grid)} nodes ({s_grid[0]:.2f} → {s_grid[-1]:.2f} m, "
               f"dx={args.dx} m, bin_half_width={bhw} m, bed={mode})")

    # ------------------------------------------------------------ STEP 6
    csv1, geo_arrays = export_csvs(
        outdir, perfil_id, progresiva_m, s_grid, depth_grid,
        axis_origin, theta_section_final, utm_crs, posgar_crs,
        water_surface_elev=args.water_surface_elev, brazo=brazo, ws_elev=ws_elev,
    )
    raw_csv = export_raw_points(outdir, perfil_id, cloud, s_pts_final, utm_crs, posgar_crs)
    log.append(f"[ok ] {csv1.name}")
    log.append(f"[ok ] {raw_csv.name}")

    geo = dict(
        s_m=s_grid, depth_m=depth_grid, bed_elev=geo_arrays[6],
        E=geo_arrays[0], N=geo_arrays[1], lon=geo_arrays[2], lat=geo_arrays[3],
        xp=geo_arrays[4], yp=geo_arrays[5],
    )

    # ------------------------------------------------------------ STEP 7
    plot1 = plot_plan_view(outdir, perfil_id, progresiva_m, utm, cloud,
                           axis_origin, theta_section_final,
                           (float(s_grid[0]), float(s_grid[-1])),
                           theta_flow=theta_flow, brazo=brazo)
    plot2 = plot_cross_section(outdir, perfil_id, progresiva_m,
                               s_pts_final, cloud["depth"], cloud["beam_index"],
                               s_grid, depth_grid,
                               brazo=brazo, ws_elev=ws_elev, datum_name=datum_name)
    log.append(f"[ok ] {plot1.name}")
    log.append(f"[ok ] {plot2.name}")

    # ------------------------------------------------------------ STEP 8
    if HAS_GPD:
        ax_path = export_axis_shp(outdir, perfil_id, progresiva_m, axis_origin,
                                  theta_section_final,
                                  (float(s_grid[0]), float(s_grid[-1])), utm_crs,
                                  brazo=brazo)
        if ax_path: log.append(f"[ok ] {ax_path.name}")
        pp = export_profile_points_shp(outdir, perfil_id, progresiva_m, geo,
                                       utm_crs, posgar_crs, brazo=brazo, ws_elev=ws_elev)
        if pp: log.append(f"[ok ] {pp.name}")
        rb = export_raw_beam_points_shp(outdir, perfil_id, cloud, s_pts_final,
                                        utm_crs, posgar_crs, brazo=brazo)
        if rb: log.append(f"[ok ] {rb.name}")
    else:
        log.append("[warn] geopandas missing — shapefiles not written")

    # ------------------------------------------------------------ STEP 9  (QA)
    acceptance_tests(theta_flow, theta_section_final, s_grid, depth_grid,
                     edge_left=edge_l, log=log)
    if route is not None and np.isfinite(dist_axis):
        width = float(s_grid[-1] - s_grid[0])
        ok_axis = dist_axis <= max(50.0, 0.75 * width)
        log.append(f"[QA] transect centroid close to river axis? "
                   f"{'PASS' if ok_axis else 'WARN'} (dist={dist_axis:.1f} m, "
                   f"section width={width:.1f} m)")
    if ws_note:
        log.append(f"[QA] water-surface interpolation: WARN — {ws_note}")

    # write the per-transect log to disk
    try:
        (outdir / "process_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    except Exception:
        pass
    if not quiet:
        print()
        for line in log:
            print(line)
        print()
        print(f"[done] outputs in: {outdir}")

    # bed elevation stats for the survey aggregates
    bed = geo["bed_elev"]
    interior = depth_grid > 1e-6
    thalweg_elev = float(np.nanmin(bed[interior])) if interior.any() else float(np.nanmin(bed))
    max_depth = float(np.nanmax(depth_grid))

    return dict(
        perfil=perfil_id, outdir=outdir,
        progresiva=progresiva_m, brazo=brazo, dist_axis=dist_axis,
        ws_elev=ws_elev, ws_note=ws_note,
        s_grid=s_grid, depth_grid=depth_grid, bed_elev=bed,
        xp=geo["xp"], yp=geo["yp"],
        width_m=float(s_grid[-1] - s_grid[0]),
        thalweg_elev=thalweg_elev, max_depth=max_depth,
        theta_flow=theta_flow, n_samples=int(len(data["vb_depth"])),
        utm_crs=utm_crs, posgar_crs=posgar_crs, log=log,
    )


# ============================================================================ #
#  SURVEY-LEVEL AGGREGATES (multi-file mode)
# ============================================================================ #

def write_survey_index(results, resumen_dir: Path, datum_name: str):
    """One row per transect: progresiva, brazo, cotas, ancho, etc."""
    path = resumen_dir / "survey_index.csv"
    rs = sorted(results, key=lambda r: (float("inf") if r["progresiva"] is None
                                        else r["progresiva"]))
    with open(path, "w", encoding="utf-8") as f:
        f.write("perfil,progresiva_m,progresiva,brazo,ws_elev_m,thalweg_elev_m,"
                "max_depth_m,width_m,n_samples,dist_axis_m,theta_flow_deg,ws_note\n")
        for r in rs:
            prog = "" if r["progresiva"] is None else f"{r['progresiva']:.3f}"
            progf = "" if r["progresiva"] is None else format_progresiva(r["progresiva"])
            ws = "" if r["ws_elev"] is None else f"{r['ws_elev']:.3f}"
            th = "" if r["ws_elev"] is None else f"{r['thalweg_elev']:.3f}"
            da = "" if not np.isfinite(r["dist_axis"]) else f"{r['dist_axis']:.1f}"
            f.write(f"{r['perfil']},{prog},{progf},{r['brazo']},{ws},{th},"
                    f"{r['max_depth']:.3f},{r['width_m']:.2f},{r['n_samples']},"
                    f"{da},{math.degrees(r['theta_flow']):.1f},{r['ws_note']}\n")
    return path


def write_survey_profiles(results, resumen_dir: Path):
    """Every gridded profile point of every transect (long, tidy CSV)."""
    path = resumen_dir / "survey_profiles_all.csv"
    rs = sorted(results, key=lambda r: (float("inf") if r["progresiva"] is None
                                        else r["progresiva"]))
    with open(path, "w", encoding="utf-8") as f:
        f.write("perfil,progresiva_m,brazo,s_m,depth_m,ws_elev_m,bed_elev_m,"
                "x_posgar07,y_posgar07\n")
        for r in rs:
            prog = "" if r["progresiva"] is None else f"{r['progresiva']:.3f}"
            ws = "" if r["ws_elev"] is None else f"{r['ws_elev']:.3f}"
            sg, dg, be = r["s_grid"], r["depth_grid"], r["bed_elev"]
            xp, yp = r["xp"], r["yp"]
            for i in range(len(sg)):
                f.write(f"{r['perfil']},{prog},{r['brazo']},{sg[i]:.4f},{dg[i]:.4f},"
                        f"{ws},{be[i]:.4f},{xp[i]:.4f},{yp[i]:.4f}\n")
    return path


_BRAZO_COLOR = {"": "0.35", "Brazo Norte": "tab:blue", "Brazo Sur": "tab:red"}


def plot_survey_planview(results, route, resumen_dir: Path):
    """Centerline + every section axis, coloured by brazo, labelled by progresiva."""
    fig, ax = plt.subplots(figsize=(11, 12))
    # centerline (route is in its own CRS = EPSG:5344, same as xp/yp)
    if route is not None:
        xs, ys = route.main.xy
        ax.plot(xs, ys, "-", color="0.6", lw=1.4, label="Eje del río (tronco)")
        for b in route.branches:
            bx, by = b["geom"].xy
            ax.plot(bx, by, color=_BRAZO_COLOR.get(b["label"], "0.4"),
                    lw=1.4, ls="--", label=b["label"] + " (eje)")
    seen = set()
    for r in sorted(results, key=lambda r: (float("inf") if r["progresiva"] is None
                                            else r["progresiva"])):
        xp, yp = r["xp"], r["yp"]
        col = _BRAZO_COLOR.get(r["brazo"], "0.35")
        lbl = f"Secciones ({r['brazo'] or 'cauce único'})"
        ax.plot([xp[0], xp[-1]], [yp[0], yp[-1]], "-", color=col, lw=2.2,
                label=lbl if lbl not in seen else None)
        seen.add(lbl)
        mx, my = 0.5 * (xp[0] + xp[-1]), 0.5 * (yp[0] + yp[-1])
        if r["progresiva"] is not None:
            ax.annotate(format_progresiva(r["progresiva"]).split("+")[0] + "k",
                        (mx, my), fontsize=7, color=col,
                        xytext=(3, 3), textcoords="offset points")

    # Zoom to the surveyed reach (the full centerline is still drawn for context
    # but would make 40-m sections invisible on a multi-km river).
    all_x = np.concatenate([r["xp"] for r in results])
    all_y = np.concatenate([r["yp"] for r in results])
    if all_x.size:
        pad = 0.12 * max(float(np.ptp(all_x)), float(np.ptp(all_y)), 100.0)
        ax.set_xlim(all_x.min() - pad, all_x.max() + pad)
        ax.set_ylim(all_y.min() - pad, all_y.max() + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("Este POSGAR07 f2 [m]")
    ax.set_ylabel("Norte POSGAR07 f2 [m]")
    ax.set_title("Vista en planta del relevamiento — ejes de secciones por progresiva y brazo")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    path = resumen_dir / "survey_plan_view.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_long_profile(results, wsp, resumen_dir: Path, datum_name: str):
    """Longitudinal profile: water surface + thalweg vs progresiva."""
    rs = [r for r in results if r["progresiva"] is not None]
    rs.sort(key=lambda r: r["progresiva"])
    if not rs:
        return None
    absolute = any(r["ws_elev"] is not None for r in rs)
    fig, ax = plt.subplots(figsize=(12, 6))

    if absolute and wsp is not None:
        ax.plot(wsp.prog, wsp.elev, "-", color="navy", lw=1.6,
                label="Perfil de pelo de agua (interpolado)")
        ax.plot(wsp.raw_prog, wsp.raw_elev, ".", color="tab:cyan", ms=5,
                alpha=0.7, label="Puntos de pelo de agua medidos")

    prog = [r["progresiva"] for r in rs]
    if absolute:
        thal = [r["thalweg_elev"] for r in rs]
        wss  = [r["ws_elev"] if r["ws_elev"] is not None else np.nan for r in rs]
        for r in rs:
            if r["ws_elev"] is not None:
                ax.plot([r["progresiva"], r["progresiva"]],
                        [r["thalweg_elev"], r["ws_elev"]], "-",
                        color="0.7", lw=0.8, zorder=1)
        ax.plot(prog, thal, "o", color="#7a5230", ms=6,
                label="Thalweg (cota mínima del lecho)", zorder=3)
        ax.plot(prog, wss, "_", color="navy", ms=10, mew=2,
                label="Pelo de agua en cada sección", zorder=2)
        ax.set_ylabel(f"Cota {datum_name} [m]")
    else:
        thal = [-r["max_depth"] for r in rs]
        ax.plot(prog, thal, "o-", color="#7a5230", ms=6,
                label="Thalweg (profundidad máx., relativa)")
        ax.axhline(0, color="navy", lw=1, ls="--", alpha=0.7, label="Pelo de agua (rel.)")
        ax.set_ylabel("Cota relativa al pelo de agua [m]")

    # brazo markers on the x baseline
    for r in rs:
        if r["brazo"]:
            ax.annotate(r["brazo"].replace("Brazo ", "B."), (r["progresiva"], thal_min(thal)),
                        fontsize=7, color=_BRAZO_COLOR.get(r["brazo"], "0.4"),
                        ha="center", va="top", xytext=(0, -2), textcoords="offset points")

    ax.set_xlabel("Progresiva [m]")
    ax.set_title("Perfil longitudinal del relevamiento — pelo de agua y thalweg")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    path = resumen_dir / "survey_long_profile.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def thal_min(vals):
    v = [x for x in vals if np.isfinite(x)]
    return min(v) if v else 0.0


def export_survey_shapefiles(results, resumen_dir: Path, datum_name: str):
    """Combined survey_axes.shp + survey_profile_points.shp in EPSG:5344."""
    if not HAS_GPD:
        return None, None
    crs5344 = CRS.from_epsg(5344)
    rs = sorted(results, key=lambda r: (float("inf") if r["progresiva"] is None
                                        else r["progresiva"]))
    # axes
    axg, arows = [], []
    for r in rs:
        xp, yp = r["xp"], r["yp"]
        axg.append(LineString([(xp[0], yp[0]), (xp[-1], yp[-1])]))
        arows.append(dict(
            perfil=r["perfil"],
            progr_m=np.nan if r["progresiva"] is None else float(r["progresiva"]),
            brazo=r["brazo"],
            ws_elev_m=np.nan if r["ws_elev"] is None else float(r["ws_elev"]),
            thalweg_m=np.nan if r["ws_elev"] is None else float(r["thalweg_elev"]),
            max_dep_m=float(r["max_depth"]),
            width_m=float(r["width_m"]),
        ))
    ax_path = resumen_dir / "survey_axes.shp"
    gpd.GeoDataFrame(arows, geometry=axg, crs=crs5344).to_file(ax_path)

    # points
    pg, prows = [], []
    for r in rs:
        sg, dg, be = r["s_grid"], r["depth_grid"], r["bed_elev"]
        xp, yp = r["xp"], r["yp"]
        for i in range(len(sg)):
            pg.append(Point(xp[i], yp[i]))
            prows.append(dict(
                perfil=r["perfil"],
                progr_m=np.nan if r["progresiva"] is None else float(r["progresiva"]),
                brazo=r["brazo"], s_m=float(sg[i]), depth_m=float(dg[i]),
                ws_elev_m=np.nan if r["ws_elev"] is None else float(r["ws_elev"]),
                bed_elev_m=float(be[i]),
            ))
    pp_path = resumen_dir / "survey_profile_points.shp"
    gpd.GeoDataFrame(prows, geometry=pg, crs=crs5344).to_file(pp_path)
    return ax_path, pp_path


# ============================================================================ #
#  MAIN
# ============================================================================ #

def gather_matfiles(inputs):
    """Expand CLI inputs (files, folders, globs) into a sorted list of .mat paths."""
    out = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            out.extend(sorted(p.glob("*.mat")))
        elif any(ch in item for ch in "*?[]"):
            out.extend(sorted(Path(x) for x in glob.glob(item)))
        else:
            out.append(p)
    # de-dup preserving order
    seen, uniq = set(), []
    for p in out:
        rp = p.resolve()
        if rp not in seen and rp.suffix.lower() == ".mat":
            seen.add(rp); uniq.append(p)
    return uniq


# ============================================================================ #
#  CONFIG FILE (campaña.ini)  +  TRACEABILITY RECORD
# ============================================================================ #

# argparse dest -> type converter for values read from the INI config file.
_CONFIG_TYPES = {
    "outdir": str, "survey_name": str, "dx": float, "bin_half_width": float,
    "water_surface_elev": float, "no_edge_extrapolation": bool,
    "blend_beams": bool, "vb_min": int, "centerline": str,
    "chainage_offset": float, "chainage_reverse": bool,
    "water_surface_csv": str, "ws_east_col": str, "ws_north_col": str,
    "ws_elev_col": str, "ws_crs": str, "datum_name": str,
    "offset_scale": float, "no_offset_weighting": bool,
}
# Config keys that are file/dir paths -> resolved relative to the config file.
_CONFIG_PATHS = {"outdir", "centerline", "water_surface_csv"}


def _to_bool(s) -> bool:
    return str(s).strip().lower() in ("1", "true", "yes", "on", "si", "sí")


def _find_section(cp, names):
    low = {s.lower(): s for s in cp.sections()}
    for n in names:
        if n in low:
            return low[n]
    return None


def load_config(path):
    """
    Read an INI config file (campaña.ini). Returns (params, campania):
      params   -> dict of argparse dest -> typed value (for parser.set_defaults)
      campania -> dict of free-form metadata (for the traceability record)

    Section [procesamiento] holds parameters; [campania]/[campaña] holds metadata.
    Parameter keys may use '-' or '_'. File/dir paths are resolved relative to the
    config file's folder, so the config can live inside the campaign folder and
    use paths relative to it. Command-line options override anything set here.
    """
    cfgpath = Path(path)
    if not cfgpath.exists():
        sys.exit(f"[error] config file not found: {cfgpath}")
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    try:
        cp.read(cfgpath, encoding="utf-8")
    except Exception as e:
        sys.exit(f"[error] could not parse config file: {e}")

    base = cfgpath.resolve().parent
    params, inputs = {}, []

    proc = _find_section(cp, ["procesamiento", "proceso", "processing",
                              "parametros", "parámetros"])
    if proc:
        for raw_key, val in cp.items(proc):
            if val is None or str(val).strip() == "":
                continue
            key = raw_key.strip().lower().replace("-", "_")
            if key in ("inputs", "matfiles", "input", "entradas"):
                items = [s.strip() for line in str(val).splitlines()
                         for s in line.split(",") if s.strip()]
                for it in items:
                    p = Path(it)
                    inputs.append(str(p if p.is_absolute() else (base / p)))
                continue
            if key not in _CONFIG_TYPES:
                continue                      # ignore unknown keys silently
            typ = _CONFIG_TYPES[key]
            try:
                conv = _to_bool(val) if typ is bool else typ(str(val).strip())
            except (TypeError, ValueError):
                sys.exit(f"[error] config: bad value for '{raw_key}': {val!r}")
            if key in _CONFIG_PATHS and conv:
                pp = Path(conv)
                conv = str(pp if pp.is_absolute() else (base / pp))
            params[key] = conv
    if inputs:
        params["matfiles"] = inputs

    camp = _find_section(cp, ["campania", "campaña", "campanha",
                              "metadatos", "metadata"])
    campania = dict(cp.items(camp)) if camp else {}
    return params, campania


def _git_commit(start_path):
    """Best-effort short git commit of the script's repo (None if unavailable)."""
    try:
        base = Path(start_path).resolve().parent
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=str(base), capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None


def write_run_record(out_dir, args, campania, matfiles, extra_lines=None):
    """Write a human-readable traceability record of this run to procesamiento.txt."""
    lines = ["REGISTRO DE PROCESAMIENTO — process_adcp_bathimetric",
             f"script version : v{__version__}"]
    commit = _git_commit(sys.argv[0] if sys.argv and sys.argv[0] else __file__)
    if commit:
        lines.append(f"git commit     : {commit}")
    lines.append(f"fecha de corrida: "
                 f"{datetime.datetime.now().isoformat(timespec='seconds')}")
    if campania:
        lines += ["", "[campaña]"]
        lines += [f"  {k}: {v}" for k, v in campania.items()]
    lines += ["", "[entradas]", f"  transectas ({len(matfiles)}):"]
    lines += [f"    - {Path(m).resolve()}" for m in matfiles]
    lines += [f"  centerline        : {args.centerline}",
              f"  water-surface CSV : {args.water_surface_csv}",
              "", "[parámetros]"]
    for key in ("dx", "bin_half_width", "vb_min", "blend_beams",
                "chainage_offset", "chainage_reverse", "ws_crs", "datum_name",
                "offset_scale", "no_offset_weighting", "no_edge_extrapolation",
                "water_surface_elev", "survey_name", "outdir"):
        lines.append(f"  {key}: {getattr(args, key, None)}")
    if extra_lines:
        lines += [""] + list(extra_lines)
    path = Path(out_dir) / "procesamiento.txt"
    try:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception:
        pass
    return path


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("matfiles", nargs="*",
                   help="One or more .mat files, a folder of .mat files, or a glob. "
                        "Two or more transects triggers SURVEY mode. May instead be "
                        "supplied via 'inputs' in a --config file.")
    p.add_argument("--config", default=None,
                   help="INI config file (campaña.ini) providing inputs, parameters "
                        "and campaign metadata. Command-line options override it.")
    p.add_argument("--outdir", default=None,
                   help="Output directory (single: ./out_<basename>/; "
                        "survey: ./<survey-name>/).")
    p.add_argument("--survey-name", default=None,
                   help="Name for the survey output folder (default: 'salida_<first perfil>').")
    p.add_argument("--dx", type=float, default=0.5,
                   help="Profile grid spacing in metres (default 0.50).")
    p.add_argument("--bin-half-width", type=float, default=None,
                   help="Half-width [m] for s-bin sampling (default: 2*dx).")
    p.add_argument("--water-surface-elev", type=float, default=0.0,
                   help="Constant water-surface elevation [m] used only when no "
                        "--water-surface-csv is given (bed_elev = WSE - depth).")
    p.add_argument("--no-edge-extrapolation", action="store_true",
                   help="Disable bank extrapolation to depth=0 from Setup.Edges_*.")
    p.add_argument("--blend-beams", action="store_true",
                   help="Blend vertical + slant beams together for the bed (older "
                        "behaviour). Default: the nadir/vertical beam defines the bed "
                        "where present, slant beams only fill near-bank/edge gaps.")
    p.add_argument("--vb-min", type=int, default=1,
                   help="Min. vertical-beam points in a node window to use VB alone "
                        "(default 1; higher = fall back to slant beams sooner).")
    # --- river chainage (progresiva) ---
    p.add_argument("--centerline", default=None,
                   help="River-axis vector file (shp/gpkg/geojson). Branched "
                        "centerlines (islands) are handled automatically.")
    p.add_argument("--chainage-offset", type=float, default=0.0,
                   help="Base chainage [m] added to every progresiva "
                        "(e.g. the km of the centerline's upstream origin).")
    p.add_argument("--chainage-reverse", action="store_true",
                   help="Measure chainage from the centerline's LAST vertex "
                        "(ignored if a water-surface CSV is given, since the route "
                        "is then oriented automatically by the water-surface slope).")
    # --- water surface / absolute datum ---
    p.add_argument("--water-surface-csv", default=None,
                   help="CSV of surveyed water-surface points (X, Y, elevation). "
                        "Enables ABSOLUTE bed elevations in the vertical datum.")
    p.add_argument("--ws-east-col", default="East", help="WS CSV easting column (default 'East').")
    p.add_argument("--ws-north-col", default="North", help="WS CSV northing column (default 'North').")
    p.add_argument("--ws-elev-col", default="H_correg", help="WS CSV elevation column (default 'H_correg').")
    p.add_argument("--ws-crs", default="EPSG:5344",
                   help="CRS of the WS CSV coordinates (default EPSG:5344 = POSGAR07 f2).")
    p.add_argument("--datum-name", default=DATUM_NAME_DEFAULT,
                   help=f"Vertical datum label for plots/columns (default '{DATUM_NAME_DEFAULT}').")
    # --- perpendicular-offset penalty ---
    p.add_argument("--offset-scale", type=float, default=OFFSET_SCALE_DEFAULT,
                   help=f"Gaussian scale [m] for the perpendicular-offset penalty "
                        f"(default {OFFSET_SCALE_DEFAULT}; larger = gentler).")
    p.add_argument("--no-offset-weighting", action="store_true",
                   help="Disable the perpendicular-offset penalty entirely.")
    return p


def main():
    # Two-phase parse: read --config first so it can supply defaults that the
    # command line then overrides (precedence: CLI > config file > built-in).
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=None)
    pre_args, _ = pre.parse_known_args()

    parser = build_parser()
    campania = {}
    if pre_args.config:
        cfg_params, campania = load_config(pre_args.config)
        parser.set_defaults(**cfg_params)
    args = parser.parse_args()
    datum_name = args.datum_name

    if not args.matfiles:
        sys.exit("[error] no inputs given (pass .mat/folder/glob on the command "
                 "line, or set 'inputs' in a --config file)")
    matfiles = gather_matfiles(args.matfiles)
    if not matfiles:
        sys.exit("[error] no .mat files found in the given inputs")
    for m in matfiles:
        if not m.exists():
            sys.exit(f"[error] file not found: {m}")

    setup_log: list[str] = [f"[info] process_adcp_bathimetric v{__version__}"]

    # --- water-surface points (parsed once; used for orientation + profile) ---
    ws_points_native = None
    if args.water_surface_csv:
        try:
            ws_points_native = read_ws_csv(
                args.water_surface_csv, args.ws_east_col, args.ws_north_col, args.ws_elev_col)
        except Exception as e:
            sys.exit(f"[error] could not read water-surface CSV: {e}")

    # --- river route (branched centerline) ---
    route = None
    if args.centerline:
        if not HAS_GPD:
            sys.exit("[error] --centerline needs geopandas/shapely installed")
        # peek at the centerline CRS to reproject WS points for orientation
        try:
            route_crs = gpd.read_file(args.centerline).crs
        except Exception as e:
            sys.exit(f"[error] could not read centerline: {e}")
        ws_xy_elev = None
        if ws_points_native is not None:
            if route_crs is not None and CRS.from_user_input(args.ws_crs) != route_crs:
                tr = Transformer.from_crs(args.ws_crs, route_crs, always_xy=True)
                ws_xy_elev = [(*tr.transform(x, y), h) for (x, y, h) in ws_points_native]
            else:
                ws_xy_elev = list(ws_points_native)
        try:
            route = RiverRoute(args.centerline, crs_override=None,
                               reverse=args.chainage_reverse, offset=args.chainage_offset,
                               ws_xy_elev=ws_xy_elev, log=setup_log)
        except Exception as e:
            sys.exit(f"[error] could not build river route: {e}")

    # --- water-surface profile (needs the route to project onto) ---
    wsp = None
    if ws_points_native is not None:
        if route is None:
            setup_log.append("[warn] --water-surface-csv ignored: needs --centerline "
                             "to convert points to progresivas")
        else:
            # reproject WS points into route CRS, then build the profile
            if CRS.from_user_input(args.ws_crs) != route.crs:
                tr = Transformer.from_crs(args.ws_crs, route.crs, always_xy=True)
                pts = [(*tr.transform(x, y), h) for (x, y, h) in ws_points_native]
            else:
                pts = list(ws_points_native)
            try:
                wsp = WaterSurfaceProfile(route, pts, log=setup_log)
            except Exception as e:
                setup_log.append(f"[warn] could not build water-surface profile: {e}")
                wsp = None

    if setup_log:
        print("\n".join(setup_log))
        print()

    # =============================== single transect ======================= #
    if len(matfiles) == 1:
        r = process_one(matfiles[0], args, route=route, wsp=wsp, datum_name=datum_name)
        if r is not None:
            rec = write_run_record(r["outdir"], args, campania, matfiles)
            print(f"[ok ] {rec.name}")
        return

    # =============================== survey mode =========================== #
    survey_name = args.survey_name or f"salida_{matfiles[0].stem}"
    survey_root = Path(args.outdir) if args.outdir else Path.cwd() / survey_name
    survey_root.mkdir(parents=True, exist_ok=True)
    resumen = survey_root / "_resumen"
    resumen.mkdir(parents=True, exist_ok=True)

    print(f"[survey] {len(matfiles)} transects -> {survey_root}")
    results = []
    for i, m in enumerate(matfiles, 1):
        r = process_one(m, args, route=route, wsp=wsp, datum_name=datum_name,
                        survey_outdir=survey_root, quiet=True)
        if r is None:
            print(f"  [{i}/{len(matfiles)}] {m.name}: FAILED (see log)")
            continue
        prog = "n/a" if r["progresiva"] is None else format_progresiva(r["progresiva"])
        ws = "" if r["ws_elev"] is None else f", pelo={r['ws_elev']:.3f}m"
        th = "" if r["ws_elev"] is None else f", thalweg={r['thalweg_elev']:.3f}m"
        print(f"  [{i}/{len(matfiles)}] {r['perfil']}: progr={prog}"
              f"{(' ['+r['brazo']+']') if r['brazo'] else ''}"
              f", ancho={r['width_m']:.1f}m, prof.máx={r['max_depth']:.2f}m{ws}{th}"
              f"{'  [WARN '+r['ws_note']+']' if r['ws_note'] else ''}")
        results.append(r)

    if not results:
        sys.exit("[error] no transect processed successfully")

    # aggregates
    idx = write_survey_index(results, resumen, datum_name)
    prof = write_survey_profiles(results, resumen)
    print(f"[ok ] {idx.relative_to(survey_root)}")
    print(f"[ok ] {prof.relative_to(survey_root)}")
    pv = plot_survey_planview(results, route, resumen)
    print(f"[ok ] {pv.relative_to(survey_root)}")
    lp = plot_long_profile(results, wsp, resumen, datum_name)
    if lp: print(f"[ok ] {lp.relative_to(survey_root)}")
    if HAS_GPD:
        sax, spp = export_survey_shapefiles(results, resumen, datum_name)
        if sax: print(f"[ok ] {sax.relative_to(survey_root)}")
        if spp: print(f"[ok ] {spp.relative_to(survey_root)}")

    # traceability record for the whole survey
    summary = [f"[resumen] {len(results)} transectas procesadas:"]
    for r in sorted(results, key=lambda r: (float("inf") if r["progresiva"] is None
                                            else r["progresiva"])):
        prog = "n/a" if r["progresiva"] is None else format_progresiva(r["progresiva"])
        summary.append(f"    {r['perfil']}: progr={prog}"
                       f"{(' ['+r['brazo']+']') if r['brazo'] else ''}"
                       f", ancho={r['width_m']:.1f} m, prof.máx={r['max_depth']:.2f} m"
                       + ("" if r["ws_elev"] is None
                          else f", pelo={r['ws_elev']:.3f} m, thalweg={r['thalweg_elev']:.3f} m")
                       + (f"  [WARN {r['ws_note']}]" if r["ws_note"] else ""))
    rec = write_run_record(resumen, args, campania, matfiles, extra_lines=summary)
    print(f"[ok ] {rec.relative_to(survey_root)}")

    print(f"\n[done] survey outputs in: {survey_root}")
    print(f"       per-transect folders + consolidated results in: {resumen}")


if __name__ == "__main__":
    main()
