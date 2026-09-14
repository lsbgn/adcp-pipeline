#!/usr/bin/env python3
"""
process_adcp_bathimetric.py
================================================================================
ADCP CROSS-SECTION EXTRACTION — Bathymetry-only pipeline (v6.1)

Builds a hydraulically consistent 1D bathymetric cross-section perpendicular
to the mean flow, from a SonTek M9/S5 .mat file (RiverSurveyor Live export).

Velocities are NOT exported. Only the bed profile and supporting QA artifacts.

Changes in 6.1 (bug-fix release; numeric core untouched — build_profile,
build_beam_cloud, flow direction and depth filters are identical to 6.0):
    * NETWORK PATHS INSTEAD OF A TRUNK. v6.0 stitched the surveyed rivers into
      ONE chain and interpolated the water surface on it. A survey touching two
      tributaries of the same confluence (Neuquen + Limay) had two upstream
      heads, the chain was undefined, every GNSS point was discarded and every
      bed came out RELATIVE, silently (campaign 2026-08-25). The centerline is
      now read as a NETWORK: each river's downstream end is linked to the river
      it flows into, and each head-to-outlet PATH (Neuquen>Negro, Limay>Negro)
      carries a continuous chainage. The reach below a confluence belongs to
      every path through it, so no river has to be chosen as trunk.
    * WATER SURFACE PER RIVER, MEETING AT THE CONFLUENCE (WaterSurfaceModel).
      Each river keeps its own surface. The stage at a confluence is estimated
      from EACH tributary with data (linear interpolation along its path
      between its last value upstream and the first value downstream of the
      node: its slope continued into the common river); the estimates are
      averaged with weight 1/gap into one junction value and their spread is
      QC'd against --ws-qc-tol. Anchors are the GNSS points plus the gauges
      not bracketed by GNSS; bracketed gauges are QC only. Beyond the last
      value of a path the surface is CONTINUED WITH ITS LOCAL SLOPE (least
      squares over >= 1 km) instead of the flat v5/v6.0 clamp: a tributary
      without data takes the junction stage and the slope of the common reach.
      Extrapolations longer than --ws-extrap-tol are warnings. ws_source is one
      of survey | bridge | station | extrap | none. _resumen/water_surface_qc.csv
      lists every GNSS point, gauge and junction: river, km, distance to axis,
      alternative river, role and reason.
    * TRANSECTS LOCATED BY SECTION CROSSING. River and km are where the section
      line (principal axis of the GPS track, extended 15 %) crosses a
      centerline, as in HEC-RAS — not the axis nearest to the track centroid,
      which put a Negro section a few metres below the confluence on the Limay.
      --transect-locate centroid restores the v6.0 rule; the INI [rios] section
      forces the river of a transect; an optional 'rio' column in the WS CSV
      forces the river of a GNSS point (points with two candidate rivers are
      flagged AMBIGUOUS).
    * STATIONS. CSV registries (X/Y in --ws-crs) are read directly — pyogrio
      returns a plain DataFrame for a CSV, which crashed with "'DataFrame' object
      has no attribute 'crs'". A station is projected on its DECLARED river;
      stations off axis, beyond the digitised reach or without gauge_zero are
      excluded with the reason. Time-series readings are matched to the transect
      time (representative_time now knows the SonTek 2000-01-01 epoch).
    * NO SILENT RELATIVE OUTPUT. A transect without water surface logs a [warn]
      and gets ws_source=none; if a water-surface source is configured and no
      transect can get one, the run stops (--allow-relative to continue).
    * CONFIG. The INI key schema is derived from the CLI parser: the v6 keys
      (depth-ref, composite, bt-avg, group-tol, ...) were silently IGNORED in
      campaign INIs. Choices are validated; unknown keys and sections are
      warned; [progresivas] river names match the centerline accent- and
      case-insensitively. ';' starts a comment in every section.
    * RUN LOG. The whole console is saved to procesamiento.log; procesamiento.txt
      records every effective parameter, the network / path / water-surface
      setup, every [warn] and the QA warnings of each transect.
    * FIXES. format_progresiva printed '1+1000.00' for 1999.999 m; survey_index
      .csv is written with the csv module (notes may contain commas; the new
      columns recorrido and loc_method are appended at the end); the section-
      orientation QC uses the tangent of the transect's own river; the long
      profile draws one panel per path.
    * COMPATIBILITY. RiverNetwork.build_survey_chain / survey_chainage,
      WaterSurfaceProfile and StationWS.elev_at are kept verbatim as legacy API
      (no longer used by the pipeline). survey_chainage_m in survey_index.csv
      is the continuous chainage along the transect's display path.

Changes vs v5 (inherited in v6):
    * AFORO GROUPING. Repeated transects at the same section (the usual >=4 runs
      alternating start bank) are detected and averaged into ONE profile.
      Grouping requires same river+brazo, |dchainage| <= --group-tol (15 m,
      measured along the section-normal), section axes parallel within
      --group-angle-tol (25 deg), and TEMPORAL CONTIGUITY — an uninterrupted run
      in time order, which is what separates a real aforo from a second visit
      hours later at a different stage. [grupos] in the INI forces grouping by
      hand. Averaging MERGES THE BEAM CLOUDS and runs build_profile once, rather
      than averaging already-built profiles: every point is projected from its
      own (E, N), so runs of different width, start bank and track combine with
      no resampling. On a synthetic 4-run aforo this cut RMSE 0.046 -> 0.028 m.
      Per-transect folders stay intact; campaign aggregates see only the group
      profile. New <grupo>/ outputs plus _resumen/grupos_qc.csv (node sigma, per-run
      RMS, chainage and azimuth spread, aforo discharge).
    * DEPTH REFERENCE ported to the QRev model. --blend-beams / --vb-min give way
      to --depth-ref {auto,vb,bt} (auto honours Setup.depthReference, i.e. the
      field operator's choice), --composite {on,off}, --bt-avg {idw,simple},
      --bt-geometry {footprints,ensemble} and --edge-anchor {ref,cloud}; the old
      flags remain as aliases. Beam averaging is now INVERSE-RANGE WEIGHTED
      (DepthData.average_depth): on a sloping bed the long beam looks downslope
      and biases a plain mean deep. >=2 valid beams required. build_profile also
      returns src_grid with the per-node QRev source code (1 BT, 2 VB, 3 DS,
      4 interpolated). vb_priority generalises to 'the selected reference', so
      --depth-ref vb reproduces v5 exactly.
    * DEPTH SPIKE FILTER. v5 had NO outlier rejection. --depth-filter
      {smooth,off} (default smooth) follows DepthData.filter_smooth: residuals
      against a robust smooth, running IQR, rejection threshold = max(IQR
      criterion, 5 % of depth, 0.10 m). Each beam filtered separately.
    * SECTION AZIMUTH weighting. compute_flow_direction was already a vector sum;
      it now weights each ensemble by interval x depth, giving the direction of
      the section's TOTAL UNIT DISCHARGE VECTOR (--flow-weighting
      {discharge,density,none}, default discharge; none = v5). Edge ensembles
      (System.Step 2/4, boat stopped at the bank) are excluded from the azimuth
      as QRev does; they stay in the cloud, since near-bank bed returns are
      scarce and valuable — only their crowding weight is removed
      (--density-weighting {on,off}, default on).
    * SECTION-ORIENTATION QC. New section_orientation_qc plus
      RiverNetwork.tangent_azimuth: track-vs-axis obliquity (and the % of width
      lost, 1 - cos theta) and flow-normal vs trace-normal angle. The diagnostic
      for "my sections don't look transversal". --section-orientation
      {flow,centerline} stays on flow (v5 behaviour).
    * AFORO DISCHARGE. The last value of Summary.Total_Q is the transect's
      measured discharge (cross-checked against the sum of the five components).
      survey_index.csv gains q_m3s; grupos_qc.csv gains mean, sd, CV, range and
      per-run values; the repeatability plot title shows Q and its CV. USGS
      practice is ~5 % per-transect spread about the mean. Independent of the
      bathymetric repeatability — read the two separately.
    * VERIFIED, NO CODE CHANGE: edge mapping (Edges_0 = LEFT always, regardless
      of startEdge) and transducer draft (file depths already include it) were
      both flagged as suspect and confirmed correct in v5.

Changes vs v4 (inherited in v5):
    * MULTI-RIVER centerline (RiverNetwork, replaces RiverRoute). The centerline
      shapefile is now read as a NETWORK of channels: one feature per continuous
      river+role reach, with attributes `river` (Neuquen|Limay|Negro|...), `role`
      (main|anabranch), a stable `seg_id`, and optional `label`/`flow_dir`. A
      survey may span more than one river (e.g. lower Neuquen -> confluence ->
      Rio Negro): each transect is tagged with its `river` and reported at that
      river's OFFICIAL km. Backward compatible: a centerline with no `river`/
      `role` columns collapses to the v4 single-river-with-islands behaviour.
    * PER-RIVER chainage. Internal chainage runs along each river's main from its
      upstream end (=0); the report shows km_oficial = river_offset + internal,
      with per-river offsets from the INI ([progresivas]). A single continuous
      INTERNAL axis is stitched across the confluence ALONG THE SURVEYED PATH for
      the longitudinal profile (Neuquen->Negro is resolved automatically because
      only those two carry transects; the Limay stays drawn but inert). The
      canonical-trunk `mainstem` hook (Limay->Negro) is intentionally NOT wired.
    * ANABRANCHES for any N (not just Norte/Sur). Islands are detected by topology
      (a set of >=2 channels sharing a bifurcation and a confluence) and labelled
      left-to-right facing downstream: MI / M2.. / MD (or an explicit `label`).
      Branch chainage is the projection onto the river's main axis.
    * MULTI-STATION water surface. GNSS water-surface points stay the PRIMARY
      source; a station registry (points shp: station_id, name, river, gauge_zero)
      plus a per-campaign readings CSV (spot or time-series, auto-detected) act as
      FALLBACK (no GNSS), GAP-FILL (beyond GNSS range, replacing the v4 clamp) and
      QC (|WS_gnss - WS_station| over a tolerance -> flag). ws_elev = gauge_zero +
      nivel. Each transect records ws_source (survey|station|clamp|bridge).
    * NEW survey output: <survey>_raw_beam_points (merged point cloud of every
      bed return of every transect, PointZ, EPSG:5344, Z = bed_elev SRVN16), in
      _resumen/. Written as .shp, falling back to .gpkg past the shapefile limits.

Changes vs v3 (inherited):
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
    survey_index.csv          one row per transect / aforo group (river, km, brazo,
                              cotas, ancho, ws_source, recorrido, loc_method)
    survey_profiles_all.csv   every gridded profile point of every transect
    survey_plan_view.png      centerline + all section axes, coloured by brazo
    survey_long_profile.png   water surface + thalweg, one panel per network path
    survey_axes.shp / survey_profile_points.shp / survey_raw_beam_points.shp
    grupos_qc.csv             repeatability and discharge of every aforo group
    water_surface_qc.csv      every GNSS point, gauge and confluence stage: river,
                              km, role (anchor / qc / dropped / excluded) and why
    procesamiento.txt         traceability record (parameters, setup, warnings)
    procesamiento.log         full console transcript of the run

Usage:
    # whole campaign from its INI (recommended; see campanha_ejemplo_v6.ini)
    python process_adcp_bathimetric.py --config campanha.ini

    # single transect
    python process_adcp_bathimetric.py <input.mat> [--outdir <dir>]
        [--centerline <river_network.shp>] [--water-surface-csv <ws.csv>]
        [--stations <estaciones.csv> --readings <niveles.csv>]

    # whole survey (folder of .mat, or several files)
    python process_adcp_bathimetric.py <survey_dir>/  --centerline <network.shp>
        --water-surface-csv <ws.csv> [--survey-name <name>] [--outdir <dir>]

Tested with: SonTek RiverSurveyor M9, .mat v5 export
================================================================================
"""

from __future__ import annotations

import argparse
import bisect
import configparser
import csv
import datetime
import glob
import io
import math
import re
import subprocess
import sys
import traceback
import warnings
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
__version__ = "6.1"

# 6.1: water-surface extrapolation beyond the last value of a path (slope
# continuation) is reported as a WARNING only past this distance [m]; shorter
# ones (e.g. a transect a few tens of metres up a tributary from the confluence)
# are informative notes. Overridable with --ws-extrap-tol.
WS_EXTRAP_TOL_DEFAULT = 200.0

# 6.1: the slope used to continue the water surface beyond the last value of a
# path is a least-squares fit over the values nearest to that end, taken until
# they span at least this distance [m] (and at least two values).
WS_SLOPE_WINDOW = 1000.0

# 6.1: a GNSS water-surface point whose second-nearest river axis is within this
# margin [m] of the nearest one is flagged as AMBIGUOUS in water_surface_qc.csv
# (typical at confluences); a 'rio' column in the CSV forces its river.
WS_AMBIG_MARGIN = 25.0

# 6.1: a gauge or GNSS point farther than this [m] from its river's axis (or
# projected onto an axis END, i.e. beyond the digitised reach, by more than
# STATION_END_TOL) is not used for the water surface: its chainage would be wrong.
STATION_AXIS_TOL = 250.0
STATION_END_TOL = 50.0

# 6.1: transect location by SECTION CROSSING — the track's principal axis is
# extended by max(MIN, FRAC x track length) on each side before intersecting the
# centerlines, so a track that stops short of the channel axis still crosses it.
TRACK_LOC_MARGIN_MIN = 15.0
TRACK_LOC_MARGIN_FRAC = 0.15

# 6.1: rise [m] between consecutive water-surface anchors, going downstream,
# tolerated before it is flagged (GNSS noise is a few cm).
WS_MONO_TOL = 0.05

# Tolerance [m] for deciding that a line's endpoint lies ON a river's main path,
# i.e. that the line is a side branch (island split) rather than a disjoint reach.
BRANCH_SNAP_TOL = 15.0

# Tolerance [m] for treating two centerline endpoints as the SAME network node
# (river-to-river connection at a confluence, branch split/rejoin). The linework
# must be snapped in QGIS; this only absorbs sub-metre rounding.
NODE_SNAP_TOL = 2.0

# River name used when the centerline has no `river` attribute (v4 compatibility:
# the whole file is one river with optional islands).
DEFAULT_RIVER_NAME = "río"

# Attribute names recognised on the centerline (case-insensitive). Everything is
# optional; missing attributes are derived (see RiverNetwork).
CL_ATTR_RIVER = ("river", "rio", "río")
CL_ATTR_ROLE  = ("role", "rol", "tipo")
CL_ATTR_SEGID = ("seg_id", "segid", "id")
CL_ATTR_LABEL = ("label", "etiqueta", "nombre", "name")
CL_ATTR_FLOW  = ("flow_dir", "flowdir", "sentido")

# Roles.
ROLE_MAIN = "main"
ROLE_ANAB = "anabranch"


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

    # v6: RiverSurveyor tags every ensemble with a Step code. QRev uses
    # Step == 3 as the "in-transect" set; 2 and 4 are the stationary EDGE
    # ensembles measured at each bank before/after the crossing. Those are
    # spatially piled up (measured 3-9x denser than the crossing) and their
    # water velocity is bank-influenced, so they must not drive the section
    # azimuth. Missing Step -> everything is treated as in-transect.
    _step = get_field(sysd, "Step", default=None)
    step = (np.asarray(_step, dtype=float) if _step is not None
            else np.full(len(vb_depth), 3.0))

    # v6: Summary.Depth is the reference depth RiverSurveyor itself used.
    summary_depth = get_field(summ, "Depth")
    if summary_depth is not None:
        summary_depth = np.asarray(summary_depth, dtype=float)
    bt_depth = get_field(bt, "BT_Depth")          # SonTek's own 4-beam average
    if bt_depth is not None:
        bt_depth = np.asarray(bt_depth, dtype=float)

    mean_vel   = get_field(summ, "Mean_Vel")     # (n, 2) east, north
    if mean_vel is not None:
        mean_vel = np.asarray(mean_vel, dtype=float)

    # Edge distances from Setup (used for bank extrapolation).
    # VERIFIED against QRev (Classes/TransectData.py, SonTek reader): QRev binds
    # Edges_0__DistanceToBank to the LEFT edge object and Edges_1 to the RIGHT
    # one UNCONDITIONALLY — startEdge is used only to decide which System.Step
    # code (2 or 4) counts as each bank's edge ensembles, never to swap the
    # distances. The mapping below is therefore correct as-is.
    # `or 0.0` is safe here because a 0.0 m edge is meaningless anyway; but for
    # startEdge a legitimate 0 (= Left) must NOT be turned into 1.
    edge_left  = float(get_field(setp, "Edges_0__DistanceToBank", default=0.0) or 0.0)
    edge_right = float(get_field(setp, "Edges_1__DistanceToBank", default=0.0) or 0.0)
    _se        = get_field(setp, "startEdge", default=1)
    start_edge = int(_se) if _se is not None else 1   # 0 = Left, 1 = Right

    # v6: transducer draft. VERIFIED against QRev (Classes/DepthData.py): the
    # depths stored in the SonTek .mat ALREADY INCLUDE the draft — QRev keeps
    # sensorDepth only so a user can later change the draft by difference. So
    # this value must NOT be added to VB_Depth / BT_Beam_Depth. It is carried
    # here for the IDW 4-beam average (which needs the range from the
    # transducer) and for the run record.
    _sd = get_field(setp, "sensorDepth", default=0.0)
    sensor_depth = float(_sd) if _sd is not None else 0.0

    # v6: the depth reference the field operator selected in RiverSurveyor.
    # QRev: depthReference < 0.5 -> vertical beam, else bottom-track beams.
    _dr = get_field(setp, "depthReference", default=None)
    depth_reference = ("vb" if (_dr is None or float(_dr) < 0.5) else "bt")

    return dict(
        vb_depth=vb_depth,
        bt_beams=bt_beams,
        bt_freq=bt_freq,
        bt_depth=bt_depth,
        summary_depth=summary_depth,
        utm=np.asarray(utm, dtype=float) if utm is not None else None,
        lat=lat, lon=lon,
        heading=heading, pitch=pitch, roll=roll,
        time=time, step=step,
        mean_vel=mean_vel,
        edge_left=edge_left, edge_right=edge_right, start_edge=start_edge,
        sensor_depth=sensor_depth, depth_reference=depth_reference,
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

def compute_flow_direction(mean_vel: np.ndarray,
                           weights: np.ndarray | None = None,
                           subset: np.ndarray | None = None
                           ) -> tuple[float, float]:
    """
    Return (theta_flow, theta_section) in radians, math convention
    (CCW from +East, +X axis).  theta_section = theta_flow + 90 deg.

    The direction is the argument of the RESULTANT vector, i.e. a vector sum,
    so faster ensembles already weigh more than slow ones (a mean of unit
    vectors would not). `weights` adds two further corrections on top:

        * DENSITY — ensembles are not equally spaced along the crossing, so a
          plain mean over-represents wherever the boat slowed down.
        * DEPTH   — an ensemble over the thalweg carries far more flow than one
          over a shallow shelf at the same speed.

    With weights = interval * depth the result is the direction of the section's
    TOTAL UNIT-DISCHARGE VECTOR, which is exactly the quantity a cross-section
    should be perpendicular to. Pass weights=None for the plain v5 behaviour.

    `subset` is an optional boolean mask (v6 uses Step == 3, the in-transect
    ensembles) so stationary edge ensembles cannot drag the azimuth.
    """
    if mean_vel is None or mean_vel.size == 0:
        raise ValueError("Cannot compute flow direction: no Mean_Vel available.")
    ve = np.asarray(mean_vel[:, 0], dtype=float)
    vn = np.asarray(mean_vel[:, 1], dtype=float)
    mask = np.isfinite(ve) & np.isfinite(vn) & ((ve != 0) | (vn != 0))
    if subset is not None:
        sub = np.asarray(subset, dtype=bool)
        if sub.shape == mask.shape and np.any(mask & sub):
            mask = mask & sub          # only narrow if something survives
    if not np.any(mask):
        raise ValueError("Cannot compute flow direction: all Mean_Vel are zero/NaN.")

    if weights is None:
        w = np.ones(int(np.count_nonzero(mask)), dtype=float)
    else:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if not np.any(w > 0):
            w = np.ones(int(np.count_nonzero(mask)), dtype=float)

    theta_flow    = math.atan2(float(np.sum(w * vn[mask])),
                               float(np.sum(w * ve[mask])))
    theta_section = theta_flow + math.pi / 2.0
    return theta_flow, theta_section


# ============================================================================ #
#  v6 — REPEATED-TRANSECT (AFORO) GROUPING
# ============================================================================ #

# SonTek System.Time is seconds since 2000-01-01 (verified against the file
# name timestamp of a real transect: 829565394.36 s -> 2026-04-15 10:49:54).
SONTEK_EPOCH = datetime.datetime(2000, 1, 1)


def sontek_time(data: dict):
    """Absolute start time of a transect, or None if unavailable."""
    t = np.asarray(data.get("time"), dtype=float)
    if t.size == 0 or not np.any(np.isfinite(t)):
        return None
    t0 = float(np.nanmin(t))
    # A plain sample counter (extract_data's fallback) is not a timestamp.
    if t0 < 1e6:
        return None
    try:
        return SONTEK_EPOCH + datetime.timedelta(seconds=t0)
    except (OverflowError, ValueError):
        return None


def _ang_diff_180(a: float, b: float) -> float:
    """Absolute difference between two AXIS directions, modulo 180 deg."""
    d = math.degrees(a - b)
    return abs((d + 90.0) % 180.0 - 90.0)


def group_transects(results, args, log: list | None = None):
    """
    Cluster repeated transects of a discharge measurement (aforo) into groups.

    A discharge measurement is an OCCUPATION: the boat parks at the section and
    repeats the crossing back to back (SonTek recommends an even number of
    passes, at least 4). Those repetitions must be averaged into ONE profile
    rather than reported as four.

    Four conditions must hold simultaneously:

      1. same river and same anabranch;
      2. |Delta chainage| <= --group-tol  (default 15 m);
      3. section axes parallel to within --group-angle-tol (mod 180 deg);
      4. CONTIGUITY IN TIME — the group must be an unbroken run in the
         time-ordered sequence of transects.

    Condition 4 replaces the absolute time window of earlier drafts, which could
    not discriminate: consecutive sections 50 m apart are often less than 5 min
    apart, so any window loose enough to hold a real aforo also holds them.
    Contiguity needs no parameter and catches the case chainage alone cannot —
    passing the same spot again hours later at a different stage, which must NOT
    be merged because the water surface has moved.

    Returns a list of lists of result dicts (singletons included).
    """
    tol      = float(getattr(args, "group_tol", 15.0) or 15.0)
    ang_tol  = float(getattr(args, "group_angle_tol", 25.0) or 25.0)
    max_span = float(getattr(args, "group_max_span", 0.0) or 0.0)   # min, 0 = off

    if getattr(args, "no_group_repeats", False):
        return [[r] for r in results]

    # Manual override from the INI's [grupos] section takes precedence: any
    # transect named there is grouped exactly as written, and only what is left
    # over goes through the automatic criteria.
    manual = dict(getattr(args, "manual_groups", {}) or {})
    n_total = len(results)
    by_id = {r["perfil"]: r for r in results}
    forced, claimed = [], set()
    for gname, ids in manual.items():
        members = [by_id[i] for i in ids if i in by_id]
        missing = [i for i in ids if i not in by_id]
        if log is not None and missing:
            log.append(f"[warn] [grupos] {gname}: sin resultados para "
                       f"{', '.join(missing)}")
        if members:
            forced.append(members)
            claimed.update(m["perfil"] for m in members)
    if forced and log is not None:
        log.append(f"[info] [grupos]: {len(forced)} grupos forzados desde el INI")
    results = [r for r in results if r["perfil"] not in claimed]

    # Order by acquisition time; fall back to file name (a timestamp anyway).
    def _tkey(r):
        return (r.get("t_start") or datetime.datetime.max, str(r["perfil"]))
    ordered = sorted(results, key=_tkey)

    groups, current = [], []
    for r in ordered:
        if not current:
            current = [r]
            continue
        ref = current[0]
        same_river = (r.get("river", "") == ref.get("river", "")
                      and r.get("brazo", "") == ref.get("brazo", ""))
        ka, kb = r.get("group_chainage"), ref.get("group_chainage")
        if ka is not None and kb is not None:
            near = abs(ka - kb) <= tol
        else:
            # No centerline -> no chainage. Fall back to the ALONG-RIVER
            # component of the centroid offset, i.e. its projection on the
            # normal to the section axis. Using the straight distance would be
            # wrong: a repetition that covered a slightly different span of the
            # SAME section line has an along-section offset which must not
            # break the group, while a few metres upstream or downstream must.
            ca, cb = r.get("centroid"), ref.get("centroid")
            if ca is None or cb is None:
                near = False
            else:
                th = ref.get("theta_section", 0.0)
                nx, ny = -math.sin(th), math.cos(th)
                near = abs((ca[0] - cb[0]) * nx + (ca[1] - cb[1]) * ny) <= tol
        par = _ang_diff_180(r.get("theta_section", 0.0),
                            ref.get("theta_section", 0.0)) <= ang_tol
        span_ok = True
        if max_span > 0 and r.get("t_start") and current[0].get("t_start"):
            span_ok = ((r["t_start"] - current[0]["t_start"]).total_seconds()
                       <= max_span * 60.0)
        if same_river and near and par and span_ok:
            current.append(r)          # contiguous run continues
        else:
            groups.append(current)     # break: a different section intervened
            current = [r]
    if current:
        groups.append(current)
    groups = forced + groups

    if log is not None:
        multi = [g for g in groups if len(g) > 1]
        log.append(f"[info] agrupación de aforos: {len(groups)} grupos "
                   f"({len(multi)} con repeticiones) sobre {n_total} transectas "
                   f"[tol={tol:.1f} m, ang={ang_tol:.0f}°, contigüidad temporal]")
        for g in multi:
            ks = [x["group_chainage"] for x in g if x["group_chainage"] is not None]
            spread = (max(ks) - min(ks)) if ks else float("nan")
            log.append(f"        {group_id(g)}: {len(g)} repeticiones "
                       f"({', '.join(x['perfil'] for x in g)}), "
                       f"dispersión de progresiva = {spread:.1f} m")
    return groups


def group_id(group) -> str:
    """Stable identifier for a group: first transect id plus the repeat count."""
    base = sorted(g["perfil"] for g in group)[0]
    return base if len(group) == 1 else f"{base}_x{len(group)}"


def merge_group_clouds(group):
    """
    Merge the beam clouds of a group's repetitions into one cloud.

    Merging the CLOUDS and running build_profile once is preferred over
    averaging the finished profiles: every point is projected from its own real
    (E, N), so repetitions with different widths, different start edges and
    different boat tracks combine without any resampling or alignment choice.
    On a synthetic 4-repetition aforo this cut the RMSE against the known bed
    from 0.046 m (mean of the individual profiles) to 0.028 m.
    """
    keys = ("E", "N", "depth", "weight", "beam_index", "beam_freq_khz",
            "ens", "src", "primary")
    out = {}
    for k in keys:
        parts = [np.asarray(r["cloud"][k]) for r in group if r.get("cloud")]
        out[k] = np.concatenate(parts) if parts else np.array([])
    # Tag each point with the repetition it came from, for the QC plot.
    rep = []
    for i, r in enumerate(group):
        if r.get("cloud"):
            rep.append(np.full(len(r["cloud"]["E"]), i, dtype=int))
    out["rep"] = np.concatenate(rep) if rep else np.array([], dtype=int)
    return out


def group_flow_direction(group, args, log: list | None = None):
    """
    Section azimuth for a GROUP: the direction of the pooled unit-discharge
    vector of every in-transect ensemble of every repetition.

    Pooling (rather than averaging one angle per transect) is deliberate. Once
    the stationary edge ensembles are excluded, a repetition with more ensembles
    contributed more real crossing over the same flow field, and that is
    information worth keeping — averaging per transect would discard it. The
    pathological case (a repetition that wandered into a different flow field)
    is handled by DETECTION below, not by blanket down-weighting.
    """
    mvs, ws = [], []
    for r in group:
        d = r.get("data_ref")
        if d is None or d.get("mean_vel") is None:
            continue
        mv = np.asarray(d["mean_vel"], dtype=float)
        step = np.asarray(d.get("step", np.full(mv.shape[0], 3.0)), dtype=float)
        w = unit_discharge_weights(d, r["utm_ref"])
        keep = (step == 3) if np.any(step == 3) else np.ones(mv.shape[0], bool)
        mvs.append(mv[keep]); ws.append(w[keep])
    if not mvs:
        return group[0]["theta_flow"], group[0]["theta_section"]

    mv_all = np.vstack(mvs); w_all = np.concatenate(ws)
    theta_flow, theta_section = compute_flow_direction(mv_all, weights=w_all)

    if log is not None and len(group) > 1:
        devs = [_ang_diff_180(r["theta_section"], theta_section) for r in group]
        worst = max(devs)
        tol = float(getattr(args, "group_angle_tol", 25.0) or 25.0)
        tag = "PASS" if worst <= tol else "WARN"
        log.append(f"[QA] dispersión de θ entre repeticiones: máx "
                   f"{worst:.1f}° (tol {tol:.0f}°) {tag}")
        for r, d in zip(group, devs):
            if d > tol:
                log.append(f"       [WARN] {r['perfil']} se aparta {d:.1f}° "
                           f"del eje del grupo — ¿divagó a otro campo de flujo?")
    return theta_flow, theta_section


def group_repeat_stats(group, s_grid, depth_grid, centroid=None,
                       theta_section=None, s_shift=0.0):
    """
    Repeatability of the group: each repetition's own profile placed on the
    MERGED profile's grid, then the node-by-node spread. This is the number an
    aforo should report and that v5 could not produce.

    Alignment is GEOMETRIC: every repetition's profile nodes are projected onto
    the merged section axis using their real (E, N), exactly as the merged cloud
    was. A node that sits at a given point on the ground lands at the same s in
    both, so the comparison is like-for-like.

    An earlier version normalised by fraction of section width instead, which
    stretched each repetition to the merged width. With widths of 44.0, 43.5,
    42.5 and 45.5 m against a merged 46.0 m, that displaced features by up to
    1.15 m horizontally, and on a steep bank slope such a shift becomes a large
    apparent vertical difference — enough to make the merged profile plot
    OUTSIDE the envelope of the four repetitions it came from. The merged
    profile was right; the comparison was not.

    Returns (stack, sigma, rms_per_rep) with stack shape (n_reps, n_nodes).
    """
    rows = []
    have_axis = (centroid is not None and theta_section is not None)
    ux = math.cos(theta_section) if have_axis else 0.0
    uy = math.sin(theta_section) if have_axis else 0.0
    for r in group:
        sg = np.asarray(r["s_grid"], dtype=float)
        dg = np.asarray(r["depth_grid"], dtype=float)
        if sg.size < 2:
            rows.append(np.full_like(s_grid, np.nan))
            continue
        pe, pn = r.get("prof_E"), r.get("prof_N")
        if have_axis and pe is not None and pn is not None and len(pe) == len(dg):
            # s of this repetition's nodes IN THE MERGED FRAME
            s_rep = -((np.asarray(pe, float) - centroid[0]) * ux
                      + (np.asarray(pn, float) - centroid[1]) * uy) + s_shift
        else:
            # fallback: absolute distance from the left bank (never normalised)
            s_rep = sg
        o = np.argsort(s_rep)
        rows.append(np.interp(s_grid, s_rep[o], dg[o], left=np.nan, right=np.nan))
    stack = np.vstack(rows) if rows else np.empty((0, s_grid.size))
    # Nodes beyond the reach of every repetition (the merged bank ramps) are
    # all-NaN columns by construction; numpy warns about them through the
    # `warnings` module, which np.errstate does not silence.
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sigma = np.nanstd(stack, axis=0) if stack.size else np.full_like(s_grid, np.nan)
        rms = [float(np.sqrt(np.nanmean((row - depth_grid) ** 2)))
               if np.any(np.isfinite(row)) else float("nan") for row in stack]
    # Nodes where EVERY repetition has data. Outside that common zone the spread
    # is dominated by how far each pass got toward the bank, not by disagreement
    # about the bed: the merged cloud reaches further than any single pass, so
    # its bank ramp starts where some repetitions have nothing to compare with.
    common = (np.all(np.isfinite(stack), axis=0) if stack.size
              else np.zeros(s_grid.size, dtype=bool))
    return stack, sigma, rms, common


def transect_discharge(mat: dict) -> dict:
    """
    Discharge of one transect as RiverSurveyor computed it in the field.

    Summary.Total_Q is CUMULATIVE along the crossing, so the transect discharge
    is its last value; it equals the sum of the five component columns
    (Top / Middle / Bottom / Left / Right) at the same index, which is used here
    as a consistency check. RiverSurveyor's own edge estimates are kept — this
    is the measured discharge, not a recomputation.
    """
    summ = mat.get("Summary")
    out = dict(q_total=float("nan"), q_top=float("nan"), q_middle=float("nan"),
               q_bottom=float("nan"), q_left=float("nan"), q_right=float("nan"))
    if summ is None:
        return out
    tq = get_field(summ, "Total_Q")
    if tq is None:
        return out
    tq = np.asarray(tq, dtype=float)
    if tq.size == 0 or not np.any(np.isfinite(tq)):
        return out
    fin = tq[np.isfinite(tq)]
    out["q_total"] = float(fin[-1])
    for key, name in (("q_top", "Top_Q"), ("q_middle", "Middle_Q"),
                      ("q_bottom", "Bottom_Q"), ("q_left", "Left_Q"),
                      ("q_right", "Right_Q")):
        v = get_field(summ, name)
        if v is None:
            continue
        v = np.asarray(v, dtype=float)
        if v.size and np.any(np.isfinite(v)):
            out[key] = float(v[np.isfinite(v)][-1])
    return out


def group_discharge_stats(group) -> dict:
    """
    Discharge statistics of an aforo: the mean of the repetitions with its
    spread. This is the deliverable of a discharge measurement.

    The coefficient of variation is the figure to watch. USGS practice is that
    individual transects should fall within about 5 % of the mean; beyond that
    the measurement usually needs more repetitions or has a problem (moving bed,
    unstable track reference, transects too short).

    Note this is INDEPENDENT of the bathymetric repeatability: a stable bed with
    a noisy discharge is perfectly possible, and the two numbers should be read
    separately.
    """
    qs = [r.get("q_total") for r in group]
    qs = [q for q in qs if q is not None and np.isfinite(q)]
    out = dict(q_mean=float("nan"), q_std=float("nan"), q_cov_pct=float("nan"),
               q_min=float("nan"), q_max=float("nan"), q_values=[])
    if not qs:
        return out
    a = np.asarray(qs, dtype=float)
    out["q_values"] = [float(x) for x in a]
    out["q_mean"] = float(np.mean(a))
    out["q_min"] = float(np.min(a))
    out["q_max"] = float(np.max(a))
    if a.size > 1:
        out["q_std"] = float(np.std(a, ddof=1))
        if out["q_mean"] != 0:
            out["q_cov_pct"] = 100.0 * out["q_std"] / abs(out["q_mean"])
    return out


def resolve_depth_reference(data: dict, args, log: list | None = None):
    """
    Resolve (depth_ref, bt_geometry, composite_on) from the CLI and the file.

    `--depth-ref auto` honours what the field operator selected in
    RiverSurveyor (Setup.depthReference), which QRev also reads. The deprecated
    v5 flag --blend-beams still maps to composite processing over the whole
    cloud so old command lines keep working.
    """
    ref = str(getattr(args, "depth_ref", "auto") or "auto").lower()
    if ref == "auto":
        ref = str(data.get("depth_reference", "vb") or "vb")
        src = "Setup.depthReference"
    else:
        src = "CLI/config"
    if ref not in ("vb", "bt"):
        ref = "vb"
    geom = str(getattr(args, "bt_geometry", "footprints") or "footprints").lower()
    comp = str(getattr(args, "composite", "on") or "on").lower() != "off"
    if getattr(args, "blend_beams", False):
        comp = True
        geom = "footprints"
    if log is not None:
        log.append(f"[info] depth reference: {ref.upper()} (from {src}), "
                   f"bt_geometry={geom}, bt_avg={getattr(args,'bt_avg','idw')}, "
                   f"composite={'on' if comp else 'off'}")
    return ref, geom, comp


def unit_discharge_weights(data: dict, utm: np.ndarray) -> np.ndarray:
    """
    Per-ensemble weight = (along-track interval) x (depth), i.e. the unit
    discharge each ensemble represents once multiplied by its velocity.
    """
    w = ensemble_interval_weights(utm)
    d = data.get("summary_depth")
    if d is None:
        d = data.get("vb_depth")
    d = np.asarray(d, dtype=float)
    if d.shape[0] != w.shape[0]:
        return w
    d = np.where(np.isfinite(d) & (d > 0), d, np.nan)
    med = float(np.nanmedian(d)) if np.any(np.isfinite(d)) else 1.0
    d = np.where(np.isfinite(d), d, med)
    return w * d


def flow_direction_for(data: dict, utm: np.ndarray, args, log: list | None = None):
    """
    Section azimuth for one transect, with the v6 corrections and a log line
    quantifying what they changed relative to the plain v5 mean.
    """
    mv = data["mean_vel"]
    mode = str(getattr(args, "flow_weighting", "discharge") or "discharge").lower()
    step = np.asarray(data.get("step", np.full(mv.shape[0], 3.0)), dtype=float)
    in_transect = step == 3
    use_step = (getattr(args, "use_edge_ensembles", False) is False
                and int(np.count_nonzero(in_transect)) >= 5)

    th_plain, _ = compute_flow_direction(mv)          # v5 reference value
    if mode == "none":
        w = None
    elif mode == "density":
        w = ensemble_interval_weights(utm)
    else:
        w = unit_discharge_weights(data, utm)

    theta_flow, theta_section = compute_flow_direction(
        mv, weights=w, subset=(in_transect if use_step else None))

    if log is not None:
        n_edge = int(np.count_nonzero(~in_transect))
        log.append(f"[info] theta_flow    = {math.degrees(theta_flow):+.1f}°  "
                   f"(math, CCW from East; weighting={mode}"
                   + (f", {n_edge} edge ensembles excluded" if use_step and n_edge else "")
                   + ")")
        log.append(f"[info] theta_section = {math.degrees(theta_section):+.1f}°  "
                   f"(perpendicular)   [v5 plain mean would give "
                   f"{math.degrees(th_plain) + 90.0:+.1f}°, "
                   f"Δ={math.degrees(theta_flow - th_plain):+.2f}°]")
    return theta_flow, theta_section


def section_orientation_qc(theta_section: float, utm: np.ndarray,
                           network=None, utm_crs=None,
                           centroid=None, log: list | None = None,
                           river=None) -> dict:
    """
    Diagnose HOW TRANSVERSE the section really is.

    Reports two angles that answer different questions:

      * TRACK vs AXIS — how obliquely the boat crossed relative to the section
        plane. The projected width is the traversed length times cos(angle), so
        this is a direct measure of how much section width is being lost.
      * CENTERLINE vs AXIS — the angle between the flow-perpendicular section
        and the perpendicular to the river centerline at this chainage. This is
        the diagnostic for "my sections don't look transverse to the river":
        if the deviation is systematic and single-signed across a campaign,
        suspect the flow computation or the digitised centerline; if it
        scatters and tracks bends, islands and bifurcations, the flow really is
        oblique there and the sections are right.
    """
    out = dict(track_angle_deg=float("nan"), width_loss_pct=float("nan"),
               centerline_angle_deg=float("nan"))
    P = np.asarray(utm, dtype=float)
    ok = np.isfinite(P[:, 0]) & np.isfinite(P[:, 1])
    if int(np.count_nonzero(ok)) >= 2:
        Q = P[ok]
        th_track = math.atan2(Q[-1, 1] - Q[0, 1], Q[-1, 0] - Q[0, 0])
        d = math.degrees(theta_section - th_track)
        d = (d + 90.0) % 180.0 - 90.0
        out["track_angle_deg"] = d
        out["width_loss_pct"] = 100.0 * (1.0 - abs(math.cos(math.radians(d))))

    if network is not None and utm_crs is not None and centroid is not None:
        try:
            th_cl = network.tangent_azimuth(centroid, utm_crs, river=river)
            if th_cl is not None and np.isfinite(th_cl):
                # centerline-perpendicular vs the flow-perpendicular section
                d2 = math.degrees(theta_section - (th_cl + math.pi / 2.0))
                out["centerline_angle_deg"] = (d2 + 90.0) % 180.0 - 90.0
        except Exception:
            pass

    if log is not None:
        if np.isfinite(out["track_angle_deg"]):
            log.append(f"[QA] orientación: recorrido vs eje = "
                       f"{out['track_angle_deg']:+.1f}°  "
                       f"(ancho proyectado pierde {out['width_loss_pct']:.1f}%)")
        if np.isfinite(out["centerline_angle_deg"]):
            tag = "PASS" if abs(out["centerline_angle_deg"]) <= 25.0 else "WARN"
            log.append(f"[QA] orientación: sección ⟂flujo vs ⟂traza = "
                       f"{out['centerline_angle_deg']:+.1f}° {tag}")
    return out


def ensemble_interval_weights(utm: np.ndarray) -> np.ndarray:
    """
    Per-ensemble weight equal to the along-track interval each ensemble
    represents (half the distance to each neighbour). Makes any ensemble
    statistic invariant to boat speed.

    Motivation: the stationary EDGE ensembles pack 3-9x more samples per metre
    than the crossing (measured: 6.2 and 16.0 ens/m at the banks vs 1.8 ens/m
    in-transect), and a field delay at a bank — a stuck motor, a late start —
    can pile up hundreds of near-identical samples at one spot. Un-weighted,
    that pile dominates the weighted median of its bin AND anchors the bank
    ramp. Weighted by the interval it represents, it cannot.

    Returns weights normalised to mean 1.0 so it composes with the beam weights
    and the perpendicular-offset penalty without changing their scale.
    """
    P = np.asarray(utm, dtype=float)
    n = P.shape[0]
    if n < 2:
        return np.ones(max(n, 1), dtype=float)
    seg = np.hypot(np.diff(P[:, 0]), np.diff(P[:, 1]))
    seg = np.where(np.isfinite(seg), seg, 0.0)
    w = np.empty(n, dtype=float)
    w[0]    = seg[0] / 2.0
    w[-1]   = seg[-1] / 2.0
    w[1:-1] = (seg[:-1] + seg[1:]) / 2.0
    good = np.isfinite(w) & (w > 0)
    if not np.any(good):
        return np.ones(n, dtype=float)
    # Clamp wild intervals (GPS jumps) to the 95th percentile before normalising
    cap = float(np.percentile(w[good], 95)) * 3.0
    w = np.clip(w, 0.0, cap if cap > 0 else None)
    m = float(np.mean(w[w > 0])) if np.any(w > 0) else 1.0
    w = np.where(w > 0, w / m, 0.0)
    # A zero-interval ensemble still measured the bed; give it a small floor so
    # it contributes without being able to dominate.
    return np.where(w > 0, w, 0.05)


def average_bt_depth(bt_beams: np.ndarray, draft: float,
                     method: str = "idw") -> np.ndarray:
    """
    Collapse the 4 bottom-track beam depths of each ensemble to one depth.

    Ported from QRev (Classes/DepthData.average_depth). QRev's default is IDW,
    NOT a simple mean: the beams are weighted by the inverse of their range from
    the transducer, so the shorter (shallower) returns weigh more. On a sloping
    bed the long beam looks down the slope and biases the mean deep; IDW pulls
    it back.

        rng = depth - draft
        w   = 1 - rng / sum(rng)
        avg = draft + sum(rng*w) / sum(w)

    NOTE the depths already include the draft (see extract_data), so it is
    subtracted here to get the true range and added back afterwards.
    """
    b = np.asarray(bt_beams, dtype=float)
    if b.ndim == 1:
        b = b.reshape(-1, 1)
    b = np.where(np.isfinite(b) & (b > 0), b, np.nan)
    if method == "simple":
        with np.errstate(invalid="ignore"):
            return np.nanmean(b, axis=1)
    rng = b - float(draft)
    with np.errstate(invalid="ignore", divide="ignore"):
        tot = np.nansum(rng, axis=1, keepdims=True)
        w = 1.0 - np.divide(rng, np.where(tot == 0, np.nan, tot))
        wsum = np.nansum(w, axis=1)
        avg = float(draft) + np.divide(np.nansum(rng * w, axis=1),
                                       np.where(wsum == 0, np.nan, wsum))
    avg = np.where(np.isfinite(avg), avg, np.nan)
    avg[avg == float(draft)] = np.nan
    return avg


def bt_beams_valid(bt_beams: np.ndarray, min_beams: int = 2) -> np.ndarray:
    """
    QRev's validity rule for a 4-beam average: with 3 valid beams the average is
    accepted; with fewer than `min_beams` it is rejected and an alternative
    reference is used instead (see DepthStructure.composite_depths docstring).
    """
    b = np.asarray(bt_beams, dtype=float)
    if b.ndim == 1:
        b = b.reshape(-1, 1)
    return np.sum(np.isfinite(b) & (b > 0), axis=1) >= int(min_beams)


def _running_iqr(data: np.ndarray, half_width: int) -> np.ndarray:
    """Running inter-quartile range, excluding the target point (QRev run_iqr)."""
    x = np.asarray(data, dtype=float)
    n = x.size
    out = np.full(n, np.nan)
    for i in range(n):
        lo = max(0, i - half_width)
        hi = min(n, i + half_width + 1)
        w = np.concatenate([x[lo:i], x[i + 1:hi]])
        w = w[np.isfinite(w)]
        if w.size >= 4:
            out[i] = float(np.percentile(w, 75) - np.percentile(w, 25))
    return out


def filter_depth_spikes(depth: np.ndarray, half_width: int = 10,
                        multiplier: float = 15.0,
                        cycles: int = 3) -> np.ndarray:
    """
    Flag unnatural spikes in a beam's depth series. Returns a boolean mask of
    VALID samples.

    Ported in spirit from QRev's DepthData.filter_smooth: residuals against a
    smoothed version of the series are tested against a running IQR, and the
    rejection threshold is the MAXIMUM of the IQR criterion, 5 % of the measured
    depth, and 0.10 m — so it never rejects on noise smaller than the instrument
    can resolve. QRev uses a robust Loess; a median filter is used here to avoid
    pulling in statsmodels for one function, which is equivalent for spike
    detection and much faster.

    v5 had NO depth outlier rejection at all: a single bad bottom detection went
    straight into the cloud and shifted its bin.
    """
    d = np.asarray(depth, dtype=float).copy()
    n = d.size
    valid = np.isfinite(d) & (d > 0)
    if n < 5 or int(np.count_nonzero(valid)) < 5:
        return valid
    work = d.copy()
    for _ in range(max(1, int(cycles))):
        filled = work.copy()
        idx = np.arange(n)
        ok = np.isfinite(filled)
        if int(np.count_nonzero(ok)) < 3:
            break
        filled = np.interp(idx, idx[ok], filled[ok])
        k = max(3, 2 * (half_width // 2) + 1)
        pad = k // 2
        padded = np.pad(filled, pad, mode="edge")
        smooth = np.array([np.median(padded[i:i + k]) for i in range(n)])
        resid = filled - smooth
        iqr = _running_iqr(resid, half_width)
        crit = np.fmax(np.fmax(iqr * multiplier, 0.05 * np.abs(filled)), 0.10)
        bad = np.isfinite(crit) & (np.abs(resid) > crit)
        if not np.any(bad & valid):
            break
        valid = valid & ~bad
        work = np.where(valid, d, np.nan)
    return valid


def filter_cloud_depths(data: dict, enabled: bool = True,
                        log: list | None = None) -> dict:
    """
    Apply filter_depth_spikes independently to the vertical beam and to each of
    the 4 bottom-track beams, and return a copy of `data` with the rejected
    samples set to NaN. Filtering each beam separately is what QRev does — a bad
    return on one beam should not discard the other three.
    """
    out = dict(data)
    if not enabled:
        return out
    vb = np.asarray(data["vb_depth"], dtype=float).copy()
    vb_ok = filter_depth_spikes(vb)
    n_vb = int(np.count_nonzero(np.isfinite(vb) & (vb > 0) & ~vb_ok))
    vb[~vb_ok] = np.nan
    out["vb_depth"] = vb

    bt = np.asarray(data["bt_beams"], dtype=float).copy()
    n_bt = 0
    if bt.ndim == 2:
        for k in range(bt.shape[1]):
            col = bt[:, k]
            ok = filter_depth_spikes(col)
            n_bt += int(np.count_nonzero(np.isfinite(col) & (col > 0) & ~ok))
            col[~ok] = np.nan
            bt[:, k] = col
    out["bt_beams"] = bt
    if log is not None:
        log.append(f"[info] depth spike filter: rejected {n_vb} VB and {n_bt} "
                   f"slant-beam samples")
    return out


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
# v6 source codes, matching QRev's DepthStructure.composite_depths comp_source
SRC_BT = 1        # bottom-track (slant) beams
SRC_VB = 2        # vertical beam (nadir)
SRC_DS = 3        # depth sounder (not present in these exports)
SRC_INT = 4       # interpolated / gap-filled
BEAM_BT_ENS = -1  # beam_index of a collapsed ensemble-averaged BT point
BEAM_DS     = -2  # beam_index of a depth-sounder point

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


def build_beam_cloud(data: dict, utm: np.ndarray,
                     depth_ref: str = "vb",
                     bt_geometry: str = "footprints",
                     bt_avg: str = "idw",
                     density_weighting: bool = True) -> dict:
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
    draft    = float(data.get("sensor_depth", 0.0) or 0.0)

    # --- v6: depth-reference model (QRev vocabulary) ------------------------
    # depth_ref     : 'vb' | 'bt'  -> which source is PRIMARY
    # bt_geometry   : 'footprints' -> the 4 slant beams keep their real
    #                                 footprints (better near-bank coverage)
    #                 'ensemble'   -> collapsed to ONE averaged point at the
    #                                 boat, which is what QRev does
    # bt_avg        : 'idw' | 'simple' averaging for the collapsed point
    # Every point carries `src` (1=BT, 2=VB, matching QRev's comp_source codes)
    # and `primary` (True if it belongs to the selected primary reference), so
    # build_profile can apply the composite rule to any reference, and the
    # exports can record where each bed point actually came from.
    depth_ref   = str(depth_ref or "vb").lower()
    bt_geometry = str(bt_geometry or "footprints").lower()
    bt_avg      = str(bt_avg or "idw").lower()

    bt_ens_depth = None
    bt_ens_ok = None
    if bt_geometry == "ensemble" or depth_ref == "bt":
        bt_ens_depth = average_bt_depth(bt_beams, draft, bt_avg)
        bt_ens_ok = bt_beams_valid(bt_beams, min_beams=2)

    if density_weighting:
        w_den = ensemble_interval_weights(utm)
    else:
        w_den = np.ones(n, dtype=float)

    pts_e, pts_n, pts_z, pts_w, pts_b, pts_f = [], [], [], [], [], []
    pts_ens, pts_src, pts_pri = [], [], []

    tilt = math.radians(TILT_DEG)

    for i in range(n):                         # loop over SAMPLES
        E0, N0 = utm[i]
        if not (np.isfinite(E0) and np.isfinite(N0)):
            continue
        hdg = heading[i] if np.isfinite(heading[i]) else 0.0
        hdg_rad = math.radians(hdg)            # ADP +CW from true north
        wd = float(w_den[i]) if np.isfinite(w_den[i]) else 1.0

        # --- vertical beam ---
        if np.isfinite(vb_depth[i]) and vb_depth[i] > 0:
            pts_e.append(E0); pts_n.append(N0); pts_z.append(float(vb_depth[i]))
            pts_w.append(W_VB * wd); pts_b.append(0); pts_f.append(FREQ_VB)
            pts_ens.append(i); pts_src.append(SRC_VB)
            pts_pri.append(depth_ref == "vb")

        f = float(bt_freq[i]) if np.isfinite(bt_freq[i]) else FREQ_HIGH

        if bt_geometry == "ensemble":
            # --- ONE averaged bottom-track point at the boat (QRev style) ---
            d = bt_ens_depth[i] if bt_ens_depth is not None else np.nan
            if np.isfinite(d) and d > 0 and bool(bt_ens_ok[i]):
                pts_e.append(E0); pts_n.append(N0); pts_z.append(float(d))
                pts_w.append(W_VB * wd)        # one point/ensemble, same as VB
                pts_b.append(BEAM_BT_ENS); pts_f.append(f)
                pts_ens.append(i); pts_src.append(SRC_BT)
                pts_pri.append(depth_ref == "bt")
        else:
            # --- 4 slant beams at their real footprints ---
            az_set = AZ_3000 if f >= FREQ_SPLIT else AZ_1000
            for k in range(4):
                d = bt_beams[i, k]
                if not np.isfinite(d) or d <= 0:
                    continue
                az_inst = math.radians(az_set[k])     # instrument frame, +CW from forward
                # World azimuth (true north, +CW) = heading + instrument azimuth
                az_world = hdg_rad + az_inst
                r = d * math.tan(tilt)                # horizontal radius of footprint
                # Compass azimuth (CW from N): dE = r*sin(az), dN = r*cos(az)
                dE = r * math.sin(az_world)
                dN = r * math.cos(az_world)
                pts_e.append(E0 + dE)
                pts_n.append(N0 + dN)
                pts_z.append(float(d))
                pts_w.append(beam_weight(az_set[k], f) * wd)
                pts_b.append(k + 1)
                pts_f.append(f)
                pts_ens.append(i); pts_src.append(SRC_BT)
                pts_pri.append(depth_ref == "bt")

    return dict(
        E=np.asarray(pts_e),
        N=np.asarray(pts_n),
        depth=np.asarray(pts_z),
        weight=np.asarray(pts_w),
        beam_index=np.asarray(pts_b, dtype=int),
        beam_freq_khz=np.asarray(pts_f),
        ens=np.asarray(pts_ens, dtype=int),
        src=np.asarray(pts_src, dtype=int),
        primary=np.asarray(pts_pri, dtype=bool),
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
    """Argentine chainage notation: metres -> 'k+mmm.mm' (e.g. 12340.5 -> '12+340.50').
    6.1: rounded to the centimetre BEFORE splitting — 1999.999 printed '1+1000.00'."""
    p = round(abs(float(p_m)), 2)
    sign = "-" if p_m < 0 and p > 0 else ""
    km = int(p // 1000)
    m = round(p - km * 1000.0, 2)
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


ISLAND_S_TOL = 25.0   # [m] tolerance grouping anabranches that share a split/rejoin


def _mkline(coords):
    return LineString(list(coords))


class RiverAxis:
    """One river's channel network: a MAIN line (carrying the river's internal
    chainage, measured from its upstream end = 0) plus zero or more ANABRANCHES
    (island channels) that leave and rejoin the main.

    The chainage of ANY point is its projection onto the main, so the
    longitudinal axis stays single-valued across islands; being on an anabranch
    only changes the point's `brazo` label, not its chainage (agreed convention:
    branch chainage = projection onto the main axis)."""

    def __init__(self, river: str, main, log=None):
        self.river = river
        self.main = main                       # LineString oriented US -> DS
        self.length = float(main.length)
        self.anabranches = []                  # dicts: geom, s0, s1, label
        self.islands = []                      # dicts: s0, s1, main_label
        self._log = log if log is not None else []

    @property
    def us_xy(self):
        c = self.main.coords[0]
        return (float(c[0]), float(c[1]))

    @property
    def ds_xy(self):
        c = self.main.coords[-1]
        return (float(c[0]), float(c[1]))

    def add_anabranch(self, g, label=None):
        s0 = float(self.main.project(Point(g.coords[0])))
        s1 = float(self.main.project(Point(g.coords[-1])))
        if s1 < s0:                            # digitised upstream -> flip
            g = _mkline(list(g.coords)[::-1])
            s0, s1 = s1, s0
        self.anabranches.append(dict(geom=g, s0=s0, s1=s1, label=label))

    def _flow_left_normal(self, s_mid):
        """Unit left-normal (facing downstream) of the main near chainage s_mid."""
        a = max(0.0, s_mid - 1.0)
        b = min(self.length, s_mid + 1.0)
        pa = self.main.interpolate(a)
        pb = self.main.interpolate(b)
        fx, fy = pb.x - pa.x, pb.y - pa.y
        n = math.hypot(fx, fy)
        if n < 1e-9:
            return (0.0, 0.0)
        return (-fy / n, fx / n)               # +90 deg CCW = left bank facing DS

    def finalize_islands(self):
        """Group anabranches that share a split/rejoin into islands, and label
        every channel of each island left-to-right facing downstream
        (MI, M2.., MD). The main sub-reach inside an island is one channel too,
        so on a 2-channel island the main gets the side opposite the anabranch."""
        if not self.anabranches:
            return
        order = sorted(range(len(self.anabranches)),
                       key=lambda i: (self.anabranches[i]["s0"],
                                      self.anabranches[i]["s1"]))
        groups = []
        for i in order:
            ab = self.anabranches[i]
            placed = False
            for grp in groups:
                r0, r1 = grp["rep"]
                if abs(ab["s0"] - r0) <= ISLAND_S_TOL and abs(ab["s1"] - r1) <= ISLAND_S_TOL:
                    grp["members"].append(i)
                    placed = True
                    break
            if not placed:
                groups.append(dict(rep=(ab["s0"], ab["s1"]), members=[i]))

        for grp in groups:
            s0 = float(np.mean([self.anabranches[i]["s0"] for i in grp["members"]]))
            s1 = float(np.mean([self.anabranches[i]["s1"] for i in grp["members"]]))
            s_mid = 0.5 * (s0 + s1)
            lnx, lny = self._flow_left_normal(s_mid)
            center = self.main.interpolate(s_mid)

            channels = [dict(kind="main", idx=None)]
            for i in grp["members"]:
                channels.append(dict(kind="anab", idx=i))

            def left_score(ch):
                if ch["kind"] == "main":
                    return 0.0
                g = self.anabranches[ch["idx"]]["geom"]
                m = g.interpolate(0.5 * g.length)
                return (m.x - center.x) * lnx + (m.y - center.y) * lny

            channels.sort(key=left_score, reverse=True)   # most-left first
            labels = self._island_labels(len(channels))
            for lab, ch in zip(labels, channels):
                if ch["kind"] == "main":
                    self.islands.append(dict(s0=s0, s1=s1, main_label=lab))
                else:
                    ab = self.anabranches[ch["idx"]]
                    if ab["label"] is None or ab["label"] == "":
                        ab["label"] = lab

    @staticmethod
    def _island_labels(n):
        if n <= 1:
            return [""]
        if n == 2:
            return ["MI", "MD"]
        return ["MI"] + [f"M{k}" for k in range(2, n)] + ["MD"]

    def locate_local(self, x, y):
        """(km_internal, brazo, role, dist_to_nearest_channel) for a point."""
        p = Point(x, y)
        best_kind, best_d, best_i = "main", self.main.distance(p), None
        for i, ab in enumerate(self.anabranches):
            d = ab["geom"].distance(p)
            if d < best_d:
                best_kind, best_d, best_i = "anab", d, i
        km_internal = float(self.main.project(p))     # always project onto main
        if best_kind == "main":
            brazo = ""
            for isl in self.islands:
                if isl["s0"] - ISLAND_S_TOL <= km_internal <= isl["s1"] + ISLAND_S_TOL:
                    brazo = isl["main_label"]
                    break
            role = ROLE_MAIN
        else:
            brazo = self.anabranches[best_i]["label"] or ""
            role = ROLE_ANAB
        return km_internal, brazo, role, float(best_d)


class RiverNetwork:
    """Multi-river centerline. Reads a vector file of channel reaches tagged by
    `river` and `role`, builds one RiverAxis per river, keeps per-river OFFICIAL
    chainage offsets, and stitches a continuous INTERNAL chainage across the
    surveyed rivers for the longitudinal profile.

    Backward compatible: if the file has no `river`/`role` columns it is treated
    as ONE river (DEFAULT_RIVER_NAME) whose longest merged line is the main and
    whose island parts (both endpoints on the main) are anabranches — i.e. exactly
    the v4 single-river behaviour.

    Chainage is computed in `crs` (the centerline's own projected CRS unless
    overridden); callers pass points already reprojected to `crs`.
    """

    def __init__(self, centerline_path, crs_override=None, reverse: bool = False,
                 river_offsets=None, ws_xy_elev=None,
                 default_river: str = DEFAULT_RIVER_NAME, log=None):
        if not HAS_GPD:
            raise RuntimeError("RiverNetwork needs geopandas/shapely")
        self._log = log if log is not None else []
        self._reverse = bool(reverse)
        self.river_offsets = {str(k): float(v) for k, v in (river_offsets or {}).items()}

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

        cols = {c.lower(): c for c in gdf.columns}
        col_river = self._resolve(cols, CL_ATTR_RIVER)
        col_role = self._resolve(cols, CL_ATTR_ROLE)
        col_label = self._resolve(cols, CL_ATTR_LABEL)
        col_flow = self._resolve(cols, CL_ATTR_FLOW)
        self._has_schema = col_river is not None
        self._log.append(
            "[info] centerline attributes: "
            f"river={col_river or '-'}, role={col_role or '-'}, "
            f"label={col_label or '-'}, flow_dir={col_flow or '-'}"
            + ("" if self._has_schema else "  -> single-river (v4) mode"))

        feats = []
        for _, row in gdf.iterrows():
            g = row.geometry
            if g is None or g.is_empty:
                continue
            river = str(row[col_river]).strip() if col_river else default_river
            if not river or river.lower() in ("nan", "none"):
                river = default_river
            role = (str(row[col_role]).strip().lower()
                    if col_role and row[col_role] is not None else "")
            label = (str(row[col_label]).strip()
                     if col_label and row[col_label] is not None else "")
            if label.lower() in ("nan", "none"):
                label = ""
            flow = (str(row[col_flow]).strip().upper()
                    if col_flow and row[col_flow] is not None else "")
            geoms = ([g] if g.geom_type == "LineString"
                     else (list(g.geoms) if g.geom_type == "MultiLineString" else []))
            for gg in geoms:
                feats.append(dict(river=river, role=role, label=label,
                                  flow=flow, geom=gg))
        if not feats:
            raise ValueError("centerline file has no usable line geometry")

        rivers = {}
        for f in feats:
            rivers.setdefault(f["river"], []).append(f)

        self.axes = []
        for river, fl in rivers.items():
            axis = self._build_axis(river, fl, ws_xy_elev, col_role is not None)
            if axis is not None:
                self.axes.append(axis)
        if not self.axes:
            raise ValueError("centerline produced no usable river axis")

        # 6.1: [progresivas] keys matched accent/case-insensitively ('Neuquén' =
        # 'Neuquen'); a key matching no river used to be ignored silently
        by_norm = {_norm_name(ax.river): ax.river for ax in self.axes}
        offs = {}
        for k, v in self.river_offsets.items():
            r = k if k in by_norm.values() else by_norm.get(_norm_name(k))
            if r is None:
                self._log.append(f"[warn] [progresivas] '{k}' matches no river of the "
                                 f"centerline {sorted(by_norm.values())} — ignored")
                continue
            offs[r] = v
        self.river_offsets = offs
        for ax in self.axes:
            self._log.append(
                f"[info] río '{ax.river}': main {ax.length:.0f} m, "
                f"{len(ax.anabranches)} anabranch(es), {len(ax.islands)} island(s); "
                f"official offset {self.river_offsets.get(ax.river, 0.0):.1f} m")

        self.survey_offsets = None             # legacy v5/v6.0 single chain
        self.survey_ok = False
        self._by_river = {ax.river: ax for ax in self.axes}
        self._build_topology()                 # 6.1: links + head->outlet paths

    # ------------------------------------------------------------------ #
    @staticmethod
    def _resolve(cols_lower, names):
        for n in names:
            if n in cols_lower:
                return cols_lower[n]
        return None

    def _build_axis(self, river, feats, ws_xy_elev, has_role):
        if has_role:
            main_geoms = [f["geom"] for f in feats if f["role"] == ROLE_MAIN]
            anab_feats = [f for f in feats if f["role"] == ROLE_ANAB]
            other = [f for f in feats if f["role"] not in (ROLE_MAIN, ROLE_ANAB)]
            if other:
                main_geoms += [f["geom"] for f in other]
                self._log.append(
                    f"[warn] río '{river}': {len(other)} feature(s) with unknown role "
                    "treated as main")
            if not main_geoms:
                self._log.append(f"[warn] río '{river}': no role=main feature — skipped")
                return None
        else:
            main_geoms = [f["geom"] for f in feats]
            anab_feats = []

        merged = linemerge(main_geoms) if len(main_geoms) > 1 else main_geoms[0]
        parts = ([merged] if merged.geom_type == "LineString" else list(merged.geoms))
        parts.sort(key=lambda g: g.length, reverse=True)
        main_raw = parts[0]

        flow_attr = ""
        for f in feats:
            if f.get("flow"):
                flow_attr = f["flow"]
                break
        main = self._orient(main_raw, ws_xy_elev, flow_attr, river)
        axis = RiverAxis(river, main, log=self._log)

        extra = parts[1:]
        if has_role:
            for f in anab_feats:
                axis.add_anabranch(f["geom"], label=(f["label"] or None))
            for g in extra:
                self._log.append(
                    f"[warn] río '{river}': extra role=main part of {g.length:.0f} m "
                    "not merged into the main — ignored for chainage")
        else:
            for g in extra:
                d0 = main.distance(Point(g.coords[0]))
                d1 = main.distance(Point(g.coords[-1]))
                if d0 <= BRANCH_SNAP_TOL and d1 <= BRANCH_SNAP_TOL:
                    axis.add_anabranch(g)
                else:
                    self._log.append(
                        f"[warn] río '{river}': part of {g.length:.0f} m not connected "
                        f"to the main (gaps {d0:.0f}/{d1:.0f} m) — ignored for chainage")
        axis.finalize_islands()
        return axis

    def _orient(self, line, ws_xy_elev, flow_attr, river):
        """Return `line` oriented so chainage increases DOWNSTREAM (from km0 at
        the upstream end). Priority: explicit flow_dir attribute > water-surface
        slope of the points on THIS river > --chainage-reverse > digitised order."""
        if flow_attr == "DS2US":
            self._log.append(f"[info] río '{river}': oriented by flow_dir=DS2US (reversed)")
            return _mkline(list(line.coords)[::-1])
        if flow_attr == "US2DS":
            self._log.append(f"[info] río '{river}': oriented by flow_dir=US2DS")
            return line
        if ws_xy_elev:
            xy = np.asarray([(p[0], p[1]) for p in ws_xy_elev], float)
            zs = np.asarray([p[2] for p in ws_xy_elev], float)
            s = np.asarray([line.project(Point(x, y)) for x, y in xy])
            d = np.asarray([line.distance(Point(x, y)) for x, y in xy])
            on = np.isfinite(s) & np.isfinite(zs) & (d <= max(BRANCH_SNAP_TOL, 50.0))
            if int(on.sum()) >= 3 and float(np.ptp(s[on])) > 1.0:
                slope = float(np.polyfit(s[on], zs[on], 1)[0])
                if slope > 0:
                    self._log.append(
                        f"[info] río '{river}': oriented from water surface — raw dir "
                        f"was UPSTREAM, reversed (WS slope +{slope*1000:.2f} m/km)")
                    return _mkline(list(line.coords)[::-1])
                self._log.append(
                    f"[info] río '{river}': orientation confirmed from water surface "
                    f"(WS slope {slope*1000:.2f} m/km, decreasing downstream)")
                return line
        if self._reverse:
            self._log.append(f"[info] río '{river}': reversed by --chainage-reverse")
            return _mkline(list(line.coords)[::-1])
        return line

    # ------------------------------------------------------------------ #
    def locate(self, x, y):
        """Locate a point. Returns dict(river, role, brazo, km_internal,
        km_oficial, dist). km_oficial = river_offset + km_internal."""
        best = None
        for ax in self.axes:
            km_i, brazo, role, dist = ax.locate_local(x, y)
            if best is None or dist < best["dist"]:
                best = dict(river=ax.river, role=role, brazo=brazo,
                            km_internal=km_i, dist=dist)
        best["km_oficial"] = self.river_offsets.get(best["river"], 0.0) + best["km_internal"]
        return best

    def locate_on(self, river, x, y):
        """Like locate(), but on the axis of a GIVEN river (6.1). Used for gauges,
        whose river is declared: near a confluence the nearest axis can belong to
        another river, and its chainage would then be meaningless. Returns None
        if the river is not in the network. `at_end` flags a projection clamped
        onto an end of the main (point beyond the digitised reach)."""
        ax = self._by_river.get(river)
        if ax is None:
            return None
        km_i, brazo, role, dist = ax.locate_local(x, y)
        return dict(river=river, role=role, brazo=brazo, km_internal=km_i,
                    dist=dist,
                    km_oficial=self.river_offsets.get(river, 0.0) + km_i,
                    at_end=(km_i <= 1e-6 or km_i >= ax.length - 1e-6))

    def tangent_azimuth(self, centroid, from_crs=None, span: float = 50.0,
                        river=None):
        """
        Local direction of the river centerline at the point nearest `centroid`,
        as a math-convention angle (CCW from +East) IN THE COORDINATES OF
        `from_crs` — v6, used by the section-orientation QC.

        `centroid` is given in `from_crs` (the transect's UTM); it is projected
        into the network CRS to find the nearest point, and the tangent is then
        transformed back so the angle is directly comparable with the section
        azimuth. The tangent is a chord over +/- `span` metres of chainage,
        which smooths the vertex-to-vertex noise of a digitised centerline.
        """
        cx, cy = float(centroid[0]), float(centroid[1])
        tr_fwd = tr_back = None
        if from_crs is not None and CRS.from_user_input(from_crs) != self.crs:
            tr_fwd = Transformer.from_crs(from_crs, self.crs, always_xy=True)
            tr_back = Transformer.from_crs(self.crs, from_crs, always_xy=True)
            cx, cy = tr_fwd.transform(cx, cy)

        p = Point(cx, cy)
        best_line, best_d = None, float("inf")
        # 6.1: on the transect's own river when known — near a confluence the
        # nearest line can belong to another river
        for ax in self.axes:
            if river and river in self._by_river and ax.river != river:
                continue
            for g in [ax.main] + [ab["geom"] for ab in ax.anabranches]:
                d = g.distance(p)
                if d < best_d:
                    best_d, best_line = d, g
        if best_line is None:
            return None
        s0 = best_line.project(p)
        a = best_line.interpolate(max(0.0, s0 - span))
        b = best_line.interpolate(min(best_line.length, s0 + span))
        ax_, ay_, bx_, by_ = a.x, a.y, b.x, b.y
        if tr_back is not None:
            ax_, ay_ = tr_back.transform(ax_, ay_)
            bx_, by_ = tr_back.transform(bx_, by_)
        if not (np.isfinite(bx_ - ax_) and np.isfinite(by_ - ay_)):
            return None
        if abs(bx_ - ax_) < 1e-9 and abs(by_ - ay_) < 1e-9:
            return None
        return math.atan2(by_ - ay_, bx_ - ax_)

    def snapped_xy(self, x, y):
        """Point on the network nearest to (x, y), in `crs` — for QA plots."""
        p = Point(x, y)
        best_line, best_d = None, float("inf")
        for ax in self.axes:
            for g in [ax.main] + [ab["geom"] for ab in ax.anabranches]:
                d = g.distance(p)
                if d < best_d:
                    best_d, best_line = d, g
        q = best_line.interpolate(best_line.project(p))
        return q.x, q.y

    # --------------- continuous internal chainage across the survey --------- #
    def _node_key(self, xy):
        return (round(xy[0] / NODE_SNAP_TOL), round(xy[1] / NODE_SNAP_TOL))

    def build_survey_chain(self, rivers_present):
        """LEGACY (v5/v6.0), kept verbatim for API compatibility; the 6.1 pipeline
        uses the head-to-outlet PATHS instead (see _build_topology), which need
        no trunk choice. Order the surveyed rivers into a single downstream
        chain and set the cumulative INTERNAL offsets used for the longitudinal
        profile.

        Returns the ordered list of river names, or None when the surveyed rivers
        do not form a single unambiguous chain — then the profile falls back to
        per-river official km and a warning is logged. That ambiguous case (e.g.
        surveying BOTH the Neuquén and the Limay above the confluence) is what the
        deferred `mainstem` declaration will resolve; it is intentionally not
        wired here.
        """
        present = list(dict.fromkeys(r for r in rivers_present if r in self._by_river))
        if not present:
            self.survey_offsets, self.survey_ok = None, False
            return None
        if len(present) == 1:
            self.survey_offsets = {present[0]: 0.0}
            self.survey_ok = True
            return list(present)

        us = {r: self._node_key(self._by_river[r].us_xy) for r in present}
        ds = {r: self._node_key(self._by_river[r].ds_xy) for r in present}
        ds_nodes = set(ds.values())

        heads = [r for r in present if us[r] not in ds_nodes]
        if len(heads) != 1:
            self._log.append(
                f"[warn] surveyed rivers {present} have {len(heads)} upstream head(s), "
                "not 1 — longitudinal profile stays per-river (declare a `mainstem` "
                "when the Limay is surveyed)")
            self.survey_offsets, self.survey_ok = None, False
            return None

        chain, current, seen = [], heads[0], set()
        while current is not None and current not in seen:
            chain.append(current)
            seen.add(current)
            nxt = [r for r in present if us[r] == ds[current] and r != current]
            if len(nxt) == 1:
                current = nxt[0]
            elif len(nxt) == 0:
                current = None
            else:
                self._log.append(
                    f"[warn] ambiguous downstream continuation after '{current}' "
                    f"({nxt}) — profile stays per-river (declare a `mainstem`)")
                self.survey_offsets, self.survey_ok = None, False
                return None

        if set(chain) != set(present):
            self._log.append(
                f"[warn] surveyed rivers {present} are not one connected chain "
                f"(chained: {chain}) — profile stays per-river")
            self.survey_offsets, self.survey_ok = None, False
            return None

        offsets, acc = {}, 0.0
        for r in chain:
            offsets[r] = acc
            acc += self._by_river[r].length
        self.survey_offsets, self.survey_ok = offsets, True
        self._log.append("[info] survey chain (downstream): "
                         + " -> ".join(f"{r}(+{offsets[r]:.0f} m)" for r in chain))
        return chain

    def survey_chainage(self, river, km_internal):
        """Continuous internal chainage along the surveyed chain, or None if the
        chain is undefined (fall back to per-river official km)."""
        if not self.survey_ok or river not in (self.survey_offsets or {}):
            return None
        return self.survey_offsets[river] + km_internal

    # ------------------------------------------------------------------ #
    #  6.1 — river-to-river topology and head-to-outlet PATHS
    # ------------------------------------------------------------------ #
    def _build_topology(self):
        """Link every river's downstream end to the river it flows into and
        enumerate the head-to-outlet PATHS (6.1).

        A river r flows into d when r's downstream end lies on d's main line
        (within NODE_SNAP_TOL) at a chainage s_join that d CONTINUES past — a
        node where two rivers END (Neuquen and Limay) is a shared inflow, not a
        link between them. s_join is 0 when d starts there (Negro) and > 0 for a
        tributary joining mid-reach.

        A PATH starts at a head (a river no other river flows into at its km 0)
        and follows the links to the outlet: Neuquen>Negro, Limay>Negro. Along a
        path the chainage p is continuous. Rivers downstream of a confluence
        belong to every path through it: that shared reach is where the water
        surfaces of the tributaries meet, so no trunk has to be chosen."""
        self.links, self.paths, self.junctions = {}, [], []
        for ax in self.axes:
            pt = Point(ax.ds_xy)
            cands = []
            for other in self.axes:
                if other is ax:
                    continue
                d = float(other.main.distance(pt))
                if d > NODE_SNAP_TOL:
                    continue
                s_join = float(other.main.project(pt))
                if s_join >= other.length - NODE_SNAP_TOL:
                    continue                    # other river also ends here
                cands.append((d, other.river, s_join))
            if len(cands) == 1:
                self.links[ax.river] = dict(down=cands[0][1], s_join=cands[0][2])
            elif len(cands) > 1:
                self._log.append(
                    f"[warn] río '{ax.river}': its downstream end touches "
                    f"{[c[1] for c in cands]} — ambiguous outflow, treated as an outlet "
                    "(check the centerline linework / flow_dir)")

        into_start = {r for r, lk in self.links.items() if lk["s_join"] <= NODE_SNAP_TOL}
        starts_fed = {self.links[r]["down"] for r in into_start}
        heads = [ax.river for ax in self.axes if ax.river not in starts_fed]
        for h in heads:
            segs, p, cur, seen = [], 0.0, h, set()
            s0 = 0.0
            while cur is not None and cur not in seen:
                seen.add(cur)
                L = self._by_river[cur].length
                segs.append(dict(river=cur, s0=s0, s1=L, p0=p))
                p += L - s0
                lk = self.links.get(cur)
                if lk is None:
                    break
                if lk["down"] in seen:
                    self._log.append(f"[warn] centerline loop at '{lk['down']}' — path "
                                     f"from '{h}' stopped there")
                    break
                cur, s0 = lk["down"], lk["s_join"]
            self.paths.append(dict(name=">".join(sg["river"] for sg in segs),
                                   head=h, segs=segs, length=p))
        self.path_order = list(range(len(self.paths)))

        # junctions: (downstream river, s_join) with the rivers flowing in there
        jmap = {}
        for r, lk in self.links.items():
            key = None
            for k in jmap:
                if k[0] == lk["down"] and abs(k[1] - lk["s_join"]) <= NODE_SNAP_TOL:
                    key = k
                    break
            key = key or (lk["down"], lk["s_join"])
            jmap.setdefault(key, []).append(r)
        for (d, sj), ups in jmap.items():
            # distance from the junction to the outlet along its downstream path
            rest = 0.0
            cur, s_from = d, sj
            seen = set()
            while cur is not None and cur not in seen:
                seen.add(cur)
                rest += self._by_river[cur].length - s_from
                lk = self.links.get(cur)
                cur, s_from = (lk["down"], lk["s_join"]) if lk else (None, 0.0)
            self.junctions.append(dict(river=d, s=float(sj), ups=sorted(ups),
                                       mid_reach=sj > NODE_SNAP_TOL, to_outlet=rest))
        self.junctions.sort(key=lambda j: -j["to_outlet"])     # upstream first

        for ax in self.axes:
            lk = self.links.get(ax.river)
            self._log.append(
                f"[info] topology: '{ax.river}' -> "
                + (f"'{lk['down']}' at its km {lk['s_join']:.0f} m" if lk else "outlet"))
        for pa in self.paths:
            self._log.append(f"[info] path {pa['name']}: {pa['length'] / 1000.0:.2f} km")

    def path_positions(self, river, km_internal, tol: float = 1.0):
        """[(path index, p)] of every path containing (river, km_internal)."""
        out = []
        if river is None or km_internal is None:
            return out
        s = float(km_internal)
        for i, pa in enumerate(self.paths):
            for sg in pa["segs"]:
                if sg["river"] == river and sg["s0"] - tol <= s <= sg["s1"] + tol:
                    out.append((i, sg["p0"] + (s - sg["s0"])))
                    break
        return out

    def rank_paths(self, located, log=None):
        """Order the paths for display (long profile, index sort) by how many of
        the `located` transects [(river, km_internal)] they hold. Only the
        display depends on this order — the water surface uses every path."""
        cnt = [0] * len(self.paths)
        for river, km in located:
            for i, _p in self.path_positions(river, km):
                cnt[i] += 1
        self.path_order = sorted(range(len(self.paths)),
                                 key=lambda i: (-cnt[i], self.paths[i]["name"]))
        if log is not None and self.paths:
            log.append("[info] paths by transect count: "
                       + ", ".join(f"{self.paths[i]['name']} ({cnt[i]})"
                                   for i in self.path_order))
        return cnt

    def display_position(self, river, km_internal):
        """(rank, path name, p) on the first path, in display order, holding
        (river, km_internal); (None, '', None) if none."""
        pos = dict(self.path_positions(river, km_internal))
        for rank, i in enumerate(self.path_order):
            if i in pos:
                return rank, self.paths[i]["name"], float(pos[i])
        return None, "", None

    def locate_all(self, x, y):
        """locate_local() of (x, y) on EVERY river, nearest first (6.1)."""
        out = []
        for ax in self.axes:
            km_i, brazo, role, dist = ax.locate_local(x, y)
            out.append(dict(river=ax.river, role=role, brazo=brazo, km_internal=km_i,
                            dist=dist, km_oficial=self.river_offsets.get(ax.river, 0.0) + km_i,
                            at_end=(km_i <= 1e-6 or km_i >= ax.length - 1e-6)))
        out.sort(key=lambda d: d["dist"])
        return out

    def locate_track(self, xs, ys, force_river=None, method="crossing"):
        """Locate a TRANSECT from its track (6.1): river, brazo and chainage where
        the SECTION LINE crosses a channel centerline (the HEC-RAS convention),
        not the channel nearest to the track centroid.

        Near a confluence the nearest-centroid rule is unreliable: a section
        across the Negro a few metres below the node, with a track covering
        mostly one bank, has its centroid closer to the end of the Limay's line
        and was reported as 'Limay'. The section line is the principal axis of
        the track, extended by max(TRACK_LOC_MARGIN_MIN, FRAC x length) on each
        side; if it crosses several channels the crossing nearest to the track
        centroid wins. If it crosses none, or for method='centroid' (the v6.0
        rule), the nearest axis to the centroid is used. `force_river` (the INI
        [rios] override) restricts the search to that river.

        Returns dict(river, role, brazo, km_internal, km_oficial, dist, method,
        note); `dist` is the distance from the track centroid to the crossing."""
        x = np.asarray(xs, dtype=float).ravel()
        y = np.asarray(ys, dtype=float).ravel()
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        if x.size == 0:
            raise ValueError("track has no finite position")
        cx, cy = float(x.mean()), float(y.mean())
        if force_river:
            nearest = self.locate_on(force_river, cx, cy)
        else:
            nearest = self.locate(cx, cy)
        tag = "manual" if force_river else None

        def _fallback(meth, note):
            out = dict(nearest)
            out.pop("at_end", None)
            out.update(method=tag or meth, note=note)
            return out

        if method == "centroid":
            return _fallback("centroid", "")
        if x.size < 3:
            return _fallback("nearest", "track too short for a section line")
        P = np.column_stack([x - cx, y - cy])
        w, v = np.linalg.eigh(P.T @ P)
        u = v[:, int(np.argmax(w))]
        t = P @ u
        lo, hi = float(t.min()), float(t.max())
        if hi - lo < 1.0:
            return _fallback("nearest", "track too short for a section line")
        m = max(TRACK_LOC_MARGIN_MIN, TRACK_LOC_MARGIN_FRAC * (hi - lo))
        seg = LineString([(cx + (lo - m) * u[0], cy + (lo - m) * u[1]),
                          (cx + (hi + m) * u[0], cy + (hi + m) * u[1])])

        def _points(g):
            if g is None or g.is_empty:
                return []
            if g.geom_type == "Point":
                return [g]
            if g.geom_type in ("MultiPoint", "GeometryCollection", "MultiLineString"):
                return [q for gg in g.geoms for q in _points(gg)]
            if g.geom_type == "LineString":          # collinear overlap
                return [g.interpolate(0.5, normalized=True)]
            return []

        cands = []
        for ax in self.axes:
            if force_river and ax.river != force_river:
                continue
            for g in [ax.main] + [ab["geom"] for ab in ax.anabranches]:
                for q in _points(seg.intersection(g)):
                    tq = (q.x - cx) * u[0] + (q.y - cy) * u[1]
                    cands.append((abs(tq), ax.river, float(q.x), float(q.y)))
        if not cands:
            if force_river:
                return _fallback("manual", f"forced to {force_river} ([rios]); the section "
                                           "does not cross its axis, km from the track "
                                           "centroid projected on it")
            return _fallback("nearest", "the section line crosses no centerline — "
                                        "nearest axis to the track centroid used")
        cands.sort(key=lambda c: c[0])
        d0, river, qx, qy = cands[0]
        km_i, brazo, role, _d = self._by_river[river].locate_local(qx, qy)
        notes = []
        others = sorted({c[1] for c in cands if c[1] != river})
        if others:
            notes.append("section also crosses " + ", ".join(others))
        if not force_river and nearest["river"] != river:
            notes.append(f"nearest axis to the track centroid is {nearest['river']} "
                         f"({nearest['dist']:.0f} m) — the crossing wins")
        return dict(river=river, role=role, brazo=brazo, km_internal=float(km_i),
                    km_oficial=self.river_offsets.get(river, 0.0) + float(km_i),
                    dist=float(d0), method=tag or "crossing", note="; ".join(notes))


def read_ws_csv(csv_path, east_col="East", north_col="North", elev_col="H_correg"):
    """Read a water-surface CSV -> list of (x, y, elevation) in the file's CRS."""
    pts, _skipped, _cols = read_ws_table(csv_path, east_col, north_col, elev_col)
    return [(p["x"], p["y"], p["h"]) for p in pts]


def read_ws_table(csv_path, east_col="East", north_col="North", elev_col="H_correg"):
    """Read a water-surface CSV keeping each point's identifier (6.1).

    Returns (points, skipped, columns): points = [dict(id, x, y, h)] in the
    file's CRS; skipped = [(row_label, reason)] for rows that could not be
    parsed, so they are REPORTED instead of dropped silently. The id comes from
    a punto/point/id/name column if present, else the 1-based data row number.
    A UTF-8 BOM (typical of Excel exports) is handled."""
    pts, skipped = [], []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rd = csv.DictReader(f)
        cols = [c.strip() for c in (rd.fieldnames or [])]
        rd.fieldnames = cols
        ec = _pick_col(cols, east_col, ["east", "x", "este", "x_utm", "x_posgar07"])
        nc = _pick_col(cols, north_col, ["north", "y", "norte", "y_utm", "y_posgar07"])
        hc = _pick_col(cols, elev_col, ["h_correg", "z", "elev", "cota", "h", "altura"])
        ic = _pick_col(cols, None, ["punto", "point", "id", "name", "nombre", "pto", "pt"])
        rc = _pick_col(cols, None, ["rio", "río", "river"])      # 6.1: optional override
        if not (ec and nc and hc):
            raise ValueError(f"WS CSV: could not resolve E/N/elev columns among {cols}")
        for k, r in enumerate(rd, 1):
            label = (str(r.get(ic) or "").strip() if ic else "") or f"row{k}"
            try:
                x, y, h = float(r[ec]), float(r[nc]), float(r[hc])
            except (TypeError, ValueError):
                skipped.append((label, f"non-numeric {ec}/{nc}/{hc}"))
                continue
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(h)):
                skipped.append((label, "non-finite value"))
                continue
            river = (str(r.get(rc) or "").strip() if rc else "")
            pts.append(dict(id=label, x=x, y=y, h=h, river=river))
    if len(pts) < 2:
        raise ValueError("WS CSV: need at least 2 valid points")
    return pts, skipped, dict(east=ec, north=nc, elev=hc, id=ic, river=rc)


class WaterSurfaceProfile:
    """Longitudinal water-surface profile: elevation as a function of the survey's
    continuous INTERNAL chainage.

    Built from GNSS water-surface points already reduced to (chainage, elevation)
    pairs (the caller projects each point onto the network and applies the survey
    chain). This is the PRIMARY water surface. Interpolating in chainage — not in
    2-D space — honours the longitudinal slope and meanders, exactly as in v4.
    """

    def __init__(self, pairs, log=None):
        """`pairs` = iterable of (chainage_m, elevation_m[, dist_m])."""
        self._log = log if log is not None else []
        rows = [(float(p[0]), float(p[1]),
                 float(p[2]) if len(p) > 2 and p[2] is not None else float("nan"))
                for p in pairs
                if p[0] is not None and np.isfinite(p[0]) and np.isfinite(p[1])]
        if len(rows) < 2:
            raise ValueError("WS: need at least 2 valid (chainage, elev) points")
        prog = np.asarray([r[0] for r in rows])
        elev = np.asarray([r[1] for r in rows])
        dist = np.asarray([r[2] for r in rows])
        order = np.argsort(prog)
        prog, elev, dist = prog[order], elev[order], dist[order]

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
        self.max_dist = float(np.nanmax(dist)) if np.isfinite(dist).any() else float("nan")
        self.slope_m_per_km = float(np.polyfit(self.prog, self.elev, 1)[0] * 1000.0)

        self._log.append(
            f"[info] water surface (GNSS): {len(self.prog)} pts, chainage "
            f"{self.prog.min():.0f}–{self.prog.max():.0f} m, elev "
            f"{self.elev.min():.3f}–{self.elev.max():.3f} m, slope "
            f"{self.slope_m_per_km:.2f} m/km"
            + ("" if not np.isfinite(self.max_dist)
               else f", max offset from axis {self.max_dist:.0f} m"))

    def covers(self, chainage):
        return float(self.prog[0]) <= chainage <= float(self.prog[-1])

    def elev_at(self, chainage):
        """(ws_elev, note). `note` flags extrapolation beyond the surveyed range."""
        lo, hi = float(self.prog[0]), float(self.prog[-1])
        if chainage < lo:
            return float(self.elev[0]), f"EXTRAP {chainage - lo:.0f} m (below WS range)"
        if chainage > hi:
            return float(self.elev[-1]), f"EXTRAP +{chainage - hi:.0f} m (above WS range)"
        return float(np.interp(chainage, self.prog, self.elev)), ""


# =============================================================================== #
#  Hydrometric stations (secondary water-surface source: fallback / gap-fill / QC)
# =============================================================================== #

STN_ATTR_ID    = ("station_id", "stationid", "id", "codigo", "código")
STN_ATTR_NAME  = ("name", "nombre")
STN_ATTR_RIVER = ("river", "rio", "río")
STN_ATTR_ZERO  = ("gauge_zero", "gaugezero", "cero", "cero_escala", "z0")
STN_ATTR_CHAIN = ("chainage", "progresiva", "km")


def _resolve_col(cols_lower, names):
    for n in names:
        if n in cols_lower:
            return cols_lower[n]
    return None


def _parse_dt(s):
    """Parse a reading timestamp into np.datetime64, or None."""
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M",
                "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%Y-%m-%d"):
        try:
            return np.datetime64(datetime.datetime.strptime(str(s).strip(), fmt))
        except ValueError:
            continue
    return None


STN_ATTR_X = ("x", "este", "east", "e", "x_posgar07", "x_utm", "coord_x")
STN_ATTR_Y = ("y", "norte", "north", "n", "y_posgar07", "y_utm", "coord_y")


def _norm_name(s):
    """Accent- and case-insensitive key for river names ('Neuquén' == 'neuquen')."""
    import unicodedata
    t = unicodedata.normalize("NFKD", str(s or "")).encode("ascii", "ignore").decode()
    return t.strip().lower()


def _read_station_table(path, crs_default):
    """Read a station registry from a CSV/TXT (X/Y columns, coordinates in
    `crs_default`) or from any vector file geopandas can open (6.1).

    geopandas/pyogrio return a plain DataFrame — no geometry, no .crs, every
    field a string — for a CSV; that is what used to crash the registry.
    Returns (rows, crs): rows = [dict(attrs={lower_col: value}, x, y)]."""
    path = Path(path)
    if path.suffix.lower() in (".csv", ".txt"):
        with open(path, newline="", encoding="utf-8-sig") as f:
            rd = csv.DictReader(f)
            recs = [{str(k).strip().lower(): v for k, v in r.items() if k is not None}
                    for r in rd]
        crs = CRS.from_user_input(crs_default)
        geom = None
    else:
        gdf = gpd.read_file(path)
        recs = [{str(k).strip().lower(): v for k, v in r.items() if k != "geometry"}
                for r in gdf.to_dict("records")]
        has_geom = (isinstance(gdf, gpd.GeoDataFrame) and "geometry" in gdf.columns
                    and gdf.geometry.notna().any())
        geom = list(gdf.geometry) if has_geom else None
        crs = (gdf.crs if has_geom and gdf.crs is not None
               else CRS.from_user_input(crs_default))
    rows = []
    cols = set(recs[0].keys()) if recs else set()
    cx = next((c for c in STN_ATTR_X if c in cols), None)
    cy = next((c for c in STN_ATTR_Y if c in cols), None)
    for i, a in enumerate(recs):
        x = y = None
        if geom is not None and geom[i] is not None and not geom[i].is_empty:
            g = geom[i] if geom[i].geom_type == "Point" else geom[i].representative_point()
            x, y = float(g.x), float(g.y)
        elif cx and cy:
            try:
                x, y = float(a[cx]), float(a[cy])
            except (TypeError, ValueError):
                x = y = None
        rows.append(dict(attrs=a, x=x, y=y))
    if geom is None and not (cx and cy):
        raise ValueError(f"stations: no geometry and no X/Y columns among {sorted(cols)}")
    return rows, crs


class StationRegistry:
    """Hydrometric stations (points): station_id, name, river, gauge_zero
    (cero de escala, cota SRVN16), optional chainage. Each station is placed on
    the network (river + internal chainage) so its readings interpolate in the
    same chainage frame as the transects. ws_elev of a reading = gauge_zero +
    nivel. The registry is built once and reused between campaigns.

    6.1: CSV registries (X/Y columns in --ws-crs) are read directly; a station
    is projected on its DECLARED river's axis (not the nearest axis of any
    river); stations off their axis, beyond the digitised reach, on a river the
    centerline does not have, or without gauge_zero are kept out of `stations`
    and listed in `excluded` with the reason."""

    def __init__(self, stations_path, network, crs_override="EPSG:5344", log=None):
        if not HAS_GPD:
            raise RuntimeError("StationRegistry needs geopandas/shapely")
        self._log = log if log is not None else []
        rows, src_crs = _read_station_table(stations_path, crs_override)
        tr = None
        if network.crs is not None and CRS.from_user_input(src_crs) != network.crs:
            tr = Transformer.from_crs(src_crs, network.crs, always_xy=True)

        cols = {c: c for c in (rows[0]["attrs"].keys() if rows else [])}
        c_id = _resolve_col(cols, STN_ATTR_ID)
        c_nm = _resolve_col(cols, STN_ATTR_NAME)
        c_rv = _resolve_col(cols, STN_ATTR_RIVER)
        c_z0 = _resolve_col(cols, STN_ATTR_ZERO)
        c_ch = _resolve_col(cols, STN_ATTR_CHAIN)
        if not (c_id and c_z0):
            raise ValueError("stations: need at least station_id and gauge_zero "
                             f"columns among {sorted(cols)}")
        by_norm = {_norm_name(ax.river): ax.river for ax in network.axes}

        def _txt(v):
            t = "" if v is None else str(v).strip()
            return "" if t.lower() in ("nan", "none") else t

        self.stations, self.excluded = {}, {}
        for row in rows:
            a = row["attrs"]
            sid = _txt(a.get(c_id))
            if not sid:
                continue
            name = _txt(a.get(c_nm)) if c_nm else ""
            name = name or sid
            rec = dict(station_id=sid, name=name, river="", gauge_zero=None,
                       x=row["x"], y=row["y"], km_internal=None, dist_axis=float("nan"),
                       survey_chainage=None, reason="")
            try:
                gz = float(_txt(a.get(c_z0)))
                rec["gauge_zero"] = gz if np.isfinite(gz) else None
            except (TypeError, ValueError):
                rec["gauge_zero"] = None
            if rec["gauge_zero"] is None:
                rec["reason"] = "no gauge_zero (not levelled)"
            if row["x"] is None or row["y"] is None:
                rec["reason"] = rec["reason"] or "no coordinates"
            else:
                x, y = row["x"], row["y"]
                if tr is not None:
                    x, y = tr.transform(x, y)
                rec["x"], rec["y"] = float(x), float(y)
                declared = _txt(a.get(c_rv)) if c_rv else ""
                river = by_norm.get(_norm_name(declared)) if declared else None
                if declared and river is None:
                    rec["river"] = declared
                    rec["reason"] = rec["reason"] or (
                        f"river '{declared}' is not in the centerline")
                else:
                    if river is None:            # undeclared -> nearest axis
                        river = network.locate(x, y)["river"]
                    rec["river"] = river
                    loc = network.locate_on(river, x, y)
                    rec["dist_axis"] = float(loc["dist"])
                    km_explicit = None
                    if c_ch and _txt(a.get(c_ch)):
                        try:
                            km_explicit = float(_txt(a.get(c_ch)))
                        except (TypeError, ValueError):
                            km_explicit = None
                    if km_explicit is not None:
                        rec["km_internal"] = km_explicit
                    else:
                        rec["km_internal"] = float(loc["km_internal"])
                        if loc["at_end"] and loc["dist"] > STATION_END_TOL:
                            rec["reason"] = rec["reason"] or (
                                f"beyond the digitised {river} reach (projects onto "
                                f"the axis end, {loc['dist']:.0f} m away)")
                        elif loc["dist"] > STATION_AXIS_TOL:
                            rec["reason"] = rec["reason"] or (
                                f"{loc['dist']:.0f} m from the {river} axis "
                                f"(> {STATION_AXIS_TOL:.0f} m)")
            if rec["reason"]:
                self.excluded[sid] = rec
            else:
                self.stations[sid] = rec

        for sid, rec in self.excluded.items():
            self._log.append(f"[info] station {sid} ({rec['name']}): {rec['reason']} "
                             "— not used for the water surface")
        if not self.stations:
            raise ValueError("stations file has no usable station (see the [info] "
                             "lines above for why each one was excluded)")
        self._log.append(
            f"[info] stations: {len(self.stations)} usable of "
            f"{len(self.stations) + len(self.excluded)} — "
            + ", ".join(f"{s['station_id']}[{s['river']} km {s['km_internal']:.0f}, "
                        f"{s['dist_axis']:.0f} m off axis]"
                        for s in self.stations.values()))

    def attach_survey_chainage(self, network):
        """Set survey_chainage on each station (after network.build_survey_chain)."""
        for s in self.stations.values():
            s["survey_chainage"] = network.survey_chainage(s["river"], s["km_internal"])


def read_level_readings(csv_path):
    """Read a per-campaign level-readings CSV. Auto-detects the mode:
        spot   -> columns: station_id, nivel               (one reading per station)
        series -> columns: station_id, datetime, nivel     (several per station)
    Returns (mode, data):
        spot   : {station_id: nivel}
        series : {station_id: [(np.datetime64|None, nivel), ...] sorted by time}
    `nivel` is the staff-gauge reading [m]; ws_elev = gauge_zero + nivel is applied
    later against the station registry."""
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rd = csv.DictReader(f)
        cols = {c.lower().strip(): c for c in (rd.fieldnames or [])}
        c_id = _resolve_col(cols, ("station_id", "stationid", "id", "estacion", "estación"))
        c_lv = _resolve_col(cols, ("nivel", "stage", "lectura", "level", "h"))
        c_dt = _resolve_col(cols, ("datetime", "fecha_hora", "fecha", "timestamp", "time"))
        if not (c_id and c_lv):
            raise ValueError("readings CSV: need station_id and nivel columns among "
                             f"{list(rd.fieldnames or [])}")
        rows = []
        for r in rd:
            sid = (r.get(c_id) or "").strip()
            if not sid:
                continue
            try:
                lv = float(r[c_lv])
            except (TypeError, ValueError):
                continue
            dt = _parse_dt(r[c_dt]) if (c_dt and (r.get(c_dt) or "").strip()) else None
            rows.append((sid, dt, lv))

    if not rows:
        raise ValueError("readings CSV: no valid rows")
    per = {}
    for sid, dt, lv in rows:
        per.setdefault(sid, []).append((dt, lv))

    has_dt = (c_dt is not None) and any(dt is not None for lst in per.values() for dt, _ in lst)
    multi = any(len(v) > 1 for v in per.values())
    if has_dt and multi:
        series = {}
        for sid, lst in per.items():
            ordered = sorted((t for t in lst if t[0] is not None), key=lambda t: t[0])
            series[sid] = ordered if ordered else [(None, lst[0][1])]
        return "series", series
    spot = {sid: float(np.mean([lv for _, lv in lst])) for sid, lst in per.items()}
    return "spot", spot


class StationWS:
    """Water surface from hydrometric stations — SECONDARY source (fallback when
    there is no GNSS, gap-fill beyond the GNSS range instead of the v4 clamp, and
    QC against the GNSS surface). Linearly interpolates ws_elev in survey chainage
    between the two bracketing stations of the SAME river; clamps outside their
    span with a note. Downstream monotonicity is logged, not enforced."""

    def __init__(self, registry, mode, data, log=None):
        self._log = log if log is not None else []
        self.registry = registry
        self.mode = mode
        self.data = data
        # 6.1: a station needs a reading, not a survey chainage — the water-surface
        # model places it on its river and checks monotonicity per network path
        usable = [s for s in registry.stations.values() if s["station_id"] in data]
        self._log.append(
            f"[info] gauge readings: mode={mode}, {len(usable)}/{len(registry.stations)} "
            "levelled station(s) have a reading"
            + (": " + ", ".join(f"{s['station_id']} = {s['gauge_zero']:.3f} + "
                                f"{self._stage_at(s['station_id'], None):.3f} m"
                                for s in usable
                                if self._stage_at(s['station_id'], None) is not None)
               if usable else ""))

    def _stage_at(self, sid, when):
        if self.mode == "spot":
            return float(self.data[sid]) if sid in self.data else None
        seq = self.data.get(sid, [])
        if not seq:
            return None
        ts = [t for t, _ in seq if t is not None]
        vs = [v for t, v in seq if t is not None]
        if not ts:
            return float(np.mean([v for _, v in seq]))
        if when is None:
            return float(np.mean(vs))
        tt = np.asarray([t.astype("datetime64[s]").astype("float64") for t in ts])
        w = float(np.datetime64(when).astype("datetime64[s]").astype("float64"))
        return float(np.interp(w, tt, np.asarray(vs, float)))

    def elev_at(self, survey_chainage, river, when=None):
        """(ws_elev, note) from the stations of `river`; (None, note) if it cannot
        be evaluated (no station on that river / no survey chainage)."""
        if survey_chainage is None:
            return None, "no survey chainage"
        pts = []
        for s in self.registry.stations.values():
            if river and s["river"] and s["river"] != river:
                continue
            if s.get("survey_chainage") is None:
                continue
            stage = self._stage_at(s["station_id"], when)
            if stage is None:
                continue
            pts.append((float(s["survey_chainage"]), float(s["gauge_zero"] + stage)))
        if not pts:
            return None, "no station for this river"
        pts.sort(key=lambda p: p[0])
        xs = np.asarray([p[0] for p in pts])
        ys = np.asarray([p[1] for p in pts])
        if len(pts) == 1:
            return float(ys[0]), f"single station (offset {survey_chainage - xs[0]:+.0f} m)"
        lo, hi = float(xs[0]), float(xs[-1])
        if survey_chainage < lo:
            return float(ys[0]), f"EXTRAP {survey_chainage - lo:.0f} m (below stations)"
        if survey_chainage > hi:
            return float(ys[-1]), f"EXTRAP +{survey_chainage - hi:.0f} m (above stations)"
        return float(np.interp(survey_chainage, xs, ys)), ""


class WaterSurfaceModel:
    """Water surface on the river NETWORK from GNSS points (primary) and gauge
    readings (secondary) — 6.1.

    Each river keeps its own water surface on its own chainage. Rivers are joined
    along the head-to-outlet PATHS of the network (Neuquen>Negro, Limay>Negro):
    the reach below a confluence belongs to every path through it, so it is the
    COMMON downstream river where the tributaries' surfaces meet, and no trunk
    has to be chosen.

      1. ANCHORS — GNSS points (on the nearest river axis, or on the river named
         in an optional 'rio' column of the CSV) and the gauges NOT bracketed by
         GNSS points along some path; a gauge between GNSS points is QC only.
      2. CONFLUENCES — the stage at each junction is estimated from EACH incoming
         branch that has values: linear interpolation along that branch's path
         between its last value upstream and the first value downstream of the
         node, i.e. each tributary's slope continued into the common river. The
         estimates are averaged (weight 1/interpolation gap) into ONE junction
         anchor: the common reach gets a single surface and every tributary
         meets it continuously. Their spread is QC'd against --ws-qc-tol. With
         nothing downstream, the tributaries extrapolate down to the node; with
         nothing upstream, the downstream surface is continued up to it.
      3. LOOKUP, on every path through the point (averaged if several):
           between two GNSS values             -> interpolation  source=survey
           mixing gauge / junction values      -> interpolation  source=bridge
           between two gauge values            -> interpolation  source=station
           beyond the first/last value         -> SLOPE continuation  =extrap
         The extrapolation continues the least-squares slope of the values
         nearest to that end (spanning >= WS_SLOPE_WINDOW): a tributary without
         values of its own takes the junction stage and continues upstream with
         the slope of the common reach. It is a warning beyond --ws-extrap-tol.
         v6.0 clamped a FLAT value beyond the last point (0.9 m/km x distance of
         error here).

    `rows` holds one record per GNSS point, gauge and junction for
    water_surface_qc.csv (river, km, distance to axis, alternative river,
    role, reason)."""

    def __init__(self, network, points=None, station_ws=None, qc_tol=0.10,
                 skipped=None, log=None, extrap_tol=WS_EXTRAP_TOL_DEFAULT):
        self._log = log if log is not None else []
        self.network = network
        self.station_ws = station_ws
        self.qc_tol = float(qc_tol)
        self.extrap_tol = float(extrap_tol)
        self.rows = []
        self.gnss, self.gauges, self.qc_gauges = [], [], []
        self._cache = {}
        by_norm = {_norm_name(ax.river): ax.river for ax in network.axes}
        nan = float("nan")

        # ---- GNSS points --------------------------------------------------
        raw, ambiguous = [], []
        for pt in points or []:
            forced = None
            if pt.get("river"):
                forced = by_norm.get(_norm_name(pt["river"]))
                if forced is None:
                    self._log.append(f"[warn] WS point {pt['id']}: river '{pt['river']}' "
                                     "('rio' column) is not in the centerline — located "
                                     "geometrically")
            cands = network.locate_all(pt["x"], pt["y"])
            loc = next(c for c in cands if c["river"] == forced) if forced else cands[0]
            alt = next((c for c in cands if c["river"] != loc["river"]), None)
            row = dict(kind="gnss", id=pt["id"], x=pt["x"], y=pt["y"], elev=pt["h"],
                       river=loc["river"], km_internal=loc["km_internal"],
                       km_oficial=loc["km_oficial"], dist_axis=loc["dist"],
                       alt_river=alt["river"] if alt else "",
                       alt_dist=alt["dist"] if alt else nan,
                       role="anchor", note="river from the 'rio' column" if forced else "")
            if loc["at_end"] and loc["dist"] > STATION_END_TOL:
                row["role"] = "dropped"
                row["note"] = (f"beyond the digitised {loc['river']} reach (projects onto "
                               f"the axis end, {loc['dist']:.0f} m away)")
            elif loc["dist"] > STATION_AXIS_TOL:
                row["role"] = "dropped"
                row["note"] = (f"{loc['dist']:.0f} m from the {loc['river']} axis "
                               f"(> {STATION_AXIS_TOL:.0f} m)")
            else:
                if not forced and alt is not None \
                        and alt["dist"] - loc["dist"] < WS_AMBIG_MARGIN:
                    row["note"] = (f"AMBIGUOUS: {alt['river']} axis at {alt['dist']:.0f} m "
                                   f"vs {loc['dist']:.0f} m — fill the 'rio' column to force")
                    ambiguous.append(f"{pt['id']} ({loc['river']} {loc['dist']:.0f} m / "
                                     f"{alt['river']} {alt['dist']:.0f} m)")
                raw.append(dict(id=str(pt["id"]), river=loc["river"],
                                s=float(loc["km_internal"]), h=float(pt["h"]),
                                dist=float(loc["dist"])))
            if row["role"] == "dropped":
                self._log.append(f"[warn] WS point {pt['id']}: {row['note']} — not used")
            self.rows.append(row)
        for label, why in (skipped or []):
            self.rows.append(dict(kind="gnss", id=label, x=None, y=None, elev=None,
                                  river="", km_internal=None, km_oficial=None,
                                  dist_axis=nan, alt_river="", alt_dist=nan,
                                  role="dropped", note=why))
            self._log.append(f"[warn] WS point {label}: {why} — ignored")
        if ambiguous:
            self._log.append("[warn] WS point(s) near a confluence with two candidate "
                             "rivers (see water_surface_qc.csv; a 'rio' column forces "
                             "the river): " + ", ".join(ambiguous))

        for river in dict.fromkeys(r["river"] for r in raw):      # merge coincident
            rs = sorted((r for r in raw if r["river"] == river), key=lambda r: r["s"])
            i = 0
            while i < len(rs):
                j = i
                while j + 1 < len(rs) and rs[j + 1]["s"] - rs[i]["s"] < 0.01:
                    j += 1
                grp = rs[i:j + 1]
                self.gnss.append(dict(kind="gnss", id="+".join(r["id"] for r in grp),
                                      river=river,
                                      s=float(np.mean([r["s"] for r in grp])),
                                      h=float(np.mean([r["h"] for r in grp])),
                                      deps=frozenset({"gnss"})))
                i = j + 1
            ss = np.asarray([r["s"] for r in rs])
            hh = np.asarray([r["h"] for r in rs])
            slope = (float(np.polyfit(ss, hh, 1)[0] * 1000.0)
                     if len(rs) >= 2 and np.ptp(ss) > 1.0 else nan)
            self._log.append(
                f"[info] water surface (GNSS), río '{river}': {len(rs)} pt(s), km "
                f"{ss.min():.0f}–{ss.max():.0f} m, elev {hh.min():.3f}–{hh.max():.3f} m"
                + ("" if not np.isfinite(slope) else f", slope {slope:.2f} m/km")
                + f", max offset from axis {max(r['dist'] for r in rs):.0f} m")

        # ---- gauges ---------------------------------------------------------
        if station_ws is not None:
            reg, data = station_ws.registry, station_ws.data
            for sid, st in reg.stations.items():
                row = dict(kind="gauge", id=sid, x=st["x"], y=st["y"], elev=None,
                           river=st["river"], km_internal=st["km_internal"],
                           km_oficial=(network.river_offsets.get(st["river"], 0.0)
                                       + st["km_internal"]),
                           dist_axis=st["dist_axis"], alt_river="", alt_dist=nan,
                           role="", note="")
                self.rows.append(row)
                stage = station_ws._stage_at(sid, None)
                if stage is None:
                    row["role"] = "no reading"
                    continue
                elev = float(st["gauge_zero"] + stage)
                row["elev"] = elev
                g = dict(kind="gauge", id=sid, river=st["river"],
                         s=float(st["km_internal"]), h=elev,
                         zero=float(st["gauge_zero"]), name=st["name"],
                         deps=frozenset({"gauge"}), row=row)
                if self._bracketed(g["river"], g["s"]):
                    row["role"] = "qc"
                    self.qc_gauges.append(g)
                else:
                    row["role"] = "anchor"
                    self.gauges.append(g)
            for sid, rec in reg.excluded.items():
                self.rows.append(dict(kind="gauge", id=sid, x=rec["x"], y=rec["y"],
                                      elev=None, river=rec["river"],
                                      km_internal=rec["km_internal"], km_oficial=None,
                                      dist_axis=rec["dist_axis"], alt_river="",
                                      alt_dist=nan, role="excluded", note=rec["reason"]))
                if sid in data:
                    self._log.append(f"[warn] gauge {sid} has a reading in this campaign "
                                     f"but is excluded: {rec['reason']}")
            for sid in data:
                if sid not in reg.stations and sid not in reg.excluded:
                    self._log.append(f"[warn] reading for station '{sid}', which is not "
                                     "in the station registry — ignored")
            if self.gauges:
                self._log.append(
                    "[info] gauge anchor(s) (not bracketed by GNSS): " + ", ".join(
                        f"{g['id']} {g['river']} km {g['s']:.0f} = {g['h']:.3f} m"
                        for g in self.gauges))

        # ---- solve (spot / mean stages): junctions, QC, monotonicity --------
        sol = self._solve(None)
        for jr in sol["junctions"]:
            J = jr["J"]
            parts = [f"{e['branch']}: {e['value']:.3f} ({e['kind']} {e['frm']}"
                     + (f"->{e['to']}" if e["to"] else "") + f" over {e['gap']:.0f} m)"
                     for e in jr["ests"]]
            note = "; ".join(parts) + f"; adopted {jr['h']:.3f} ({jr['how']})"
            self.rows.append(dict(kind="junction", id=jr["anchor"]["id"], x=None, y=None,
                                  elev=jr["h"], river=J["river"], km_internal=J["s"],
                                  km_oficial=network.river_offsets.get(J["river"], 0.0) + J["s"],
                                  dist_axis=nan, alt_river="", alt_dist=nan,
                                  role="junction", note=note))
            self._log.append(f"[info] confluence into '{J['river']}' at km {J['s']:.0f} "
                             f"({' + '.join(J['ups'])}): stage {jr['h']:.3f} m — {note}")
            ni = sum(1 for e in jr["ests"] if e["kind"] == "interp")
            if ni >= 2:
                tag = "PASS" if jr["spread"] <= self.qc_tol else "WARN"
                self._log.append(
                    f"[QA] confluence into '{J['river']}': the tributaries' estimates of the "
                    f"junction stage differ by {jr['spread']:.3f} m (tol {self.qc_tol:.2f}) "
                    f"{tag}")
            elif jr["how"] != "interp":
                self._log.append(f"[warn] confluence into '{J['river']}': stage "
                                 f"EXTRAPOLATED ({jr['how']}) — no water-surface value on "
                                 "one side of the node")
        for g in self.qc_gauges:
            res = self.resolve(g["river"], g["s"])
            if res["ws_elev"] is None:
                continue
            dif = g["h"] - res["ws_elev"]
            tag = "PASS" if abs(dif) <= self.qc_tol else "WARN"
            g["row"]["note"] = (f"GNSS {res['ws_elev']:.3f} m; gauge-GNSS {dif:+.3f} m "
                                f"(tol {self.qc_tol:.2f}) {tag}")
            self._log.append(
                f"[QA] gauge {g['id']} ({g['name']}) between GNSS points: {g['h']:.3f} m vs "
                f"GNSS {res['ws_elev']:.3f} m, Δ={dif:+.3f} m (tol {self.qc_tol:.2f}) {tag}")
        seen = set()
        for pi, seq in enumerate(sol["per_path"]):
            for (pa, a), (pb, b) in zip(seq, seq[1:]):
                if b["h"] - a["h"] > WS_MONO_TOL and (a["id"], b["id"]) not in seen:
                    seen.add((a["id"], b["id"]))
                    self._log.append(
                        f"[warn] water surface RISES downstream on "
                        f"{network.paths[pi]['name']}: {a['id']} {a['h']:.3f} m -> "
                        f"{b['id']} {b['h']:.3f} m over {pb - pa:.0f} m "
                        f"(+{b['h'] - a['h']:.3f} m) — misplaced point, wrong gauge zero or "
                        "reading?")
            if seq:
                kinds = {}
                for _p, a in seq:
                    kinds[a["kind"]] = kinds.get(a["kind"], 0) + 1
                self._log.append(
                    f"[info] path {network.paths[pi]['name']}: "
                    + ", ".join(f"{v} {k}" for k, v in kinds.items())
                    + f" value(s), p {seq[0][0] / 1000.0:.3f}–{seq[-1][0] / 1000.0:.3f} km")
            else:
                self._log.append(f"[info] path {network.paths[pi]['name']}: no "
                                 "water-surface value")

    # ------------------------------------------------------------------ #
    def has_anchors(self):
        return bool(self.gnss or self.gauges)

    def _on_path(self, pi, anchors):
        """[(p, anchor)] of the anchors lying on path `pi`, sorted by p."""
        out = []
        for a in anchors:
            for i, p in self.network.path_positions(a["river"], a["s"]):
                if i == pi:
                    out.append((float(p), a))
                    break
        out.sort(key=lambda t: t[0])
        return out

    def _bracketed(self, river, s):
        """True when GNSS points lie both upstream and downstream of (river, s)
        along some path — a gauge there is a QC point, not an anchor."""
        for pi, p in self.network.path_positions(river, s):
            on = self._on_path(pi, self.gnss)
            if any(q < p - 0.5 for q, _ in on) and any(q > p + 0.5 for q, _ in on):
                return True
        return False

    @staticmethod
    def _end_slope(seq, end):
        """Least-squares slope dh/dp [m/m] over the values nearest to one end of
        `seq` ([(p, anchor)] ascending; end=0 upstream, -1 downstream), taken
        until they span WS_SLOPE_WINDOW. A slope rising downstream is replaced
        by the whole-path fit, else by 0. Returns (slope, description)."""
        if len(seq) < 2:
            return 0.0, "flat: single value"
        order = range(len(seq)) if end == 0 else range(len(seq) - 1, -1, -1)
        pts = []
        for i in order:
            pts.append(seq[i])
            if len(pts) >= 2 and abs(pts[-1][0] - pts[0][0]) >= WS_SLOPE_WINDOW:
                break
        ps = np.asarray([q[0] for q in pts])
        hs = np.asarray([q[1]["h"] for q in pts])
        if np.ptp(ps) < 1.0:
            return 0.0, "flat: coincident values"
        slope = float(np.polyfit(ps, hs, 1)[0])
        why = f"{len(pts)} values over {np.ptp(ps):.0f} m"
        if slope > 0:
            ps2 = np.asarray([q[0] for q in seq])
            hs2 = np.asarray([q[1]["h"] for q in seq])
            s2 = float(np.polyfit(ps2, hs2, 1)[0]) if np.ptp(ps2) >= 1.0 else 0.0
            if s2 <= 0:
                return s2, f"whole path, {len(seq)} values (end window rose)"
            return 0.0, "flat: the surface rises downstream here"
        return slope, why

    def _junction(self, J, anchors):
        """Stage at junction J from every incoming branch (see class doc)."""
        net = self.network
        d, sj = J["river"], J["s"]
        branches = {}
        for pi, pj in net.path_positions(d, sj):
            segs = net.paths[pi]["segs"]
            k = next(i for i, sg in enumerate(segs) if sg["river"] == d)
            br = d if (segs[k]["s0"] < sj - NODE_SNAP_TOL or k == 0) else segs[k - 1]["river"]
            branches.setdefault(br, []).append((pi, pj))
        if not branches:
            return None
        pi0, pj0 = next(iter(branches.values()))[0]
        on0 = self._on_path(pi0, anchors)
        if any(abs(p - pj0) <= 0.5 and a["kind"] != "junction" for p, a in on0):
            return None                              # a measured value sits on the node
        down = [(p - pj0, p, a) for p, a in on0 if p > pj0 + 0.5]
        ests = []
        for br, plist in branches.items():
            best = None
            for pi, pj in plist:
                on = self._on_path(pi, anchors)
                ups = [(pj - p, p, a) for p, a in on if p < pj - 0.5]
                if ups:
                    du, pu, au = min(ups, key=lambda t: t[0])
                    if best is None or du < best[0]:
                        best = (du, pu, au, [(p, a) for p, a in on if p < pj - 0.5])
            if best is None:
                continue                             # no value on this branch
            du, pu, au, upseq = best
            if down:
                dd, _pd, ad = down[0]
                gap = du + dd
                v = au["h"] + (ad["h"] - au["h"]) * du / gap
                ests.append(dict(branch=br, value=float(v), kind="interp", gap=float(gap),
                                 frm=au["id"], to=ad["id"], deps=au["deps"] | ad["deps"]))
            else:
                slope, why = self._end_slope(upseq, -1)
                ests.append(dict(branch=br, value=float(au["h"] + slope * du),
                                 kind="extrap", gap=float(du), frm=au["id"], to="",
                                 deps=au["deps"]))
        interp = [e for e in ests if e["kind"] == "interp"]
        use = interp or ests
        if use:
            w = np.asarray([1.0 / max(e["gap"], 1.0) for e in use])
            h = float(np.sum(w * np.asarray([e["value"] for e in use])) / np.sum(w))
            deps = frozenset().union(*[e["deps"] for e in use])
            how = "interp" if interp else "extrap-down"
            spread = float(np.ptp([e["value"] for e in use])) if len(use) > 1 else 0.0
        elif down:
            dd, _pd, ad = down[0]
            slope, why = self._end_slope([(p, a) for _dd, p, a in down], 0)
            h = float(ad["h"] - slope * dd)
            deps, how, spread = ad["deps"], "extrap-up", 0.0
            ests.append(dict(branch=d, value=h, kind="extrap-up", gap=float(dd),
                             frm=ad["id"], to="", deps=deps))
        else:
            return None
        anchor = dict(kind="junction", id=f"J:{d}@{sj:.0f}", river=d, s=float(sj), h=h,
                      deps=deps)
        return dict(anchor=anchor, J=J, ests=ests, how=how, spread=spread, h=h)

    def _solve(self, when):
        """Anchors (with gauge stages at `when` for time-series readings) plus
        the junction anchors, per path. Cached per minute."""
        key = None
        if (when is not None and self.station_ws is not None
                and self.station_ws.mode == "series" and self.gauges):
            key = str(np.datetime64(when, "m"))
        if key in self._cache:
            return self._cache[key]
        anchors = list(self.gnss)
        for g in self.gauges:
            h = g["h"]
            if key is not None:
                st = self.station_ws._stage_at(g["id"], when)
                if st is not None:
                    h = g["zero"] + float(st)
            anchors.append(dict(g, h=float(h)))
        junctions = []
        for J in self.network.junctions:            # upstream junctions first
            jr = self._junction(J, anchors)
            if jr is not None:
                anchors.append(jr["anchor"])
                junctions.append(jr)
        per_path = [self._on_path(i, anchors) for i in range(len(self.network.paths))]
        sol = dict(anchors=anchors, per_path=per_path, junctions=junctions)
        self._cache[key] = sol
        return sol

    def _eval_on_path(self, seq, p):
        """Water surface at path chainage p from the path's values `seq`."""
        if not seq:
            return None
        ps = [q[0] for q in seq]
        if len(seq) == 1 or p < ps[0] - 0.5 or p > ps[-1] + 0.5:
            end = 0 if (len(seq) == 1 or p < ps[0]) else -1
            pe, ae = seq[end]
            slope, why = self._end_slope(seq, end)
            dist = p - pe
            tag = "EXTRAP" if abs(dist) > self.extrap_tol else "extrap"
            side = "upstream of" if dist < 0 else "downstream of"
            note = (f"{tag} {abs(dist):.0f} m {side} {ae['id']} at "
                    f"{slope * 1000.0:+.2f} m/km ({why})")
            return dict(value=float(ae["h"] + slope * dist), source="extrap", note=note,
                        extrap=abs(float(dist)), deps=ae["deps"])
        k = bisect.bisect_right(ps, p) - 1
        k = min(max(k, 0), len(seq) - 2)
        (pa, a), (pb, b) = seq[k], seq[k + 1]
        t = 0.0 if pb - pa < 1e-9 else min(max((p - pa) / (pb - pa), 0.0), 1.0)
        val = a["h"] + t * (b["h"] - a["h"])
        deps = a["deps"] | b["deps"]
        if deps == {"gnss"}:
            src, note = "survey", ""
        else:
            src = "station" if deps == {"gauge"} else "bridge"
            note = (f"{src} {a['id']} ({a['h']:.3f}) -> {b['id']} ({b['h']:.3f}) "
                    f"over {pb - pa:.0f} m")
        return dict(value=float(val), source=src, note=note, extrap=0.0, deps=deps)

    def resolve(self, river, km_internal, when=None):
        """Water surface at (river, km_internal). Returns dict(ws_elev, source,
        note, path, p, qc, extrap_m); ws_elev None (source='none') when no path
        through the river carries any value."""
        out = dict(ws_elev=None, source="none", note="", path="", p=None, qc=None,
                   extrap_m=0.0)
        if not river or km_internal is None:
            out["note"] = "transect not located on the network"
            return out
        pos = self.network.path_positions(river, km_internal)
        if not pos:
            out["note"] = f"river '{river}' lies on no network path"
            return out
        sol = self._solve(when)
        evals = []
        for pi, p in pos:
            e = self._eval_on_path(sol["per_path"][pi], p)
            if e is not None:
                evals.append((pi, p, e))
        if not evals:
            out["note"] = f"no GNSS point, gauge or junction value on any path through {river}"
            return out
        interp = [t for t in evals if t[2]["source"] != "extrap"]
        use = interp or evals
        vals = np.asarray([t[2]["value"] for t in use])
        pi, p, e = min(use, key=lambda t: t[2]["extrap"])
        note = e["note"]
        if len(use) > 1 and float(np.ptp(vals)) > 0.0005:
            note = ((note + "; ") if note else "") + \
                f"mean of {len(use)} paths (spread {float(np.ptp(vals)):.3f} m)"
        out.update(ws_elev=float(vals.mean()), source=e["source"], note=note,
                   path=self.network.paths[pi]["name"], p=float(p),
                   extrap_m=float(e["extrap"]))
        return out

    # -------------------------------------------------------------- plots
    def path_series(self, pi, p_values):
        """Water surface along path `pi` at chainages `p_values` (mean stages)."""
        seq = self._solve(None)["per_path"][pi]
        out = []
        for p in p_values:
            e = self._eval_on_path(seq, float(p))
            out.append(np.nan if e is None else e["value"])
        return np.asarray(out, dtype=float)

    def path_markers(self, pi):
        """{'gnss'|'gauge'|'junction'|'qc': [(p, elev, id)]} on path `pi`."""
        seq = self._solve(None)["per_path"][pi]
        out = dict(gnss=[], gauge=[], junction=[], qc=[])
        for p, a in seq:
            out[a["kind"]].append((p, a["h"], a["id"]))
        for g in self.qc_gauges:
            for i, p in self.network.path_positions(g["river"], g["s"]):
                if i == pi:
                    out["qc"].append((p, g["h"], g["id"]))
        return out

    def write_csv(self, path):
        """water_surface_qc.csv: one row per GNSS point, gauge and junction."""
        def f(v, fmt):
            return ("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                    else fmt.format(v))
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["kind", "id", "river", "km_internal_m", "km_oficial_m",
                        "dist_axis_m", "alt_river", "alt_dist_m", "elev_m", "role",
                        "note", "x", "y"])
            for r in self.rows:
                w.writerow([r["kind"], r["id"], r["river"], f(r["km_internal"], "{:.2f}"),
                            f(r["km_oficial"], "{:.2f}"), f(r["dist_axis"], "{:.1f}"),
                            r.get("alt_river", ""), f(r.get("alt_dist"), "{:.1f}"),
                            f(r["elev"], "{:.3f}"), r["role"], r["note"],
                            f(r["x"], "{:.3f}"), f(r["y"], "{:.3f}")])
        return Path(path)


def _ws_is_warn(source, note) -> bool:
    """Whether a transect's water-surface note is a warning (6.1): no surface,
    a legacy clamp, an extrapolation beyond --ws-extrap-tol (upper-case
    'EXTRAP'), or a failed QC. Interpolations (survey, bridge, station) and a
    short slope continuation ('extrap', lower case) are informative."""
    n = str(note or "")
    return (source in ("clamp", "none") or "EXTRAP" in n or "QC dif" in n
            or n.startswith("single"))


def representative_time(data):
    """Best-effort representative datetime for a transect from System.Time.
    Returns np.datetime64 or None.

    SonTek time units vary by export/config (MATLAB datenum, Unix s or ms), and
    the timezone may be UTC or local — so this only handles the common encodings
    heuristically and returns None when unsure, in which case station time-series
    matching falls back to the per-station series mean (see StationWS). UNTESTED
    against real .mat time: confirm the encoding/timezone per instrument before
    relying on time-series station matching."""
    t = data.get("time")
    if t is None:
        return None
    t = np.asarray(t, float)
    t = t[np.isfinite(t)]
    if t.size == 0:
        return None
    v = float(np.median(t))
    # 6.1: RiverSurveyor System.Time is seconds since 2000-01-01, verified in
    # sontek_time() against the file-name timestamp. 2026 is ~8.3e8 s, which fell
    # between the datenum and Unix branches below, so this returned None and
    # time-series gauges silently used their mean. Same epoch as sontek_time.
    if 1e8 < v < 2e9:
        try:
            return np.datetime64(SONTEK_EPOCH + datetime.timedelta(seconds=v))
        except (OverflowError, ValueError):
            return None
    if 6e5 < v < 8e5:                          # MATLAB datenum (days), year ~2000+
        try:
            dt = (datetime.datetime.fromordinal(int(v) - 366)
                  + datetime.timedelta(days=v % 1.0))
            return np.datetime64(dt)
        except Exception:
            return None
    if 9e8 < v < 4e9:                          # Unix epoch seconds
        return np.datetime64(int(round(v)), "s")
    if 9e11 < v < 4e12:                        # Unix epoch milliseconds
        return np.datetime64(int(round(v / 1000.0)), "s")
    return None


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
    primary: np.ndarray | None = None,
    src: np.ndarray | None = None,
    composite: bool = True,
    primary_min: int = 1,
    edge_anchor: str = "ref",
    vb_priority: bool | None = None,
    vb_min: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build the single-valued bed profile depth = f(s) on a regular grid.

    `bin_half_width` is the half-width [m] of the s-window used to sample the
    beam cloud at each grid node. When None it defaults to 2*dx, i.e. every node
    is informed by points within ±2 grid steps — wide enough to always capture
    neighbours at dx = 0.5 m, narrow enough not to over-smooth real bedforms.

    v6 — COMPOSITE DEPTHS (QRev model, generalised from v5's VB-priority).
    `primary` marks the points belonging to the selected primary depth reference
    (see build_beam_cloud). At each node:

        * where the window holds at least `primary_min` PRIMARY points, the bed
          is the weighted median of those points ALONE;
        * where it does not, and `composite` is on, the remaining points take
          over — this is the near-bank / edge zone the primary reference never
          reached, and is exactly QRev's composite fill;
        * with `composite` off the node is left empty and gap-filled later.

    This is v5's vb_priority rule with 'VB' replaced by 'the selected reference',
    so it works identically for depth_ref='vb' (the v5 default, byte-for-byte
    equivalent results) and for depth_ref='bt'. The point is to stop the more
    numerous secondary points from out-voting the reference and flattening the
    bank slope.

    `vb_priority` / `vb_min` are accepted as deprecated aliases so v5 call sites
    and configs keep working.

    Also returns `src_grid`: the QRev source code actually used at each node
    (1=BT, 2=VB, 3=DS, 4=interpolated), so every exported bed point is traceable
    to where it came from.
    """
    if bin_half_width is None:
        bin_half_width = 2 * dx

    # --- deprecated v5 aliases ---------------------------------------------
    if vb_priority is not None:
        composite = bool(vb_priority)
        if primary is None and beam_index is not None:
            primary = np.asarray(beam_index) == 0
    if vb_min is not None:
        primary_min = int(vb_min)

    _pri = None if primary is None else np.asarray(primary, dtype=bool)

    # Reference the bank ramps to the PRIMARY reference extent, not the full
    # cloud. The secondary points reach further toward the banks but are less
    # reliable on the bank slope, so they must NOT define the bank: the bed is
    # measured only across the reference span, and from the outermost reference
    # point on each side a straight ramp descends to depth 0 at the bank
    # (edge_left/right metres beyond it). With edge_anchor='cloud' the full
    # cloud sets the extent instead.
    ref_min, ref_max = s_data_min, s_data_max
    if edge_anchor == "ref" and _pri is not None:
        _ok = (np.isfinite(s_pts) & np.isfinite(depth_pts)
               & (depth_pts > 0) & _pri)
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

    src_grid = np.full_like(s_grid, SRC_INT, dtype=float)

    finite = np.isfinite(sp_shift) & np.isfinite(depth_pts) & (depth_pts > 0)
    sp = sp_shift[finite]; dp = depth_pts[finite]; wp = w_pts[finite]
    is_pri = None if _pri is None else _pri[finite]
    sp_src = None if src is None else np.asarray(src)[finite]

    def _mode_src(m):
        """Dominant source code among the points actually used at a node."""
        if sp_src is None or not np.any(m):
            return SRC_INT
        vals, cnts = np.unique(sp_src[m], return_counts=True)
        return float(vals[int(np.argmax(cnts))])

    # Bin-and-weighted-median in the data zone
    for i, sc in enumerate(s_grid):
        if sc < s_data_min_new - dx or sc > s_data_max_new + dx:
            continue
        mask = (sp >= sc - bin_half_width) & (sp <= sc + bin_half_width)
        if not np.any(mask):
            continue
        if is_pri is not None:
            pri_mask = mask & is_pri
            if int(np.count_nonzero(pri_mask)) >= primary_min:
                depth_grid[i] = weighted_median(dp[pri_mask], wp[pri_mask])
                src_grid[i] = _mode_src(pri_mask)
                continue
            if not composite:
                continue          # leave empty; gap-fill handles it below
        depth_grid[i] = weighted_median(dp[mask], wp[mask])
        src_grid[i] = _mode_src(mask)

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

    return s_grid, depth_grid, shift, src_grid


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

def _manual_river(args, network, perfil_id, log=None):
    """River forced for a transect by the INI [rios] section (6.1), matched
    accent- and case-insensitively against the centerline river names."""
    want = (getattr(args, "manual_rivers", None) or {}).get(perfil_id)
    if not want or network is None:
        return None
    river = {_norm_name(ax.river): ax.river for ax in network.axes}.get(_norm_name(want))
    if river is None and log is not None:
        log.append(f"[warn] [rios] {perfil_id} = {want}: that river is not in the "
                   "centerline — override ignored")
    return river


def process_one(matpath: Path, args, network=None, wsp=None, station_ws=None,
                datum_name: str = DATUM_NAME_DEFAULT,
                survey_outdir: Path | None = None, quiet: bool = False):
    """
    Full pipeline for ONE .mat transect. Returns a result dict (or None on
    failure). When `survey_outdir` is given, outputs go to <survey_outdir>/<perfil>/.

    v5: `network` is a RiverNetwork (multi-river). Each transect is tagged with
    its river and reported at that river's official km; a continuous internal
    survey chainage feeds the longitudinal profile and the water-surface lookup.
    `wsp` is the PRIMARY (GNSS) water surface; `station_ws` (optional) is the
    SECONDARY station-based surface used for fallback / gap-fill / QC.
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
        qinfo = transect_discharge(mat)
    except Exception as e:
        log.append(f"[error] could not read/extract {matpath.name}: {e}")
        try:
            (outdir / "process_log.txt").write_text("\n".join(log) + "\n",
                                                    encoding="utf-8")
        except Exception:
            pass
        print("\n".join(log) if not quiet else f"    {log[-1]}")
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

    # --------------------------------------------------- v6 depth reference
    depth_ref, bt_geometry, composite_on = resolve_depth_reference(data, args, log)

    # ------------------------------------------------- v6 depth spike filter
    data = filter_cloud_depths(data, enabled=(args.depth_filter != "off"), log=log)

    t_start = sontek_time(data)
    if t_start is not None:
        log.append(f"[info] inicio de transecta: {t_start:%Y-%m-%d %H:%M:%S} "
                   f"(System.Time, época 2000-01-01)")

    # ------------------------------------------------------------ STEP 2
    # v6: the section azimuth is the direction of the section's total unit-
    # discharge vector, restricted to the in-transect ensembles. See
    # compute_flow_direction() for why both corrections matter.
    theta_flow, theta_section = flow_direction_for(data, utm, args, log)

    # ------------------------------------------------------------ STEP 3
    cloud = build_beam_cloud(data, utm, depth_ref=depth_ref,
                             bt_geometry=bt_geometry, bt_avg=args.bt_avg,
                             density_weighting=(args.density_weighting != "off"))
    log.append(f"[info] beam cloud points: {len(cloud['E'])} "
               f"({int(cloud['primary'].sum())} primary [{depth_ref}])")

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
    # River, brazo and chainage from the network (if given). The report uses the
    # river's OFFICIAL km (progresiva_m); a continuous INTERNAL survey chainage
    # (across the surveyed rivers) feeds the longitudinal profile and WS lookup.
    river = ""
    brazo = ""
    role = ""
    progresiva_m = None            # official km of the transect's river (reported)
    km_internal = None             # per-river internal chainage
    survey_chainage = None         # 6.1: continuous chainage along its display path
    recorrido, path_rank = "", None
    loc_method, loc_note = "", ""
    dist_axis = float("nan")
    if network is not None:
        # 6.1: where the SECTION crosses a centerline, not the axis nearest to
        # the centroid (which put a Negro section below the confluence on the Limay)
        tr_c = Transformer.from_crs(utm_crs, network.crs, always_xy=True)
        tx, ty = tr_c.transform(utm[:, 0], utm[:, 1])
        force = _manual_river(args, network, perfil_id, log)
        loc = network.locate_track(tx, ty, force_river=force,
                                   method=getattr(args, "transect_locate", "crossing"))
        river = loc["river"]; role = loc["role"]; brazo = loc["brazo"]
        km_internal = loc["km_internal"]; progresiva_m = loc["km_oficial"]
        dist_axis = loc["dist"]; loc_method = loc["method"]; loc_note = loc["note"]
        path_rank, recorrido, survey_chainage = network.display_position(river, km_internal)
        log.append(
            f"[info] río={river} [{role}]"
            + (f" ({brazo})" if brazo else "")
            + f"   km oficial {progresiva_m:.2f} m ({format_progresiva(progresiva_m)})"
            + f"   located by {loc_method}, centroid-to-axis {dist_axis:.1f} m")
        if loc_note:
            log.append(("[warn] location: " if loc_method == "nearest"
                        else "[info] location: ") + loc_note)
        if survey_chainage is not None:
            log.append(f"[info] path {recorrido}: continuous chainage "
                       f"{survey_chainage:.2f} m")
    else:
        log.append("[info] chainage: not computed (no --centerline)")

    # ------------------------------------------------------------ STEP 4c
    # Water-surface elevation -> absolute datum. Precedence:
    #   1) GNSS water-surface profile where it covers the transect (ws_source=survey)
    #   2) stations, as gap-fill beyond the GNSS range or as sole source (=station)
    #   3) clamp to the nearest GNSS endpoint (=clamp) when neither of the above
    # If both GNSS and a station are available, |WS_gnss - WS_station| is checked
    # against the QC tolerance and flagged.
    when = representative_time(data)
    ws_elev = None
    ws_note = ""
    ws_source = ""
    x_ws = (network.survey_chainage(river, km_internal)       # legacy single chain
            if network is not None and km_internal is not None else None)
    qc_dif = None
    if wsp is not None and hasattr(wsp, "resolve"):
        # 6.1: WaterSurfaceModel — frames, GNSS/gauge bridge, per-river fallback
        res = wsp.resolve(river, km_internal, when=when)
        if res["ws_elev"] is not None:
            ws_elev, ws_note, ws_source = res["ws_elev"], res["note"], res["source"]
            if res.get("qc"):
                qc_dif = res["qc"]["dif"]
        else:
            ws_note = res["note"]
    else:
        # legacy path: a bare WaterSurfaceProfile (+ StationWS), as in v5/v6.0
        gnss_val, gnss_note, gnss_covers = None, "", False
        if wsp is not None and x_ws is not None:
            gnss_val, gnss_note = wsp.elev_at(x_ws)
            gnss_covers = wsp.covers(x_ws)
        stn_val, stn_note = None, ""
        if station_ws is not None and x_ws is not None:
            stn_val, stn_note = station_ws.elev_at(x_ws, river, when=when)

        if gnss_val is not None and gnss_covers:
            ws_elev, ws_note, ws_source = gnss_val, gnss_note, "survey"
        elif gnss_val is not None and stn_val is not None:
            ws_elev, ws_note, ws_source = stn_val, stn_note, "station"   # gap-fill
        elif gnss_val is not None:
            ws_elev, ws_note, ws_source = gnss_val, gnss_note, "clamp"
        elif stn_val is not None:
            ws_elev, ws_note, ws_source = stn_val, stn_note, "station"
        if gnss_val is not None and gnss_covers and stn_val is not None:
            qc_dif = abs(gnss_val - stn_val)
    if ws_elev is None and args.water_surface_elev:
        ws_elev, ws_note, ws_source = float(args.water_surface_elev), "", "constant"

    if ws_elev is not None:
        msg = (f"[info] water surface = {ws_elev:.3f} m {datum_name} "
               f"[source={ws_source}]  -> bed elevations ABSOLUTE ({datum_name})")
        if ws_note:
            msg += (f"   [WARN {ws_note}]" if _ws_is_warn(ws_source, ws_note)
                    else f"   [{ws_note}]")
        log.append(msg)
        # QC cross-check GNSS vs station where both exist
        if qc_dif is not None:
            tol = float(getattr(args, "ws_qc_tol", 0.10) or 0.10)
            tag = "PASS" if qc_dif <= tol else "WARN"
            log.append(f"[QA] WS GNSS vs estación: |Δ|={qc_dif:.3f} m (tol {tol:.3f}) {tag}")
            if qc_dif > tol:
                ws_note = (ws_note + "; " if ws_note else "") + f"QC dif={qc_dif:.3f} m"
    else:
        # 6.1: never silent — this transect's bed is RELATIVE (-depth)
        ws_source = "none"
        why = ws_note or ("no water-surface source configured" if wsp is None
                          and station_ws is None else "no water surface at this chainage")
        ws_note = f"SIN pelo de agua: {why}"
        log.append(f"[warn] NO water surface for this transect ({why}) — bed "
                   f"elevations are RELATIVE (-depth), not {datum_name}")

    # --------------------------------------------- v6 section-orientation QC
    orient_qc = section_orientation_qc(theta_section_final, utm,
                                       network=network, utm_crs=utm_crs,
                                       centroid=centroid, log=log, river=river)

    # ------------------------------------------------------------ STEP 5
    s_grid, depth_grid, s_shift, src_grid = build_profile(
        s_pts, cloud["depth"], eff_weight,
        s_data_min, s_data_max, dx=args.dx,
        bin_half_width=args.bin_half_width,
        edge_left_dist=edge_l, edge_right_dist=edge_r,
        beam_index=cloud["beam_index"],
        primary=cloud["primary"], src=cloud["src"],
        composite=composite_on, primary_min=args.vb_min,
        edge_anchor=args.edge_anchor,
    )
    bhw = args.bin_half_width if args.bin_half_width is not None else 2 * args.dx
    s_pts_final = s_pts + s_shift
    # Geographic origin of the profile (s = 0 = left bank), consistent with the
    # shift build_profile actually applied (reference-anchored banks).
    axis_origin = (Ec - s_shift * ux_LR, Nc - s_shift * uy_LR)
    mode = (f"ref={depth_ref}/{bt_geometry}"
            f", composite={'on' if composite_on else 'off'}"
            f", primary_min={args.vb_min}")
    log.append(f"[info] grid: {len(s_grid)} nodes ({s_grid[0]:.2f} → {s_grid[-1]:.2f} m, "
               f"dx={args.dx} m, bin_half_width={bhw} m, bed={mode})")
    if np.isfinite(qinfo["q_total"]):
        log.append(f"[info] caudal medido (RiverSurveyor) = {qinfo['q_total']:.3f} m³/s "
                   f"(margen izq {qinfo['q_left']:.3f}, der {qinfo['q_right']:.3f})")
    _n_int = int(np.count_nonzero(src_grid == SRC_INT))
    log.append(f"[info] bed source: {int(np.count_nonzero(src_grid == SRC_VB))} nodes VB, "
               f"{int(np.count_nonzero(src_grid == SRC_BT))} nodes BT, "
               f"{_n_int} nodes interpolated")

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
    if network is not None and np.isfinite(dist_axis):
        width = float(s_grid[-1] - s_grid[0])
        ok_axis = dist_axis <= max(50.0, 0.75 * width)
        log.append(f"[QA] transect centroid close to river axis? "
                   f"{'PASS' if ok_axis else 'WARN'} (dist={dist_axis:.1f} m, "
                   f"section width={width:.1f} m, located by {loc_method})")
    if ws_note and _ws_is_warn(ws_source, ws_note):
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

    # Raw beam cloud in POSGAR07 (for the merged survey point cloud). Z is the
    # ABSOLUTE bed elevation (ws_elev - depth) when a water surface is available,
    # otherwise a relative depth (-depth) flagged as such.
    tr_pg = Transformer.from_crs(utm_crs, posgar_crs, always_xy=True)
    cxp, cyp = tr_pg.transform(cloud["E"], cloud["N"])
    cloud_rel = ws_elev is None
    if ws_elev is None:
        cloud_z = -cloud["depth"]
    else:
        cloud_z = float(ws_elev) - cloud["depth"]

    return dict(
        perfil=perfil_id, outdir=outdir,
        river=river, role=role, brazo=brazo, dist_axis=dist_axis,
        km_internal=km_internal, km_oficial=progresiva_m,
        survey_chainage=survey_chainage, recorrido=recorrido, path_rank=path_rank,
        loc_method=loc_method, loc_note=loc_note,
        ws_elev=ws_elev, ws_note=ws_note, ws_source=ws_source,
        s_grid=s_grid, depth_grid=depth_grid, bed_elev=bed,
        xp=geo["xp"], yp=geo["yp"],
        width_m=float(s_grid[-1] - s_grid[0]),
        thalweg_elev=thalweg_elev, max_depth=max_depth,
        theta_flow=theta_flow, n_samples=int(len(data["vb_depth"])),
        utm_crs=utm_crs, posgar_crs=posgar_crs, log=log,
        # merged-cloud stash (one entry per (sample, beam) footprint)
        cloud_xp=np.asarray(cxp), cloud_yp=np.asarray(cyp), cloud_z=cloud_z,
        cloud_depth=cloud["depth"], cloud_beam_index=cloud["beam_index"],
        cloud_ens=cloud["ens"], cloud_rel=cloud_rel,
        # ---- v6: everything the aforo grouping engine needs ----
        theta_section=theta_section_final,
        t_start=t_start,
        prof_E=geo["E"], prof_N=geo["N"],
        q_total=qinfo["q_total"], q_left=qinfo["q_left"], q_right=qinfo["q_right"],
        group_chainage=progresiva_m,        # 6.1: same river is required anyway
        cloud=cloud, data_ref=data, utm_ref=utm, centroid=centroid,
        src_grid=src_grid, edge_l=edge_l, edge_r=edge_r,
        depth_ref=depth_ref, bt_geometry=bt_geometry, composite=composite_on,
        orient_qc=orient_qc, is_group=False, n_reps=1,
        members=[perfil_id],
    )


# ============================================================================ #
#  v6 — GROUP PROCESSING (averaged aforo profile)
# ============================================================================ #

def plot_group_qc(outdir: Path, gid: str, s_grid, depth_grid, stack, sigma,
                  group, ws_elev=None, datum_name: str = DATUM_NAME_DEFAULT,
                  qstat=None):
    """
    Repeatability plot of an aforo: every repetition over the merged profile,
    plus the +/-1 sigma band. This is the figure that shows whether the
    measurement is reproducible, and v5 could not produce it.
    """
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]})

    bed = (np.asarray(ws_elev) - depth_grid) if ws_elev is not None else -depth_grid
    for i, row in enumerate(stack):
        b = (ws_elev - row) if ws_elev is not None else -row
        ax1.plot(s_grid, b, lw=0.9, alpha=0.55,
                 label=f"rep {i+1}: {group[i]['perfil']}")
    ok = np.isfinite(sigma)
    if np.any(ok):
        ax1.fill_between(s_grid[ok], (bed - sigma)[ok], (bed + sigma)[ok],
                         color="0.5", alpha=0.30, lw=0, label="±1σ entre repeticiones")
    ax1.plot(s_grid, bed, "k-", lw=2.0, label="perfil promediado (nube fusionada)")
    if ws_elev is not None:
        ax1.axhline(ws_elev, color="tab:blue", lw=1.2, ls="--",
                    label=f"pelo de agua {ws_elev:.3f} m")
        ax1.set_ylabel(f"Cota [m {datum_name}]")
    else:
        ax1.set_ylabel("Cota relativa [m]")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8, loc="lower right", ncol=2)
    title = f"Aforo {gid} — {len(group)} repeticiones promediadas"
    if qstat and np.isfinite(qstat.get("q_mean", float("nan"))):
        title += (f"   |   Q = {qstat['q_mean']:.1f} m³/s"
                  f"  (CV {qstat['q_cov_pct']:.1f}%)")
    ax1.set_title(title)

    ax2.plot(s_grid, sigma, color="tab:red", lw=1.2)
    ax2.set_ylabel("σ [m]")
    ax2.set_xlabel("Progresiva transversal s [m]  (margen izquierda → derecha)")
    ax2.grid(alpha=0.3)
    if np.any(ok):
        ax2.axhline(float(np.nanmean(sigma)), color="0.4", ls=":", lw=1.0,
                    label=f"σ media = {np.nanmean(sigma):.3f} m")
        ax2.legend(fontsize=8)

    fig.tight_layout()
    path = outdir / "grupo_repetibilidad.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_group_qc_csv(groups_out, resumen_dir: Path):
    """One row per aforo group: repeatability statistics."""
    path = resumen_dir / "grupos_qc.csv"
    with open(path, "w", encoding="utf-8") as f:
        f.write("grupo,n_repeticiones,miembros,river,brazo,km_oficial_m,"
                "dispersion_progresiva_m,ancho_m,sigma_media_m,sigma_max_m,"
                "rms_por_repeticion_m,theta_section_deg,dispersion_theta_deg,"
                "sigma_media_comun_m,sigma_max_comun_m,"
                "q_media_m3s,q_desvio_m3s,q_cv_pct,q_min_m3s,q_max_m3s,"
                "q_por_repeticion_m3s\n")
        for g in groups_out:
            if not g.get("is_group"):
                continue
            rms = ";".join(f"{v:.4f}" if np.isfinite(v) else "nan"
                           for v in g.get("rep_rms", []))
            km = g.get("km_oficial")
            km_s = "" if km is None else f"{km:.2f}"
            miembros = " ".join(g.get("members", []))
            q = g.get("q_stats", {}) or {}
            qvals = ";".join(f"{v:.3f}" for v in q.get("q_values", []))
            f.write(f"{g['perfil']},{g['n_reps']},\"{miembros}\","
                    f"{g.get('river', '')},{g.get('brazo', '')},{km_s},"
                    f"{g.get('chainage_spread', float('nan')):.2f},"
                    f"{g.get('width_m', float('nan')):.2f},"
                    f"{g.get('sigma_mean', float('nan')):.4f},"
                    f"{g.get('sigma_max', float('nan')):.4f},"
                    f"\"{rms}\","
                    f"{math.degrees(g.get('theta_section', float('nan'))):.2f},"
                    f"{g.get('theta_spread_deg', float('nan')):.2f},"
                    f"{g.get('sigma_mean_common', float('nan')):.4f},"
                    f"{g.get('sigma_max_common', float('nan')):.4f},"
                    f"{q.get('q_mean', float('nan')):.3f},"
                    f"{q.get('q_std', float('nan')):.3f},"
                    f"{q.get('q_cov_pct', float('nan')):.2f},"
                    f"{q.get('q_min', float('nan')):.3f},"
                    f"{q.get('q_max', float('nan')):.3f},"
                    f"\"{qvals}\"\n")
    return path


def process_group(group, args, network=None, wsp=None, station_ws=None,
                  datum_name: str = DATUM_NAME_DEFAULT,
                  survey_outdir: Path | None = None):
    """
    Build ONE averaged profile from a group of repeated transects (an aforo).

    Strategy: merge the beam clouds of every repetition, project them onto a
    single group axis, and run build_profile ONCE. See merge_group_clouds() for
    why this beats averaging the finished profiles.

    Singleton groups are returned unchanged — no averaging, no extra folder.
    """
    if len(group) == 1:
        return group[0]

    gid = group_id(group)
    outdir = (survey_outdir / gid) if survey_outdir is not None else Path.cwd() / gid
    outdir.mkdir(parents=True, exist_ok=True)

    ref = group[0]
    log: list[str] = [f"[info] GRUPO {gid}: {len(group)} repeticiones promediadas",
                      f"[info] miembros: {', '.join(r['perfil'] for r in group)}",
                      f"[info] outdir : {outdir}"]

    # ---- group axis: pooled unit-discharge direction of all in-transect ens.
    theta_flow, theta_section = group_flow_direction(group, args, log)
    ux_LR, uy_LR = -math.cos(theta_section), -math.sin(theta_section)
    theta_section_final = math.atan2(uy_LR, ux_LR)

    # ---- merged cloud on one common centroid ------------------------------
    cloud = merge_group_clouds(group)
    centroid = (float(np.nanmean(cloud["E"])), float(np.nanmean(cloud["N"])))
    s_TL, d_off = project_to_axis(cloud["E"], cloud["N"], centroid, theta_section)
    s_pts = -s_TL
    log.append(f"[info] nube fusionada: {len(cloud['E'])} puntos "
               f"({int(cloud['primary'].sum())} de la referencia primaria)")

    # perpendicular-offset penalty, as for a single transect
    if args.no_offset_weighting or args.offset_scale <= 0:
        offset_penalty = np.ones_like(cloud["weight"])
    else:
        offset_penalty = np.exp(-0.5 * (np.abs(d_off) / args.offset_scale) ** 2)
    eff_weight = cloud["weight"] * offset_penalty

    # edge distances: median across repetitions (they alternate start edge but
    # Edges_0/1 are bound to left/right, so no swap is needed)
    edge_l = float(np.nanmedian([r["edge_l"] for r in group]))
    edge_r = float(np.nanmedian([r["edge_r"] for r in group]))
    log.append(f"[info] márgenes (mediana del grupo): izq={edge_l:.2f} m, der={edge_r:.2f} m")

    s_grid, depth_grid, s_shift, src_grid = build_profile(
        s_pts, cloud["depth"], eff_weight,
        float(np.nanmin(s_pts)), float(np.nanmax(s_pts)), dx=args.dx,
        bin_half_width=args.bin_half_width,
        edge_left_dist=edge_l, edge_right_dist=edge_r,
        beam_index=cloud["beam_index"], primary=cloud["primary"], src=cloud["src"],
        composite=ref.get("composite", True), primary_min=args.vb_min,
        edge_anchor=args.edge_anchor,
    )
    Ec, Nc = centroid
    axis_origin = (Ec - s_shift * ux_LR, Nc - s_shift * uy_LR)

    # ---- chainage / water surface: recomputed from the MERGED centroid -----
    river = ref.get("river", ""); brazo = ref.get("brazo", ""); role = ref.get("role", "")
    km_internal = ref.get("km_internal"); progresiva_m = ref.get("km_oficial")
    survey_chainage = ref.get("survey_chainage"); dist_axis = ref.get("dist_axis", float("nan"))
    recorrido, path_rank = ref.get("recorrido", ""), ref.get("path_rank")
    loc_method, loc_note = ref.get("loc_method", ""), ""
    utm_crs = ref["utm_crs"]; posgar_crs = ref["posgar_crs"]
    if network is not None:
        # 6.1: section crossing of the pooled tracks, on the members' river
        txs, tys = [], []
        for r in group:
            u = r.get("utm_ref")
            if u is None:
                continue
            trm = Transformer.from_crs(r["utm_crs"], network.crs, always_xy=True)
            a, b = trm.transform(np.asarray(u)[:, 0], np.asarray(u)[:, 1])
            txs.append(np.asarray(a)); tys.append(np.asarray(b))
        if txs:
            loc = network.locate_track(np.concatenate(txs), np.concatenate(tys),
                                       force_river=(river or None),
                                       method=getattr(args, "transect_locate", "crossing"))
        else:
            tr_c = Transformer.from_crs(utm_crs, network.crs, always_xy=True)
            loc = network.locate(*tr_c.transform(centroid[0], centroid[1]))
            loc.update(method="centroid", note="")
        river, role, brazo = loc["river"], loc["role"], loc["brazo"]
        km_internal, progresiva_m, dist_axis = loc["km_internal"], loc["km_oficial"], loc["dist"]
        loc_method = ref.get("loc_method") or loc["method"]
        loc_note = loc["note"] if loc["method"] == "nearest" else ""
        path_rank, recorrido, survey_chainage = network.display_position(river, km_internal)
        log.append(f"[info] río={river} [{role}]" + (f" ({brazo})" if brazo else "")
                   + f"   km oficial {progresiva_m:.2f} m "
                     f"({format_progresiva(progresiva_m)})  dist eje={dist_axis:.1f} m"
                   + (f"   recorrido {recorrido}" if recorrido else ""))

    # water surface: mean of the repetitions (same occupation, same stage)
    ws_vals = [r["ws_elev"] for r in group if r.get("ws_elev") is not None]
    ws_elev = float(np.mean(ws_vals)) if ws_vals else None
    ws_source = ref.get("ws_source", "")
    ws_note = "; ".join(sorted({r["ws_note"] for r in group if r.get("ws_note")}))
    if ws_elev is None:
        ws_source = "none"
        log.append(f"[warn] NO water surface for this group — bed elevations are "
                   f"RELATIVE (-depth), not {datum_name}")
    if ws_vals and len(ws_vals) > 1:
        spread = float(np.max(ws_vals) - np.min(ws_vals))
        log.append(f"[info] pelo de agua del grupo = {ws_elev:.3f} m {datum_name} "
                   f"(rango entre repeticiones {spread:.3f} m)")
        if spread > float(getattr(args, "ws_qc_tol", 0.10) or 0.10):
            log.append(f"[QA] WARN: el pelo de agua varió {spread:.3f} m dentro del "
                       f"grupo — ¿son realmente la misma ocupación?")

    # ---- repeatability ----------------------------------------------------
    stack, sigma, rep_rms, common = group_repeat_stats(
        group, s_grid, depth_grid,
        centroid=centroid, theta_section=theta_section, s_shift=s_shift)
    qstat = group_discharge_stats(group)
    if np.isfinite(qstat["q_mean"]):
        log.append(f"[info] caudal del aforo = {qstat['q_mean']:.3f} m³/s "
                   f"(n={len(qstat['q_values'])}, σ={qstat['q_std']:.3f}, "
                   f"CV={qstat['q_cov_pct']:.1f}%, rango "
                   f"{qstat['q_min']:.3f}–{qstat['q_max']:.3f})")
        log.append("[info] caudales por repetición: "
                   + ", ".join(f"{v:.3f}" for v in qstat["q_values"]) + " m³/s")
        if np.isfinite(qstat["q_cov_pct"]):
            tag = "PASS" if qstat["q_cov_pct"] <= 5.0 else "WARN"
            log.append(f"[QA] CV del caudal entre repeticiones: "
                       f"{qstat['q_cov_pct']:.1f}% (referencia USGS ≤5%) {tag}")
    sigma_mean = float(np.nanmean(sigma)) if np.any(np.isfinite(sigma)) else float("nan")
    sigma_max  = float(np.nanmax(sigma))  if np.any(np.isfinite(sigma)) else float("nan")
    if np.any(common):
        sc = sigma[common]
        sigma_mean_c = float(np.nanmean(sc)); sigma_max_c = float(np.nanmax(sc))
        s_lo, s_hi = float(s_grid[common][0]), float(s_grid[common][-1])
    else:
        sigma_mean_c = sigma_max_c = float("nan"); s_lo = s_hi = float("nan")
    log.append(f"[QA] repetibilidad del lecho medido (s={s_lo:.1f}–{s_hi:.1f} m, "
               f"donde las {len(group)} repeticiones se solapan): "
               f"σ media = {sigma_mean_c:.3f} m, σ máx = {sigma_max_c:.3f} m")
    log.append(f"[QA] repetibilidad incluyendo las rampas de margen: "
               f"σ media = {sigma_mean:.3f} m, σ máx = {sigma_max:.3f} m "
               f"(dominada por cuánto se acercó cada pasada al banco)")
    log.append(f"[QA] RMS de cada repetición contra el promedio: "
               + ", ".join(f"{v:.3f}" for v in rep_rms) + " m")

    # ---- exports ----------------------------------------------------------
    csv1, geo_arrays = export_csvs(
        outdir, gid, progresiva_m, s_grid, depth_grid,
        axis_origin, theta_section_final, utm_crs, posgar_crs,
        water_surface_elev=args.water_surface_elev, brazo=brazo, ws_elev=ws_elev)
    geo = dict(s_m=s_grid, depth_m=depth_grid, bed_elev=geo_arrays[6],
               E=geo_arrays[0], N=geo_arrays[1], lon=geo_arrays[2], lat=geo_arrays[3],
               xp=geo_arrays[4], yp=geo_arrays[5])
    log.append(f"[ok ] {csv1.name}")

    p2 = plot_cross_section(outdir, gid, progresiva_m,
                            s_pts + s_shift, cloud["depth"], cloud["beam_index"],
                            s_grid, depth_grid, brazo=brazo, ws_elev=ws_elev,
                            datum_name=datum_name)
    log.append(f"[ok ] {p2.name}")
    pqc = plot_group_qc(outdir, gid, s_grid, depth_grid, stack, sigma, group,
                        ws_elev=ws_elev, datum_name=datum_name, qstat=qstat)
    log.append(f"[ok ] {pqc.name}")

    if HAS_GPD:
        ax_path = export_axis_shp(outdir, gid, progresiva_m, axis_origin,
                                  theta_section_final,
                                  (float(s_grid[0]), float(s_grid[-1])), utm_crs,
                                  brazo=brazo)
        if ax_path: log.append(f"[ok ] {ax_path.name}")
        pp = export_profile_points_shp(outdir, gid, progresiva_m, geo,
                                       utm_crs, posgar_crs, brazo=brazo, ws_elev=ws_elev)
        if pp: log.append(f"[ok ] {pp.name}")

    acceptance_tests(theta_flow, theta_section_final, s_grid, depth_grid,
                     edge_left=edge_l, log=log)
    try:
        (outdir / "process_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    except Exception:
        pass

    bed = geo["bed_elev"]
    interior = depth_grid > 1e-6
    thalweg_elev = float(np.nanmin(bed[interior])) if interior.any() else float(np.nanmin(bed))

    ks = [r["group_chainage"] for r in group if r.get("group_chainage") is not None]
    devs = [_ang_diff_180(r["theta_section"], theta_section_final) for r in group]

    # POSGAR cloud for the merged survey point cloud
    tr_p = Transformer.from_crs(utm_crs, posgar_crs, always_xy=True)
    cxp, cyp = tr_p.transform(cloud["E"], cloud["N"])
    cloud_z = ((ws_elev - cloud["depth"]) if ws_elev is not None
               else -cloud["depth"])

    return dict(
        perfil=gid, outdir=outdir,
        river=river, role=role, brazo=brazo, dist_axis=dist_axis,
        km_internal=km_internal, km_oficial=progresiva_m,
        survey_chainage=survey_chainage, recorrido=recorrido, path_rank=path_rank,
        loc_method=loc_method, loc_note=loc_note,
        ws_elev=ws_elev, ws_note=ws_note, ws_source=ws_source,
        s_grid=s_grid, depth_grid=depth_grid, bed_elev=bed,
        xp=geo["xp"], yp=geo["yp"],
        width_m=float(s_grid[-1] - s_grid[0]),
        thalweg_elev=thalweg_elev, max_depth=float(np.nanmax(depth_grid)),
        theta_flow=theta_flow, theta_section=theta_section_final,
        n_samples=int(sum(r["n_samples"] for r in group)),
        utm_crs=utm_crs, posgar_crs=posgar_crs, log=log,
        cloud_xp=np.asarray(cxp), cloud_yp=np.asarray(cyp), cloud_z=cloud_z,
        cloud_depth=cloud["depth"], cloud_beam_index=cloud["beam_index"],
        cloud_ens=cloud["ens"], cloud_rel=(ws_elev is None),
        src_grid=src_grid, edge_l=edge_l, edge_r=edge_r,
        t_start=min((r["t_start"] for r in group if r.get("t_start")), default=None),
        group_chainage=progresiva_m,
        is_group=True, n_reps=len(group),
        members=[r["perfil"] for r in group],
        rep_stack=stack, rep_sigma=sigma, rep_rms=rep_rms,
        sigma_mean=sigma_mean, sigma_max=sigma_max,
        sigma_mean_common=sigma_mean_c, sigma_max_common=sigma_max_c,
        common_s=(s_lo, s_hi),
        prof_E=geo["E"], prof_N=geo["N"],
        q_total=qstat["q_mean"], q_stats=qstat,
        chainage_spread=(float(max(ks) - min(ks)) if ks else float("nan")),
        theta_spread_deg=(float(max(devs)) if devs else float("nan")),
        depth_ref=ref.get("depth_ref", ""), bt_geometry=ref.get("bt_geometry", ""),
        composite=ref.get("composite", True), orient_qc=ref.get("orient_qc", {}),
    )


# ============================================================================ #
#  SURVEY-LEVEL AGGREGATES (multi-file mode)
# ============================================================================ #

def _prof_x(r):
    """X of a transect on the longitudinal profile: the continuous chainage
    along its display path (6.1) when available, else the river's official km."""
    x = r.get("survey_chainage")
    if x is None:
        x = r.get("km_oficial")
    return x


def _sort_key(r):
    """Display order: by path (most transects first), then along it (6.1)."""
    x = _prof_x(r)
    rank = r.get("path_rank")
    return (float("inf") if rank is None else float(rank),
            float("inf") if x is None else float(x))


def write_survey_index(results, resumen_dir: Path, datum_name: str):
    """One row per transect: river, official km, brazo, cotas, ancho, WS source.
    6.1: written with the csv module (notes may contain commas); `recorrido`
    and `loc_method` appended at the end so existing column positions hold."""
    path = resumen_dir / "survey_index.csv"
    rs = sorted(results, key=_sort_key)

    def _f(v, fmt):
        return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) \
            else fmt.format(v)

    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["perfil", "river", "km_oficial_m", "km_oficial", "brazo",
                    "survey_chainage_m", "ws_elev_m", "ws_source", "thalweg_elev_m",
                    "max_depth_m", "width_m", "n_samples", "dist_axis_m",
                    "theta_flow_deg", "ws_note", "q_m3s", "recorrido", "loc_method"])
        for r in rs:
            kmo = r.get("km_oficial")
            absolute = r["ws_elev"] is not None
            w.writerow([
                r["perfil"], r.get("river", ""), _f(kmo, "{:.3f}"),
                "" if kmo is None else format_progresiva(kmo), r["brazo"],
                _f(r.get("survey_chainage"), "{:.3f}"),
                _f(r["ws_elev"], "{:.3f}"), r.get("ws_source", ""),
                _f(r["thalweg_elev"], "{:.3f}") if absolute else "",
                f"{r['max_depth']:.3f}", f"{r['width_m']:.2f}", r["n_samples"],
                _f(r["dist_axis"], "{:.1f}"), f"{math.degrees(r['theta_flow']):.1f}",
                r["ws_note"], _f(r.get("q_total"), "{:.3f}"),
                r.get("recorrido", ""), r.get("loc_method", "")])
    return path


def write_survey_profiles(results, resumen_dir: Path):
    """Every gridded profile point of every transect (long, tidy CSV)."""
    path = resumen_dir / "survey_profiles_all.csv"
    rs = sorted(results, key=_sort_key)
    with open(path, "w", encoding="utf-8") as f:
        f.write("perfil,river,km_oficial_m,brazo,s_m,depth_m,ws_elev_m,bed_elev_m,"
                "x_posgar07,y_posgar07\n")
        for r in rs:
            kmo = r.get("km_oficial")
            prog = "" if kmo is None else f"{kmo:.3f}"
            ws = "" if r["ws_elev"] is None else f"{r['ws_elev']:.3f}"
            sg, dg, be = r["s_grid"], r["depth_grid"], r["bed_elev"]
            xp, yp = r["xp"], r["yp"]
            for i in range(len(sg)):
                f.write(f"{r['perfil']},{r.get('river','')},{prog},{r['brazo']},"
                        f"{sg[i]:.4f},{dg[i]:.4f},{ws},{be[i]:.4f},"
                        f"{xp[i]:.4f},{yp[i]:.4f}\n")
    return path


_BRAZO_PALETTE = ["tab:blue", "tab:red", "tab:green", "tab:purple",
                  "tab:orange", "tab:brown", "tab:pink", "tab:olive"]


def _brazo_color(label):
    """Stable colour for an arbitrary brazo label (MI / MD / M2 ... or ''). The
    empty label (main channel, no island) is neutral grey."""
    if not label:
        return "0.35"
    return _BRAZO_PALETTE[(hash(label) % len(_BRAZO_PALETTE))]


def plot_survey_planview(results, network, resumen_dir: Path):
    """Centerlines (per river) + every section axis, coloured by brazo, labelled
    by the river's official km."""
    fig, ax = plt.subplots(figsize=(11, 12))
    # centerlines (network is in its own CRS = EPSG:5344, same as xp/yp)
    if network is not None:
        seenr = set()
        for axr in network.axes:
            xs, ys = axr.main.xy
            lbl = f"Eje {axr.river} (cauce principal)"
            ax.plot(xs, ys, "-", color="0.6", lw=1.4,
                    label=lbl if lbl not in seenr else None)
            seenr.add(lbl)
            for b in axr.anabranches:
                bx, by = b["geom"].xy
                ax.plot(bx, by, color=_brazo_color(b.get("label")),
                        lw=1.4, ls="--",
                        label=f"{axr.river} {b.get('label') or 'brazo'} (eje)")
    seen = set()
    for r in sorted(results, key=_sort_key):
        xp, yp = r["xp"], r["yp"]
        col = _brazo_color(r["brazo"])
        lbl = f"Secciones ({r.get('river','')}{' '+r['brazo'] if r['brazo'] else ''})"
        ax.plot([xp[0], xp[-1]], [yp[0], yp[-1]], "-", color=col, lw=2.2,
                label=lbl if lbl not in seen else None)
        seen.add(lbl)
        mx, my = 0.5 * (xp[0] + xp[-1]), 0.5 * (yp[0] + yp[-1])
        kmo = r.get("km_oficial")
        if kmo is not None:
            ax.annotate(format_progresiva(kmo).split("+")[0] + "k",
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


def _legend_unique(ax, **kw):
    h, l = ax.get_legend_handles_labels()
    seen, hh, ll = set(), [], []
    for hi, li in zip(h, l):
        if li and li not in seen:
            seen.add(li); hh.append(hi); ll.append(li)
    if hh:
        ax.legend(hh, ll, **kw)


def plot_long_profile(results, wsp, resumen_dir: Path, datum_name: str, network=None):
    """Longitudinal profile(s): water surface + thalweg.

    6.1: one panel per head-to-outlet PATH holding transects not drawn in an
    earlier panel (paths in display order), on that path's continuous chainage;
    the reach below a confluence is common to the paths through it. The water
    surface is the network model evaluated along the path, with its GNSS points,
    gauges (anchor / control) and junction stages. Without a network, the
    v6.0 single-axis plot."""
    net = network if network is not None else getattr(wsp, "network", None)
    located = [r for r in results
               if r.get("river") and r.get("km_internal") is not None]
    if net is None or not getattr(net, "paths", None) or not located:
        return _plot_long_profile_simple(results, wsp, resumen_dir, datum_name)
    panels, drawn = [], set()
    for pi in net.path_order:
        on = []
        for r in located:
            pos = dict(net.path_positions(r["river"], r["km_internal"]))
            if pi in pos:
                on.append((float(pos[pi]), r))
        if not on or all(r["perfil"] in drawn for _p, r in on):
            continue
        on.sort(key=lambda t: t[0])
        panels.append((pi, on))
        drawn.update(r["perfil"] for _p, r in on)
    if not panels:
        return None
    absolute = any(r["ws_elev"] is not None for r in located)
    model = wsp if hasattr(wsp, "path_series") else None
    fig, axes = plt.subplots(len(panels), 1, figsize=(12, 5.8 * len(panels)),
                             squeeze=False)
    for (pi, on), ax in zip(panels, axes[:, 0]):
        pa = net.paths[pi]
        nombre = pa["name"].replace(">", " → ")
        ps = np.asarray([p for p, _r in on])
        pad = max(250.0, 0.04 * float(ps.max() - ps.min()))
        x0, x1 = float(ps.min()) - pad, float(ps.max()) + pad
        if absolute and model is not None:
            pg = np.linspace(x0, x1, 800)
            ax.plot(pg / 1000.0, model.path_series(pi, pg), "-", color="navy", lw=1.5,
                    label="Pelo de agua (modelo de la red)", zorder=2)
            mk = model.path_markers(pi)

            def inr(t, lo=x0, hi=x1):
                return lo <= t[0] <= hi
            g = [t for t in mk["gnss"] if inr(t)]
            if g:
                ax.plot([t[0] / 1000.0 for t in g], [t[1] for t in g], ".",
                        color="tab:cyan", ms=6, label="Puntos GNSS de pelo de agua", zorder=4)
            for key, lab, fill in (("gauge", "Escala usada como ancla", "tab:red"),
                                   ("qc", "Escala de control (entre puntos GNSS)", "none")):
                lst = [t for t in mk[key] if inr(t)]
                if lst:
                    ax.plot([t[0] / 1000.0 for t in lst], [t[1] for t in lst], "s",
                            color="tab:red", mfc=fill, ms=7, mew=1.5, label=lab, zorder=5)
                    for t in lst:
                        ax.annotate(t[2], (t[0] / 1000.0, t[1]), fontsize=7,
                                    color="tab:red", xytext=(4, 4),
                                    textcoords="offset points")
            for t in [t for t in mk["junction"] if inr(t)]:
                ax.plot(t[0] / 1000.0, t[1], "D", color="k", ms=6,
                        label="Cota en la confluencia (estimada)", zorder=5)
        for sg in pa["segs"][1:]:
            if x0 <= sg["p0"] <= x1:
                ax.axvline(sg["p0"] / 1000.0, color="0.45", ls="--", lw=1.0, zorder=1)
                kmj = net.river_offsets.get(sg["river"], 0.0) + sg["s0"]
                ax.annotate(f"confluencia → río {sg['river']} (km {format_progresiva(kmj)})",
                            (sg["p0"] / 1000.0, 1.0),
                            xycoords=("data", "axes fraction"), rotation=90, fontsize=8,
                            color="0.3", ha="right", va="top", xytext=(-3, -4),
                            textcoords="offset points")
        xs = ps / 1000.0
        if absolute:
            for p, r in on:
                if r["ws_elev"] is not None:
                    ax.plot([p / 1000.0] * 2, [r["thalweg_elev"], r["ws_elev"]], "-",
                            color="0.7", lw=0.8, zorder=1)
            thal = [r["thalweg_elev"] if r["ws_elev"] is not None else np.nan
                    for _p, r in on]
            wss = [r["ws_elev"] if r["ws_elev"] is not None else np.nan for _p, r in on]
            ax.plot(xs, thal, "o", color="#7a5230", ms=6,
                    label="Thalweg (cota mínima del lecho)", zorder=3)
            ax.plot(xs, wss, "_", color="navy", ms=10, mew=2,
                    label="Pelo de agua en cada sección", zorder=3)
            n_rel = sum(1 for _p, r in on if r["ws_elev"] is None)
            if n_rel:
                ax.text(0.01, 0.02, f"{n_rel} sección(es) sin pelo de agua (cotas "
                                    "relativas) no graficadas", transform=ax.transAxes,
                        fontsize=8, color="tab:red")
            ax.set_ylabel(f"Cota {datum_name} [m]")
        else:
            thal = [-r["max_depth"] for _p, r in on]
            ax.plot(xs, thal, "o-", color="#7a5230", ms=6,
                    label="Thalweg (profundidad máx., relativa)")
            ax.axhline(0, color="navy", lw=1, ls="--", alpha=0.7,
                       label="Pelo de agua (rel.)")
            ax.set_ylabel("Cota relativa al pelo de agua [m]")
        for (p, r), th in zip(on, thal):
            if r["brazo"]:
                ax.annotate(r["brazo"], (p / 1000.0, thal_min(thal)), fontsize=7,
                            color=_brazo_color(r["brazo"]), ha="center", va="top",
                            xytext=(0, -2), textcoords="offset points")
        ax.set_xlim(x0 / 1000.0, x1 / 1000.0)
        ax.set_xlabel(f"Progresiva continua del recorrido {nombre} [km]")
        ax.set_title(f"Perfil longitudinal — recorrido {nombre}: pelo de agua y thalweg")
        ax.grid(True, alpha=0.3)
        _legend_unique(ax, loc="best", fontsize=8)
    fig.tight_layout()
    path = resumen_dir / "survey_long_profile.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _plot_long_profile_simple(results, wsp, resumen_dir: Path, datum_name: str):
    """v6.0 single-axis longitudinal profile (no network / legacy profile)."""
    rs = [r for r in results if _prof_x(r) is not None]
    rs.sort(key=_sort_key)
    if not rs:
        return None
    absolute = any(r["ws_elev"] is not None for r in rs)
    fig, ax = plt.subplots(figsize=(12, 6))
    if absolute and wsp is not None and getattr(wsp, "prog", None) is not None:
        ax.plot(wsp.prog, wsp.elev, "-", color="navy", lw=1.6,
                label="Perfil de pelo de agua (interpolado)")
        if getattr(wsp, "raw_prog", None) is not None:
            ax.plot(wsp.raw_prog, wsp.raw_elev, ".", color="tab:cyan", ms=5,
                    alpha=0.7, label="Puntos de pelo de agua medidos")
    prog = [_prof_x(r) for r in rs]
    if absolute:
        thal = [r["thalweg_elev"] for r in rs]
        wss = [r["ws_elev"] if r["ws_elev"] is not None else np.nan for r in rs]
        for r in rs:
            if r["ws_elev"] is not None:
                ax.plot([_prof_x(r), _prof_x(r)], [r["thalweg_elev"], r["ws_elev"]],
                        "-", color="0.7", lw=0.8, zorder=1)
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
    for r in rs:
        if r["brazo"]:
            ax.annotate(r["brazo"], (_prof_x(r), thal_min(thal)), fontsize=7,
                        color=_brazo_color(r["brazo"]), ha="center", va="top",
                        xytext=(0, -2), textcoords="offset points")
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
    rs = sorted(results, key=_sort_key)

    def _kmo(r):
        return np.nan if r.get("km_oficial") is None else float(r["km_oficial"])

    # axes
    axg, arows = [], []
    for r in rs:
        xp, yp = r["xp"], r["yp"]
        axg.append(LineString([(xp[0], yp[0]), (xp[-1], yp[-1])]))
        arows.append(dict(
            perfil=r["perfil"], river=r.get("river", ""),
            km_ofic_m=_kmo(r), brazo=r["brazo"], ws_src=r.get("ws_source", ""),
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
        kmo = _kmo(r)
        for i in range(len(sg)):
            pg.append(Point(xp[i], yp[i]))
            prows.append(dict(
                perfil=r["perfil"], river=r.get("river", ""),
                km_ofic_m=kmo, brazo=r["brazo"], s_m=float(sg[i]),
                depth_m=float(dg[i]),
                ws_elev_m=np.nan if r["ws_elev"] is None else float(r["ws_elev"]),
                bed_elev_m=float(be[i]),
            ))
    pp_path = resumen_dir / "survey_profile_points.shp"
    gpd.GeoDataFrame(prows, geometry=pg, crs=crs5344).to_file(pp_path)
    return ax_path, pp_path


# Fall back from .shp to .gpkg past ~this many features (shapefile 2 GB / field
# limits) or on any write error.
CLOUD_SHP_MAX = 1_500_000


def export_survey_beam_cloud(results, resumen_dir: Path):
    """Merged point cloud of EVERY bed return of EVERY transect, as PointZ in
    EPSG:5344 with Z = bed_elev (SRVN16). Written as a shapefile, falling back to
    GeoPackage past the shapefile limits. Attributes (km deliberately excluded):
        pt_id, transect, ens, beam_id, beam_type (slant|vert),
        bed_elev, depth, river, ws_elev, flag
    `flag` carries the water-surface source of the parent transect (survey /
    bridge / station / extrap / constant) plus 'REL' when Z is a relative depth
    (no WS)."""
    if not HAS_GPD:
        return None
    crs5344 = CRS.from_epsg(5344)
    rs = sorted(results, key=_sort_key)

    n_total = int(sum(len(r.get("cloud_z", [])) for r in rs))
    if n_total == 0:
        return None

    pid = np.empty(n_total, dtype=np.int64)
    transect = np.empty(n_total, dtype=object)
    ens = np.empty(n_total, dtype=np.int64)
    beam_id = np.empty(n_total, dtype=object)
    beam_type = np.empty(n_total, dtype=object)
    bed_elev = np.empty(n_total, dtype=float)
    depth = np.empty(n_total, dtype=float)
    river = np.empty(n_total, dtype=object)
    ws_elev = np.empty(n_total, dtype=float)
    flag = np.empty(n_total, dtype=object)
    geom = np.empty(n_total, dtype=object)

    k = 0
    for r in rs:
        z = r.get("cloud_z")
        if z is None or len(z) == 0:
            continue
        xp = r["cloud_xp"]; yp = r["cloud_yp"]
        dep = r["cloud_depth"]; bidx = r["cloud_beam_index"]; ce = r["cloud_ens"]
        wsv = r["ws_elev"]
        base_flag = r.get("ws_source", "") or ""
        if r.get("cloud_rel"):
            base_flag = (base_flag + ";REL") if base_flag else "REL"
        m = len(z)
        for i in range(m):
            pid[k] = k
            transect[k] = r["perfil"]
            ens[k] = int(ce[i])
            b = int(bidx[i])
            beam_id[k] = "VB" if b == 0 else str(b)
            beam_type[k] = "vert" if b == 0 else "slant"
            bed_elev[k] = float(z[i])
            depth[k] = float(dep[i])
            river[k] = r.get("river", "")
            ws_elev[k] = np.nan if wsv is None else float(wsv)
            flag[k] = base_flag
            geom[k] = Point(float(xp[i]), float(yp[i]), float(z[i]))
            k += 1

    gdf = gpd.GeoDataFrame({
        "pt_id": pid[:k], "transect": transect[:k], "ens": ens[:k],
        "beam_id": beam_id[:k], "beam_type": beam_type[:k],
        "bed_elev": bed_elev[:k], "depth": depth[:k], "river": river[:k],
        "ws_elev": ws_elev[:k], "flag": flag[:k],
    }, geometry=list(geom[:k]), crs=crs5344)

    use_gpkg = k > CLOUD_SHP_MAX
    if not use_gpkg:
        try:
            path = resumen_dir / "survey_raw_beam_points.shp"
            gdf.to_file(path)
            return path
        except Exception:
            use_gpkg = True
    path = resumen_dir / "survey_raw_beam_points.gpkg"
    gdf.to_file(path, driver="GPKG", layer="raw_beam_points")
    return path


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


def prescan_transects(matfiles, network, args=None, log: list | None = None):
    """Locate every transect on the network from its GPS alone (6.1).

    Loads only the GPS struct of each .mat (scipy variable_names), rebuilds the
    track exactly as process_one does (GPS.UTM when finite, else from Lat/Lon in
    the auto UTM zone) and locates it with the SAME rule (section crossing, INI
    [rios] override). Used BEFORE processing to rank the paths for display and
    to check that every transect will get a water surface. Returns
    [dict(perfil, river, km_internal, km_oficial, dist, method, note)]; files
    that fail here are reported by process_one."""
    out = []
    method = getattr(args, "transect_locate", "crossing") if args is not None else "crossing"
    for m in matfiles:
        try:
            mat = sio.loadmat(m, struct_as_record=False, squeeze_me=True,
                              variable_names=["GPS"])
            gps = mat.get("GPS")
            lat = np.asarray(get_field(gps, "Latitude"), dtype=float).ravel()
            lon = np.asarray(get_field(gps, "Longitude"), dtype=float).ravel()
            utm = get_field(gps, "UTM")
            utm = np.asarray(utm, dtype=float) if utm is not None else None
            utm_crs = auto_utm_crs(lat, lon)
            if utm is None or utm.ndim != 2 or not np.isfinite(utm).all():
                tr = Transformer.from_crs(CRS.from_epsg(4326), utm_crs, always_xy=True)
                ex, ny = tr.transform(lon, lat)
                utm = np.column_stack([ex, ny])
            tr_c = Transformer.from_crs(utm_crs, network.crs, always_xy=True)
            tx, ty = tr_c.transform(utm[:, 0], utm[:, 1])
            stem = Path(m).stem
            force = _manual_river(args, network, stem, log) if args is not None else None
            loc = network.locate_track(tx, ty, force_river=force, method=method)
            out.append(dict(perfil=stem, river=loc["river"],
                            km_internal=loc["km_internal"], km_oficial=loc["km_oficial"],
                            dist=loc["dist"], method=loc["method"], note=loc["note"]))
        except Exception as e:
            if log is not None:
                log.append(f"[warn] pre-scan: could not locate {Path(m).name} from its "
                           f"GPS ({e})")
    if log is not None and out:
        cnt = {}
        for t in out:
            cnt[t["river"]] = cnt.get(t["river"], 0) + 1
        log.append("[info] transects per river (GPS pre-scan, located by "
                   f"{method}): " + ", ".join(f"{k}: {v}" for k, v in cnt.items()))
        for t in out:
            if t["note"] or t["method"] not in ("crossing", "centroid"):
                tag = "[warn]" if t["method"] == "nearest" else "[info]"
                log.append(f"{tag} {t['perfil']}: {t['river']} km "
                           f"{format_progresiva(t['km_oficial'])} ({t['method']})"
                           + (f" — {t['note']}" if t["note"] else ""))
    return out


# ============================================================================ #
#  CONFIG FILE (campaña.ini)  +  TRACEABILITY RECORD
# ============================================================================ #

# 6.1: the INI key schema is derived from the CLI parser (_config_schema); these
# dests are not configuration values and are skipped there and in the record.
_CONFIG_SKIP_DESTS = {"help", "version", "config"}
# [procesamiento] keys that list the input .mat files / folders / globs.
_CONFIG_INPUT_KEYS = ("inputs", "matfiles", "input", "entradas")
# Config keys that are file/dir paths -> resolved relative to the config file.
_CONFIG_PATHS = {"outdir", "centerline", "water_surface_csv", "stations", "readings"}

# v6: keys of the [grupos] section are group names, values are space/comma
# separated transect ids (file stems) that must be averaged together. Use it
# when the automatic criteria get an occupation wrong.
_GROUP_SECTIONS = ["grupos", "grupo", "aforos", "groups"]


def _to_bool(s) -> bool:
    return str(s).strip().lower() in ("1", "true", "yes", "on", "si", "sí")


def _find_section(cp, names):
    low = {s.lower(): s for s in cp.sections()}
    for n in names:
        if n in low:
            return low[n]
    return None


def _config_schema(parser):
    """{dest: (kind, choices)} for every option of the CLI parser (6.1).

    Derived from the parser itself so a new CLI option is automatically valid in
    the INI. v6.0 kept a hand-written table that nobody updated for the v6 flags:
    depth-ref, composite, group-tol, ... were silently IGNORED in campaign INIs."""
    schema = {}
    for a in parser._actions:
        if not a.option_strings or a.dest in _CONFIG_SKIP_DESTS:
            continue
        if a.nargs == 0:                              # store_true switches
            kind = bool
        else:
            kind = a.type if callable(a.type) else str
        schema[a.dest] = (kind, tuple(a.choices) if a.choices else None)
    return schema


def load_config(path, parser=None, log: list | None = None):
    """
    Read an INI config file (campaña.ini). Returns (params, campania,
    river_offsets, manual_groups):
      params   -> dict of argparse dest -> typed value (for parser.set_defaults)
      campania -> dict of free-form metadata (for the traceability record)

    Section [procesamiento] holds parameters; [campania]/[campaña] holds metadata.
    Parameter keys may use '-' or '_'. File/dir paths are resolved relative to the
    config file's folder, so the config can live inside the campaign folder and
    use paths relative to it. Command-line options override anything set here.

    6.1: every CLI option is a valid key (schema from the parser); values of
    choice options are validated (argparse does not validate defaults); unknown
    keys are reported in `log` instead of being dropped silently. In-line ';'
    starts a comment in EVERY section (metadata included), as in v6.0.
    """
    cfgpath = Path(path)
    if not cfgpath.exists():
        sys.exit(f"[error] config file not found: {cfgpath}")
    cp = configparser.ConfigParser(interpolation=None,
                                   inline_comment_prefixes=(";",))
    cp.optionxform = str
    try:
        cp.read(cfgpath, encoding="utf-8")
    except Exception as e:
        sys.exit(f"[error] could not parse config file: {e}")
    if parser is None:
        parser = build_parser()
    schema = _config_schema(parser)
    log = log if log is not None else []

    base = cfgpath.resolve().parent
    params, inputs = {}, []

    proc = _find_section(cp, ["procesamiento", "proceso", "processing",
                              "parametros", "parámetros"])
    if proc:
        for raw_key, val in cp.items(proc):
            if val is None or str(val).strip() == "":
                continue
            key = raw_key.strip().lower().replace("-", "_")
            if key in _CONFIG_INPUT_KEYS:
                items = [s.strip() for line in str(val).splitlines()
                         for s in line.split(",") if s.strip()]
                for it in items:
                    p = Path(it)
                    inputs.append(str(p if p.is_absolute() else (base / p)))
                continue
            if key not in schema:
                log.append(f"[warn] config [{proc}]: unknown key '{raw_key}' — ignored "
                           "(typo? run with --help for the valid options)")
                continue
            typ, choices = schema[key]
            sval = str(val).strip()
            try:
                if typ is bool:
                    conv = _to_bool(sval)
                elif choices:
                    conv = sval.lower()
                    if conv not in [str(c).lower() for c in choices]:
                        sys.exit(f"[error] config: '{raw_key} = {sval}' — valid values: "
                                 f"{', '.join(str(c) for c in choices)}")
                else:
                    conv = typ(sval)
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
    campania = {k: (v or "").strip() for k, v in cp.items(camp)} if camp else {}

    # [progresivas]: per-river official chainage offsets {river: offset_m}.
    # km_oficial = river_offset + internal_chainage. Rivers not listed use 0
    # (or --chainage-offset as a global fallback). River names must match the
    # centerline `river` attribute (ASCII: Neuquen, Limay, Negro).
    river_offsets = {}
    prog = _find_section(cp, ["progresivas", "progresiva", "offsets",
                              "chainage", "km"])
    if prog:
        for raw_key, val in cp.items(prog):
            if val is None or str(val).strip() == "":
                continue
            try:
                river_offsets[str(raw_key).strip()] = float(str(val).strip())
            except (TypeError, ValueError):
                sys.exit(f"[error] config [progresivas]: bad offset for "
                         f"'{raw_key}': {val!r}")

    # v6: [grupos] — manual override of the aforo grouping.
    #   <group name> = <transect id> <transect id> ...
    manual_groups = {}
    gsec = _find_section(cp, _GROUP_SECTIONS)
    if gsec:
        for gname, val in cp.items(gsec):
            if val is None or str(val).strip() == "":
                continue
            ids = [t for t in re.split(r"[,\s]+", str(val).strip()) if t]
            if ids:
                manual_groups[str(gname).strip()] = ids

    known = [proc, camp, prog, gsec, _find_section(cp, _RIVER_SECTIONS)]
    for sec in cp.sections():
        if sec not in known:
            log.append(f"[warn] config: unknown section [{sec}] — ignored (valid: "
                       "[campanha], [procesamiento], [progresivas], [grupos], [rios])")

    return params, campania, river_offsets, manual_groups


# 6.1: [rios] — force the river of a transect (file stem = river name), for the
# rare section the crossing rule still gets wrong.
_RIVER_SECTIONS = ["rios", "ríos", "rio", "río", "asignacion_rio", "rivers"]


def load_river_overrides(path):
    """{transect id: river} from the INI [rios] section (6.1). Kept separate
    from load_config() so its return signature stays as in v6.0."""
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=(";",))
    cp.optionxform = str
    try:
        cp.read(Path(path), encoding="utf-8")
    except Exception:
        return {}
    sec = _find_section(cp, _RIVER_SECTIONS)
    if not sec:
        return {}
    return {str(k).strip(): str(v).strip() for k, v in cp.items(sec)
            if v is not None and str(v).strip()}


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


def write_run_record(out_dir, args, campania, matfiles, extra_lines=None,
                     river_offsets=None, setup_lines=None, warn_lines=None,
                     transect_warn_lines=None, log_name=None):
    """Write a human-readable traceability record of this run to procesamiento.txt.

    6.1: records EVERY effective parameter (from the parser, so new options are
    never missing), the network / paths / water-surface setup, every [warn] of
    the run and the QA warnings of each transect; `log_name` points to the full
    console transcript saved next to it."""
    lines = ["REGISTRO DE PROCESAMIENTO — process_adcp_bathimetric",
             f"script version : v{__version__}"]
    commit = _git_commit(sys.argv[0] if sys.argv and sys.argv[0] else __file__)
    if commit:
        lines.append(f"git commit     : {commit}")
    lines.append(f"fecha de corrida: "
                 f"{datetime.datetime.now().isoformat(timespec='seconds')}")
    if log_name:
        lines.append(f"consola completa: {log_name}")
    if campania:
        lines += ["", "[campaña]"]
        lines += [f"  {k}: {v}" for k, v in campania.items()]
    lines += ["", "[entradas]", f"  transectas ({len(matfiles)}):"]
    lines += [f"    - {Path(m).resolve()}" for m in matfiles]
    lines += [f"  centerline        : {args.centerline}",
              f"  water-surface CSV : {args.water_surface_csv}",
              f"  stations          : {getattr(args, 'stations', None)}",
              f"  readings          : {getattr(args, 'readings', None)}",
              "", "[parámetros] (valores efectivos: CLI > INI > default)"]
    seen = set()
    for a in build_parser()._actions:
        if not a.option_strings or a.dest in _CONFIG_SKIP_DESTS or a.dest in seen:
            continue
        seen.add(a.dest)
        lines.append(f"  {a.dest}: {getattr(args, a.dest, None)}")
    mg = getattr(args, "manual_groups", None)
    if mg:
        lines += ["", "[grupos] (forzados desde el INI)"]
        lines += [f"  {k}: {' '.join(v)}" for k, v in mg.items()]
    mr = getattr(args, "manual_rivers", None)
    if mr:
        lines += ["", "[rios] (río forzado por transecta desde el INI)"]
        lines += [f"  {k}: {v}" for k, v in mr.items()]
    if river_offsets:
        lines += ["", "[progresivas] (offsets oficiales por río, m)"]
        lines += [f"  {k}: {v}" for k, v in river_offsets.items()]
    if setup_lines:
        lines += ["", "[red, recorridos y pelo de agua]"]
        lines += [f"  {x}" for x in setup_lines]
    lines += ["", f"[advertencias de la corrida] ({len(warn_lines or [])})"]
    lines += [f"  {x}" for x in (warn_lines or [])] or ["  (ninguna)"]
    if transect_warn_lines is not None:
        lines += ["", f"[advertencias QA por transecta] ({len(transect_warn_lines)})"]
        lines += [f"  {x}" for x in transect_warn_lines] or ["  (ninguna)"]
    if extra_lines:
        lines += [""] + list(extra_lines)
    path = Path(out_dir) / "procesamiento.txt"
    try:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception:
        pass
    return path


def collect_transect_warnings(results):
    """Every WARN / [warn] / [error] line of every transect or group log (6.1),
    prefixed with its profile id, for the campaign-level record."""
    out = []
    for r in sorted(results, key=_sort_key):
        for line in r.get("log") or []:
            t = str(line).strip()
            if t.startswith("[warn]") or t.startswith("[error]") or " WARN" in t \
                    or t.startswith("[WARN"):
                out.append(f"{r['perfil']}: {t}")
    return out


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version",
                   version=f"%(prog)s {__version__}")
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
                   help="DEPRECATED v5 flag, kept so old command lines run: forces "
                        "composite=on and bt-geometry=footprints. Use --depth-ref / "
                        "--composite / --bt-geometry instead.")
    p.add_argument("--vb-min", type=int, default=1,
                   help="Min. PRIMARY-reference points in a node window to use the "
                        "reference alone (default 1; higher = fall back to the "
                        "secondary points sooner).")
    # --- v6: depth reference (QRev model) ---
    p.add_argument("--depth-ref", choices=["auto", "vb", "bt"], default="auto",
                   help="Primary depth reference: 'vb' vertical beam, 'bt' bottom-track "
                        "beams, 'auto' = whatever the operator selected in the field "
                        "(Setup.depthReference). Default auto.")
    p.add_argument("--composite", choices=["on", "off"], default="on",
                   help="QRev composite depths: fill nodes the primary reference did "
                        "not reach with the other source (default on).")
    p.add_argument("--bt-avg", choices=["idw", "simple"], default="idw",
                   help="How the 4 bottom-track beams collapse to one depth. 'idw' is "
                        "QRev's default (shallower beams weigh more); 'simple' is a "
                        "plain mean. Only used with --bt-geometry ensemble or "
                        "--depth-ref bt.")
    p.add_argument("--bt-geometry", choices=["footprints", "ensemble"],
                   default="footprints",
                   help="'footprints' keeps the 4 slant beams at their real ground "
                        "positions (better near-bank coverage, v5 behaviour); "
                        "'ensemble' collapses them to one averaged point at the boat, "
                        "as QRev does. Default footprints.")
    p.add_argument("--edge-anchor", choices=["ref", "cloud"], default="ref",
                   help="Anchor the bank ramps to the primary reference extent (default) "
                        "or to the whole cloud.")
    p.add_argument("--depth-filter", choices=["smooth", "off"], default="smooth",
                   help="Reject depth spikes per beam before building the cloud "
                        "(running IQR on residuals from a robust smooth, QRev-style). "
                        "Default smooth; v5 had no filter at all.")
    # --- v6: weighting and section azimuth ---
    p.add_argument("--density-weighting", choices=["on", "off"], default="on",
                   help="Weight each ensemble by the along-track interval it represents "
                        "so a stalled boat cannot pile up points and dominate a bin "
                        "(default on).")
    p.add_argument("--flow-weighting", choices=["discharge", "density", "none"],
                   default="discharge",
                   help="Weighting of the section azimuth: 'discharge' = interval x depth "
                        "x velocity (direction of the total unit-discharge vector, "
                        "default), 'density' = interval only, 'none' = v5 plain mean.")
    p.add_argument("--use-edge-ensembles", action="store_true",
                   help="Include the stationary bank ensembles (System.Step 2/4) in the "
                        "flow-direction estimate. Off by default, as in QRev.")
    p.add_argument("--section-orientation", choices=["flow", "centerline"],
                   default="flow",
                   help="Orient the section perpendicular to the mean flow (default, v5 "
                        "behaviour) or to the river centerline at that chainage. The "
                        "angle between the two is always reported as QA.")
    # --- v6: repeated-transect (aforo) grouping ---
    p.add_argument("--no-group-repeats", action="store_true",
                   help="Disable aforo grouping and report every repetition as its own "
                        "profile (v5 behaviour).")
    p.add_argument("--group-tol", type=float, default=15.0,
                   help="Max chainage difference [m] between repetitions of the same "
                        "aforo (default 15).")
    p.add_argument("--group-angle-tol", type=float, default=25.0,
                   help="Max angle [deg] between the section axes of two repetitions, "
                        "and QA tolerance for per-repetition azimuth spread (default 25).")
    p.add_argument("--group-max-span", type=float, default=0.0,
                   help="Optional cap [min] on the total duration of one occupation. "
                        "0 = off (default); time contiguity is the primary criterion.")
    # --- river chainage (progresiva) ---
    p.add_argument("--centerline", default=None,
                   help="River-network vector file (shp/gpkg/geojson). One feature "
                        "per river+role reach, attributes river/role[/seg_id/label/"
                        "flow_dir]. Islands and multi-river surveys are handled "
                        "automatically. A file with no river/role columns collapses "
                        "to v4 single-river-with-islands behaviour.")
    p.add_argument("--chainage-offset", type=float, default=0.0,
                   help="Global fallback official-km offset [m] for rivers not "
                        "listed in the INI [progresivas] section (per-river offsets "
                        "there take precedence). Default 0.")
    p.add_argument("--transect-locate", choices=["crossing", "centroid"],
                   default="crossing",
                   help="6.1: how a transect is placed on the network. 'crossing' "
                        "(default): river and km where the section line (principal "
                        "axis of the track) crosses a centerline, HEC-RAS style — "
                        "robust at confluences. 'centroid': nearest axis to the track "
                        "centroid (v6.0 rule, reproduces old chainages). The INI "
                        "[rios] section forces the river of a given transect.")
    p.add_argument("--chainage-reverse", action="store_true",
                   help="Measure chainage from the centerline's LAST vertex "
                        "(ignored where a water-surface CSV is given, since each "
                        "river is then oriented by its water-surface slope).")
    # --- water surface / absolute datum ---
    p.add_argument("--water-surface-csv", default=None,
                   help="CSV of surveyed GNSS water-surface points (X, Y, elevation) "
                        "— the PRIMARY water surface. Enables ABSOLUTE bed "
                        "elevations in the vertical datum.")
    p.add_argument("--ws-east-col", default="East", help="WS CSV easting column (default 'East').")
    p.add_argument("--ws-north-col", default="North", help="WS CSV northing column (default 'North').")
    p.add_argument("--ws-elev-col", default="H_correg", help="WS CSV elevation column (default 'H_correg').")
    p.add_argument("--ws-crs", default="EPSG:5344",
                   help="CRS of the WS CSV coordinates (default EPSG:5344 = POSGAR07 f2).")
    p.add_argument("--datum-name", default=DATUM_NAME_DEFAULT,
                   help=f"Vertical datum label for plots/columns (default '{DATUM_NAME_DEFAULT}').")
    # --- hydrometric stations (secondary water surface: fallback/gap-fill/QC) ---
    p.add_argument("--stations", default=None,
                   help="Station registry: CSV with X/Y columns in --ws-crs, or a "
                        "points vector file (station_id, name, river, gauge_zero "
                        "SRVN16, optional chainage). SECONDARY water surface: bridge "
                        "beyond the GNSS range, sole source without GNSS, and QC.")
    p.add_argument("--readings", default=None,
                   help="Per-campaign level readings CSV. Spot (station_id, nivel) "
                        "or time-series (station_id, datetime, nivel), auto-detected. "
                        "ws_elev = gauge_zero + nivel.")
    p.add_argument("--ws-qc-tol", type=float, default=0.10,
                   help="Tolerance [m] for the water-surface QC: gauges between GNSS "
                        "points vs the GNSS surface, and the spread of the tributaries' "
                        "estimates of a confluence stage (default 0.10).")
    p.add_argument("--ws-extrap-tol", type=float, default=WS_EXTRAP_TOL_DEFAULT,
                   help="6.1: beyond the first/last water-surface value of a network "
                        "path the surface is continued with its local slope; that "
                        "extrapolation is flagged as a WARNING past this distance [m] "
                        f"(default {WS_EXTRAP_TOL_DEFAULT:.0f}).")
    p.add_argument("--allow-relative", action="store_true",
                   help="6.1: when a water-surface source is configured but no "
                        "transect can get a water surface, continue with RELATIVE "
                        "bed elevations instead of stopping.")
    # --- perpendicular-offset penalty ---
    p.add_argument("--offset-scale", type=float, default=OFFSET_SCALE_DEFAULT,
                   help=f"Gaussian scale [m] for the perpendicular-offset penalty "
                        f"(default {OFFSET_SCALE_DEFAULT}; larger = gentler).")
    p.add_argument("--no-offset-weighting", action="store_true",
                   help="Disable the perpendicular-offset penalty entirely.")
    return p


def _line_label(r):
    """One-line description of a transect result for console/record."""
    kmo = r.get("km_oficial")
    prog = "n/a" if kmo is None else format_progresiva(kmo)
    tag = f" ({r['brazo']})" if r.get("brazo") else ""
    ws = "" if r["ws_elev"] is None else f", pelo={r['ws_elev']:.3f}m[{r.get('ws_source','')}]"
    th = "" if r["ws_elev"] is None else f", thalweg={r['thalweg_elev']:.3f}m"
    note = r.get("ws_note") or ""
    warn = ("" if not note else
            f"  [WARN {note}]" if _ws_is_warn(r.get("ws_source", ""), note) else f"  [{note}]")
    return (f"{r['perfil']}: {r.get('river','')} km {prog}{tag}"
            f", ancho={r['width_m']:.1f}m, prof.máx={r['max_depth']:.2f}m{ws}{th}{warn}")


class _Tee:
    """File-like object writing to a console stream AND a capture buffer."""

    def __init__(self, primary, buf):
        self._primary, self._buf = primary, buf

    def write(self, s):
        try:
            self._buf.write(s)
        except Exception:
            pass
        return self._primary.write(s)

    def flush(self):
        try:
            self._primary.flush()
        except Exception:
            pass

    def __getattr__(self, name):                 # encoding, isatty, fileno, ...
        return getattr(self._primary, name)


class RunLog:
    """Transcript of everything the run prints (6.1).

    stdout and stderr are tee'd into memory from the first line of main(), so
    the file also holds Python warnings (numpy, geopandas, ...), the reason of
    an early stop and, on a crash, the traceback. It is written as
    procesamiento.log in `target` (the _resumen folder, or the single-transect
    output folder) when the run ends — whether it ends well or not. Before 6.1
    the setup warnings only ever reached the console."""

    NAME = "procesamiento.log"

    def __init__(self):
        self._buf = io.StringIO()
        self._saved = None
        self.target = None

    def start(self):
        if self._saved is None:
            self._saved = (sys.stdout, sys.stderr)
            sys.stdout = _Tee(sys.stdout, self._buf)
            sys.stderr = _Tee(sys.stderr, self._buf)
        return self

    def stop(self):
        if self._saved is not None:
            for st in (sys.stdout, sys.stderr):
                try:
                    st.flush()
                except Exception:
                    pass
            sys.stdout, sys.stderr = self._saved
            self._saved = None

    def note(self, text):
        """Append to the transcript only (the console already shows it)."""
        self._buf.write(text if str(text).endswith("\n") else f"{text}\n")

    def text(self):
        return self._buf.getvalue()

    def warnings(self):
        """Run-level warnings: [warn]/[error] lines, failed [QA] checks and
        Python warnings, in the order they appeared."""
        out = []
        for line in self.text().splitlines():
            t = line.strip()
            if (t.startswith("[warn]") or t.startswith("[error]")
                    or (t.startswith("[QA]") and "WARN" in t)
                    or re.search(r"\b\w*Warning: ", t)):
                out.append(t)
        return out

    def path(self):
        return None if self.target is None else Path(self.target) / self.NAME

    def dump(self):
        p = self.path()
        if p is None:
            return None
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            head = (f"CONSOLA — process_adcp_bathimetric v{__version__}   "
                    f"{datetime.datetime.now().isoformat(timespec='seconds')}\n"
                    f"comando: {' '.join(sys.argv)}\n\n")
            p.write_text(head + self.text(), encoding="utf-8")
            return p
        except Exception:
            return None


def main():
    runlog = RunLog().start()
    try:
        _run(runlog)
    except SystemExit as e:
        if e.code not in (None, 0) and not isinstance(e.code, int):
            runlog.note(str(e.code))
        raise
    except KeyboardInterrupt:
        runlog.note("[error] interrupted by the user (Ctrl+C)")
        raise
    except BaseException as e:
        runlog.note(f"[error] unhandled exception: {e!r} — traceback follows")
        runlog.note(traceback.format_exc())
        raise
    finally:
        if runlog.path() is not None:
            print(f"[log] console transcript: {runlog.path()}")
        runlog.stop()
        runlog.dump()


def _run(runlog):
    # Two-phase parse: read --config first so it can supply defaults that the
    # command line then overrides (precedence: CLI > config file > built-in).
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=None)
    pre_args, _ = pre.parse_known_args()

    parser = build_parser()
    campania, river_offsets, manual_groups = {}, {}, {}
    config_log: list[str] = []
    if pre_args.config:
        cfg_params, campania, river_offsets, manual_groups = load_config(
            pre_args.config, parser=parser, log=config_log)
        parser.set_defaults(**cfg_params)
    args = parser.parse_args()
    args.manual_groups = manual_groups        # v6: [grupos] override
    args.manual_rivers = (load_river_overrides(pre_args.config)   # 6.1: [rios]
                          if pre_args.config else {})
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

    # 6.1: output locations are known from here on, so the console transcript
    # is saved even if the run stops during setup.
    survey_mode = len(matfiles) > 1
    if survey_mode:
        survey_name = args.survey_name or f"salida_{matfiles[0].stem}"
        survey_root = Path(args.outdir) if args.outdir else Path.cwd() / survey_name
        resumen = survey_root / "_resumen"
        runlog.target = resumen
    else:
        survey_root = resumen = None
        runlog.target = (Path(args.outdir) if args.outdir
                         else Path.cwd() / f"out_{matfiles[0].stem}")

    setup_log: list[str] = [f"[info] process_adcp_bathimetric v{__version__}"] + config_log
    ws_requested = bool(args.water_surface_csv) or bool(args.stations and args.readings)

    # --- GNSS water-surface points (parsed once; orientation + primary surface) --
    ws_table, ws_skipped = None, []
    if args.water_surface_csv:
        try:
            ws_table, ws_skipped, ws_cols = read_ws_table(
                args.water_surface_csv, args.ws_east_col, args.ws_north_col,
                args.ws_elev_col)
        except Exception as e:
            sys.exit(f"[error] could not read water-surface CSV: {e}")
        setup_log.append(
            f"[info] WS CSV: {len(ws_table)} point(s) read (E={ws_cols['east']}, "
            f"N={ws_cols['north']}, H={ws_cols['elev']}, id={ws_cols['id'] or 'row#'}, "
            f"CRS {args.ws_crs})"
            + (f"; {len(ws_skipped)} unreadable row(s)" if ws_skipped else ""))

    # --- river network (multi-river centerline) --------------------------- #
    network = None
    ws_pts_net = None
    if args.centerline:
        if not HAS_GPD:
            sys.exit("[error] --centerline needs geopandas/shapely installed")
        try:
            net_crs = gpd.read_file(args.centerline).crs
        except Exception as e:
            sys.exit(f"[error] could not read centerline: {e}")
        if ws_table is not None:
            if net_crs is not None and CRS.from_user_input(args.ws_crs) != net_crs:
                tr = Transformer.from_crs(args.ws_crs, net_crs, always_xy=True)
                ws_pts_net = []
                for p in ws_table:
                    x, y = tr.transform(p["x"], p["y"])
                    ws_pts_net.append(dict(p, x=float(x), y=float(y)))
            else:
                ws_pts_net = [dict(p) for p in ws_table]
        ws_xy_elev = ([(p["x"], p["y"], p["h"]) for p in ws_pts_net]
                      if ws_pts_net else None)
        try:
            network = RiverNetwork(
                args.centerline, crs_override=None, reverse=args.chainage_reverse,
                river_offsets=river_offsets, ws_xy_elev=ws_xy_elev, log=setup_log)
        except Exception as e:
            sys.exit(f"[error] could not build river network: {e}")
        # --chainage-offset is the global fallback for rivers with no INI offset
        for ax in network.axes:
            network.river_offsets.setdefault(ax.river, float(args.chainage_offset))

    # --- 6.1: transects located from their GPS (section crossing), paths ranked
    # for display only — the water surface needs no trunk.
    prescan = []
    if network is not None:
        if args.manual_rivers:
            setup_log.append("[info] [rios] overrides: " + ", ".join(
                f"{k} -> {v}" for k, v in args.manual_rivers.items()))
            stems = {m.stem for m in matfiles}
            for k in args.manual_rivers:
                if k not in stems:
                    setup_log.append(f"[warn] [rios] {k}: no transect with that id among "
                                     "the inputs — override unused")
        prescan = prescan_transects(matfiles, network, args=args, log=setup_log)
        network.rank_paths([(t["river"], t["km_internal"]) for t in prescan],
                           log=setup_log)

    # --- hydrometric stations (secondary water surface) -------------------- #
    station_ws = None
    if args.stations:
        if network is None:
            setup_log.append("[warn] --stations ignored: needs --centerline")
        elif not args.readings:
            setup_log.append("[warn] --stations given without --readings: no level "
                             "readings to build a water surface — stations ignored")
        else:
            try:
                registry = StationRegistry(args.stations, network,
                                           crs_override=args.ws_crs, log=setup_log)
                mode, rdata = read_level_readings(args.readings)
                station_ws = StationWS(registry, mode, rdata, log=setup_log)
            except Exception as e:
                setup_log.append(f"[warn] could not build station water surface: {e}")
                station_ws = None

    # --- water-surface model (GNSS + gauges on the network paths) ---------- #
    wsm = None
    if network is not None and (ws_pts_net or station_ws is not None):
        wsm = WaterSurfaceModel(network, points=ws_pts_net, station_ws=station_ws,
                                qc_tol=args.ws_qc_tol, skipped=ws_skipped,
                                log=setup_log, extrap_tol=args.ws_extrap_tol)
    elif ws_table is not None and network is None:
        setup_log.append("[warn] --water-surface-csv ignored: needs --centerline")

    # --- 6.1: coverage check BEFORE processing (no silent relative output) - #
    no_ws = []
    for t in prescan:
        res = (wsm.resolve(t["river"], t["km_internal"]) if wsm is not None
               else dict(ws_elev=None, note="no water-surface model"))
        if res["ws_elev"] is None:
            no_ws.append(f"{t['perfil']} ({t['river']} km "
                         f"{format_progresiva(t['km_oficial'])}: {res['note']})")
    fatal = None
    if ws_requested and not args.water_surface_elev:
        if wsm is None or not wsm.has_anchors():
            fatal = "no usable GNSS point or gauge reading"
        elif prescan and len(no_ws) == len(prescan):
            fatal = (f"none of the {len(prescan)} transects lies on a network path "
                     "carrying a water-surface value")
    if no_ws and not fatal:
        setup_log.append(f"[warn] {len(no_ws)} of {len(prescan)} transect(s) will get NO "
                         "water surface (RELATIVE bed):")
        setup_log += [f"[warn]   - {x}" for x in no_ws]

    if setup_log:
        print("\n".join(setup_log))
        print()
    if wsm is not None:
        try:
            Path(runlog.target).mkdir(parents=True, exist_ok=True)
            qcp = wsm.write_csv(Path(runlog.target) / "water_surface_qc.csv")
            print(f"[ok ] {qcp}")
        except Exception as e:
            print(f"[warn] could not write water_surface_qc.csv: {e}")
    if fatal:
        msg = (f"[error] a water-surface source is configured but {fatal}: every bed "
               f"elevation would be RELATIVE, not {datum_name}. See the [warn] lines "
               "above" + (" and water_surface_qc.csv" if wsm is not None else "")
               + "; fix the inputs, or set 'allow-relative = true' (--allow-relative) "
               "to continue anyway.")
        if args.allow_relative:
            print(msg.replace("[error]", "[warn]", 1) + "  -> continuing (allow-relative)")
        else:
            sys.exit(msg)

    # =============================== single transect ======================= #
    if not survey_mode:
        r = process_one(matfiles[0], args, network=network, wsp=wsm,
                        station_ws=station_ws, datum_name=datum_name)
        if r is not None:
            rec = write_run_record(r["outdir"], args, campania, matfiles,
                                   river_offsets=river_offsets, setup_lines=setup_log,
                                   warn_lines=runlog.warnings(),
                                   transect_warn_lines=collect_transect_warnings([r]),
                                   log_name=RunLog.NAME)
            print(f"[ok ] {rec.name}")
        return

    # =============================== survey mode =========================== #
    survey_root.mkdir(parents=True, exist_ok=True)
    resumen.mkdir(parents=True, exist_ok=True)

    print(f"[survey] {len(matfiles)} transects -> {survey_root}")
    results = []
    for i, m in enumerate(matfiles, 1):
        r = process_one(m, args, network=network, wsp=wsm, station_ws=station_ws,
                        datum_name=datum_name, survey_outdir=survey_root, quiet=True)
        if r is None:
            print(f"  [{i}/{len(matfiles)}] {m.name}: FAILED "
                  f"(see {Path(m.stem) / 'process_log.txt'})")
            continue
        print(f"  [{i}/{len(matfiles)}] {_line_label(r)}")
        results.append(r)

    if not results:
        sys.exit("[error] no transect processed successfully")

    # ------------------------------------------------------------------ v6
    # AFORO GROUPING. Repetitions of one occupation collapse into a single
    # averaged profile. The per-transect folders written above are kept
    # untouched for traceability, but from here on every aggregate — index,
    # longitudinal profile, plan view, shapefiles, merged point cloud — sees
    # ONLY the group profiles, so a 4-repetition aforo appears once.
    all_transects = list(results)
    group_log: list[str] = []
    groups = group_transects(results, args, log=group_log)
    for line in group_log:
        print(line)
    aggregated = []
    for g in groups:
        if len(g) == 1:
            aggregated.append(g[0])
            continue
        gr = process_group(g, args, network=network, wsp=wsm,
                           station_ws=station_ws, datum_name=datum_name,
                           survey_outdir=survey_root)
        print(f"  [grupo] {gr['perfil']}: {gr['n_reps']} repeticiones -> "
              f"ancho={gr['width_m']:.1f} m, prof.máx={gr['max_depth']:.2f} m, "
              f"σ lecho={gr['sigma_mean_common']:.3f} m"
              + ("" if not np.isfinite((gr.get('q_stats') or {}).get('q_mean', float('nan')))
                 else f", Q={gr['q_stats']['q_mean']:.1f} m³/s "
                      f"(CV {gr['q_stats']['q_cov_pct']:.1f}%)"))
        aggregated.append(gr)

    n_grp = sum(1 for r in aggregated if r.get("is_group"))
    if n_grp:
        qc = write_group_qc_csv(aggregated, resumen)
        print(f"[ok ] {qc.relative_to(survey_root)}")
    results = aggregated

    # 6.1: water-surface source tally — the first thing to check in a campaign
    tally = {}
    for r in results:
        k = r.get("ws_source") or "none"
        tally[k] = tally.get(k, 0) + 1
    ws_line = ("[info] pelo de agua por perfil: "
               + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(ws_line)
    if tally.get("none"):
        print(f"[warn] {tally['none']} perfil(es) SIN pelo de agua: cotas RELATIVAS")

    # aggregates
    idx = write_survey_index(results, resumen, datum_name)
    prof = write_survey_profiles(results, resumen)
    print(f"[ok ] {idx.relative_to(survey_root)}")
    print(f"[ok ] {prof.relative_to(survey_root)}")
    pv = plot_survey_planview(results, network, resumen)
    print(f"[ok ] {pv.relative_to(survey_root)}")
    lp = plot_long_profile(results, wsm, resumen, datum_name, network=network)
    if lp:
        print(f"[ok ] {lp.relative_to(survey_root)}")
    if HAS_GPD:
        sax, spp = export_survey_shapefiles(results, resumen, datum_name)
        if sax:
            print(f"[ok ] {sax.relative_to(survey_root)}")
        if spp:
            print(f"[ok ] {spp.relative_to(survey_root)}")
        cloud = export_survey_beam_cloud(results, resumen)
        if cloud:
            print(f"[ok ] {cloud.relative_to(survey_root)}")

    # traceability record for the whole survey
    summary = [ws_line, "", f"[resumen] {len(results)} transectas procesadas:"]
    for r in sorted(results, key=_sort_key):
        summary.append("    " + _line_label(r))
    twarn = collect_transect_warnings(
        all_transects + [r for r in results if r.get("is_group")])
    rec = write_run_record(resumen, args, campania, matfiles,
                           extra_lines=summary, river_offsets=river_offsets,
                           setup_lines=setup_log, warn_lines=runlog.warnings(),
                           transect_warn_lines=twarn, log_name=RunLog.NAME)
    print(f"[ok ] {rec.relative_to(survey_root)}")

    print(f"\n[done] survey outputs in: {survey_root}")
    print(f"       per-transect folders + consolidated results in: {resumen}")


if __name__ == "__main__":
    main()
