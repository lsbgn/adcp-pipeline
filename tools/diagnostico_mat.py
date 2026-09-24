#!/usr/bin/env python3
"""
diagnostico_mat.py  (v2.1, pipeline 6.5)
================================================================================
Inspect RiverSurveyor Live .mat exports (SonTek M9/S5) before trusting the
beam geometry used by process_adcp_bathimetric.py.

Per file:
  [1] every Setup field — edge shape (Edges_0/1__Method), depth reference,
      start edge and draft flagged;
  [2] pitch / roll distribution and the footprint shift the tilt implies at
      the median depth;
  [3] beam NUMBERING, two independent ways:
        - the file's Transformation_Matrices: RSL defines XYZ right-handed with
          +X toward beam 1 and +Z up, so beam 2 at +90 deg of beam 1 means beam
          2 to PORT ('ccw');
        - the data: slant depths compared with the vertical-beam profile at the
          footprint each of 8 hypotheses predicts (beam 1 at 0/90/180/270 deg,
          clockwise or counter-clockwise);
  [4] beam directions per frequency from the matrix (pointing down, so the
      old 180-degree ambiguity is gone).

Campaign level (every file given on the command line):
  [5] pitch / roll SIGN and GNSS antenna height, from repeated passes over the
      same section in OPPOSITE directions (aforos). A tilt moves every
      footprint by (depth + antenna height) * tan(tilt) along the boat's bow;
      when the boat reverses, that shift reverses along the section, so the
      right sign is the one that makes opposite passes agree. A GNSS time lag
      (shift = speed * lag) is fitted jointly, because bow-up pitch grows with
      speed and would otherwise be confused with it.
      v2.1: the POWER of the test is estimated first (bed slope x footprint
      shift against the same-direction noise): a flat aforo section, or little
      pitch, cannot decide anything. Pitch and roll signs are decided
      separately (roll-sign none when only the pitch is clear), and a
      per-occupation table is printed.

Not tested any more: whether RSL compensated the depths for tilt. The RSL
manual states it does, for the vertical and the slant beams ("including the
Transducer Depth and compensation for tilt", BottomTrack structure), and QRev
relies on it. The v1 regression that said otherwise was confounded by the
mirrored beam numbering.

Usage (PowerShell, venv active):
    python tools\\diagnostico_mat.py <carpeta de .mat | archivos> [...] > diag_mat.txt
================================================================================
"""
from __future__ import annotations

import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import scipy.io as sio

TILT_DEG = 25.0
FREQ_SPLIT = 2500.0
FLAGGED = ("depthReference", "startEdge", "sensorDepth", "coordinateSystem",
           "magneticDeclination", "headingSource", "offsetX", "offsetY")

# attitude test
OCC_DIST = 20.0          # [m] centroids of two passes of the same occupation
OCC_ANGLE = 25.0         # [deg] their track axes
OCC_GAP_MIN = 30.0       # [min] max time between consecutive passes
ANT_GRID = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5)   # [m] antenna height candidates
LAG_GRID = np.round(np.arange(-2.0, 2.0001, 0.1), 2)   # [s] GNSS lag candidates
PITCH_HYP = (("bow-up", 1.0), ("bow-down", -1.0))
ROLL_HYP = (("stbd-down", 1.0), ("port-down", -1.0))


# --------------------------------------------------------------------------- #
#  reading
# --------------------------------------------------------------------------- #

def get(st, name):
    if st is None or name not in getattr(st, "_fieldnames", []):
        return None
    return getattr(st, name)


def per_sample(v, n):
    """One value per sample (RSL lists Compass.Pitch/Roll as NS x 3)."""
    if v is None:
        return None
    a = np.asarray(v, dtype=float)
    if a.ndim == 2:
        a = a[:, 0] if a.shape[0] == n else a[0, :] if a.shape[1] == n else a.ravel()
    a = a.ravel()
    return a if a.size == n else None


def read(path: Path) -> dict:
    mat = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
    setp, comp, sysd, bt, gps = (mat.get(k) for k in
                                 ("Setup", "Compass", "System", "BottomTrack", "GPS"))
    vb = get(bt, "VB_Depth")
    vb = None if vb is None else np.asarray(vb, dtype=float).ravel()
    n = 0 if vb is None else vb.size
    if vb is not None:
        vb = np.where(vb > 0, vb, np.nan)
    B = get(bt, "BT_Beam_Depth")
    if B is not None:
        B = np.asarray(B, dtype=float)
        if B.ndim == 2 and B.shape[1] != 4 and B.shape[0] == 4:
            B = B.T
        B = np.where(np.isfinite(B) & (B > 0), B, np.nan) if B.shape == (n, 4) else None
    pitch = per_sample(get(comp, "Pitch"), n)
    roll = per_sample(get(comp, "Roll"), n)
    src = "Compass"
    if pitch is None or roll is None:
        pitch, roll, src = per_sample(get(sysd, "Pitch"), n), per_sample(get(sysd, "Roll"), n), "System"
    hdg = get(sysd, "True_North_ADP_Heading")
    if hdg is None:
        hdg = get(sysd, "Heading")
    hdg = per_sample(hdg, n)
    freq = per_sample(get(bt, "BT_Frequency"), n)
    t = per_sample(get(sysd, "Time"), n)
    U = get(gps, "UTM")
    E = N = None
    if U is not None:
        U = np.asarray(U, dtype=float)
        if U.ndim == 2 and U.shape[0] == n and U.shape[1] >= 2:
            E, N = U[:, 0].copy(), U[:, 1].copy()
            bad = ~np.isfinite(E) | ~np.isfinite(N) | ((E == 0) & (N == 0))
            E[bad] = np.nan
            N[bad] = np.nan
    return dict(mat=mat, setp=setp, n=n, vb=vb, B=B, pitch=pitch, roll=roll,
                att_src=src, hdg=hdg, freq=freq, t=t, E=E, N=N,
                tm=mat.get("Transformation_Matrices"))


# --------------------------------------------------------------------------- #
#  geometry
# --------------------------------------------------------------------------- #

def dist(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return "sin datos"
    p1, p5, p50, p95, p99 = np.percentile(x, [1, 5, 50, 95, 99])
    a = np.abs(x)
    return (f"n={x.size:5d}  p1={p1:+6.2f} p5={p5:+6.2f} med={p50:+6.2f} "
            f"p95={p95:+6.2f} p99={p99:+6.2f}   |x|>2: {100 * np.mean(a > 2):5.1f}%"
            f"  >3: {100 * np.mean(a > 3):5.1f}%  >5: {100 * np.mean(a > 5):5.1f}%")


def matrices(tm):
    """[(freq, 4x3 beam directions in XYZ with Z up, pointing DOWN)]."""
    M, F = get(tm, "Matrix"), get(tm, "Frequency")
    if M is None:
        return []
    M = np.asarray(M, dtype=float)
    if M.ndim == 2:
        M = M[:, :, None]
    elif M.ndim == 3 and M.shape[0] != 4 and M.shape[-1] == 4:
        M = np.transpose(M, (1, 2, 0))
    F = np.asarray(F, dtype=float).ravel() if F is not None else np.array([])
    out = []
    for k in range(M.shape[2]):
        T = M[:, :, k]
        if T.shape != (4, 4) or not np.all(np.isfinite(T)) or abs(np.linalg.det(T)) < 1e-9:
            continue
        D = np.linalg.inv(T)[:, :3]
        nrm = np.linalg.norm(D, axis=1, keepdims=True)
        if np.any(nrm < 1e-9):
            continue
        D = D / nrm
        D = np.where(D[:, 2:3] > 0, -D, D)          # beams point down (Z up)
        incl = np.degrees(np.arccos(np.clip(np.abs(D[:, 2]), 0, 1)))
        if not np.all((incl > 15) & (incl < 35)):   # the vertical beam's matrix
            continue
        out.append((F[k] if k < F.size else float("nan"), D))
    return out


def layout_of(D):
    z = D[0, 0] * D[1, 1] - D[0, 1] * D[1, 0]
    return "ccw" if z > 0 else "cw"


def track_frame(E, N, ok):
    P = np.column_stack([E, N])
    c = P[ok].mean(axis=0)
    u = np.linalg.svd(P[ok] - c, full_matrices=False)[2][0]
    return c, u


def vb_profile(s, d):
    m = np.isfinite(s) & np.isfinite(d)
    if m.sum() < 10:
        return None
    o = np.argsort(s[m])
    sv, dv = s[m][o], d[m][o]
    pad = np.concatenate([dv[:1], dv, dv[-1:]])            # running median of 3
    med = np.median(np.stack([pad[:-2], pad[1:-1], pad[2:]]), axis=0)
    return sv, med


def numbering_test(f):
    """RMS (slant depth - VB profile at the footprint) for 8 layouts."""
    E, N, hdg, vb, B, freq = f["E"], f["N"], f["hdg"], f["vb"], f["B"], f["freq"]
    if any(v is None for v in (E, N, hdg, vb, B)):
        return None
    ok = np.isfinite(E) & np.isfinite(N) & np.isfinite(hdg)
    if ok.sum() < 20:
        return None
    c, u = track_frame(E, N, ok)
    s_boat = (np.column_stack([E, N]) - c) @ u
    prof = vb_profile(np.where(ok, s_boat, np.nan), vb)
    if prof is None:
        return None
    sv, med = prof
    base = np.where((freq if freq is not None else np.full(E.size, 3000.0)) >= FREQ_SPLIT, 0.0, 45.0)
    tan_t = math.tan(math.radians(TILT_DEG))
    out = {}
    for sense in (+1, -1):
        for off in (0.0, 90.0, 180.0, 270.0):
            res = []
            for k in range(4):
                az = np.radians(hdg + off + sense * (base + 90.0 * k))
                d = B[:, k]
                sk = s_boat + d * tan_t * (np.sin(az) * u[0] + np.cos(az) * u[1])
                m = ok & np.isfinite(d) & (sk > sv[0]) & (sk < sv[-1])
                res.append(d[m] - np.interp(sk[m], sv, med))
            r = np.concatenate(res)
            if r.size >= 30:
                out[(sense, off)] = (float(np.sqrt(np.mean(r ** 2))), int(r.size))
    return out or None


# --------------------------------------------------------------------------- #
#  per-file report
# --------------------------------------------------------------------------- #

def report(path: Path, f: dict) -> dict:
    out = {"file": path.stem}
    print("=" * 100)
    print(path.name)

    print("\n[1] Setup")
    setp = f["setp"]
    for fld in getattr(setp, "_fieldnames", []):
        a = np.asarray(getattr(setp, fld))
        if a.dtype.kind in "biuf" and a.size <= 6:
            txt = np.array2string(a.ravel(), precision=4)
        elif a.dtype.kind in "US":
            txt = str(getattr(setp, fld))
        else:
            txt = f"<{a.dtype} {a.shape}>"
        flag = "   <--" if (fld.lower().startswith("edges_") and fld.endswith(("Method", "DistanceToBank"))
                             or fld in FLAGGED) else ""
        print(f"    {fld:38s} {txt}{flag}")
        if fld in ("depthReference", "Edges_0__Method", "Edges_1__Method") \
                and a.dtype.kind in "biuf" and a.size == 1:
            out[fld] = float(a.ravel()[0])
    shape = {2.0: "triangular", 1.0: "rectangular", 0.0: "Q de usuario"}
    for side, key in (("izquierda (Edges_0)", "Edges_0__Method"),
                      ("derecha   (Edges_1)", "Edges_1__Method")):
        if key in out:
            print(f"    -> margen {side}: {shape.get(out[key], 'código ' + str(out[key]))}")
    if "depthReference" in out:
        print(f"    -> referencia de profundidad en RSL: "
              f"{'VB' if out['depthReference'] < 0.5 else 'BT'}")

    print(f"\n[2] Actitud [grados] ({f['att_src']})")
    pitch, roll, vb = f["pitch"], f["roll"], f["vb"]
    if pitch is not None and roll is not None:
        print(f"    {'pitch':16s} {dist(pitch)}")
        print(f"    {'roll':16s} {dist(roll)}")
        tilt = np.degrees(np.arccos(np.clip(np.cos(np.radians(pitch))
                                            * np.cos(np.radians(roll)), -1.0, 1.0)))
        print(f"    {'inclinación':16s} {dist(tilt)}")
        out.update(pitch_med=float(np.nanmedian(pitch)),
                   pitch95=float(np.nanpercentile(np.abs(pitch), 95)),
                   roll95=float(np.nanpercentile(np.abs(roll), 95)),
                   tilt95=float(np.nanpercentile(tilt, 95)))
        if vb is not None and np.any(np.isfinite(vb)):
            dmed = float(np.nanmedian(vb))
            dmax = float(np.nanpercentile(vb, 95))
            out["shift95"] = dmax * math.tan(math.radians(out["tilt95"]))
            print(f"    -> con inclinación p95 ({out['tilt95']:.1f}°) la huella del VB se "
                  f"corre {dmed * math.tan(math.radians(out['tilt95'])):.2f} m a la "
                  f"profundidad mediana ({dmed:.2f} m) y {out['shift95']:.2f} m a la p95 "
                  f"({dmax:.2f} m), más altura de antena x sen(inclinación)")
    else:
        print("    sin pitch/roll en el archivo")

    print("\n[3] Numeración de haces")
    mats = matrices(f["tm"])
    lays = sorted({layout_of(D) for _fq, D in mats})
    if lays:
        out["matrix"] = "/".join(lays)
        print(f"    matriz: {'antihorario (haz 2 a BABOR)' if lays == ['ccw'] else 'horario (haz 2 a estribor)' if lays == ['cw'] else 'contradictoria'}"
              "   [RSL: +X hacia el haz 1, +Z arriba]")
    sc = numbering_test(f)
    if sc:
        best = min(sc, key=lambda h: sc[h][0])
        print("    datos: RMS (haz inclinado - perfil VB en la huella), mejor primero")
        for (sense, off), (rms, n) in sorted(sc.items(), key=lambda kv: kv[1][0]):
            tag = ""
            if (sense, off) == (-1, 0.0):
                tag = "   <- M9 (ccw), pipeline 6.5"
            elif (sense, off) == (1, 0.0):
                tag = "   <- regla v3-v6.4 (cw)"
            print(f"      haz 1 a {off:5.0f}°  {'horario    ' if sense > 0 else 'antihorario'}"
                  f"   RMS {rms:.3f} m  (n={n}){tag}")
        if (-1, 0.0) in sc and (1, 0.0) in sc:
            out["rms_ccw"], out["rms_cw"] = sc[(-1, 0.0)][0], sc[(1, 0.0)][0]
        out["numbering"] = ("ccw" if best == (-1, 0.0) else "cw" if best == (1, 0.0)
                            else f"{best[1]:.0f}°/{'h' if best[0] > 0 else 'ah'}")
        print(f"    -> datos: mejor {out['numbering']}")

    if mats:
        print("\n[4] Direcciones de haz según Transformation_Matrices "
              "(acimut desde el haz 1, sentido horario visto desde arriba)")
        for fq, D in mats:
            cells = []
            for i, w in enumerate(D):
                az = math.degrees(math.atan2(-w[1], w[0])) % 360   # +Y = port
                cells.append(f"h{i + 1} {az:5.1f}° incl {math.degrees(math.acos(min(1.0, abs(w[2])))):4.1f}°")
            print(f"    {fq:7.0f} kHz   " + "   ".join(cells))
    return out


# --------------------------------------------------------------------------- #
#  [5] attitude sign from opposite passes
# --------------------------------------------------------------------------- #

def _rot_down(hdg, pitch, roll):
    """Down axis of the instrument in (north, east, down) — aerospace R[:, 2]."""
    ps, cs = np.sin(np.radians(hdg)), np.cos(np.radians(hdg))
    pt, ct = np.sin(np.radians(pitch)), np.cos(np.radians(pitch))
    pr, cr = np.sin(np.radians(roll)), np.cos(np.radians(roll))
    return (cs * pt * cr + ps * pr, ps * pt * cr - cs * pr, ct * cr)


def _stem_time(stem):
    m = re.match(r"(\d{14})", stem)
    if not m:
        return None
    try:
        import datetime
        return datetime.datetime.strptime(m.group(1), "%Y%m%d%H%M%S")
    except ValueError:
        return None


def occupations(files):
    """Consecutive crossings of the same section: close centroids, parallel
    axes, contiguous in time. Only groups with passes in BOTH directions."""
    info = []
    for stem, f in files:
        E, N = f["E"], f["N"]
        if E is None or f["vb"] is None or f["hdg"] is None or f["pitch"] is None:
            continue
        ok = np.isfinite(E) & np.isfinite(N)
        if ok.sum() < 20:
            continue
        c, u = track_frame(E, N, ok)
        s = (np.column_stack([E, N])[ok] - c) @ u
        if np.ptp(s) < 10.0:
            continue
        info.append(dict(stem=stem, f=f, c=c, u=u, t=_stem_time(stem), span=float(np.ptp(s))))
    info.sort(key=lambda r: (r["t"] or 0, r["stem"]))
    groups, cur = [], []
    for r in info:
        if cur:
            a = cur[0]
            ang = math.degrees(math.acos(min(1.0, abs(float(a["u"] @ r["u"])))))
            near = float(np.hypot(*(r["c"] - a["c"]))) <= OCC_DIST
            gap_ok = (r["t"] is None or cur[-1]["t"] is None
                      or (r["t"] - cur[-1]["t"]).total_seconds() <= OCC_GAP_MIN * 60)
            if near and ang <= OCC_ANGLE and gap_ok:
                cur.append(r)
                continue
            groups.append(cur)
        cur = [r]
    if cur:
        groups.append(cur)
    out = []
    for g in groups:
        if len(g) < 2:
            continue
        u = g[0]["u"]
        dirs = []
        for r in g:
            E, N = r["f"]["E"], r["f"]["N"]
            ok = np.isfinite(E) & np.isfinite(N)
            s = np.column_stack([E, N])[ok] @ u
            dirs.append(1 if s[-1] > s[0] else -1)
        if len(set(dirs)) == 2:
            for r, d in zip(g, dirs):
                r["dir"] = d
            out.append(g)
    return out


def _pass_points(r, c, u, ps, rs, H, lag):
    """VB footprint s on the group axis under one hypothesis."""
    f = r["f"]
    E, N, hdg, vb = f["E"], f["N"], f["hdg"], f["vb"]
    p = np.nan_to_num(f["pitch"]) * ps if ps else np.zeros(E.size)
    q = np.nan_to_num(f["roll"]) * rs if rs else np.zeros(E.size)
    dn, de, dz = _rot_down(hdg, p, q)
    dz = np.where(dz > 0.2, dz, np.nan)
    # the antenna sits H above the transducer along the instrument axis, so the
    # transducer is antenna + H*down; the VB footprint is a further vb/dz along
    # the same axis (vb is already the VERTICAL depth)
    k = H + vb / dz
    ex = E + k * de
    nx = N + k * dn
    if lag:
        t = f["t"] if f["t"] is not None else np.arange(E.size, dtype=float)
        vE = np.gradient(E, t)
        vN = np.gradient(N, t)
        ex = ex + vE * lag
        nx = nx + vN * lag
    s = (np.column_stack([ex, nx]) - c) @ u
    ok = np.isfinite(s) & np.isfinite(vb)
    return s[ok], vb[ok]


def _pair_rms(a, b):
    sa, da = a
    prof = vb_profile(*b)
    if prof is None:
        return None
    sv, med = prof
    m = (sa > sv[0]) & (sa < sv[-1])
    if m.sum() < 10:
        return None
    r = da[m] - np.interp(sa[m], sv, med)
    return float(np.sum(r ** 2)), int(m.sum())


def group_score(g, ps, rs, H, lag):
    E = np.concatenate([r["f"]["E"] for r in g])
    N = np.concatenate([r["f"]["N"] for r in g])
    ok = np.isfinite(E) & np.isfinite(N)
    c, u = track_frame(E, N, ok)
    pts = [_pass_points(r, c, u, ps, rs, H, lag) for r in g]
    ss, nn, ss_same, nn_same = 0.0, 0, 0.0, 0
    for i in range(len(g)):
        for j in range(len(g)):
            if i == j:
                continue
            q = _pair_rms(pts[i], pts[j])
            if q is None:
                continue
            if g[i]["dir"] != g[j]["dir"]:
                ss += q[0]; nn += q[1]
            else:
                ss_same += q[0]; nn_same += q[1]
    return ss, nn, ss_same, nn_same


def group_power(g):
    """What a real tilt WOULD do to this occupation, before testing it.
    Each footprint moves depth*tan(pitch) toward the bow (projected on the
    section); on an opposite pair the shifts add along the section, and a shift
    only shows up as a depth disagreement where the bed slopes. Point by point:
    expected = RMS[ slope(s) x (shift_a(s) + median shift of the other pass) ].
    Antenna height is left out, so this is a LOWER bound."""
    E = np.concatenate([r["f"]["E"] for r in g])
    N = np.concatenate([r["f"]["N"] for r in g])
    ok = np.isfinite(E) & np.isfinite(N)
    c, u = track_frame(E, N, ok)
    sp, dp = [], []
    for r in g:
        f = r["f"]
        r["s_lvl"] = (np.column_stack([f["E"], f["N"]]) - c) @ u
        proj = np.abs(np.sin(np.radians(f["hdg"])) * u[0] + np.cos(np.radians(f["hdg"])) * u[1])
        r["shift_pts"] = f["vb"] * np.tan(np.radians(np.abs(f["pitch"]))) * proj
        r["shift"] = float(np.nanmedian(r["shift_pts"]))
        r["pitch_med"] = float(np.nanmedian(f["pitch"]))
        sp.append(r["s_lvl"]); dp.append(f["vb"])
    s = np.concatenate(sp); d = np.concatenate(dp)
    m = np.isfinite(s) & np.isfinite(d)
    if m.sum() < 20:
        return float("nan"), float("nan")
    edges = np.arange(np.floor(s[m].min()), np.ceil(s[m].max()) + 1.0, 1.0)
    idx = np.digitize(s[m], edges)
    xb, zb = [], []
    for k in np.unique(idx):
        q = idx == k
        if q.sum() >= 2:
            xb.append(float(np.median(s[m][q]))); zb.append(float(np.median(d[m][q])))
    if len(xb) < 5:
        return float("nan"), float("nan")
    xb, zb = np.asarray(xb), np.asarray(zb)
    slope = np.gradient(zb, xb)
    slope_rms = float(np.sqrt(np.mean(slope ** 2)))
    ss, nn = 0.0, 0
    for a in g:
        for b in g:
            if a["dir"] == b["dir"]:
                continue
            q = np.isfinite(a["s_lvl"]) & np.isfinite(a["shift_pts"])
            e = np.interp(a["s_lvl"][q], xb, slope) * (a["shift_pts"][q] + b["shift"])
            ss += float(np.sum(e ** 2)); nn += int(e.size)
    return slope_rms, (math.sqrt(ss / nn) if nn else float("nan"))


def attitude_test(files):
    groups = occupations(files)
    print("\n" + "=" * 100)
    print("[5] SIGNO DE PITCH / ROLL Y ALTURA DE ANTENA — pasadas opuestas sobre la misma sección")
    if not groups:
        print("    No hay ocupaciones con pasadas en los dos sentidos entre los archivos dados.")
        print("    Pasá la carpeta con los aforos (repeticiones alternando margen de inicio).")
        return None
    exp_all, floor_g = [], []
    for g in groups:
        slope_rms, expected = group_power(g)
        _ss, _nn, s0, n0 = group_score(g, 0.0, 0.0, 0.0, 0.0)
        floor = math.sqrt(s0 / n0) if n0 else float("nan")
        exp_all.append(expected); floor_g.append(floor)
        print(f"    ocupación: {', '.join(r['stem'] for r in g)}  "
              f"(sentidos {''.join('+' if r['dir'] > 0 else '-' for r in g)})")
        pm = ", ".join(f"{r['pitch_med']:+.1f}" for r in g)
        shv = ", ".join(f"{r['shift']:.2f}" for r in g)
        print(f"        pitch mediano por pasada {pm}°; corrimiento de huella proyectado "
              f"{shv} m; pendiente RMS del lecho {slope_rms:.2f}")
        print(f"        desacuerdo que produciría una inclinación real: ~{expected:.3f} m, "
              f"contra un piso de {floor:.3f} m entre pasadas del mismo sentido")
    hyps = [("nivelado (off)", 0.0, 0.0)] + [(f"pitch {pn} / roll {rn}", pv, rv)
                                             for pn, pv in PITCH_HYP for rn, rv in ROLL_HYP]

    def _total(ps, rs, H, lag):
        tot = [group_score(g, ps, rs, H, lag) for g in groups]
        nn = sum(t[1] for t in tot); n0 = sum(t[3] for t in tot)
        return (math.sqrt(sum(t[0] for t in tot) / nn) if nn else float("nan"),
                math.sqrt(sum(t[2] for t in tot) / n0) if n0 else float("nan"))

    rows = []
    for name, ps, rs in hyps:
        Hs = (0.0,) if ps == 0 else ANT_GRID
        # without GNSS lag: the physical model alone
        nolag = min(((_total(ps, rs, H, 0.0)[0], H) for H in Hs), key=lambda q: q[0])
        # with the lag as a free nuisance parameter
        wlag = min(((_total(ps, rs, H, lag)[0], H, float(lag)) for H in Hs for lag in LAG_GRID),
                   key=lambda q: q[0])
        rows.append(dict(name=name, ps=ps, rs=rs, rms0=nolag[0], H0=nolag[1],
                         rms=wlag[0], H=wlag[1], lag=wlag[2]))
    same = _total(0.0, 0.0, 0.0, 0.0)[1]
    rows.sort(key=lambda r: r["rms0"])
    print("\n    RMS entre pasadas OPUESTAS del perfil VB (menor = mejor), con la altura de "
          "antena ajustada; a la derecha, dejando libre además un retardo GNSS")
    print(f"    {'hipótesis':34s} {'RMS':>7s} {'antena':>7s}   |  {'RMS':>7s} {'antena':>7s} "
          f"{'retardo':>8s}")
    for r in rows:
        print(f"    {r['name']:34s} {r['rms0']:7.3f} {r['H0']:6.2f}m   |  {r['rms']:7.3f} "
              f"{r['H']:6.2f}m {r['lag']:+7.1f}s")
    print(f"    (piso de ruido: RMS entre pasadas del MISMO sentido = {same:.3f} m)")
    print("\n    por ocupación, sin retardo ni antena: mismo sentido / opuestas nivelado / "
          "mejor hipótesis con signo")
    for g in groups:
        s_off, n_off, s0, n0 = group_score(g, 0.0, 0.0, 0.0, 0.0)
        cand = []
        for pn, pv in PITCH_HYP:
            for rn, rv in ROLL_HYP:
                q = group_score(g, pv, rv, 0.0, 0.0)
                if q[1]:
                    cand.append((math.sqrt(q[0] / q[1]), f"{pn}/{rn}"))
        cand.sort()
        f0 = math.sqrt(s0 / n0) if n0 else float("nan")
        fo = math.sqrt(s_off / n_off) if n_off else float("nan")
        print(f"      {g[0]['stem']}..{g[-1]['stem'][-6:]} ({len(g)} pasadas): "
              f"{f0:.3f} / {fo:.3f} / {cand[0][0]:.3f} {cand[0][1]}"
              + ("   <- las opuestas discrepan más que las del mismo sentido"
                 if np.isfinite(f0) and fo > 1.15 * f0 else ""))

    off = next(r for r in rows if r["ps"] == 0.0)
    signed = [r for r in rows if r["ps"] != 0.0]
    first = signed[0]
    gain = off["rms0"] - first["rms0"]
    other_pitch = min(r["rms0"] for r in signed if r["ps"] != first["ps"])
    other_roll = min(r["rms0"] for r in signed if r["ps"] == first["ps"] and r["rs"] != first["rs"])
    margin_p = other_pitch - first["rms0"]
    margin_r = other_roll - first["rms0"]
    exp_ok = [e for e in exp_all if np.isfinite(e)]
    fl_ok = [f for f in floor_g if np.isfinite(f)]
    power = (max(exp_ok) / float(np.median(fl_ok))) if exp_ok and fl_ok else 0.0
    print(f"\n    -> mejor hipótesis con signo: {first['name']}, antena {first['H0']:.2f} m")
    print(f"       mejora contra nivelado {gain:.3f} m; ventaja sobre el otro signo de "
          f"pitch {margin_p:.3f} m, sobre el otro signo de roll {margin_r:.3f} m")
    print(f"       poder de la prueba: el desacuerdo que produciría una inclinación real es "
          f"{power:.2f} veces el piso de ruido (hace falta ~0,5 o más para decidir)")
    if power < 0.5:
        verdict = ("SIN PODER DE RESOLUCIÓN: en estas ocupaciones una inclinación real "
                   "movería los perfiles menos que el ruido entre pasadas (lecho plano o "
                   "poco cabeceo). El signo no se puede decidir con estos datos: dejar "
                   "pitch-roll = off o determinar el signo en banco")
    elif gain < 0.005:
        verdict = ("la actitud no mejora la coincidencia aunque la prueba tenía poder: "
                   "dejar pitch-roll = off y revisar si el pitch reportado es la "
                   "inclinación real del transductor")
    elif margin_p < 0.003:
        verdict = ("el signo del pitch no queda determinado (los dos signos casi "
                   "empatan): dejar pitch-roll = off")
    else:
        pn = "bow-up" if first["ps"] > 0 else "bow-down"
        rn = ("stbd-down" if first["rs"] > 0 else "port-down") if margin_r >= 0.003 else "none"
        # the lag is taken only if it really improves on the physical model
        use_lag = first["rms"] < first["rms0"] - 0.005
        H = first["H"] if use_lag else first["H0"]
        verdict = (f"pitch-roll = on, pitch-sign = {pn}, roll-sign = {rn}, "
                   f"antenna-height = {H:.2f}"
                   + (f", gnss-lag = {first['lag']:+.1f}" if use_lag else "")
                   + "  (medí la altura real de la antena)")
        if rn == "none":
            verdict += ("\n       el signo del roll no queda determinado (efecto chico): "
                        "roll-sign = none aplica sólo el pitch")
        if off["rms"] <= min(first["rms0"], first["rms"]) + 0.005:
            verdict += (f"\n       ojo: un retardo GNSS de {off['lag']:+.1f} s sin actitud explica "
                        f"casi lo mismo ({off['rms']:.3f} m); confirmá el signo en banco")
    print(f"    -> {verdict}")
    return verdict


# --------------------------------------------------------------------------- #

def main(argv):
    paths = []
    for a in argv:
        p = Path(a)
        paths += sorted(p.glob("*.mat")) if p.is_dir() else [p]
    if not paths:
        sys.exit("uso: python diagnostico_mat.py <archivo.mat | carpeta> [...]")
    rows, files = [], []
    for p in paths:
        try:
            f = read(p)
            rows.append(report(p, f))
            files.append((p.stem, f))
        except Exception as e:                       # keep going on a bad file
            print(f"[error] {p.name}: {e}")
    if len(rows) > 1:
        def fmt(v, spec):
            return "" if v is None else format(v, spec)
        w = max(16, max(len(r["file"]) for r in rows))
        print("\n" + "=" * 100 + "\nRESUMEN POR ARCHIVO")
        print(f"  {'archivo':{w}s} {'ref':>4s} {'Edg0':>5s} {'Edg1':>5s} {'pitch med':>9s} "
              f"{'|pitch|95':>9s} {'|roll|95':>8s} {'huella95':>8s} {'matriz':>7s} "
              f"{'RMS ccw':>8s} {'RMS cw':>7s}  datos")
        for r in rows:
            dr = r.get("depthReference")
            ref = "" if dr is None else ("VB" if dr < 0.5 else "BT")
            print(f"  {r['file']:{w}s} {ref:>4s} {fmt(r.get('Edges_0__Method'), '5.0f'):>5s} "
                  f"{fmt(r.get('Edges_1__Method'), '5.0f'):>5s} "
                  f"{fmt(r.get('pitch_med'), '+9.2f'):>9s} {fmt(r.get('pitch95'), '9.2f'):>9s} "
                  f"{fmt(r.get('roll95'), '8.2f'):>8s} {fmt(r.get('shift95'), '8.2f'):>8s} "
                  f"{r.get('matrix', '-'):>7s} {fmt(r.get('rms_ccw'), '8.3f'):>8s} "
                  f"{fmt(r.get('rms_cw'), '7.3f'):>7s}  {r.get('numbering', '-')}")
        both = [r for r in rows if "rms_ccw" in r]
        if both:
            k = sum(r["rms_ccw"] < r["rms_cw"] for r in both)
            print(f"\n  numeración: ccw mejor que cw en {k} de {len(both)} archivos; "
                  f"matrices: {sorted({r.get('matrix', '-') for r in rows})}")
    attitude_test(files)


if __name__ == "__main__":
    main(sys.argv[1:])
