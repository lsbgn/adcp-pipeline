#!/usr/bin/env python3
"""
tools/comparar_corridas.py
================================================================================
Compare two output folders of process_adcp_bathimetric.py (run A = reference,
run B = new) and say, per profile, whether it changed, how much and where.

    python tools/comparar_corridas.py salida_v66 salida_v67
    python tools/comparar_corridas.py salida_v66 salida_v67 --plots all
    python tools/comparar_corridas.py salida_v66 salida_off --check   # regression

What it does
    1. IDENTITY. SHA-256 of every bathymetric_profile.csv and raw_bed_points.csv:
       'igual' means byte for byte the same file. For the raw cloud only its
       first 12 columns are hashed (v6.8 appended `ens` and `descarte`), and
       the points dropped by [descartes] are counted (descartes_a/_b).
    2. MAGNITUDE. Profiles that are not identical are compared on the GROUND,
       not on s: the nodes of B are projected (x_utm, y_utm) onto the section
       line of A, so a moved s = 0 (a wider or narrower section) or a slightly
       rotated axis does not show up as a false bed change. On the common
       stretch the bed of B is interpolated at the nodes of A and
       dz = z_B - z_A gives RMS, mean (bias), max |dz| and the s where it
       happens. z is bed_elev_m when both runs have a water surface, else
       -depth_m. Width, max depth, thalweg, water surface and the axis shift
       (mean distance of B's nodes to A's line, angle between axes) are
       reported separately.
    3. WHERE. One PNG per changed profile: beds A and B overlaid, dz(s) below.
    4. PARAMETERS. The [parámetros] blocks of both procesamiento.txt are
       diffed, so a change can be traced to a key (or to the script version).

Outputs (default <B>/_comparacion/):
    comparacion.csv       one row per profile, sorted by RMS (largest first)
    parametros_diff.txt   parameters that differ between the two runs
    <perfil>.png          overlay for each changed profile (--plots)

Status per profile
    igual       bathymetric_profile.csv identical
    igual_num   different file, but max |dz| <= --tol-z and same width
                (+-1 cm): formatting, float noise
    cambia      a real change
    solo_A / solo_B   the profile exists in one run only
    sin_solape  the two sections do not overlap on the ground

--check exits with code 1 when any profile is 'cambia', 'solo_A' or 'solo_B'
(for regression tests: e.g. position-fill = off against v6.6).
================================================================================
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import sys
from pathlib import Path

import numpy as np

__version__ = "1.0"

PROFILE_CSV = "bathymetric_profile.csv"
RAW_CSV = "raw_bed_points.csv"
RESUMEN = "_resumen"
RECORD = "procesamiento.txt"
INDEX = "survey_index.csv"
# parameters that always differ between two runs and say nothing
PARAM_IGNORE = {"outdir", "survey_name", "matfiles", "config"}
# survey_index.csv columns of run B copied to the comparison (6.7 screening)
POS_COLS = ["gga_q", "pos_hdop", "pos_jump", "pos_rebuilt", "pos_shift_max_m",
            "s_fold_m", "edge_drift_m"]


# ---------------------------------------------------------------- discovery

def find_profiles(root: Path) -> dict[str, Path]:
    """{profile id: folder} for every folder under `root` holding a profile.
    `root` itself counts when it is a single-transect output folder."""
    out = {}
    if (root / PROFILE_CSV).exists():
        out[root.name.removeprefix("out_")] = root
    for p in sorted(root.rglob(PROFILE_CSV)):
        d = p.parent
        if RESUMEN in d.parts or d == root:
            continue
        out.setdefault(d.name, d)
    return out


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


RAW_BASE_COLS = 12      # columns of raw_bed_points.csv up to v6.7


def raw_digest(path: Path):
    """(hash of the measured cloud, points discarded by [descartes]).

    v6.8 appended `ens` and `descarte` to raw_bed_points.csv. The hash covers
    only the first 12 columns (the cloud itself, as in v6.7), so a v6.7 and a
    v6.8 run of the same data compare equal; discarded points stay in the
    file and are counted apart."""
    if not path.exists():
        return None, None
    h = hashlib.sha256()
    n_disc, i_disc = 0, None
    with open(path, encoding="utf-8") as f:
        for k, line in enumerate(f):
            parts = line.rstrip("\n").split(",")
            if k == 0 and "descarte" in parts:
                i_disc = parts.index("descarte")
            elif i_disc is not None and len(parts) > i_disc and parts[i_disc].strip():
                n_disc += 1
            h.update((",".join(parts[:RAW_BASE_COLS]) + "\n").encode("utf-8"))
    return h.hexdigest(), (n_disc if i_disc is not None else None)


# ---------------------------------------------------------------- reading

def _num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def read_profile(path: Path) -> dict:
    """Columns of bathymetric_profile.csv as float arrays (NaN when empty)."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    cols = {}
    for k in ("s_m", "depth_m", "ws_elev_m", "bed_elev_m", "x_utm", "y_utm",
              "progresiva_m"):
        cols[k] = np.array([_num(r.get(k)) for r in rows], dtype=float)
    cols["brazo"] = rows[0].get("brazo", "") if rows else ""
    return cols


def read_index(root: Path) -> dict[str, dict]:
    path = root / RESUMEN / INDEX
    if not path.exists():
        return {}
    with open(path, encoding="utf-8", newline="") as f:
        return {r["perfil"]: r for r in csv.DictReader(f)}


def read_params(root: Path) -> tuple[dict, str]:
    """({key: value}, version line) from the run record, or ({}, '')."""
    path = root / RESUMEN / RECORD
    if not path.exists():
        return {}, ""
    params, version, inside = {}, "", False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("script version"):
            version = line.split(":", 1)[1].strip()
        if line.startswith("[parámetros]"):
            inside = True
            continue
        if inside:
            if not line.startswith("  ") or ":" not in line:
                if line.strip() == "" or line.startswith("["):
                    inside = False
                continue
            k, v = line.strip().split(":", 1)
            params[k.strip()] = v.strip()
    return params, version


_RE_POS = re.compile(r"posiciones: (\d+)/(\d+) ensambles reconstruidos.*?"
                     r"(\d+) por HDOP, (\d+) por salto.*?corrimiento máx ([\d.]+) m")


def pos_from_log(path: Path) -> dict:
    """6.7 screening counts of one transect from its process_log.txt, for the
    profiles that have no survey_index.csv row (members of an aforo group)."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _RE_POS.search(line)
        if m:
            out.update(pos_rebuilt=m.group(1), pos_hdop=m.group(3),
                       pos_jump=m.group(4), pos_shift_max_m=m.group(5))
        elif "posiciones: GGA válido en todos" in line:
            out.update(pos_rebuilt="0", pos_hdop="0", pos_jump="0",
                       pos_shift_max_m="0.00")
        elif "calidad GGA" in line:
            q = re.search(r"calidad GGA ([\d: ]+?) \(", line)
            if q:
                out["gga_q"] = q.group(1).strip()
        elif line.startswith("[QA] retroceso de s"):
            f = re.search(r"cruce: ([\d.]+) m", line)
            if f:
                out["s_fold_m"] = f.group(1)
        elif line.startswith("[QA] deriva en el borde"):
            d = [float(x) for x in re.findall(r"([\d.]+) m (?:transversal|a lo largo)", line)]
            if d:
                prev = float(out.get("edge_drift_m") or 0.0)
                out["edge_drift_m"] = f"{max([prev] + d):.1f}"
    return out


# ---------------------------------------------------------------- comparison

def _z(p: dict, absolute: bool) -> np.ndarray:
    return p["bed_elev_m"] if absolute else -p["depth_m"]


def _nanmed(x):
    x = np.asarray(x, dtype=float)
    return float(np.median(x[np.isfinite(x)])) if np.isfinite(x).any() else float("nan")


def _thalweg(z, depth):
    interior = depth > 1e-6
    return float(np.nanmin(z[interior])) if interior.any() else float(np.nanmin(z))


def compare_profiles(a: dict, b: dict) -> dict:
    """Ground-aligned comparison of two profiles of the same transect."""
    absolute = bool(np.isfinite(a["ws_elev_m"]).any() and np.isfinite(b["ws_elev_m"]).any())
    za, zb = _z(a, absolute), _z(b, absolute)
    out = dict(
        z_kind="cota" if absolute else "-prof",
        width_a=float(np.nanmax(a["s_m"]) - np.nanmin(a["s_m"])),
        width_b=float(np.nanmax(b["s_m"]) - np.nanmin(b["s_m"])),
        maxd_a=float(np.nanmax(a["depth_m"])), maxd_b=float(np.nanmax(b["depth_m"])),
        thal_a=_thalweg(za, a["depth_m"]), thal_b=_thalweg(zb, b["depth_m"]),
        ws_a=_nanmed(a["ws_elev_m"]) if absolute else float("nan"),
        ws_b=_nanmed(b["ws_elev_m"]) if absolute else float("nan"),
        km_a=_nanmed(a["progresiva_m"]), km_b=_nanmed(b["progresiva_m"]),
    )
    geo = (np.isfinite(a["x_utm"]).sum() >= 2 and np.isfinite(b["x_utm"]).sum() >= 2)
    if geo:
        # A's section line: first to last node (the profile is a straight line)
        P = np.column_stack([a["x_utm"], a["y_utm"]])
        Q = np.column_stack([b["x_utm"], b["y_utm"]])
        okp = np.isfinite(P).all(axis=1)
        p0, p1 = P[okp][0], P[okp][-1]
        u = (p1 - p0) / max(np.linalg.norm(p1 - p0), 1e-9)
        n = np.array([-u[1], u[0]])
        s_b = (Q - p0) @ u + a["s_m"][okp][0]
        off_b = (Q - p0) @ n
        q0, q1 = Q[np.isfinite(Q).all(axis=1)][[0, -1]]
        ub = (q1 - q0) / max(np.linalg.norm(q1 - q0), 1e-9)
        ang = math.degrees(math.acos(min(1.0, abs(float(u @ ub)))))
        out.update(align="geo", axis_offset=float(np.nanmean(np.abs(off_b))), axis_angle=ang)
    else:
        s_b = b["s_m"]
        out.update(align="s", axis_offset=float("nan"), axis_angle=float("nan"))
    s_a = a["s_m"]
    order = np.argsort(s_b)
    s_b, zb_s = s_b[order], zb[order]
    ok = np.isfinite(s_b) & np.isfinite(zb_s)
    s_b, zb_s = s_b[ok], zb_s[ok]
    lo, hi = max(np.nanmin(s_a), s_b.min()), min(np.nanmax(s_a), s_b.max())
    common = (s_a >= lo) & (s_a <= hi) & np.isfinite(za)
    out["s_b_on_a"] = s_b
    out["zb_on_a"] = zb_s
    out["za"] = za
    if hi <= lo or common.sum() < 2:
        out.update(n_common=0, rms=float("nan"), bias=float("nan"),
                   max_abs=float("nan"), s_at_max=float("nan"), dz=None)
        return out
    zb_i = np.interp(s_a[common], s_b, zb_s)
    dz = zb_i - za[common]
    k = int(np.argmax(np.abs(dz)))
    out.update(n_common=int(common.sum()), overlap=(float(lo), float(hi)),
               rms=float(np.sqrt(np.mean(dz ** 2))), bias=float(np.mean(dz)),
               max_abs=float(np.abs(dz[k])), s_at_max=float(s_a[common][k]),
               dz=(s_a[common], dz))
    return out


# ---------------------------------------------------------------- plotting

def plot_overlay(path: Path, pid: str, a: dict, b: dict, c: dict,
                 label_a: str, label_b: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid = "#1f1f1f", "#6b6b6b", "#e6e6e6"
    col_a, col_b, col_dz = "#9a9a9a", "#1f1f1f", "#2b6cb0"
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True,
                                  gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.08))
    za = c["za"]
    ax.plot(a["s_m"], za, color=col_a, lw=2.0, ls=(0, (5, 3)), label=f"A · {label_a}")
    ax.plot(c["s_b_on_a"], c["zb_on_a"], color=col_b, lw=1.6, label=f"B · {label_b}")
    if c["z_kind"] == "cota" and np.isfinite(c["ws_a"]):
        ax.axhline(c["ws_a"], color="#3b5b92", lw=1.0, ls="--", alpha=0.7)
    ax.set_ylabel("Cota [m]" if c["z_kind"] == "cota" else "−Profundidad [m]", color=ink)
    ax.legend(loc="lower right", frameon=False, fontsize=9, labelcolor=ink)
    stats = (f"RMS Δz {c['rms']:.3f} m · sesgo {c['bias']:+.3f} m · máx |Δz| "
             f"{c['max_abs']:.3f} m en s = {c['s_at_max']:.1f} m · Δancho "
             f"{c['width_b'] - c['width_a']:+.2f} m")
    if c["align"] == "geo":
        stats += f" · eje B a {c['axis_offset']:.2f} m / {c['axis_angle']:.1f}° del de A"
    fig.text(0.08, 0.955, f"Comparación de corridas — perfil {pid}",
             fontsize=13, fontweight="bold", color=ink, va="bottom")
    fig.text(0.08, 0.925, stats, fontsize=8.5, color=muted, va="bottom")
    if c.get("dz") is not None:
        s, dz = c["dz"]
        axd.axhline(0, color=muted, lw=0.8)
        axd.fill_between(s, 0, dz, color=col_dz, alpha=0.18, lw=0)
        axd.plot(s, dz, color=col_dz, lw=1.6)
        m = max(0.05, float(np.nanmax(np.abs(dz))) * 1.15)
        axd.set_ylim(-m, m)
        axd.plot([c["s_at_max"]], [dz[int(np.argmax(np.abs(dz)))]], "o", ms=8,
                 color=col_dz, mec="white", mew=2)
    axd.set_ylabel("Δz = B − A [m]", color=ink)
    axd.set_xlabel("Distancia transversal en el eje de A, s [m]", color=ink)
    for x in (ax, axd):
        x.grid(True, color=grid, lw=0.8)
        x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            x.spines[sp].set_color(muted)
        x.tick_params(colors=muted, labelsize=9)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.1)
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- main

def _f(v, fmt="{:.3f}"):
    return "" if v is None or (isinstance(v, float) and not math.isfinite(v)) else fmt.format(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a", help="reference output folder (A)")
    ap.add_argument("run_b", help="new output folder (B)")
    ap.add_argument("--outdir", default=None, help="default <B>/_comparacion")
    ap.add_argument("--tol-z", type=float, default=0.005,
                    help="max |dz| [m] still counted as 'igual_num' (default 0.005)")
    ap.add_argument("--plots", choices=["cambia", "all", "none"], default="cambia",
                    help="overlay PNGs for changed profiles (default), all, or none")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any profile changed or is missing on one side")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    A, B = Path(args.run_a), Path(args.run_b)
    for d in (A, B):
        if not d.is_dir():
            sys.exit(f"[error] no existe la carpeta {d}")
    out = Path(args.outdir) if args.outdir else B / "_comparacion"
    out.mkdir(parents=True, exist_ok=True)
    la, lb = A.name, B.name

    pa, pb = find_profiles(A), find_profiles(B)
    idx_b = read_index(B)
    ids = sorted(set(pa) | set(pb))
    if not ids:
        sys.exit("[error] no hay bathymetric_profile.csv en ninguna de las dos carpetas")

    rows = []
    for pid in ids:
        r = dict(perfil=pid)
        if pid not in pb:
            r["estado"] = "solo_A"; rows.append(r); continue
        if pid not in pa:
            r["estado"] = "solo_B"; rows.append(r); continue
        fa, fb = pa[pid] / PROFILE_CSV, pb[pid] / PROFILE_CSV
        same = sha256(fa) == sha256(fb)
        (ha, da), (hb, db) = raw_digest(pa[pid] / RAW_CSV), raw_digest(pb[pid] / RAW_CSV)
        r["nube_igual"] = "" if ha is None or hb is None else ("si" if ha == hb else "no")
        r["descartes_a"] = "" if da is None else da          # 6.8
        r["descartes_b"] = "" if db is None else db
        a, b = read_profile(fa), read_profile(fb)
        c = compare_profiles(a, b)
        dw = c["width_b"] - c["width_a"]
        if same:
            estado = "igual"
        elif c["n_common"] == 0:
            estado = "sin_solape"
        elif c["max_abs"] <= args.tol_z and abs(dw) <= 0.01:
            estado = "igual_num"
        else:
            estado = "cambia"
        r.update(estado=estado, z=c["z_kind"], alineacion=c["align"],
                 rms_dz_m=c["rms"], sesgo_dz_m=c["bias"], max_abs_dz_m=c["max_abs"],
                 s_max_dz_m=c["s_at_max"],
                 ancho_a_m=c["width_a"], ancho_b_m=c["width_b"], d_ancho_m=dw,
                 prof_max_a_m=c["maxd_a"], prof_max_b_m=c["maxd_b"],
                 d_prof_max_m=c["maxd_b"] - c["maxd_a"],
                 thalweg_a_m=c["thal_a"], thalweg_b_m=c["thal_b"],
                 d_thalweg_m=c["thal_b"] - c["thal_a"],
                 d_pelo_agua_m=c["ws_b"] - c["ws_a"],
                 d_progresiva_m=c["km_b"] - c["km_a"],
                 eje_offset_m=c["axis_offset"], eje_angulo_deg=c["axis_angle"])
        if pid in idx_b:
            for k in POS_COLS:
                r[k] = idx_b[pid].get(k, "")
            r["en_indice"] = "si"
        else:
            # a member of an aforo group has its own folder but no row in
            # survey_index.csv (only the group has one): read its own log
            r.update(pos_from_log(pb[pid] / "process_log.txt"))
            r["en_indice"] = "no (miembro de grupo)" if idx_b else ""
        if args.plots == "all" or (args.plots == "cambia" and estado == "cambia"):
            try:
                plot_overlay(out / f"{pid}.png", pid, a, b, c, la, lb)
                r["figura"] = f"{pid}.png"
            except Exception as e:                      # never fatal
                r["figura"] = f"error: {e}"
        rows.append(r)

    # ---- comparacion.csv, largest change first
    order = {"cambia": 0, "sin_solape": 1, "solo_A": 2, "solo_B": 2, "igual_num": 3, "igual": 4}
    rows.sort(key=lambda r: (order.get(r["estado"], 9),
                             -(r.get("rms_dz_m") if isinstance(r.get("rms_dz_m"), float)
                               and math.isfinite(r["rms_dz_m"]) else -1.0)))
    cols = ["perfil", "estado", "rms_dz_m", "sesgo_dz_m", "max_abs_dz_m", "s_max_dz_m",
            "ancho_a_m", "ancho_b_m", "d_ancho_m", "prof_max_a_m", "prof_max_b_m",
            "d_prof_max_m", "thalweg_a_m", "thalweg_b_m", "d_thalweg_m",
            "d_pelo_agua_m", "d_progresiva_m", "eje_offset_m", "eje_angulo_deg",
            "z", "alineacion", "nube_igual", "descartes_a", "descartes_b",
            *POS_COLS, "en_indice", "figura"]
    with open(out / "comparacion.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([_f(r[k], "{:.2f}" if k in ("eje_angulo_deg", "d_progresiva_m")
                           else "{:.3f}") if isinstance(r.get(k), float) else r.get(k, "")
                        for k in cols])

    # ---- parameters
    par_a, ver_a = read_params(A)
    par_b, ver_b = read_params(B)
    diff = [(k, par_a.get(k, "(no existe)"), par_b.get(k, "(no existe)"))
            for k in sorted(set(par_a) | set(par_b))
            if k not in PARAM_IGNORE and par_a.get(k) != par_b.get(k)]
    lines = [f"A: {A}   (script {ver_a or '?'})", f"B: {B}   (script {ver_b or '?'})", ""]
    if not par_a or not par_b:
        lines.append("falta procesamiento.txt en alguna de las dos corridas")
    elif diff:
        w1 = max(len(k) for k, _, _ in diff)
        lines += [f"{k:<{w1}}  A = {va}   ->   B = {vb}" for k, va, vb in diff]
    else:
        lines.append("mismos parámetros")
    (out / "parametros_diff.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- console summary
    cnt = {}
    for r in rows:
        cnt[r["estado"]] = cnt.get(r["estado"], 0) + 1
    print(f"\nComparación  A = {A}  ({ver_a or '?'})   B = {B}  ({ver_b or '?'})")
    print("  " + ", ".join(f"{v} {k}" for k, v in sorted(cnt.items(), key=lambda t: order.get(t[0], 9))))
    ch = [r for r in rows if r["estado"] in ("cambia", "sin_solape")]
    if ch:
        print(f"\n  {'perfil':<22}{'RMS Δz':>9}{'máx|Δz|':>9}{'en s':>8}{'Δancho':>9}"
              f"{'Δprof.máx':>11}  pos_rebuilt")
        for r in ch[:15]:
            print(f"  {r['perfil']:<22}{_f(r['rms_dz_m']):>9}{_f(r['max_abs_dz_m']):>9}"
                  f"{_f(r['s_max_dz_m'], '{:.1f}'):>8}{_f(r['d_ancho_m'], '{:+.2f}'):>9}"
                  f"{_f(r['d_prof_max_m'], '{:+.3f}'):>11}  {r.get('pos_rebuilt', '')}")
        if len(ch) > 15:
            print(f"  ... y {len(ch) - 15} más (ver comparacion.csv)")
    for r in rows:
        if r["estado"] in ("solo_A", "solo_B"):
            print(f"  {r['perfil']}: {r['estado']}")
    print(f"\n  parámetros distintos: {len(diff)}"
          + ("" if not diff else "  (" + ", ".join(k for k, _, _ in diff[:6])
             + (", ..." if len(diff) > 6 else "") + ")"))
    print(f"  salidas en: {out}\n")

    if args.check and any(r["estado"] in ("cambia", "solo_A", "solo_B", "sin_solape")
                          for r in rows):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
