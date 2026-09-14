#!/usr/bin/env python3
"""
compara_corridas.py — compare two runs of process_adcp_bathimetric.py.

    python compara_corridas.py <procesado_viejo> <procesado_nuevo>

Checks three things:

1. DEPTHS. `depth_m` must be IDENTICAL between v6.0 and v6.1: the numeric core
   did not change, so any difference is a regression, not a new feature.
   Depth does not depend on the location or on the water surface.

2. COVERAGE. Folders present in only one run are listed. The first version of
   this script walked the OLD run only, so a group that exists solely in the
   NEW run (grouping can change when the location changes) was silently left
   out of the comparison.

3. WHAT IS EXPECTED TO CHANGE. From survey_index.csv: river, official km and
   absolute bed elevation. Reported, not flagged — 6.1 changes the location
   rule (section crossing) and the water-surface model.

Needs pandas. `q_m3s` and the water surface are only compared where both runs
have the same profile id.
"""
import sys
from pathlib import Path

import pandas as pd


def profiles(run: Path):
    return {f.parent.name: f for f in sorted(run.glob("*/bathymetric_profile.csv"))}


def main(old: Path, new: Path):
    a, b = profiles(old), profiles(new)
    common = sorted(set(a) & set(b))
    only_old, only_new = sorted(set(a) - set(b)), sorted(set(b) - set(a))

    print(f"perfiles: {len(a)} en '{old}', {len(b)} en '{new}', "
          f"{len(common)} comparables")
    if only_old:
        print(f"  sólo en '{old}': {', '.join(only_old)}")
    if only_new:
        print(f"  sólo en '{new}': {', '.join(only_new)}")
    if only_old or only_new:
        print("  (normal si cambió la agrupación de aforos: las carpetas "
              "individuales sí se comparan)")

    print("\n--- profundidades (deben ser idénticas) ---")
    worst, worst_id = 0.0, ""
    for pid in common:
        da, db = pd.read_csv(a[pid]), pd.read_csv(b[pid])
        if len(da) != len(db):
            print(f"{pid:22s} DISTINTA CANTIDAD DE NODOS: {len(da)} vs {len(db)}")
            worst, worst_id = float("inf"), pid
            continue
        d = (da["depth_m"] - db["depth_m"]).abs().max()
        if d > worst:
            worst, worst_id = d, pid
        if d > 0:
            print(f"{pid:22s} max|Δdepth| = {d:.4f} m")
    print(f"PEOR: {worst:.4f} m" + (f" ({worst_id})" if worst else "")
          + " -> " + ("IDÉNTICAS" if worst == 0 else "REVISAR"))

    ia, ib = old / "_resumen" / "survey_index.csv", new / "_resumen" / "survey_index.csv"
    if not (ia.exists() and ib.exists()):
        return
    ja = pd.read_csv(ia).set_index("perfil")
    jb = pd.read_csv(ib).set_index("perfil")
    idx = ja.index.intersection(jb.index)

    print("\n--- cambios esperados (ubicación y pelo de agua) ---")
    moved = [(i, ja.at[i, "river"], jb.at[i, "river"]) for i in idx
             if ja.at[i, "river"] != jb.at[i, "river"]]
    if moved:
        print(f"cambiaron de río ({len(moved)}):")
        for i, r0, r1 in moved:
            print(f"  {i:22s} {r0} -> {r1}")
    else:
        print("ningún perfil cambió de río")

    def _report(label, d, nd):
        d = d.dropna()
        if d.empty:
            print(f"{label}: sin pares comparables")
            return
        print(f"{label}: mediana {d.median():.{nd}f} m, máx {d.max():.{nd}f} m "
              f"({d.idxmax()})")

    same = [i for i in idx if ja.at[i, "river"] == jb.at[i, "river"]]
    _report("|Δ km oficial| (mismo río)",
            (jb.loc[same, "km_oficial_m"] - ja.loc[same, "km_oficial_m"]).abs(), 2)

    if "ws_elev_m" in ja.columns and "ws_elev_m" in jb.columns:
        za, zb = ja.loc[idx, "ws_elev_m"], jb.loc[idx, "ws_elev_m"]
        _report("|Δ pelo de agua|", (zb - za).abs(), 3)
        rel = int((za.isna() & zb.notna()).sum())
        if rel:
            print(f"  {rel} perfil(es) pasaron de cota RELATIVA a absoluta")
        rel_new = int((zb.isna() & za.notna()).sum())
        if rel_new:
            print(f"  ATENCIÓN: {rel_new} perfil(es) pasaron de absoluta a RELATIVA")
    if "ws_source" in jb.columns:
        print("fuentes del pelo de agua (nueva corrida): "
              + ", ".join(f"{k} {v}" for k, v in
                          jb.loc[idx, "ws_source"].value_counts().items()))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(Path(sys.argv[1]), Path(sys.argv[2]))
