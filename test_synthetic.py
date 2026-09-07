#!/usr/bin/env python3
"""Synthetic test for the GIS layer (RiverNetwork / stations / water surface).

No .mat and no real multi-river shapefile are needed: we fabricate a small
network in EPSG:5344-like planar metres —

    Neuquen main:  (0,0) --> (10000,0)            km0 at El Chanar (0,0)
      island:      anabranch (4000,0)->(5000,300)->(6000,0)  [north = left = MI]
    Negro   main:  (10000,0) --> (20000,0)        km0 at confluence (10000,0)
    Limay   main:  (10000,-8000) --> (10000,0)    drawn but NOT surveyed

the confluence node (10000,0) is shared exactly by the three mains. Water-surface
points and stations are placed only along Neuquen+Negro (the surveyed path), so
the survey chain must resolve to [Neuquen, Negro] and ignore the Limay.

Run:  python test_v5_synthetic.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Point

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_adcp_bathimetric as P

CRS = "EPSG:5344"
OK = True


def check(name, cond, detail=""):
    global OK
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))
    OK = OK and bool(cond)


def approx(a, b, tol=1.0):
    return abs(float(a) - float(b)) <= tol


def build_files(d: Path):
    # --- centerline (v5 schema) ---
    rows = [
        dict(seg_id=1, river="Neuquen", role="main",      label="",
             geometry=LineString([(0, 0), (10000, 0)])),
        dict(seg_id=2, river="Neuquen", role="anabranch", label="",
             geometry=LineString([(4000, 0), (5000, 300), (6000, 0)])),
        dict(seg_id=3, river="Negro",   role="main",      label="",
             geometry=LineString([(10000, 0), (20000, 0)])),
        dict(seg_id=4, river="Limay",   role="main",      label="",
             geometry=LineString([(10000, -8000), (10000, 0)])),
    ]
    cl = d / "centerline.shp"
    gpd.GeoDataFrame(rows, geometry="geometry", crs=CRS).to_file(cl)

    # --- stations (spot readings applied later) ---
    st = [
        dict(station_id="NEU-01", name="A", river="Neuquen", gauge_zero=290.0,
             geometry=Point(2000, 0)),
        dict(station_id="NEU-02", name="B", river="Neuquen", gauge_zero=285.0,
             geometry=Point(8000, 0)),
        dict(station_id="NEG-01", name="C", river="Negro",   gauge_zero=280.0,
             geometry=Point(15000, 0)),
    ]
    stp = d / "stations.shp"
    gpd.GeoDataFrame(st, geometry="geometry", crs=CRS).to_file(stp)

    # --- readings (spot) ---
    rd = d / "readings_spot.csv"
    rd.write_text("station_id,nivel\nNEU-01,1.00\nNEU-02,1.00\nNEG-01,1.00\n",
                  encoding="utf-8")

    # --- centerline with NO schema (v4 fallback) ---
    cl0 = d / "centerline_noschema.shp"
    gpd.GeoDataFrame([dict(id=7, geometry=LineString([(0, 0), (5000, 0)]))],
                     geometry="geometry", crs=CRS).to_file(cl0)
    return cl, stp, rd, cl0


def main():

    print(f"=== test sintético — pipeline v{P.__version__} ===")
    tmp = Path(tempfile.mkdtemp())
    cl, stp, rd, cl0 = build_files(tmp)

    # water-surface points along the surveyed path (planar x == survey chainage)
    ws_xy_elev = [(float(x), 0.0, 300.0 - 0.001 * x) for x in range(0, 20001, 1000)]

    log = []
    net = P.RiverNetwork(cl, crs_override=None, reverse=False,
                         river_offsets={"Neuquen": 0.0, "Negro": 0.0, "Limay": 0.0},
                         ws_xy_elev=ws_xy_elev, log=log)

    print("\n=== RiverNetwork ===")
    rivers = sorted(ax.river for ax in net.axes)
    check("3 river axes (Limay, Negro, Neuquen)", rivers == ["Limay", "Negro", "Neuquen"], str(rivers))

    nqn = net._by_river["Neuquen"]
    check("Neuquen main length ~10000 m", approx(nqn.length, 10000, 5), f"{nqn.length:.1f}")
    check("Neuquen has 1 anabranch", len(nqn.anabranches) == 1, str(len(nqn.anabranches)))
    check("Neuquen has 1 island", len(nqn.islands) == 1, str(len(nqn.islands)))
    ab_label = nqn.anabranches[0].get("label")
    main_lab = nqn.islands[0]["main_label"] if nqn.islands else None
    check("north anabranch labelled MI (left facing downstream)", ab_label == "MI", str(ab_label))
    check("main sub-reach in island labelled MD", main_lab == "MD", str(main_lab))

    print("\n=== locate() ===")
    L = net.locate(2000, 0)
    check("point on Neuquen main -> river Neuquen", L["river"] == "Neuquen", L["river"])
    check("  km_internal ~2000", approx(L["km_internal"], 2000, 2), f"{L['km_internal']:.1f}")
    check("  km_oficial ~2000 (offset 0)", approx(L["km_oficial"], 2000, 2), f"{L['km_oficial']:.1f}")
    check("  brazo empty on main outside island", L["brazo"] == "", repr(L["brazo"]))

    La = net.locate(5000, 300)
    check("point on north anabranch -> brazo MI", La["brazo"] == "MI", La["brazo"])
    check("  anabranch chainage projected onto main ~5000",
          approx(La["km_internal"], 5000, 60), f"{La['km_internal']:.1f}")

    Lm = net.locate(5000, 0)   # main sub-reach inside the island span
    check("point on main inside island -> brazo MD", Lm["brazo"] == "MD", Lm["brazo"])

    Ln = net.locate(15000, 0)
    check("point past confluence -> river Negro", Ln["river"] == "Negro", Ln["river"])
    check("  Negro km_internal ~5000", approx(Ln["km_internal"], 5000, 2), f"{Ln['km_internal']:.1f}")

    print("\n=== survey chain ===")
    chain = net.build_survey_chain(["Neuquen", "Negro"])
    check("chain = [Neuquen, Negro]", chain == ["Neuquen", "Negro"], str(chain))
    check("survey_ok True", net.survey_ok is True)
    check("Negro survey offset ~10000 (= Neuquen length)",
          approx(net.survey_offsets.get("Negro", -1), 10000, 5),
          f"{net.survey_offsets.get('Negro'):.1f}")
    sc_neg = net.survey_chainage("Negro", 5000)
    check("survey_chainage(Negro, 5000) ~15000 (continuous, monotonic)",
          approx(sc_neg, 15000, 5), f"{sc_neg:.1f}")
    sc_nqn = net.survey_chainage("Neuquen", 2000)
    check("survey_chainage(Neuquen, 2000) ~2000", approx(sc_nqn, 2000, 2), f"{sc_nqn:.1f}")
    check("continuity: Negro@0 == Neuquen@end", approx(net.survey_chainage("Negro", 0), nqn.length, 5))

    print("\n=== primary water-surface profile (continuous chainage) ===")
    pairs = []
    for (x, y, h) in ws_xy_elev:
        lo = net.locate(x, y)
        s = net.survey_chainage(lo["river"], lo["km_internal"])
        if s is not None:
            pairs.append((s, h, lo["dist"]))
    wsp = P.WaterSurfaceProfile(pairs, log=log)
    check("WSP covers survey chainage 15000", wsp.covers(15000) is True)
    e15, note = wsp.elev_at(15000)
    check("WSP elev @15000 ~285.0", approx(e15, 285.0, 0.2), f"{e15:.3f} note={note!r}")
    check("WSP flags extrapolation past the range", wsp.elev_at(25000)[1] != "",
          wsp.elev_at(25000)[1])
    check("WSP slope negative (drops downstream)", wsp.slope_m_per_km < 0,
          f"{wsp.slope_m_per_km:.2f} m/km")

    print("\n=== stations (secondary WS: interp / single / fallback) ===")
    reg = P.StationRegistry(stp, net, crs_override=CRS, log=log)
    reg.attach_survey_chainage(net)
    mode, data = P.read_level_readings(rd)
    check("readings auto-detected as spot", mode == "spot", mode)
    sws = P.StationWS(reg, mode, data, log=log)
    e_mid, n_mid = sws.elev_at(5000, "Neuquen")
    # NEU-01(sc2000, ws291) .. NEU-02(sc8000, ws286): at 5000 -> 288.5
    check("station interp on Neuquen @5000 ~288.5", approx(e_mid, 288.5, 0.05),
          f"{e_mid:.3f}")
    e_neg, n_neg = sws.elev_at(15000, "Negro")
    check("single Negro station -> ws 281.0 with note", approx(e_neg, 281.0, 0.05) and n_neg,
          f"{e_neg:.3f} note={n_neg!r}")
    e_none, n_none = sws.elev_at(5000, "Limay")
    check("no station on Limay -> None", e_none is None, f"{e_none} note={n_none!r}")

    print("\n=== QC cross-check magnitude (GNSS vs station) ===")
    gnss5000 = wsp.elev_at(5000)[0]              # ~295.0
    check("|WS_gnss - WS_station| @5000 reasonable",
          abs(gnss5000 - e_mid) < 10.0, f"gnss={gnss5000:.3f} stn={e_mid:.3f}")

    print("\n=== merged beam-cloud export ===")
    # two fake transects with tiny clouds (VB + one slant each)
    def fake_result(perfil, river, kmo, sc, ws):
        return dict(perfil=perfil, river=river, role="main", brazo="",
                    dist_axis=float("nan"),
                    km_internal=kmo, km_oficial=kmo, survey_chainage=sc,
                    ws_elev=ws, ws_note="", ws_source="survey",
                    s_grid=np.array([0.0, 1.0]), depth_grid=np.array([0.0, 2.0]),
                    bed_elev=np.array([ws, ws - 2.0]),
                    xp=np.array([0.0, 1.0]), yp=np.array([0.0, 0.0]),
                    width_m=1.0, thalweg_elev=ws - 2.0, max_depth=2.0,
                    theta_flow=0.0, n_samples=2, utm_crs=None, posgar_crs=None, log=[],
                    cloud_xp=np.array([100.0 + kmo, 101.0 + kmo]),
                    cloud_yp=np.array([0.0, 0.5]),
                    cloud_z=np.array([ws - 3.0, ws - 2.5]),
                    cloud_depth=np.array([3.0, 2.5]),
                    cloud_beam_index=np.array([0, 1]),
                    cloud_ens=np.array([0, 0]), cloud_rel=False)

    results = [fake_result("T1", "Neuquen", 2000, 2000, 291.0),
               fake_result("T2", "Negro", 5000, 15000, 285.0)]
    idx = P.write_survey_index(results, tmp, "SRVN16")
    check("survey_index.csv written", idx.exists())
    header = idx.read_text(encoding="utf-8").splitlines()[0]
    check("index has river + km_oficial + ws_source cols",
          all(c in header for c in ("river", "km_oficial_m", "ws_source")), header)
    cloud = P.export_survey_beam_cloud(results, tmp)
    check("beam cloud file written", cloud is not None and cloud.exists(),
          str(cloud.name) if cloud else "None")
    g = gpd.read_file(cloud)
    check("beam cloud has 4 points (2 transects x 2 beams)", len(g) == 4, str(len(g)))
    check("beam cloud attrs present",
          all(c in g.columns for c in ("pt_id", "transect", "ens", "beam_id",
                                       "beam_type", "bed_elev", "depth", "river",
                                       "ws_elev", "flag")),
          str(list(g.columns)))
    check("beam cloud excludes km", not any("km" in c.lower() for c in g.columns),
          str(list(g.columns)))
    types = set(g["beam_type"])
    check("beam_type has vert + slant", types == {"vert", "slant"}, str(types))
    check("Z carried as 3D geometry", bool(g.geometry.iloc[0].has_z))

    print("\n=== v4 fallback (no river/role columns) ===")
    log2 = []
    net0 = P.RiverNetwork(cl0, crs_override=None, log=log2)
    check("single default river axis", len(net0.axes) == 1, str(len(net0.axes)))
    check("default river name applied", net0.axes[0].river == P.DEFAULT_RIVER_NAME,
          net0.axes[0].river)
    L0 = net0.locate(2500, 0)
    check("v4-mode locate km_internal ~2500", approx(L0["km_internal"], 2500, 2),
          f"{L0['km_internal']:.1f}")

    print("\n" + ("ALL PASS ✅" if OK else "SOME FAILED ❌"))
    return 0 if OK else 1

if __name__ == "__main__":
    raise SystemExit(main())
