#!/usr/bin/env python3
"""Ignition-fraction forecast for a UNIFORM accretion boost.

A uniform multiplier C on the Bondi rate shifts the whole Eddington-ratio
distribution by a constant factor, so the fraction of BHs whose *raw*
(unsuppressed) ratio would cross a gate-entry threshold T after boosting is
simply

    f_ignite(C, T) = fraction of BHs with  (mdot/mdot_Edd)_unsuppressed > T / C .

This is the quantity that decides how often the LRD gate would *fire*. It is
NOT the sustained super-Eddington fraction -- once the gate fires and feedback
switches off, the runaway is set by the gated run, which cannot be forecast in
post-processing. Use this only to pick C so the gate ignites in a sensible
fraction (~1%, the LRD duty cycle) at z ~ 5-7, then RUN boost+gate.

No unit conversion is needed: the Eddington ratio is dimensionless.

Usage
-----
    python3 ignition_fraction.py LRD_base_free/snap_00*.hdf5
    python3 ignition_fraction.py --consts 10,30,100 --thresholds 1,3,10 snap_*.hdf5
"""

import argparse
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("This script needs h5py (pip install h5py).")


def load_edd(path):
    with h5py.File(path, "r") as f:
        z = float(np.atleast_1d(f["Header"].attrs.get("Redshift", 0.0))[0])
        n_bh = int(np.atleast_1d(f["Header"].attrs["NumPart_Total"])[5])
        if n_bh == 0 or "PartType5" not in f:
            return z, None
        bh = f["PartType5"]
        key = (
            "UnsuppressedEddingtonFractions"
            if "UnsuppressedEddingtonFractions" in bh
            else "EddingtonFractions"
        )
        used_unsuppressed = key == "UnsuppressedEddingtonFractions"
        edd = np.asarray(bh[key][:], dtype=float)
        edd = edd[np.isfinite(edd) & (edd > 0)]
        return z, (edd, used_unsuppressed)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("snapshots", nargs="+")
    ap.add_argument("--consts", default="10,30,100",
                    help="Uniform boost multipliers C to test (default 10,30,100).")
    ap.add_argument("--thresholds", default="1,10",
                    help="Gate-entry thresholds T in Eddington units (default 1,10).")
    args = ap.parse_args(argv)

    consts = [float(c) for c in args.consts.split(",") if c.strip()]
    thr = [float(t) for t in args.thresholds.split(",") if t.strip()]

    rows = []
    for path in args.snapshots:
        try:
            z, res = load_edd(path)
        except (OSError, KeyError) as e:
            sys.stderr.write(f"[skip] {path}: {e}\n")
            continue
        if res is None:
            continue
        edd, used = res
        if edd.size == 0:
            continue
        rows.append((z, edd, used, os.path.basename(path)))

    if not rows:
        sys.exit("No snapshots with black holes.")

    rows.sort(key=lambda r: -r[0])  # high z first

    warned = False
    for z, edd, used, name in rows:
        if not used and not warned:
            sys.stderr.write(
                "  [warn] a snapshot only has POST-CAP EddingtonFractions; "
                "ignition fractions are lower limits.\n"
            )
            warned = True
        print(f"\n=== {name}  z={z:.3f}   N_BH={edd.size}   "
              f"max(mdot/mdot_Edd)={edd.max():.3g} ===")
        header = "  C \\ T   " + "   ".join(f"T>{t:g}" for t in thr)
        print(header)
        for C in consts:
            cells = []
            for t in thr:
                frac = np.mean(edd > t / C)
                nbh = int(np.sum(edd > t / C))
                cells.append(f"{frac:7.4f}({nbh})")
            print(f"  C={C:<5g} " + "  ".join(cells))

    print(
        "\nHow to read: each cell is the fraction (count) of BHs that would "
        "cross the gate-entry threshold T after a uniform boost C.\n"
        "Pick C so that at z ~ 5-7 the ignition fraction at your chosen entry "
        "threshold is ~0.01 (the LRD duty cycle): too high over-ignites (every "
        "BH becomes super-Eddington -> BHMF over-production), too low never "
        "fires the gate. Then run boost (boost_alpha_only:1, boost_alpha=C) "
        "WITH the gate ON -- the sustained growth is the gated run's to decide."
    )


if __name__ == "__main__":
    main()
