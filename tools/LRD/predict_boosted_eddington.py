#!/usr/bin/env python3
"""Predict the *boosted* Eddington-ratio distribution from baseline snapshots.

Post-processing companion to ``measure_bh_distributions.py`` ("step C"):
before running LRD_base_boost, apply the density-dependent Booth & Schaye
(2009) boost to the *recorded* unsuppressed Eddington ratios of a baseline
(no-boost, uncapped) run:

    alpha(n_H) = 1                      for n_H <= n_H*
               = (n_H / n_H*)**beta     for n_H >  n_H*      [B&S09 eq. 4]

    mdot_boosted / mdot_Edd = alpha(n_H) * (mdot / mdot_Edd)_unsuppressed

where n_H is the *volume* hydrogen number density at the BH (from
``GasDensities`` -- the same rho_gas the Bondi rate uses), NOT the column
proxy. This tells us, for beta in e.g. {1, 1.5, 2}, whether a meaningful
sub-population (say the top 10-30% in n_H) would cross ~10x Eddington at
z = 4-8, as D. Schleicher requested ("super-Eddington in a relevant way").

Caveat: this is NOT self-consistent -- in a real boosted run the extra
feedback (or, with the LRD gate, its absence) changes rho and c_s. It is a
*starting-point selector* for beta, nothing more.

Usage
-----
    python3 predict_boosted_eddington.py snap_*.hdf5
    python3 predict_boosted_eddington.py --betas 1,1.5,2 --plot \
        --outdir figs_boost_forecast/ snap_000[5-9].hdf5

Requires: numpy, h5py; matplotlib only with --plot.
"""

import argparse
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("This script needs h5py (pip install h5py).")

M_PROTON_CGS = 1.67262192369e-24

# Thresholds (in Eddington units) at which we report flagged fractions.
EDD_THRESHOLDS = [1.0, 3.0, 10.0, 30.0, 100.0]
# n_H quantile cuts defining the "dense" sub-populations.
DENSE_QUANTILES = [0.70, 0.90]  # top 30%, top 10%


def _read_physical_cgs(dset, scale_factor):
    """Dataset -> physical CGS, using SWIFT's per-dataset attributes.

    Tries, in order:
      1. the direct physical-CGS factor (includes cosmological corrections);
      2. a comoving-CGS factor times a**(a-scale exponent), accepting the
         attribute-name variants used by different SWIFT versions.
    Raises if neither is found -- code units would poison everything.
    """
    raw = dset[:]
    attrs = dset.attrs

    key_phys = (
        "Conversion factor to physical CGS (including cosmological corrections)"
    )
    if key_phys in attrs:
        return raw * float(np.atleast_1d(attrs[key_phys])[0])

    cgs_keys = [
        "Conversion factor to CGS (not including cosmological corrections)",
        "Conversion factor to CGS (including cosmological corrections)",
    ]
    a_keys = ["a-scale exponent", "a-scale exponent of the field"]
    cgs = next(
        (float(np.atleast_1d(attrs[k])[0]) for k in cgs_keys if k in attrs),
        None,
    )
    a_exp = next(
        (float(np.atleast_1d(attrs[k])[0]) for k in a_keys if k in attrs),
        None,
    )
    if cgs is not None and a_exp is not None:
        return raw * cgs * (scale_factor ** a_exp)

    sys.stderr.write(
        f"  [ERROR] {dset.name}: no usable unit attributes found.\n"
        f"          attributes present: {sorted(map(str, attrs.keys()))}\n"
        "          Values would stay in CODE units and every downstream\n"
        "          number would be wrong.\n"
    )
    raise KeyError(f"unit attributes missing on {dset.name}")


def _get(group, *names):
    for n in names:
        if n in group:
            return group[n]
    return None


def load_snapshot(path, hydrogen_fraction):
    """Per-BH subgrid mass, unsuppressed Edd ratio and volume n_H [cm^-3]."""
    with h5py.File(path, "r") as f:
        header = f["Header"].attrs
        a = float(np.atleast_1d(header.get("Scale-factor", 1.0))[0])
        z = float(np.atleast_1d(header.get("Redshift", 0.0))[0])

        n_bh = int(np.atleast_1d(header["NumPart_Total"])[5])
        if n_bh == 0 or "PartType5" not in f:
            return None
        bh = f["PartType5"]

        msun_cgs = 1.98841e33
        m_dset = _get(bh, "SubgridMasses", "DynamicalMasses", "Masses")
        m_sub = _read_physical_cgs(m_dset, a) / msun_cgs

        edd_dset = _get(bh, "UnsuppressedEddingtonFractions")
        used_unsuppressed = edd_dset is not None
        if edd_dset is None:
            edd_dset = _get(bh, "EddingtonFractions")
            if edd_dset is None:
                return None
        edd = np.asarray(edd_dset[:], dtype=float)

        rho_dset = _get(bh, "GasDensities")
        if rho_dset is None:
            sys.stderr.write(f"  [warn] {path}: no GasDensities; skipped.\n")
            return None
        rho_cgs = _read_physical_cgs(rho_dset, a)  # g cm^-3, physical
        n_h_vol = hydrogen_fraction * rho_cgs / M_PROTON_CGS  # cm^-3

    return {
        "path": path,
        "z": z,
        "m_sub_Msun": m_sub,
        "edd": edd,
        "used_unsuppressed": used_unsuppressed,
        "n_h_vol": n_h_vol,
    }


def boost_factor(n_h_vol, beta, n_h_star, boost_max=None):
    """Density-dependent B&S09 boost (floor 1), computed in log space.

    The log-space form cannot overflow; log10(alpha) is clipped at 30 (any
    boost that large means the input units are wrong anyway).
    """
    n = np.asarray(n_h_vol, dtype=float)
    log_alpha = np.zeros_like(n)
    hi = n > n_h_star
    log_alpha[hi] = beta * np.log10(n[hi] / n_h_star)
    log_alpha = np.minimum(log_alpha, 30.0)
    alpha = 10.0 ** log_alpha
    if boost_max is not None:
        alpha = np.minimum(alpha, boost_max)
    return alpha


def _frac_line(edd_boosted, sel, label):
    x = edd_boosted[sel]
    x = x[np.isfinite(x)]
    if x.size == 0:
        return f"      {label:<22s} (empty)"
    parts = [
        f"f(>{t:g})={np.mean(x > t):.3g}({int(np.sum(x > t))})"
        for t in EDD_THRESHOLDS
    ]
    return f"      {label:<22s} n={x.size:<5d} " + "  ".join(parts)


def report_snapshot(d, betas, n_h_star, boost_max):
    print(f"\n=== {os.path.basename(d['path'])}  (z={d['z']:.3f}) ===")
    edd = d["edd"]
    n_h = d["n_h_vol"]
    good = np.isfinite(edd) & np.isfinite(n_h) & (edd > 0)
    edd, n_h = edd[good], n_h[good]
    print(f"  black holes used: {edd.size}")
    if not d["used_unsuppressed"]:
        print(
            "  [warn] using POST-CAP EddingtonFractions (no unsuppressed "
            "field); boosted values are lower limits."
        )
    if edd.size == 0:
        return None

    q = np.percentile(n_h, [50, 70, 90, 99])
    print(
        "  n_H at BHs [cm^-3]: p50={:.3g} p70={:.3g} p90={:.3g} p99={:.3g}"
        .format(*q)
    )
    if q[0] > 1e6:
        print(
            "  [ERROR] median n_H > 1e6 cm^-3: physically impossible at "
            "kernel scale.\n"
            "          Unit conversion has failed -- do NOT use these "
            "boosted numbers."
        )
        return None
    print(f"  unboosted max(mdot/mdot_Edd) = {edd.max():.3g}")

    # Dense-subset masks (defined on n_H, fixed across betas).
    subsets = [("all BHs", np.ones_like(edd, dtype=bool))]
    for qc in DENSE_QUANTILES:
        cut = np.quantile(n_h, qc)
        subsets.append((f"top {100 * (1 - qc):.0f}% n_H", n_h >= cut))

    rows = {}
    for beta in betas:
        alpha = boost_factor(n_h, beta, n_h_star, boost_max)
        eb = edd * alpha
        print(
            f"  beta={beta:g}: median boost={np.median(alpha):.3g}, "
            f"p90 boost={np.percentile(alpha, 90):.3g}, "
            f"max boosted ratio={eb.max():.3g}"
        )
        for label, sel in subsets:
            print(_frac_line(eb, sel, label))
        # store f(>10 Edd) in the top-30% subset for the summary plot
        sel30 = subsets[1][1]
        rows[beta] = {
            "f10_all": float(np.mean(eb > 10.0)),
            "f10_top30": float(np.mean(eb[sel30] > 10.0)) if sel30.any() else 0.0,
            "f1_all": float(np.mean(eb > 1.0)),
        }
    return rows


def make_plots(snapshots, results, betas, n_h_star, boost_max, outdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)

    # 1) Boosted-PDF overlay, one panel per snapshot.
    for d in snapshots:
        edd = d["edd"]
        n_h = d["n_h_vol"]
        good = np.isfinite(edd) & np.isfinite(n_h) & (edd > 0)
        edd, n_h = edd[good], n_h[good]
        if edd.size < 3:
            continue
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(
            np.log10(edd), bins=40, histtype="step", density=True,
            color="k", lw=1.5, label="no boost",
        )
        for beta in betas:
            eb = edd * boost_factor(n_h, beta, n_h_star, boost_max)
            lg = np.log10(eb)
            lg = lg[np.isfinite(lg)]
            if lg.size == 0:
                continue
            ax.hist(
                lg, bins=40, histtype="step", density=True,
                label=f"beta={beta:g}",
            )
        ax.axvline(0.0, ls=":", c="gray", lw=1)
        ax.axvline(1.0, ls="--", c="k", lw=1, label="10 x Edd")
        ax.set_xlabel(r"$\log_{10}(\dot m/\dot m_{\rm Edd})$ (forecast, boosted)")
        ax.set_ylabel("PDF")
        ax.set_title(f"z={d['z']:.2f}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        p = os.path.join(outdir, f"boost_forecast_z{d['z']:.2f}.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        print(f"  wrote {p}")

    # 2) Summary: f(>10x Edd) vs z per beta (all BHs and top-30% n_H).
    zs = [d["z"] for d in snapshots]
    fig, ax = plt.subplots(figsize=(7, 5))
    for beta in betas:
        f_all = [results[d["path"]][beta]["f10_all"] for d in snapshots]
        f_top = [results[d["path"]][beta]["f10_top30"] for d in snapshots]
        (line,) = ax.plot(zs, f_top, "o-", label=f"beta={beta:g} (top 30% n_H)")
        ax.plot(zs, f_all, "s--", color=line.get_color(), alpha=0.5,
                label=f"beta={beta:g} (all)")
    ax.set_xlabel("z")
    ax.set_ylabel(r"fraction with boosted $\dot m > 10\,\dot m_{\rm Edd}$")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, "f_above_10Edd_vs_z.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"  wrote {p}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("snapshots", nargs="+", help="SWIFT snapshot HDF5 file(s).")
    ap.add_argument(
        "--betas", default="1,1.5,2",
        help="Comma-separated boost_beta values to forecast (default 1,1.5,2).",
    )
    ap.add_argument(
        "--n-h-star", type=float, default=0.1,
        help="Boost threshold density n_H* [cm^-3] (default 0.1, B&S09).",
    )
    ap.add_argument(
        "--boost-max", type=float, default=None,
        help="Optional ceiling on the boost factor (forecast a capped boost).",
    )
    ap.add_argument(
        "--hydrogen-fraction", type=float, default=0.752,
        help="Hydrogen mass fraction X_H (default 0.752; SWIFT boost uses 0.75).",
    )
    ap.add_argument("--plot", action="store_true", help="Write PNG plots.")
    ap.add_argument("--outdir", default="figs_boost_forecast",
                    help="Output dir for --plot.")
    args = ap.parse_args(argv)

    betas = [float(b) for b in args.betas.split(",") if b.strip()]

    snapshots, results = [], {}
    for path in args.snapshots:
        try:
            d = load_snapshot(path, args.hydrogen_fraction)
        except (OSError, KeyError) as e:
            sys.stderr.write(f"[skip] {path}: {e}\n")
            continue
        if d is None:
            sys.stderr.write(f"[skip] {path}: no BHs / no usable fields.\n")
            continue
        rows = report_snapshot(d, betas, args.n_h_star, args.boost_max)
        if rows is None:
            continue
        snapshots.append(d)
        results[d["path"]] = rows

    if not snapshots:
        sys.exit("No usable snapshots.")

    if args.plot:
        make_plots(snapshots, results, betas, args.n_h_star,
                   args.boost_max, args.outdir)

    print(
        "\nGuidance: pick the smallest beta for which the top-30% n_H subset "
        "shows a meaningful f(>10 x Edd) at z=4-8 (say >~0.1) while the "
        "overall max ratio stays short of absurd (<~1e3). If even beta=2 "
        "undershoots, raise beta or lower n_H*; if beta=2 overshoots wildly, "
        "consider --boost-max. Remember: no feedback back-reaction here -- "
        "this only selects the starting beta for LRD_base_boost."
    )


if __name__ == "__main__":
    main()