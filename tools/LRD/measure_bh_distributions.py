#!/usr/bin/env python3
"""Measure the black-hole accretion-rate and column-density distributions in
SWIFT SPIN_JET snapshots, to guide the choice of the LRD ("Little Red Dot")
confined-feedback threshold.

This is the "step C" diagnostic for the LRD confined-feedback project: before
committing compute time to gated runs, we want to know, in a *baseline* run
(no gating), the range of

  * the raw, unsuppressed Bondi-based Eddington ratio  (mdot / mdot_Edd),
  * the *volume* hydrogen number density n_H at the BH (cm^-3) -- this is the
    quantity the Booth & Schaye boost threshold n_H* acts on, and
  * a proxy nuclear column density  N_H (cm^-2),

as a function of redshift and black-hole mass, and how tightly they correlate.

IMPORTANT about the Eddington ratio
-----------------------------------
The standard ``EddingtonFractions`` field written by SWIFT is the *post-cap*
value (clipped at ``max_eddington_fraction``). To see the super-Eddington tail
that sets the LRD threshold you must use either:

  * the ``UnsuppressedEddingtonFractions`` field (added for this project), or
  * a baseline run with ``max_eddington_fraction`` set very large.

Density / column proxy
----------------------
Volume density:  n_H = X_H * rho_gas_phys / m_p                    [cm^-3]
Column proxy:    N_H = X_H * rho_gas_phys / m_p * L                [cm^-2]
with L the line-of-sight length (default kernel_gamma*h). n_H is directly
comparable to the boost threshold n_H* (default 0.1 cm^-3); N_H is a
kernel-scale proxy, not a face-value column.

Usage
-----
    python3 measure_bh_distributions.py snap_0000.hdf5 [snap_0001.hdf5 ...]
    python3 measure_bh_distributions.py --plot --outdir figs/ snap_*.hdf5

Requires: numpy, h5py. matplotlib only if --plot is given.
"""

import argparse
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("This script needs h5py (pip install h5py).")

# Proton mass in cgs.
M_PROTON_CGS = 1.67262192369e-24
# Kernel-support / smoothing-length ratio for Wendland-C2 in 3D (kernel_gamma).
KERNEL_GAMMA_WENDLAND_C2_3D = 1.936492

# Candidate accretion-rate thresholds (Eddington units).
EDD_THRESHOLDS = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
# Candidate column thresholds (cm^-2); 1.5e24 is the Compton-thick boundary.
NH_THRESHOLDS = [1e22, 1e23, 1.5e24, 1e25]
# Candidate *volume* density thresholds (cm^-3) bracketing plausible n_H* values.
NHSTAR_CANDIDATES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]


def _read_physical_cgs(dset, scale_factor):
    """Dataset -> physical CGS, using SWIFT's per-dataset attributes.

    Tries, in order:
      1. the direct physical-CGS factor (already includes cosmological terms);
      2. a comoving-CGS factor times a**(a-scale exponent), accepting the
         attribute-name variants used by different SWIFT versions.
    Raises if neither is found -- code units would poison everything downstream.
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
        "          Values would stay in CODE units -- aborting this field.\n"
    )
    raise KeyError(f"unit attributes missing on {dset.name}")


def _get(group, *names):
    """Return the first dataset among ``names`` that exists, else None."""
    for n in names:
        if n in group:
            return group[n]
    return None


def load_snapshot(path, column_length, hydrogen_fraction):
    """Load the per-BH quantities of interest from one snapshot."""
    with h5py.File(path, "r") as f:
        header = f["Header"].attrs
        a = float(np.atleast_1d(header.get("Scale-factor", 1.0))[0])
        z = float(np.atleast_1d(header.get("Redshift", 0.0))[0])

        n_bh = int(np.atleast_1d(header["NumPart_Total"])[5])
        if n_bh == 0 or "PartType5" not in f:
            return None
        bh = f["PartType5"]

        msun_cgs = 1.98841e33

        m_sub_dset = _get(bh, "SubgridMasses", "DynamicalMasses", "Masses")
        m_sub = _read_physical_cgs(m_sub_dset, a) / msun_cgs

        edd_dset = _get(bh, "UnsuppressedEddingtonFractions")
        used_unsuppressed = edd_dset is not None
        if edd_dset is None:
            edd_dset = _get(bh, "EddingtonFractions")
        edd = edd_dset[:] if edd_dset is not None else None

        rho_dset = _get(bh, "GasDensities")
        h_dset = _get(bh, "SmoothingLengths")
        n_h_vol = None  # volume number density [cm^-3]
        n_h = None       # column proxy [cm^-2]
        if rho_dset is not None:
            rho_cgs = _read_physical_cgs(rho_dset, a)  # g cm^-3 physical
            n_h_vol = hydrogen_fraction * rho_cgs / M_PROTON_CGS
            if h_dset is not None:
                h_cgs = _read_physical_cgs(h_dset, a)  # cm
                if column_length == "support":
                    length = KERNEL_GAMMA_WENDLAND_C2_3D * h_cgs
                elif column_length == "h":
                    length = h_cgs
                else:  # fixed physical length in kpc
                    length = np.full_like(
                        h_cgs, float(column_length) * 3.0857e21
                    )
                n_h = n_h_vol * length

        flags_dset = _get(bh, "LRDConfinedFlags")
        flags = flags_dset[:].astype(bool) if flags_dset is not None else None

        conf_dset = _get(bh, "LRDConfinedMasses")
        conf_mass = (
            _read_physical_cgs(conf_dset, a) / msun_cgs
            if conf_dset is not None
            else None
        )

        mode_dset = _get(bh, "AccretionModes")
        modes = mode_dset[:] if mode_dset is not None else None

    return {
        "path": path,
        "a": a,
        "z": z,
        "m_sub_Msun": m_sub,
        "edd": edd,
        "used_unsuppressed": used_unsuppressed,
        "n_h_vol": n_h_vol,
        "n_h": n_h,
        "flags": flags,
        "conf_mass_Msun": conf_mass,
        "modes": modes,
    }


def _pct(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None
    q = np.percentile(x, [1, 10, 50, 90, 99])
    return q, x.min(), x.max(), x.size


def _fmt_pct(label, res, fmt="{:.3g}"):
    if res is None:
        return f"  {label:<28s} (no finite values)"
    q, lo, hi, n = res
    body = "  ".join(fmt.format(v) for v in q)
    return (
        f"  {label:<28s} n={n:<6d} min={fmt.format(lo)}  "
        f"[p1 p10 p50 p90 p99]= {body}  max={fmt.format(hi)}"
    )


def report_snapshot(d):
    print(f"\n=== {os.path.basename(d['path'])}  (z={d['z']:.3f}, a={d['a']:.4f}) ===")
    n = d["m_sub_Msun"].size
    print(f"  black holes: {n}")

    print(_fmt_pct("subgrid mass [Msun]", _pct(d["m_sub_Msun"])))

    if d["edd"] is not None:
        tag = "unsuppressed" if d["used_unsuppressed"] else "POST-CAP (fallback!)"
        print(_fmt_pct(f"Eddington ratio ({tag})", _pct(d["edd"])))
        if not d["used_unsuppressed"]:
            print(
                "  [warn] Only post-cap EddingtonFractions available; the "
                "super-Eddington tail is clipped."
            )
        edd = np.asarray(d["edd"], dtype=float)
        edd = edd[np.isfinite(edd)]
        if edd.size:
            print("  fraction above Eddington-ratio thresholds:")
            for t in EDD_THRESHOLDS:
                frac = np.mean(edd > t)
                print(f"      f(>{t:>6g} x Edd) = {frac:.4g}  ({int(np.sum(edd > t))} BHs)")

    # Volume number density -- the quantity the boost threshold n_H* acts on.
    if d["n_h_vol"] is not None:
        print(_fmt_pct("volume n_H [cm^-3]", _pct(d["n_h_vol"]), fmt="{:.3e}"))
        nhv = np.asarray(d["n_h_vol"], dtype=float)
        nhv = nhv[np.isfinite(nhv)]
        if nhv.size:
            if np.median(nhv) > 1e6:
                print(
                    "  [ERROR] median volume n_H > 1e6 cm^-3: impossible at "
                    "kernel scale -> unit conversion still wrong for this file."
                )
            print("  fraction above candidate boost thresholds n_H*:")
            for t in NHSTAR_CANDIDATES:
                frac = np.mean(nhv > t)
                print(
                    f"      f(n_H > {t:>5g} cm^-3) = {frac:.4g}  "
                    f"({int(np.sum(nhv > t))} BHs)  <- would be boosted if n_H*={t}"
                )

    if d["n_h"] is not None:
        print(_fmt_pct("proxy column N_H [cm^-2]", _pct(d["n_h"]), fmt="{:.3e}"))
        nh = np.asarray(d["n_h"], dtype=float)
        nh = nh[np.isfinite(nh)]
        if nh.size:
            print("  fraction above column thresholds:")
            for t in NH_THRESHOLDS:
                frac = np.mean(nh > t)
                print(f"      f(N_H > {t:.1e}) = {frac:.4g}  ({int(np.sum(nh > t))} BHs)")

    # Correlation between the two candidate criteria.
    if d["edd"] is not None and d["n_h"] is not None:
        edd = np.asarray(d["edd"], dtype=float)
        nh = np.asarray(d["n_h"], dtype=float)
        good = np.isfinite(edd) & np.isfinite(nh) & (edd > 0) & (nh > 0)
        if good.sum() > 2:
            r = _spearman(edd[good], nh[good])
            print(
                f"  Spearman rho(Eddington ratio, N_H) = {r:+.3f}  "
                f"(n={good.sum()})  -> high |rho| => one criterion may suffice"
            )

    if d["flags"] is not None:
        n_lrd = int(d["flags"].sum())
        print(f"  LRD-flagged BHs this snapshot: {n_lrd} / {n}  ({n_lrd / n:.3g})")
    if d["conf_mass_Msun"] is not None:
        tot = d["m_sub_Msun"]
        conf = d["conf_mass_Msun"]
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(tot > 0, conf / tot, 0.0)
        print(
            f"  confined-grown mass fraction: median={np.median(frac):.3g} "
            f"max={np.max(frac):.3g}"
        )


def _spearman(x, y):
    """Spearman rank correlation without SciPy."""
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else float("nan")


def make_plots(snapshots, outdir, n_h_star):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)

    # 1) Eddington-ratio distribution per snapshot.
    fig, ax = plt.subplots(figsize=(7, 5))
    for d in snapshots:
        if d["edd"] is None:
            continue
        edd = np.asarray(d["edd"], dtype=float)
        edd = edd[np.isfinite(edd) & (edd > 0)]
        if edd.size == 0:
            continue
        ax.hist(np.log10(edd), bins=40, histtype="step", density=True,
                label=f"z={d['z']:.2f}")
    ax.set_xlabel(r"$\log_{10}(\dot m / \dot m_{\rm Edd})$ (unsuppressed)")
    ax.set_ylabel("PDF")
    ax.axvline(0.0, ls="--", c="k", lw=1, label="Eddington")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    p = os.path.join(outdir, "eddington_ratio_pdf.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"  wrote {p}")

    # 2) Volume n_H distribution per snapshot, with the boost threshold marked.
    fig, ax = plt.subplots(figsize=(7, 5))
    plotted = False
    for d in snapshots:
        if d["n_h_vol"] is None:
            continue
        nhv = np.asarray(d["n_h_vol"], dtype=float)
        nhv = nhv[np.isfinite(nhv) & (nhv > 0)]
        if nhv.size == 0:
            continue
        ax.hist(np.log10(nhv), bins=40, histtype="step", density=True,
                label=f"z={d['z']:.2f}")
        plotted = True
    if plotted:
        ax.axvline(np.log10(n_h_star), ls="--", c="k", lw=1.5,
                   label=f"n_H* = {n_h_star:g}")
        ax.set_xlabel(r"$\log_{10}(n_{\rm H}\,/\,{\rm cm^{-3}})$ at BH (volume)")
        ax.set_ylabel("PDF")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        p = os.path.join(outdir, "nH_volume_pdf.png")
        fig.savefig(p, dpi=130)
        print(f"  wrote {p}")
    plt.close(fig)

    # 3) Fraction of BHs above n_H* vs redshift.
    zs, fracs = [], []
    for d in snapshots:
        if d["n_h_vol"] is None:
            continue
        nhv = np.asarray(d["n_h_vol"], dtype=float)
        nhv = nhv[np.isfinite(nhv)]
        if nhv.size == 0:
            continue
        zs.append(d["z"])
        fracs.append(np.mean(nhv > n_h_star))
    if zs:
        order = np.argsort(zs)
        zs = np.array(zs)[order]
        fracs = np.array(fracs)[order]
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(zs, fracs, "o-")
        ax.set_xlabel("z")
        ax.set_ylabel(fr"fraction of BHs with $n_{{\rm H}} > {n_h_star:g}$ cm$^{{-3}}$")
        ax.invert_xaxis()
        fig.tight_layout()
        p = os.path.join(outdir, "frac_above_nHstar_vs_z.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        print(f"  wrote {p}")

    # 4) Eddington ratio vs proxy column (last snapshot with both).
    for d in reversed(snapshots):
        if d["edd"] is None or d["n_h"] is None:
            continue
        edd = np.asarray(d["edd"], dtype=float)
        nh = np.asarray(d["n_h"], dtype=float)
        good = np.isfinite(edd) & np.isfinite(nh) & (edd > 0) & (nh > 0)
        if good.sum() < 3:
            continue
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        ax.scatter(np.log10(edd[good]), np.log10(nh[good]), s=8, alpha=0.5)
        ax.axhline(np.log10(1.5e24), ls="--", c="r", lw=1, label="Compton-thick")
        ax.axvline(0.0, ls="--", c="k", lw=1, label="Eddington")
        ax.set_xlabel(r"$\log_{10}(\dot m / \dot m_{\rm Edd})$")
        ax.set_ylabel(r"$\log_{10}(N_{\rm H}\,/\,{\rm cm^{-2}})$")
        ax.set_title(f"z={d['z']:.2f}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        p = os.path.join(outdir, "edd_vs_column.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        print(f"  wrote {p}")
        break


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("snapshots", nargs="+", help="SWIFT snapshot HDF5 file(s).")
    ap.add_argument(
        "--column-length", default="support",
        help="LOS length for N_H proxy: 'support' (kernel_gamma*h), 'h', "
        "or a fixed physical length in kpc.",
    )
    ap.add_argument(
        "--hydrogen-fraction", type=float, default=0.752,
        help="Hydrogen mass fraction X_H (default 0.752).",
    )
    ap.add_argument(
        "--n-h-star", type=float, default=0.1,
        help="Boost threshold density n_H* [cm^-3] to mark in plots "
        "(default 0.1, B&S09).",
    )
    ap.add_argument("--plot", action="store_true", help="Write summary PNG plots.")
    ap.add_argument("--outdir", default="lrd_figs", help="Output dir for --plot.")
    args = ap.parse_args(argv)

    snapshots = []
    for path in args.snapshots:
        try:
            d = load_snapshot(path, args.column_length, args.hydrogen_fraction)
        except (OSError, KeyError) as e:
            sys.stderr.write(f"[skip] {path}: {e}\n")
            continue
        if d is None:
            sys.stderr.write(f"[skip] {path}: no black holes.\n")
            continue
        snapshots.append(d)
        report_snapshot(d)

    if not snapshots:
        sys.exit("No snapshots with black holes were read.")

    if args.plot:
        make_plots(snapshots, args.outdir, args.n_h_star)

    print(
        "\nGuidance on n_H*: read 'fraction above candidate boost thresholds' "
        "above. If f(n_H > 0.1) is tiny at z=4-8, the boost almost never fires "
        "and the bottleneck is n_H*, not beta -- lower n_H* until a meaningful "
        "(~10-30%) dense subset sits above it, OR recalibrate the accretion "
        "normalisation instead. Caveat: lowering n_H* below the EOS/SF density "
        "means boosting gas the resolution has already smoothed -- state it."
    )


if __name__ == "__main__":
    main()