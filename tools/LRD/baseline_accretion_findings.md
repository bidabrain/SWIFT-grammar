# LRD baseline runs: the accretion never reaches super-Eddington

**Date:** 2026-09-14
**Runs analysed:** `LRD_base` (SPIN_JET_EAGLE-XL, `max_eddington_fraction = 1`) and
`LRD_base_free` (identical but `max_eddington_fraction = 999`). LRD confinement
gate **off** in both (`use_lrd_confinement: 0`), i.e. these are *calibration /
baseline* runs meant to measure the unsuppressed Eddington-ratio distribution
before turning the gate on.

---

## TL;DR

- The `UnsuppressedEddingtonFractions` field works and is genuinely captured
  **before** the Eddington cap and the accretion-efficiency suppression
  (`black_holes.h:866–876`). The analysis pipeline reads it correctly.
- **The raw, unsuppressed Bondi/Eddington ratio never exceeds ~0.86 at any
  redshift** (peak `0.863` at z = 6.79; `n(>1×Edd) = 0` in *every* snapshot from
  z = 13.8 down to z = 0). So there is **no super-Eddington accretion anywhere**.
- Because nothing reaches 1×Edd, the `max_eddington_fraction` cap **never
  binds** → `base` (cap = 1) and `base_free` (cap = 999) are effectively
  identical, and an LRD gate triggered on "unsuppressed Edd > 1" would **never
  fire**.
- This is **not a bug**. It is the expected behaviour of an EAGLE-type accretion
  model at this resolution (25 Mpc/h box), driven by (a) feedback
  self-regulation and (b) unresolved nuclear gas.

---

## The diagnostic

Scanned all snapshots of `LRD_base` for the per-BH unsuppressed and capped
Eddington ratios:

| quantity | result |
|---|---|
| redshift coverage | z = 45 → 0 (BHs appear from z ≈ 13.8, `snap_0005`) |
| `max(UnsuppressedEddingtonFractions)` over all snapshots | **0.863** (at z = 6.79) |
| number of BHs with unsuppressed Edd > 1 | **0**, at every redshift |
| `max(UnsuppressedEddingtonFractions)` vs `max(EddingtonFractions)` | equal in every snapshot (cap never engages) |

The maxima wander (0.5–0.86) just **below** 1 and do **not** pile up at exactly
1.0 → it is a *soft* ceiling, not a hard numerical cap.

The first PDF figure (`figs_base/eddington_ratio_pdf.png`) had only plotted the
late (z ≲ 0.1) snapshots, which is why it initially looked "capped at 1"; the
full scan shows the ceiling holds at every epoch.

---

## Why there is no super-Eddington

### 1. Feedback self-regulation (the soft Eddington ceiling)

The unsuppressed ratio scales as

```
mdot / mdot_Edd  ∝  M_BH · ρ / c_s³
```

If the gas conditions (ρ, c_s) were fixed this would *grow* as the BH grows
(∝ M_BH). Instead it stays pinned near ~1, which means ρ/c_s³ around the BH
*drops* as M_BH grows: whenever a BH approaches Eddington its AGN feedback
heats/expels the surrounding gas (ρ↓, c_s↑), cutting the fuel supply
(Bondi ∝ ρ/c_s³). This negative-feedback loop is the standard self-regulation
that puts BHs on the M_BH–M* relation, and it produces an order-unity ceiling
without any hard cap.

### 2. Resolution — the nuclear gas is not resolved

`use_subgrid_gas_properties: 1` (on in these runs) only corrects the **phase**
underestimate: at fixed *resolved* pressure P it computes the cold-phase sound
speed `c = sqrt(γP/ρ_sub)` (`black_holes_iact.h:93–97`), lowering c_s and hence
raising the Bondi rate (∝ 1/c_s³). In SPIN_JET it changes **only c_s** — the
density used in the Bondi rate is always the resolved SPH density.

It does **not** correct the **resolution** underestimate: ρ_sub is derived from
the resolved pressure, and at 25 Mpc/h (softening ~kpc) the nuclear
pressure/density peak (the pc–100 pc gas that drives super-Eddington feeding) is
simply not resolved. So even with subgrid on, the peak accretion is
systematically underestimated. The empirical proof is this run itself:
subgrid c_s is on and the ratio still never exceeds 0.86.

**Relevant scales:** Bondi radius ~0.1–10 pc; nuclear gas reservoir ~10–100 pc;
this run's softening ~kpc → ~2 orders of magnitude too coarse to resolve the
feeding scale.

---

## Accretion-prescription genealogy (where "boost" and "subgrid" fit)

BH accretion is computed in layers; the models differ mostly in layer ①.

| generation | raw Bondi uses | density boost | other |
|---|---|---|---|
| **Booth & Schaye 2009 / OWLS** | resolved ρ + c_s | **yes:** α = (n_H/n_H*)^β, β=2 ← origin of "boost" | — |
| **original EAGLE (Schaye+2015)** | resolved ρ + c_s (with T floor) | **no (α=1), removed** | Rosas-Guevara (2015) angular-momentum suppression + Eddington cap |
| **EAGLE-XL / SPIN_JET** | **subgrid** cold-phase ρ and/or c_s | usually off | spin-dependent disc physics (SPIN_JET) |

Notes:
- The density **boost** belongs to Booth & Schaye (2009)/OWLS, **not** to the
  EAGLE (2015) fiducial, which *removed* it and instead *suppressed* accretion
  via the angular-momentum limiter.
- The SWIFT **EAGLE** module has the subgrid option too, named
  `use_subgrid_bondi`, and it is *more complete* than SPIN_JET's: it swaps
  **both** ρ and c_s for the cold-phase values (`EAGLE/black_holes.h:756–786`),
  whereas SPIN_JET's `use_subgrid_gas_properties` swaps **only c_s**.
- The CAMELS `CV_23` EAGLE-XL config uses `use_subgrid_bondi: 0` and a
  nominal-but-unity boost (`boost_alpha_only: 1, boost_alpha: 1.0`) → bare
  resolved-density Bondi, the least dense-gas compensation of all.

### Why EAGLE removed boost and added the Eddington cap

- **Removed boost:** it was a resolution-dependent fudge (β, n_H* re-tuned per
  resolution); and because of self-regulation the accretion normalisation is
  **degenerate** with feedback — the z~0 calibration target (galaxy stellar-mass
  function) is insensitive to it — so the boost was unnecessary to hit the
  targets and was dropped for parsimony, replaced by the physically-motivated
  angular-momentum limiter.
- **Added Eddington cap:** physically motivated (classical limit for a
  radiatively efficient thin disc; keeps the subgrid feedback model valid), and
  it provides numerical control — preventing resolution-dependent runaway
  accretion spikes and keeping BH growth smooth, predictable and calibratable.

**Net:** EAGLE deliberately trades away high-z rapid/super-Eddington growth for
robustness and calibratability. That is exactly why it (and CAMELS, and this
run) does not produce super-Eddington BHs.

---

## Edd–N_H correlation (Schleicher's redundancy check)

`figs_base/edd_vs_column.png` shows a clear positive correlation between the
unsuppressed Eddington ratio and the proxy column density N_H. Both scale with
the nuclear gas density ρ (Edd ∝ M·ρ/c_s³, N_H ∝ ρ·L), so:

- **High accretion ⇔ high column (obscured)** — physically consistent with LRDs
  being red/obscured high-accretors.
- A single accretion-rate criterion (as implemented) is likely **sufficient**;
  a separate column-density gate would select a similar (not identical, ~1–2 dex
  scatter) population → supports the single-parameter LRD trigger.

**Caveats:** the plot shown is at z = 0 (should be redone at z ≈ 4–7); N_H is a
kernel-scale *proxy*, not a true line-of-sight column; and the correlation is
partly *built-in* because the proxy N_H uses the same `rho_gas` as the Bondi
rate.

---

## Tension with observations (the LRD motivation)

- Cosmological EAGLE-type models struggle to make high/super-Eddington BHs
  **by construction** (Eddington cap + self-regulation + resolution + angular-
  momentum suppression); they are calibrated to z~0, not to high-z rapid growth.
- Observationally (JWST, z ~ 4–11): abundant, often **overmassive** BHs (LRDs)
  that require **rapid early growth**. Super-Eddington accretion is a leading
  candidate mechanism (alongside heavy seeds), and there are individual
  super-Eddington candidates — **but** the Eddington ratios of LRDs specifically
  are **debated** (many are ~Eddington or sub-Eddington; their nature is
  contested). So the honest statement is "observations require rapid growth and
  *some* super-Eddington candidates", not "many confirmed super-Eddington BHs".
- The LRD confined-feedback model is in the **super-Eddington** camp: by
  confining feedback it removes the self-regulation that pins accretion at
  Eddington, letting BHs grow super-Eddington. Understanding *why* EAGLE capped
  clarifies both what we are undoing and the risk we inherit (loss of
  robustness / calibratability, resolution/parameter sensitivity).

---

## Implications & path forward

1. **As configured, the LRD gate cannot fire** (nothing exceeds ~0.86×Edd), and
   `base` ≈ `base_free`. To test the LRD model at all, BHs must first be pushed
   into the super-Eddington regime.
2. Options, in order of physical cleanliness:
   - **Higher resolution / zoom-in** — resolves the nuclear pressure peak so
     ρ_sub rises self-consistently. Cleanest, most expensive.
   - **Add a real Booth-Schaye boost** — keep `use_subgrid_gas_properties: 1`
     (do **not** turn it off — that only lowers c_s-driven accretion and works
     against the goal) and set `with_boost_factor: 1` with a genuine β (e.g. 2),
     on the `base_free` (cap = 999) config so slim-disc super-Eddington growth is
     allowed. Fast, but ad-hoc and double-counts dense-gas compensation; needs
     recalibration for a science run.
   - **Tune the subgrid accretion** within the existing framework.
3. Whichever path: this is a **methodology decision** to align with Schleicher —
   accept an ad-hoc boost vs. commit to resolution — and it inherits exactly the
   robustness/calibration issues EAGLE used the cap to avoid. Mitigations worth
   considering: a finite-but-higher cap, and a parameter scan over boost/gate
   strength to check numerical robustness.

---

## Caveats / to verify against the literature

- The self-regulation, degeneracy and prescription-genealogy statements are the
  standard modelling picture; exact parameter values and the precise stated
  rationale should be checked against Booth & Schaye (2009), Schaye et al.
  (2015) §4.4, Rosas-Guevara et al. (2015), and the Husko et al. SPIN_JET
  papers.
- High-z BH / LRD observations are a fast-moving frontier (2024–2026); the
  Eddington-ratio distribution and overmassive fraction should be checked
  against the latest LRD literature (Greene, Matthee, Maiolino, Kokorev, …).
