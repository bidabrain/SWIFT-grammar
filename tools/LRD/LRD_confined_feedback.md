# LRD confined-feedback phase — implementation notes & progress

**Branch:** `LRD`  ·  **Subgrid model:** `SPIN_JET` (`src/black_holes/SPIN_JET/`)
**Last updated:** 2026-09-07

## 1. Motivation (one paragraph)

The envelope / "black-hole star" picture of Little Red Dots (Kido et al. 2025;
Naidu et al. 2025) implies that during the LRD phase AGN feedback is not merely
weak but *gravitationally confined* and re-radiated as soft (~6000 K)
photospheric emission that does not couple to galaxy-scale gas. We test whether
a single mechanism — **feedback confinement during a super-Eddington phase** —
lets a hydro simulation realise the rapid-growth phase (helping massive BHs
reach the JWST-inferred masses at z > 6) while changing when/how AGN energy
reaches the CGM. Accretion physics is left untouched: we reuse SWIFT's
three-mode disc model (Husko et al.); its slim-disc branch already treats
super-Eddington radiative efficiency reasonably. The only addition is a
**delay/gate in front of the standard AGN feedback**.

## 2. Design decision (what we actually implemented)

Following D. Schleicher's reply, we adopt the **single-parameter
accretion-rate criterion** as the primary trigger (the column-density route is
kept only as a post-processing robustness check, see §6):

- **Trigger (entry/exit):** the BH is in the LRD confined phase whenever the
  **raw, unsuppressed Bondi-based Eddington ratio** exceeds a critical value
  `f_LRD`. "Unsuppressed" = *before* the `max_eddington_fraction` cap and
  *before* accretion-efficiency suppression — so the criterion is decoupled
  from the Eddington cap (this is the "use option 1" choice we agreed on). The
  flag is recomputed every active step, so exit is automatic when `mdot` drops.
- **Gate (what is confined):** while flagged, the AGN feedback energy is **not
  added to the coupling reservoirs** (`energy_reservoir`, `jet_reservoir`).
  I.e. the confined energy is neither stored nor injected — it is assumed
  radiated away softly and lost to CGM coupling. Accretion and BH growth
  continue normally, so the BH grows without self-regulation.

The trigger uses `mdot`; the thing gated is *feedback coupling*. Accretion is
never modified.

## 3. Code changes

All in `src/black_holes/SPIN_JET/`.

### New parameters (`black_holes_properties.h`)
- `int use_lrd_confinement` — master switch (default **0/off**, so existing
  runs are byte-for-byte unaffected).
- `float lrd_eddington_fraction` — the critical unsuppressed Eddington ratio
  `f_LRD`. Required (and must be > 0) when the switch is on; set to a large
  sentinel otherwise.
- Parsed from YAML keys:
  - `SPINJETAGN:use_lrd_confinement` (optional int, default 0)
  - `SPINJETAGN:lrd_confinement_eddington_fraction` (float, required if on)

### New per-BH fields (`black_holes_part.h`)
- `float eddington_fraction_unsuppressed` — the raw Bondi Eddington ratio (the
  decision variable), stored for output.
- `char in_lrd_phase` — 1 while confined, 0 otherwise (recomputed each step).
- `float lrd_confined_mass` — cumulative subgrid mass grown while confined
  (diagnostic for the confined-growth prediction).

### Logic (`black_holes.h`, `black_holes_prepare_feedback`-side accretion routine)
- Capture the unsuppressed ratio right before the Eddington cap:
  `eddington_fraction_unsuppressed = accr_rate / Eddington_rate` (pre-cap,
  pre-efficiency) and store it on the bpart.
- After efficiencies are set, set `in_lrd_phase = use_lrd_confinement &&
  (eddington_fraction_unsuppressed > lrd_eddington_fraction)`.
- Gate reservoir accumulation: the `jet_reservoir` / `energy_reservoir`
  increments are wrapped in `if (!bp->in_lrd_phase) { ... }`.
- Track confined growth: `if (in_lrd_phase) lrd_confined_mass += delta_m_real`.
- Initialise the three new fields in `black_holes_first_init_bpart` (ICs) and
  `black_holes_create_from_gas` (runtime seeding); sum `lrd_confined_mass` in
  `black_holes_swallow_bpart` (mergers). `in_lrd_phase` is not merged (it is
  recomputed each step).

### New snapshot output fields (`black_holes_io.h`, `num_fields` 63 → 66)
- `UnsuppressedEddingtonFractions` — raw Bondi `mdot / mdot_Edd`.
- `LRDConfinedFlags` — 0/1 current confinement state.
- `LRDConfinedMasses` — cumulative mass grown while confined.

### Example parameter file
- `examples/parameter_example.yml` (SPINJETAGN block) documents the two new
  keys with `use_lrd_confinement: 0` and
  `lrd_confinement_eddington_fraction: 100.`.

## 4. Important caveat: which Eddington ratio, and the cap

- SWIFT's standard `EddingtonFractions` output is the **post-cap** value. With
  `max_eddington_fraction = 1` its super-Eddington tail is clipped and useless
  for choosing `f_LRD`. **Use `UnsuppressedEddingtonFractions`** (added here),
  or run the baseline with a large `max_eddington_fraction`.
- The internal Eddington rate uses the spin-dependent Novikov–Thorne radiative
  efficiency (4%–40%), **not** a fixed 0.1. So `f_LRD` is defined relative to
  this spin-dependent Eddington rate; the same threshold corresponds to
  different absolute `mdot` for different spins. This is intrinsic to the model
  and consistent with how `mdot_crit_ADAF` is defined. Worth stating in the
  paper.

## 5. Step C — measuring the distributions (do this first)

`tools/LRD/measure_bh_distributions.py` reads baseline snapshots and reports,
per snapshot / redshift:
- percentiles of subgrid mass, unsuppressed Eddington ratio, and a proxy
  nuclear column `N_H = X_H · rho_gas · L / m_p` (default `L = kernel_gamma·h`);
- the flagged fraction above a grid of Eddington-ratio and column thresholds;
- the **Spearman correlation between the Eddington ratio and the column** — the
  key number for deciding whether the column criterion is redundant;
- LRD-flag counts and confined-mass fractions (if those fields are present).

Run, e.g.:
```
python3 tools/LRD/measure_bh_distributions.py --plot --outdir lrd_figs snap_*.hdf5
```
Pick `f_LRD` from the high tail (e.g. p99 / knee of the PDF). If |rho| is high,
one parameter suffices (answers Schleicher's post-processing question).

## 6. Open questions / TODO

- [ ] **Choose `f_LRD` from data** (step C on a baseline box). Advisor's
  suggested starting point: ~100 × Eddington.
- [ ] **Confinement = lost vs delayed?** We currently *drop* the confined
  energy (don't accumulate the reservoir). Alternative: keep the reservoir but
  gate *injection*, releasing it after the phase — contradicts "cannot couple",
  so we chose lost. Revisit if the CGM response looks wrong.
- [ ] **Leftover reservoir at entry.** Energy accumulated *before* entering the
  LRD phase can still inject on the entry step. Negligible in practice; add an
  injection-side gate if we want strictly zero coupling during confinement.
- [ ] **Column-density criterion as robustness check.** Implement an optional
  second gate on the proxy `N_H` (Compton-thick), then compare LRD demographics
  with the mdot-only run. Needs a boost-factor-style proxy column (advisor's
  main worry — most ad hoc part).
- [ ] **Exit hysteresis?** Currently a hard threshold; consider separate
  entry/exit thresholds to avoid rapid flickering near `f_LRD`.
- [ ] Verify a full build (`./configure ... --with-subgrid=SPIN_JET_EAGLE-XL`
  then `make`); the tree here is not yet configured.

## 7. Mapping to Schleicher's reply
- "boost factor for the column is the most ad hoc part" → we made the
  accretion-rate criterion primary; column is demoted to a check.
- "critical accretion rate, e.g. 100 × Edd, one parameter" → exactly
  `lrd_confinement_eddington_fraction`.
- "tell me the range of mdot and column in your sims" → `tools/LRD/`
  measurement script (step C).
- "check via post-processing whether you need both criteria" → the Spearman
  correlation + threshold-fraction tables in the same script.
