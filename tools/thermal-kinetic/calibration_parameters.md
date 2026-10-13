# Calibration parameters — SPIN_JET + dual-channel thermal-kinetic SN feedback

**Scope:** candidate parameters to tune when calibrating this model
(`--with-subgrid=SPIN_JET_FLAMINGO/EAGLE-XL --with-stellar-feedback-mode=thermal-kinetic`),
i.e. SPIN_JET black holes + dual-channel SN feedback on an eEOS (FLAMINGO-style)
baryonic setup. Values are *starting points* — the model is **not calibrated**;
see the caveats at the end.

Calibration = adjust these so the z=0 simulation reproduces the target
observables below.

## Calibration targets (what to match)

| target | constrains mainly |
|---|---|
| **z=0 galaxy stellar mass function (GSMF)** | SN feedback (low-mass end) + AGN (high-mass end) |
| **galaxy sizes** (stellar half-mass radius vs M*) | SN feedback distribution / ISM |
| **M_BH – M\*** (or M_BH – σ) relation | AGN feedback strength + BH seeding |
| **cluster/group gas fractions** (if relevant) | AGN feedback (ejection) |

SN feedback sets the **low-mass** end; AGN sets the **high-mass** end + BH
masses. They are **co-calibrated** — do not tune one in isolation.

---

## A. SN (stellar) feedback — `EAGLEFeedback` block

### A1. Dual-channel split & injection (the new knobs)

| parameter | meaning | start / ref |
|---|---|---|
| **`SNII_f_kinetic`** (f_kin) | Fraction of coupled SN energy injected **kinetically** (rest is thermal). 1 = pure kinetic (FLAMINGO), 0 = pure thermal (EAGLE). **Primary knob.** | COLIBRE ≈ 0.3 (cold-ISM; recalibrate for eEOS) |
| **`SNII_delta_v_km_p_s`** (Δv) | Target **kick velocity** of the kinetic channel → sets wind/turbulence strength. | COLIBRE 50; FLAMINGO 562 (very different roles!) |
| **`SNII_delta_T_K`** (ΔT) | **Heating temperature** of the thermal channel → higher = less overcooling but hotter winds. Required when f_kin < 1. | 10^7.5 ≈ 3.16e7 K |

### A2. Total coupled energy fraction f_E (shared by both channels, split by f_kin)

f_E form in this module (density/metallicity power-law sigmoid):
`f_E = f_E,min + (f_E,max − f_E,min) / (1 + (Z/Z_0)^n_Z · (nH/n_0)^(−n_n))`

| parameter | meaning | note |
|---|---|---|
| **`SNII_energy_fraction_min`** | f_E at low density / high metallicity (per-SN energy in units of 1e51 erg) | = `max` → constant f_E |
| **`SNII_energy_fraction_max`** | f_E at high density / low metallicity | > `min` → density/Z-dependent |
| **`SNII_energy_fraction_n_0_H_p_cm3`** (n_0) | **density pivot**: where f_E transitions min→max | **only matters if min ≠ max** |
| `SNII_energy_fraction_n_n` | steepness of the density dependence | secondary |
| `SNII_energy_fraction_Z_0` | metallicity pivot | secondary |
| `SNII_energy_fraction_n_Z` | steepness of the metallicity dependence | secondary |
| `SNII_energy_erg` | energy per SN | usually fixed at 1e51 |

**Note (density vs pressure):** COLIBRE makes f_E depend on *birth pressure*
(it dropped the eEOS, so density and pressure decouple). You keep the eEOS, where
`P ∝ ρ^(4/3)` locks pressure to density — so this module's **density pivot n_0 is
equivalent to a pressure pivot on the EOS**. Density parameterisation is fine for
your setup.

---

## B. AGN (black-hole) feedback — `SPINJETAGN` block

### B1. Feedback strength (primary AGN knobs)

| parameter | meaning |
|---|---|
| **`coupling_efficiency`** (ε_f) | Fraction of radiated AGN energy coupled to gas (thermal/quasar mode). **Primary AGN knob** for M_BH–M* and high-mass GSMF. |
| **`eps_f_jet`** | Jet coupling efficiency (kinetic/jet mode). |
| `AGN_delta_T_K` | AGN heating temperature (thermal mode) → sets overcooling vs burstiness. |
| `v_jet_km_p_s` | Jet velocity (if `AGN_jet_velocity_model: Constant`). |

### B2. Accretion (set how fast BHs grow → when/whether super-Eddington)

| parameter | meaning |
|---|---|
| `use_subgrid_gas_properties` | Use cold-phase c_s for Bondi rate (boosts accretion). |
| `with_boost_factor` / `boost_beta` | Booth-Schaye accretion boost (off by default; turn on to reach super-Eddington — see `../LRD/baseline_accretion_findings.md`). |
| `include_slim_disk` | Allow super-Eddington (slim-disc) growth. |
| `mdot_crit_ADAF`, `accretion_efficiency_*` | Thick/thin/slim disc transition + ADIOS wind suppression. |
| `max_eddington_fraction` | Hard Eddington cap on accretion (1 = standard; large = uncapped). |

---

## C. Black-hole seeding

| parameter | meaning |
|---|---|
| **`SPINJETAGN:subgrid_seed_mass_Msun`** | BH **seed mass** → higher = faster early BH growth, stronger AGN in low-mass/high-z galaxies. **Primary knob** (COLIBRE emulated this). |
| **`FOF:black_hole_seed_halo_mass_Msun`** | Halo mass at which BHs are seeded → controls *when/where* AGN turns on. |

---

## D. (Optional) LRD confined-feedback gate — `SPINJETAGN`

Only if you enable the LRD model (off by default):

| parameter | meaning |
|---|---|
| `use_lrd_confinement` | Master switch (default 0 = off). |
| `lrd_confinement_eddington_fraction` | Unsuppressed-Eddington threshold above which feedback is confined. |

---

## Priority tiers (tune in this order)

COLIBRE's emulator varied essentially **4 parameters**: `f_kin`, AGN pivot
pressure (≈ AGN strength), `m_BH,seed`, and a density pivot. Map that to here:

1. **`SNII_f_kinetic`** — SN thermal/kinetic split.
2. **AGN strength** — `coupling_efficiency` (and/or `eps_f_jet`, `AGN_delta_T_K`).
3. **`subgrid_seed_mass_Msun`** — BH seed mass.
4. **SN energy scale** — `SNII_energy_fraction_min/max` (f_E) and/or `Δv` / `ΔT`;
   add `n_0` only if you want a density-dependent f_E.

Secondary (fix first, revisit if needed): f_E metallicity params, jet velocity,
accretion-disc efficiencies, seeding halo mass.

## Method

- **Rigorous:** Gaussian-process **emulator** over a Latin hypercube of small-box
  runs, fit to the target observables (FLAMINGO / COLIBRE approach). A substantial
  standalone effort.
- **Pragmatic:** manual / coarse-grid scan of the 4 primary knobs on small boxes;
  check GSMF + M_BH–M* each iteration; converge by hand.

## Caveats

- **Not calibrated for an eEOS setup.** COLIBRE's values (f_kin≈0.3, Δv=50,
  ΔT=10^7.5) are for its cold-ISM model; use only as starting points.
- **Off-design regime:** in an eEOS the kinetic channel's cold-ISM-turbulence role
  is muted, so the optimal f_kin may drift toward an endpoint (0 or 1) — see
  `../LRD/eagle_vs_flamingo.md` §6c.
- **SN ↔ AGN are coupled** — calibrate them together, not separately.
- The dual-channel module is **implemented, compiled, but not yet physics-validated**
  — run the `README.md` checklist (f_kin=1 regression, energy conservation) first.

## References

- Chaikin et al. (2026), arXiv 2509.04067 — COLIBRE calibration (emulated f_kin,
  P_E,pivot, m_BH,seed, n_H,pivot; f_E vs birth pressure, eq. 2).
- Chaikin et al. (2023), arXiv 2211.04619 — thermal-kinetic SN feedback.
- Schaye et al. (2023) — FLAMINGO calibration (emulator, GSMF + cluster gas).
- Husko et al. (2022) — SPIN_JET AGN model.
