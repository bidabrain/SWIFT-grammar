# EAGLE-XL vs FLAMINGO in SWIFT — what's the same, what differs

**Date:** 2026-09-16
**Context:** notes made while deciding whether to compile the LRD project with
`SPIN_JET_EAGLE-XL` or `SPIN_JET_FLAMINGO`. "Code-verified" = read directly from
`configure.ac` / `src/`; "literature" = from the papers / my reading, flagged.

---

## 0. What "EAGLE-XL" and "FLAMINGO" mean as compile targets

`--with-subgrid=X` in `configure.ac` maps each model to a set of subgrid
modules. All three SPIN_JET variants map the **black holes** module to
`SPIN_JET`:

```
SPIN_JET_EAGLE      → black_holes = SPIN_JET
SPIN_JET_EAGLE-XL   → black_holes = SPIN_JET
SPIN_JET_FLAMINGO   → black_holes = SPIN_JET
FLAMINGO            → black_holes = EAGLE   (the published FLAMINGO fiducial)
```

So `SPIN_JET_FLAMINGO` = FLAMINGO's non-BH physics + the Husko SPIN_JET black
holes; it is **not** the published FLAMINGO AGN model (see §4).

---

## 1. Module-level diff (code-verified, from `configure.ac`)

`SPIN_JET_EAGLE-XL` vs `SPIN_JET_FLAMINGO`:

| subgrid module | EAGLE-XL | FLAMINGO | same? |
|---|---|---|---|
| cooling | PS2020 | PS2020 | ✅ |
| chemistry | EAGLE | EAGLE | ✅ |
| entropy_floor | EAGLE | EAGLE | ✅ |
| stars | EAGLE | EAGLE | ✅ |
| star_formation | EAGLE | EAGLE | ✅ |
| **black_holes** | **SPIN_JET** | **SPIN_JET** | ✅ (LRD code lives here) |
| **feedback (stellar/SN)** | **EAGLE (thermal)** | **EAGLE-kinetic** | ❌ |
| **tracers** | EAGLE | FLAMINGO | ❌ |
| **extra_io** | none | EAGLE | ❌ |
| FOF | yes | yes | ✅ |

**Only three modules differ: stellar feedback, tracers, extra_io.** Everything
else — including the black-hole model — is identical.

---

## 2. The identical parts (in detail)

### 2.1 Cooling — PS2020 (Ploeckinger & Schaye 2020)
Both use the PS2020 radiative cooling tables, which also provide the **subgrid
cold-phase density and temperature** used elsewhere (BH accretion `c_s`, the
`Subgrid` SF threshold).

### 2.2 Star formation — same law
Both use `star_formation = EAGLE`, `SF_model: PressureLaw`:

- **Formula:** the pressure-law reformulation of the Kennicutt–Schmidt law,
  **Schaye & Dalla Vecchia (2008), eq. 21** (cited verbatim in
  `star_formation.h:495`):
  ```
  ṁ_* = m_g · A · (1 M⊙/pc²)^(-n) · (γ/G · f_g · P)^((n-1)/2)
  ```
  n = `KS_exponent` (1.4), breaking to 2.0 above a high-density threshold.
- **SF threshold:** `SF_threshold: Subgrid` → `star_formation.h:268` uses the
  PS2020 subgrid T and n_H: star-forming if
  `T_sub < T1` **or** `(T_sub < T2 and n_H,sub > n_H*)`.
  (Classic EAGLE alternative = metallicity-dependent threshold, Schaye 2004,
  eq. 19 & 24.)

### 2.3 Entropy floor / equation of state — same Jeans EoS
Both use `entropy_floor = EAGLE`. The star-forming (Jeans) branch is a
**polytropic EoS with γ_eff = 4/3**:

- FLAMINGO writes it `T = 8000 K · (n_H/0.1 cm⁻³)^(1/3)`.
- The EAGLE-XL config writes it `T = 800 K · (n_H/1e-4)^(1/3)` (pivot at lower
  density). These are the **same line**: both equal `17241 · n_H^(1/3)` K
  (8000 K at n_H = 0.1, slope 1/3, overdensity > 10).
- "Pressure floor" = "temperature floor" = "entropy floor" = "EoS": the same
  polytrope written in different variables (P = n k T).

**Jeans mass.** M_J ∝ T^(3/2) ρ^(-1/2). Substituting T ∝ ρ^(γ_eff-1) gives
M_J ∝ ρ^((3γ_eff-4)/2); γ_eff = 4/3 makes the exponent zero → **M_J is
density-independent** (~10⁷ M⊙ at this normalisation). That constancy is the
whole point — it stops the fragmentation scale from dropping as gas is
compressed (anti-fragmentation), and works at **any** resolution without being
re-scaled. The value ~10⁷ M⊙ is a *calibration choice* (set by
`Jeans_temperature_norm_K`), not a physical constant, and in FLAMINGO it is
*below* the particle mass (unresolved — the EoS T should not be read as a real
temperature). Different resolutions use the *same* floor because it represents a
physical scale (warm ISM ~10⁴ K) and the unresolved cold phase is handled
separately via the PS2020 subgrid quantities.

### 2.4 Black holes — SPIN_JET, identical
Both use the Husko et al. (2022) spin/jet model (`src/black_holes/SPIN_JET/`).
**All LRD additions** (`UnsuppressedEddingtonFractions`, `LRDConfinedFlags`,
`LRDConfinedMasses`, the `use_lrd_confinement` gate, pre-cap capture) live here
and work **identically** under both compile targets.
Caveat: because the stellar feedback differs, the gas *around* BHs differs, so
accretion *results* can differ even though the BH *code* is the same.

---

## 3. The key difference — stellar (SN) feedback: thermal vs kinetic

### 3.1 EAGLE-XL → thermal (Dalla Vecchia & Schaye 2012)
Stochastic **thermal** feedback:
- Heat a *few* particles by a *large, fixed* ΔT (EAGLE: **ΔT = 10^7.5 K**) so the
  cooling time exceeds the sound-crossing time → avoids numerical overcooling.
- **Stochastic:** each eligible neighbour is heated with probability
  `p_heat ∝ 1/ΔT`, set so ⟨injected energy⟩ = f_th · E_SN. If `p_heat > 1`, ΔT is
  raised so `p_heat = 1`.
- f_th (thermal coupling fraction) is density- & metallicity-dependent in EAGLE
  (`SNII_energy_fraction_min…max` in the yaml).
- **Neighbour selection matters** (Chaikin et al. 2022, MNRAS 514, 249):
  mass-weighted selection is biased toward dense gas → more overcooling; the
  **isotropic (ray)** method they introduce is more efficient. Selection is "as
  important as changing the energy by factors of a few."
- Directly generates a **hot ISM phase** that drives winds and ejects gas.

### 3.2 FLAMINGO → kinetic (`EAGLE_kinetic`; Chaikin et al. 2023, arXiv 2211.04619)
Stochastic **kinetic** feedback (code-verified as essentially pure kinetic):
- Energy budget is `SNII_E_kinetic` only — **no separate thermal channel**.
- **Kicks gas particles in pairs, in exactly opposite directions**, conserving
  energy, linear and angular momentum, statistically isotropic.
- **Pairing (ray method):** each kick event picks a random isotropic direction →
  a "true" ray + a "mirror" ray (opposite). Each ray's target = the gas neighbour
  whose direction is **angularly closest** (min arclength) to the ray. A valid
  kick needs *both* rays to point at real neighbours **and** the particle to be
  "won" by this star (largest-id ownership, to avoid double-kicks).
- **Heat fallback (rare, not a thermal channel):** if a pair can't be kicked
  (a ray points at nothing, or the particle is claimed by a larger-id star), the
  undeliverable energy `E_kinetic_unused = 0.5·energy_per_pair` is added to the
  particle's *internal energy* purely to conserve total energy
  (`u_new = (old KE + old u + injected + E_kinetic_unused − new KE)/m`), then the
  particle is time-step-synced. No fixed ΔT; not the designed thermal channel.
- The kicked particles remain **hydrodynamically coupled** (no decoupling, unlike
  old OWLS winds) → they **shock-heat** the surrounding gas downstream. So
  kinetic feedback *does* produce heating, indirectly, via shocks — which resists
  overcooling (kinetic energy doesn't radiate until it thermalises in a shock).

### 3.3 The "hot phase" and why the channels differ
The **hot phase** = the ~10⁶–10⁷ K, low-density, SN-shock-heated component of the
multiphase ISM. It drives hot galactic winds, feeds the CGM/hot halo, and
regulates SF. Thermal feedback creates it **directly**; pure kinetic creates it
**indirectly via shocks** (may be under-resolved). This is why the full model
(§3.4) uses both channels.

### 3.4 The dual-channel "thermal-kinetic" model (COLIBRE) — NOT in this codebase
Chaikin et al. (2023, arXiv 2211.04619), *"A thermal-kinetic subgrid model…"*,
combines **two channels**: an EAGLE-style large-ΔT isotropic thermal channel
(makes the hot phase, ejects gas) **and** an OWLS-style kinetic kick channel
(drives turbulence in the neutral ISM, suppresses SF). Designed for simulations
that resolve a **cold ISM** (gas cooling to 10 K), i.e. COLIBRE.

**This SWIFT tree only offers** `--with-feedback = none | EAGLE | EAGLE-thermal |
EAGLE-kinetic | GEAR | AGORA`. There is **no dual-channel / COLIBRE feedback
module** here (`COLIBRE` in `configure.ac` is only a deprecated *cooling* alias
that errors → use PS2020). To run the dual-channel model you would need a
COLIBRE / newer SWIFT branch.

### 3.5 Which is "better"?
No universal winner — depends on goal/resolution/what you resolve. Thermal:
simple, well-calibrated, good for large-volume/no-cold-ISM (EAGLE/FLAMINGO
fiducial), but winds can be too hot at the resolution limit and it is sensitive
to neighbour selection. Kinetic: momentum-driven (resists overcooling in dense
gas), better resolution convergence, drives ISM turbulence, good for cold-ISM
sims, but doesn't directly make a hot phase. Dual-channel is best for cold-ISM
sims. Both thermal and kinetic can be calibrated to the same z~0 data; they
differ in *un-calibrated* predictions (winds, CGM, cluster gas fractions,
convergence) — cf. the FLAMINGO point that different feedback implementations
give different predictions at fixed calibration.

---

## 4. AGN feedback — three models to keep straight

| model | mechanism | who uses it |
|---|---|---|
| EAGLE thermal AGN | stochastic ΔT_AGN heating (Booth-Schaye/EAGLE) | EAGLE, EAGLE-XL |
| FLAMINGO fiducial AGN | **same** thermal mechanism, **recalibrated** to cluster gas fractions; + two jet variants | published FLAMINGO |
| SPIN_JET AGN (Husko+22) | spin-dependent 3-mode disc + jets | **your `SPIN_JET_*` builds** |

So the published FLAMINGO's fiducial AGN ≈ EAGLE thermal AGN (recalibrated); its
jet variants differ mechanically (the point of the Schaye+2023 quote about
implementation mattering at fixed calibration). But `SPIN_JET_FLAMINGO` uses
**neither** — it uses the SPIN_JET model. Bottom line: switching EAGLE-XL →
FLAMINGO in your SPIN_JET builds leaves the AGN identical; only the **stellar**
feedback changes.

---

## 5. Accretion-prescription genealogy (context for the LRD work)

Raw BH accretion is Bondi-based; the "capture unresolved dense gas" term differs:

| generation | dense-gas term |
|---|---|
| Booth & Schaye 2009 / OWLS | density **boost** α=(n_H/n_H*)^β (β=2) |
| original EAGLE (Schaye+2015) | **no boost**; angular-momentum (Rosas-Guevara 2015) suppression + Eddington cap |
| EAGLE-XL / SPIN_JET | **subgrid** cold-phase density/`c_s` (`use_subgrid_gas_properties` in SPIN_JET changes only `c_s`; EAGLE's `use_subgrid_bondi` changes both ρ and `c_s`) |

TNG: no boost, Eddington-capped Bondi (like EAGLE); original Illustris: constant
α~100 boost (Springel+2005). See `baseline_accretion_findings.md` for why none
of these reach super-Eddington at cosmological resolution.

---

## 6. Practical implications of switching `SPIN_JET_EAGLE-XL` → `SPIN_JET_FLAMINGO`

- ✅ BH module and your LRD code: unchanged, work as-is.
- ✅ `SPINJETAGN` yaml block: carries over.
- ⚠️ **`EAGLEFeedback` block must change**: EAGLE-XL uses thermal params
  (`SNII_delta_T_K`, energy fractions); FLAMINGO uses `EAGLE_kinetic` params
  (`SNII_delta_v` kick velocity, etc.). Do not copy the thermal feedback block.
- ⚠️ tracers → FLAMINGO, extra_io → EAGLE (more output fields).
- ⚠️ Science: different SN feedback → different ISM around BHs → accretion
  results may differ; re-check rather than assume equality.
- For the LRD question specifically, SN feedback is a **secondary** effect; the
  primary blocker is reaching super-Eddington (boost / resolution).

---

## 7. Key references

- **Schaye & Dalla Vecchia (2008)**, MNRAS 383, 1210 — SF pressure law (eq. 21).
- **Schaye (2004)**, ApJ 609, 667 — metallicity-dependent SF threshold.
- **Dalla Vecchia & Schaye (2012)** — stochastic thermal SN feedback (ΔT).
- **Ploeckinger & Schaye (2020)**, MNRAS 497, 4857 — PS2020 cooling + subgrid
  ISM (source of the subgrid SF threshold and BH subgrid `c_s`).
- **Chaikin et al. (2022)**, MNRAS 514, 249 — isotropic (ray) energy
  distribution for thermal feedback (neighbour selection matters).
- **Chaikin et al. (2023)**, arXiv 2211.04619 — thermal-kinetic dual-channel SN
  feedback (kinetic channel = FLAMINGO's; full model = COLIBRE).
- **Schaye et al. (2015)** — EAGLE. **Schaye et al. (2023)** — FLAMINGO.
- **Husko et al. (2022)** — SPIN_JET black holes.

## 8. Confidence / caveats

- The **module-level table (§1), SPIN_JET-BH identity, feedback options (§3.4),
  code-cited formulae** (SF eq. 21, ΔT=10^7.5 K, kinetic pair-kick + heat
  fallback, EoS equivalence 800 K@1e-4 = 8000 K@0.1) are **read directly from the
  code/config — high confidence.**
- Journal volumes/pages, the "original EAGLE = 8000 K@0.1", FLAMINGO's fiducial
  thermal AGN + jet variants, and the COLIBRE attribution of the dual-channel
  model are from memory/the papers — **verify against the cited references** for
  a write-up.
