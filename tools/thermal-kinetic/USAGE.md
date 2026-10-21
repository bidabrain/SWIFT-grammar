# Dual-channel (thermal + kinetic) SN feedback — usage guide

How to compile, configure, and run the dual-channel stellar SN feedback module
(`EAGLE_thermal_kinetic`), an opt-in model that replaces the single-channel EAGLE
feedback while keeping the rest of a `--with-subgrid` preset (SPIN_JET black
holes, cooling, EoS, star formation, …).

> **Status: implemented, compiles, NOT yet physics-validated.** Run the
> validation checklist (bottom) before any science use. The kinetic kick code is
> copied verbatim from SWIFT's `EAGLE_kinetic` (= Chaikin+2023 kinetic channel,
> FLAMINGO's); the thermal channel, the f_kin split and the COLIBRE f_E option
> were added here.

---

## 1. Compile

```bash
./autogen.sh                          # REQUIRED: configure.ac / Makefile.am changed
./configure <your usual flags> \
    --with-subgrid=SPIN_JET_FLAMINGO \        # or SPIN_JET_EAGLE-XL / FLAMINGO / EAGLE-XL
    --with-stellar-feedback-mode=thermal-kinetic
make
```

Check the configure summary shows:

```
Star feedback model  : EAGLE thermal-kinetic (dual-channel, Chaikin+23) ...
```

and **no** `unrecognized options` warning. (If you see "EAGLE kinetic" + a
warning, you forgot `./autogen.sh`.)

- The flag **only** swaps the stellar feedback. SPIN_JET black holes + cooling +
  EoS + SF are unchanged.
- Without `--with-stellar-feedback-mode=thermal-kinetic`, the build is identical
  to a normal FLAMINGO / EAGLE-XL run (the flag is opt-in).

## 2. Run

No change to the runtime command; `--flamingo` and `--eagle` are interchangeable
(identical engine policies):

```bash
./swift_mpi --threads=32 --cosmology --flamingo --verbose=1 flamingo.yml
```

---

## 3. Parameters (`EAGLEFeedback` block)

### 3.1 Dual-channel core

| parameter | required? | meaning |
|---|---|---|
| `SNII_f_kinetic` | optional, default **1.0** | Fraction of coupled SN energy injected **kinetically** (rest thermal). 1 = pure kinetic (= FLAMINGO); 0 = pure thermal (= EAGLE-thermal); 0<f<1 = dual channel. |
| `SNII_delta_v_km_p_s` | **always** | Kinetic kick velocity [km/s] (COLIBRE 50; FLAMINGO 562). |
| `SNII_delta_T_model` | optional, default **Constant** | Thermal-channel ΔT model: `Constant` or `ColibreDensity` (see 3.1b). |
| `SNII_delta_T_K` | **if f_kinetic < 1** (Constant model) | Thermal-channel heating ΔT [K] (~10^7.5 = 3.16e7). |

f_kin only splits the **already-coupled** energy f_E·N·E_SN; it does not change
the total. kinetic gets `f_kin·f_E`, thermal gets `(1−f_kin)·f_E`.

### 3.1b Thermal ΔT model — `SNII_delta_T_model`

Optional, default `Constant` (= existing behaviour; omit to keep old configs).

**`Constant`** — single heating temperature `ΔT = SNII_delta_T_K` (EAGLE REF).

**`ColibreDensity`** — density-dependent ΔT (Chaikin+2026 eq. 8):
`ΔT(n_H,SN) = ΔT_pivot · (n_H,SN / n_H,pivot)^(2/3)`, clipped to `[ΔT_min, ΔT_max]`.
`n_H,SN` is the **feedback-time** (not birth) SPH-smoothed gas density around the
star — independent of the f_E `use_birth_density` switch. The 2/3 slope is the
Dalla Vecchia & Schaye (2012) cooling argument.

| parameter | required? | meaning |
|---|---|---|
| `SNII_delta_T_pivot_nH_p_cm3` | **yes** (ColibreDensity) | Pivot density n_H,pivot [cm⁻³] — **the free calibration knob**. |
| `SNII_delta_T_pivot_K` | optional, default **10^6.5** | ΔT_pivot [K]. Degenerate with n_H,pivot (COLIBRE fixes it); usually leave default. |
| `SNII_delta_T_min_K` | optional, default **10^6.5** | lower clip [K]. |
| `SNII_delta_T_max_K` | optional, default **10^7.5** | upper clip [K]. |

> Note: `ColibrePressure` f_E (birth pressure) and `ColibreDensity` ΔT
> (feedback-time density) are **independent** — COLIBRE uses birth properties
> for f_E but feedback-time density for ΔT. You can mix any f_E model with any
> ΔT model.

### 3.2 f_E model — `SNII_energy_fraction_model`

Optional, default `PowerLaw` (= existing behaviour; omit to keep old configs
working).

**`PowerLaw`** — density/metallicity power-law sigmoid (EAGLE/FLAMINGO):
`f_E = min + (max−min)/(1 + (Z/Z_0)^n_Z·(n_H/n_0)^(−n_n))`

| parameter | meaning |
|---|---|
| `SNII_energy_fraction_min` / `max` | f_E range (constant if equal) |
| `SNII_energy_fraction_n_0_H_p_cm3` | density pivot (only matters if min≠max) |
| `SNII_energy_fraction_n_n` | density slope |
| `SNII_energy_fraction_Z_0`, `n_Z` | metallicity pivot / slope |
| `SNII_energy_fraction_use_birth_density` | use birth (1) or feedback-time (0) density |
| `SNII_energy_fraction_use_birth_metallicity` | use birth (1) or feedback-time (0) metallicity |

**`ColibrePressure`** — COLIBRE pressure-based logistic (Chaikin+2026, eq. 2):
`f_E = min + (max−min)/(1 + exp(−(log10 P_birth − log10 P_pivot)/σ_P))`
with `P_birth/k_B ≈ n_H,birth · T_birth` [K cm^-3].

| parameter | meaning |
|---|---|
| `SNII_energy_fraction_min` / `max` | f_E,min / f_E,max (COLIBRE fiducial 0.1 / 4) |
| `SNII_energy_fraction_P_pivot_K_p_cm3` | pivot birth pressure P/k_B [K cm^-3] (COLIBRE ~10^3.4; **calibrate**) |
| `SNII_energy_fraction_sigma_P` | logistic transition width (COLIBRE fiducial 0.3) |

**Important — birth override:** under `ColibrePressure`, f_E uses the stellar
**birth** pressure, so `use_birth_density_for_f_th` is **forced to 1 internally**
(the yaml value is overridden; a message is printed if you set 0). Under
`PowerLaw`, `use_birth_density` is read from the yaml as usual.

> Note (pressure proxy): P_birth/k_B is approximated as `n_H,birth · T_birth`,
> which is proportional to the true thermal pressure/k_B; the constant is absorbed
> into the calibrated `P_pivot`. On an eEOS, pressure and density are locked, so
> this is also equivalent to a density pivot.

### 3.3 f_E: which to use

- **Pure kinetic / regression:** keep the **constant** FLAMINGO f_E
  (`PowerLaw`, min=max). Required to reproduce FLAMINGO with f_kinetic=1.
- **eEOS runs:** `PowerLaw` is enough (density ≡ pressure on the EoS). Constant
  f_E (min=max), or density-dependent (min≠max + n_0) if you want it.
- **`ColibrePressure`:** only if you specifically want COLIBRE's birth-pressure
  form (e.g. to mirror COLIBRE, or if you later drop the eEOS).

---

## 4. Example `EAGLEFeedback` snippets

Full template: `EAGLEFeedback_thermal_kinetic_template.yml`.

**(A) Pure kinetic — regression vs FLAMINGO (f_kinetic = 1):**
```yaml
  SNII_f_kinetic:                       1.0
  SNII_delta_v_km_p_s:                  561.807095573
  SNII_energy_fraction_model:           PowerLaw
  SNII_energy_fraction_min:             0.23829533
  SNII_energy_fraction_max:             0.23829533
  # ... FLAMINGO n_0/n_Z/n_n, use_birth_* = 0, etc.
```

**(B) Dual channel, constant f_E (PowerLaw):**
```yaml
  SNII_f_kinetic:                       0.3
  SNII_delta_v_km_p_s:                  50.0
  SNII_delta_T_K:                       3.16228e7
  SNII_energy_fraction_model:           PowerLaw
  SNII_energy_fraction_min:             1.0
  SNII_energy_fraction_max:             1.0
```

**(C) Dual channel, COLIBRE pressure-logistic f_E:**
```yaml
  SNII_f_kinetic:                       0.3
  SNII_delta_v_km_p_s:                  50.0
  SNII_delta_T_K:                       3.16228e7
  SNII_energy_fraction_model:           ColibrePressure
  SNII_energy_fraction_min:             0.1
  SNII_energy_fraction_max:             4.0
  SNII_energy_fraction_P_pivot_K_p_cm3: 8.0e3
  SNII_energy_fraction_sigma_P:         0.3
  # use_birth_density is forced to 1 internally
```

---

## 5. Validation checklist (do before trusting it)

1. **Compiles** (`./autogen.sh && ./configure ... && make`).
2. **f_kinetic = 1 regression** → reproduces pure-kinetic (FLAMINGO). Confirms the
   kinetic path + f_kin split are intact.
3. **f_kinetic = 0 regression** → reproduces pure thermal (EAGLE-thermal).
4. **Energy conservation** (`--enable-debugging-checks`): total injected ≈
   f_E·N_SN·E_SN.
5. **PowerLaw regression**: `model: PowerLaw` should match the pre-COLIBRE-option
   behaviour exactly.
6. **ColibrePressure smoke test**: a small run starts, and f_E varies sensibly
   with birth pressure (check it's between min and max).

## 6. Caveats

- **Not calibrated** for an eEOS setup. COLIBRE values (f_kin≈0.3, Δv=50,
  ΔT=10^7.5, f_E,min/max=0.1/4, P_pivot~10^3.4, σ_P=0.3) are starting points.
- **Off-design regime:** in an eEOS the kinetic channel's cold-ISM-turbulence
  role is muted (see `../LRD/eagle_vs_flamingo.md` §6c).
- SN ↔ AGN feedback are **co-calibrated** — see `calibration_parameters.md`.
- Blind-written thermal channel + f_E option — **needs the checks above**.

## 7. Files

- Module: `src/feedback/EAGLE_thermal_kinetic/`
- Wiring: `configure.ac`, `src/feedback*.h`, `src/Makefile.am`
- This dir: `USAGE.md` (this file), `README.md`, parameter template,
  `calibration_parameters.md`.

## 8. References

- Chaikin et al. (2023), arXiv 2211.04619 — thermal-kinetic SN feedback.
- Chaikin et al. (2026), arXiv 2509.04067 — COLIBRE (f_E vs birth pressure, eq. 2).
- Dalla Vecchia & Schaye (2012) — stochastic thermal feedback (ΔT).
- Husko et al. (2022) — SPIN_JET black holes.
