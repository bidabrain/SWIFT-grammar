# Dual-channel (thermal + kinetic) SN feedback — usage

A from-scratch implementation of the Chaikin et al. (2023) thermal-kinetic
stellar SN feedback, added as an **opt-in** module that replaces the
single-channel EAGLE feedback while keeping the rest of a `--with-subgrid`
preset (BH model, cooling, EoS, SF, ...).

> **Status: NOT yet compile-tested.** Written without a local build. Compile on
> the cluster and run the validation checklist below before any science use.

## What was added (all additive; existing options unchanged)

- New configure flag `--with-stellar-feedback-mode=thermal-kinetic` (does
  nothing unless given; default `single` = original behaviour).
- New feedback module `src/feedback/EAGLE_thermal_kinetic/` (copy of
  `EAGLE_kinetic` + an `f_kin` energy split + a stochastic ΔT thermal channel).
- Wiring: `configure.ac`, `src/feedback.h`, `src/Makefile.am`.

Not in upstream SWIFT; not affecting FLAMINGO / EAGLE-XL builds when the flag is
absent.

## Compile

```bash
./autogen.sh                        # required: configure.ac / Makefile.am changed
./configure <your usual flags> \
    --with-subgrid=SPIN_JET_FLAMINGO \    # or SPIN_JET_EAGLE-XL / FLAMINGO / EAGLE-XL
    --with-stellar-feedback-mode=thermal-kinetic
make
```

The flag only swaps the stellar feedback. `SPIN_JET` black holes (and cooling,
EoS, SF, ...) are preserved.

## Parameters

See `EAGLEFeedback_thermal_kinetic_template.yml`. The key additions to a normal
EAGLEFeedback block are:

| parameter | required? | meaning |
|---|---|---|
| `SNII_f_kinetic` | optional, default 1.0 | kinetic energy fraction; `1`=pure kinetic (=FLAMINGO), `0`=pure thermal (=EAGLE-thermal), `0<f<1`=dual channel |
| `SNII_delta_v_km_p_s` | **always** | kinetic kick velocity [km/s] |
| `SNII_delta_T_K` | **when f_kinetic < 1** | thermal-channel heating ΔT [K] (~10^7.5) |

`f_E` (`SNII_energy_fraction_*`) is the total coupled fraction; kinetic gets
`f_kinetic·f_E`, thermal gets `(1-f_kinetic)·f_E`. The split does not change the
total coupled energy.

## Validation checklist (do before trusting it)

1. **Compiles** (`./autogen.sh && ./configure ... && make`); fix errors.
2. **f_kinetic = 1 regression** — should reproduce the pure-kinetic (FLAMINGO)
   result. Confirms the kinetic path wasn't broken.
3. **f_kinetic = 0 regression** — should reproduce pure thermal (EAGLE-thermal).
4. **Energy conservation** — with `--enable-debugging-checks`, verify the total
   injected SN energy ≈ f_E·N_SN·E_SN.
5. **Intermediate f_kinetic** (e.g. 0.1) — run a test galaxy; check the SFH and
   wind mass loading are sane.

## Known simplifications (see code comments)

- The thermal single ray reuses the kinetic "true" ray direction (thermal target
  aligns with the kinetic-true target). Fully independent directions would need a
  separate random draw.
- Constant ΔT (EAGLE-thermal style), not the COLIBRE resolution-dependent ΔT
  range.
- **Not calibrated** for an eEOS setup — `f_kinetic`/`Δv`/`ΔT`/`f_E` must be
  re-calibrated against z=0 galaxy properties for science. (In an eEOS model the
  kinetic channel's cold-ISM-turbulence role is muted; the intermediate mix may
  offer little over the endpoints — see `../LRD/eagle_vs_flamingo.md` §6c.)

## References

- Chaikin et al. (2023), arXiv 2211.04619 — thermal-kinetic SN feedback.
- Chaikin et al. (2022), MNRAS 514, 249 — isotropic ray energy distribution.
- Dalla Vecchia & Schaye (2012) — stochastic thermal SN feedback (ΔT).
