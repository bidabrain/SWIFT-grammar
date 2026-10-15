import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import textwrap

matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

def wrap(s, w):
    out = []
    for part in s.split("\n"):
        ww = textwrap.wrap(part, w)
        out += ww if ww else [""]
    return out

PAD = 0.55          # vertical padding (in line units) added per row
LPAD = 0.006        # horizontal text pad inside a cell (axes frac)

ax = None
cursor = [0.0]
def advance(n): cursor[0] += n

def new_canvas():
    global ax
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_axes([0.012, 0.01, 0.976, 0.98])
    ax.set_xlim(0, 1)
    ax.axis("off")
    cursor[0] = 0.0
    return fig

def finish(fig, fname):
    total = cursor[0]
    ax.set_ylim(total, 0)
    fig.set_size_inches(20, total * 0.42)
    fig.savefig(fname, dpi=150, bbox_inches="tight", facecolor="white")
    print("saved", fname, "units=", round(total, 1))
    plt.close(fig)

def band(y, h, color, x0=0.0, x1=1.0, alpha=1.0, z=0):
    ax.add_patch(Rectangle((x0, y), x1 - x0, h, facecolor=color,
                            edgecolor="none", alpha=alpha, zorder=z))

def title_bar(text, fs=18, bg="#1f3b57", fg="white"):
    h = 1.0 + PAD
    y = cursor[0]
    band(y, h, bg)
    ax.text(0.5, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=fg, fontweight="bold", zorder=3)
    advance(h)

def subtitle(text, fs=15, color="#1f3b57"):
    y = cursor[0]
    ax.text(0.012, y + 0.15, text, ha="left", va="top",
            fontsize=fs, color=color, fontweight="bold", zorder=3)
    advance(1.0 + 0.25)

def para(lines, fs=12.5, color="black", gap=0.2):
    for s in lines:
        y = cursor[0]
        ax.text(0.03, y + 0.1, s, ha="left", va="top", fontsize=fs,
                family="DejaVu Sans", color=color, zorder=3)
        advance(1.0 + gap)
    advance(gap)

def formula(s, fs=21):
    h = 1.9
    y = cursor[0]
    ax.text(0.5, y + h / 2, s, ha="center", va="center", fontsize=fs, zorder=3)
    advance(h + 0.2)

def draw_table(col_frac, col_wrap, header, rows, fs=12, hfs=13.5,
               header_bg="#2c3e50", header_fg="white",
               highlight_last=False, halign=None, mono_cols=()):
    W = 1.0
    edges = [0.0]
    for f in col_frac:
        edges.append(edges[-1] + f * W)
    ncol = len(col_frac)
    if halign is None:
        halign = ["left"] * ncol

    hl = [wrap(header[i], col_wrap[i]) for i in range(ncol)]
    hn = max(len(c) for c in hl) + PAD
    y = cursor[0]
    band(y, hn, header_bg)
    for i in range(ncol):
        cx = (edges[i] + edges[i + 1]) / 2
        ax.text(cx, y + hn / 2, "\n".join(hl[i]), ha="center", va="center",
                fontsize=hfs, color=header_fg, fontweight="bold", zorder=3)
    for e in edges:
        ax.plot([e, e], [y, y + hn], color="white", lw=0.8, zorder=4)
    advance(hn)

    for ri, row in enumerate(rows):
        rl = [wrap(row[i], col_wrap[i]) for i in range(ncol)]
        rn = max(len(c) for c in rl) + PAD
        y = cursor[0]
        bg = "#f4f6f8" if ri % 2 == 0 else "#ffffff"
        band(y, rn, bg)
        if highlight_last:
            band(y, rn, "#fff3c4", x0=edges[-2], x1=edges[-1], alpha=0.75)
        for i in range(ncol):
            fam = "DejaVu Sans Mono" if i in mono_cols else "DejaVu Sans"
            if halign[i] == "center":
                tx, ha = (edges[i] + edges[i + 1]) / 2, "center"
            else:
                tx, ha = edges[i] + LPAD, "left"
            ax.text(tx, y + 0.12, "\n".join(rl[i]), ha=ha, va="top",
                    fontsize=fs, family=fam, zorder=3,
                    fontweight=("bold" if i == 0 else "normal"))
        advance(rn)
    advance(0.6)

def note(text, fs=11):
    y = cursor[0]
    ax.text(0.012, y, text, ha="left", va="top", fontsize=fs,
            style="italic", color="#444444", zorder=3,
            wrap=False)
    advance(1.6)

# ======================================================================
# TABLE 1
# ======================================================================
fig1 = new_canvas()
title_bar("Table 1.   Model comparison:   EAGLE  vs  FLAMINGO  vs  COLIBRE  vs  my dual-channel run", fs=18)
advance(0.25)

header = ["Aspect", "EAGLE", "FLAMINGO", "COLIBRE", "My run\n(SPIN_JET_FLAMINGO +\ndual-channel flag)"]
col_frac = [0.095, 0.205, 0.205, 0.225, 0.27]
col_wrap = [14, 29, 29, 32, 38]

rows = [
 ["Star formation",
  "Pressure law (Schaye & Dalla Vecchia 2008); Z-dependent density threshold (Schaye 2004)",
  "Same as EAGLE: pressure law + Z-dep density threshold",
  "Schmidt law, eps_SF = 0.01 per free-fall time; threshold from gravitational instability of resolved cold gas (no Z-dep KS)",
  "Inherited from preset (FLAMINGO): pressure law + Z-dep threshold. UNCHANGED by the flag"],

 ["EoS / ISM floor",
  "eEOS: Jeans entropy floor, P propto rho^4/3 (gamma_eff = 4/3), + temperature floor",
  "Same eEOS + floor (identical to EAGLE)",
  "NO eEOS / no floor; cold ISM resolved via CHIMES non-equilibrium cooling down to ~10 K",
  "Same as FLAMINGO: eEOS + floor. UNCHANGED"],

 ["SN feedback",
  "Single-channel THERMAL: stochastic large-dT heating (Dalla Vecchia & Schaye 2012), dT = 10^7.5 K. f_E = density/Z power-law sigmoid",
  "Single-channel KINETIC: stochastic momentum-conserving pair kicks (Chaikin+2023 ~ mod. DV&S 2008), dv ~ 562 km/s. f_E power-law (~const)",
  "DUAL-channel THERMAL + KINETIC (Chaikin+2023): f_kin ~ 0.3, dv ~ 50 km/s. f_E = logistic in BIRTH PRESSURE (eq. 2, Chaikin+2026)",
  "DUAL-channel THERMAL + KINETIC (Chaikin+2023), opt-in via --with-stellar-feedback-mode=thermal-kinetic. f_kin, dv, dT free; f_E SELECTABLE: PowerLaw OR ColibrePressure"],

 ["BH accretion",
  "Bondi-Hoyle; Booth & Schaye (2009) density boost; Eddington-capped",
  "Bondi-Hoyle; no boost; Eddington-capped",
  "Bondi-based; spin-driven disc; super-Eddington allowed (Chaikin+2026)",
  "SPIN_JET (Husko+2022): Bondi (Krumholz), subgrid gas properties ON, slim disc super-Edd, ADIOS variable efficiency; Edd cap = 1"],

 ["AGN feedback",
  "Single-mode THERMAL: stochastic heating dT ~ 10^8.5 K (Schaye+2015)",
  "THERMAL isotropic stochastic heating (fiducial); jet variant exists",
  "Spin-driven JET (kinetic) + THERMAL; thin/thick/slim disc (Chaikin+2026)",
  "SPIN_JET (Husko+2022): spin-dependent JET (kinetic) + THERMAL/quasar + ADIOS winds; thin/thick/slim disc. UNCHANGED by the flag"],
]

draw_table(col_frac, col_wrap, header, rows, fs=12, hfs=13.5,
           header_bg="#2c3e50", highlight_last=True)

note("Notes:  The flag ONLY swaps the stellar (SN) feedback; star formation, EoS, cooling and the SPIN_JET black-hole model come unchanged from the --with-subgrid preset.")
note("            EAGLE and FLAMINGO share the SAME eEOS + floor.  SPIN_JET_FLAMINGO and SPIN_JET_EAGLE-XL differ only in output/tracers.  COLIBRE BH/AGN entries are indicative (Chaikin+2026).")
advance(0.6)
finish(fig1, "/private/tmp/claude-501/-Volumes-Elements-scratch-SWIFT-LRD/01d829b5-b735-4b1e-a243-74b74350d281/scratchpad/table1_models.png")

# ======================================================================
# TABLE 2
# ======================================================================
fig2 = new_canvas()
title_bar("Table 2.   My dual-channel SN feedback:   energy budget split,   f_E formula,   and calibration knobs", fs=18)
advance(0.3)

subtitle("A.   How the SN energy budget is split between thermal and kinetic")
para([
 "Total coupled SN energy per star:   E_SN,total  =  f_E  x  N_SN  x  E_SN        ( E_SN = 10^51 erg per SN )",
 "The kinetic fraction f_kin only splits the ALREADY-COUPLED energy  -  it does NOT change the total:",
 "     - Kinetic channel:   E_kin =  f_kin        x E_SN,total   ->  stochastic momentum-conserving pair kicks at dv   (Chaikin+2023 ~ mod. Dalla Vecchia & Schaye 2008)",
 "     - Thermal channel:   E_th  = (1 - f_kin)   x E_SN,total   ->  stochastic large-dT heating                       (Dalla Vecchia & Schaye 2012)",
 "     f_kin = 1  ->  pure kinetic  (= FLAMINGO) ;      f_kin = 0  ->  pure thermal  (= EAGLE).",
], fs=12.5, gap=0.22)

subtitle("B.   f_E formula  (COLIBRE version:  logistic in the stellar BIRTH pressure,  Chaikin+2026 eq. 2)")
formula(r"$f_E \;=\; f_{E,\min} \;+\; \dfrac{f_{E,\max}-f_{E,\min}}"
        r"{1+\exp\!\left[-\,\dfrac{\log_{10}P_{\rm birth}-\log_{10}P_{\rm pivot}}{\sigma_P}\right]}$", fs=23)
formula(r"$P_{\rm birth}/k_B \;=\; n_{H,\rm birth}\,T_{\rm birth}\,/\,(X_H\,\mu)"
        r"\;=\; \rho_{\rm birth}\,T_{\rm birth}\,/\,(\mu\,m_H)"
        r"\qquad [\,\mathrm{K\,cm^{-3}}\,],\quad \mu=\mu_{\rm ionised}\approx0.59,\ X_H\approx0.752$", fs=15)

subtitle("C1.   Calibration parameters  -  SUPERNOVA feedback  (EAGLEFeedback block)")
p_header = ["Parameter (yaml)", "Symbol", "Role / effect"]
p_frac = [0.33, 0.08, 0.59]
p_wrap = [36, 40, 82]
sn_rows = [
 ["SNII_f_kinetic", r"$f_{\rm kin}$",
  "Kinetic/thermal split of the SN budget.  PRIMARY knob: trades momentum (winds, turbulence) against heating.  COLIBRE ~ 0.3."],
 ["SNII_delta_v_km_p_s", r"$\Delta v$",
  "Kick velocity of the kinetic channel -> sets wind speed / injected turbulent energy.  COLIBRE 50, FLAMINGO 562 km/s."],
 ["SNII_delta_T_K", r"$\Delta T$",
  "Heating temperature of the thermal channel -> higher = less numerical overcooling.  EAGLE 10^7.5 K.  Required when f_kin < 1."],
 ["SNII_energy_fraction_min", r"$f_{E,\min}$",
  "Floor of the coupled-energy fraction, reached by stars born at LOW pressure."],
 ["SNII_energy_fraction_max", r"$f_{E,\max}$",
  "Ceiling of the coupled-energy fraction, reached by stars born at HIGH pressure."],
 ["SNII_energy_fraction_P_pivot_K_p_cm3", r"$P_{\rm pivot}$",
  "Birth pressure P/k_B [K cm^-3] at which f_E transitions min -> max.  COLIBRE ~ 10^3.4.  Sets WHICH stars get boosted feedback."],
 ["SNII_energy_fraction_sigma_P", r"$\sigma_P$",
  "Width / sharpness of the logistic transition, in dex of log10 P.  Small = sharp switch, large = gradual.  COLIBRE ~ 0.3."],
]
draw_table(p_frac, p_wrap, p_header, sn_rows, fs=12, hfs=13.5,
           header_bg="#34495e", halign=["left", "center", "left"], mono_cols=(0,))

subtitle("C2.   Calibration parameters  -  BLACK HOLES & AGN  (SPINJETAGN / FOF blocks)")
agn_rows = [
 ["SPINJETAGN:subgrid_seed_mass_Msun", r"$M_{\rm seed}$",
  "BH seed mass at creation.  Higher -> faster early BH growth, stronger AGN in low-mass / high-z galaxies.  PRIMARY (COLIBRE emulated this)."],
 ["FOF:black_hole_seed_halo_mass_Msun", r"$M_{\rm halo,seed}$",
  "Halo mass at which a BH is seeded -> controls WHEN / WHERE AGN switches on."],
 ["SPINJETAGN:coupling_efficiency", r"$\epsilon_f$",
  "Fraction of radiated AGN energy coupled to gas (thermal / quasar mode).  PRIMARY AGN knob: high-mass GSMF, M_BH-M*."],
 ["SPINJETAGN:AGN_delta_T_K", r"$\Delta T_{\rm AGN}$",
  "AGN heating temperature (thermal mode) -> higher = burstier, less overcooling; lower = gentler, more frequent heating."],
 ["SPINJETAGN:v_jet_km_p_s", r"$v_{\rm jet}$",
  "Jet (kinetic AGN) velocity -> sets jet feedback strength and reach.  Used when AGN_jet_velocity_model = Constant."],
 ["SPINJETAGN:eps_f_jet", r"$\epsilon_{\rm jet}$",
  "Jet coupling efficiency (kinetic mode).  Scales how much spin-extracted energy drives the jet."],
 ["SPINJETAGN:include_slim_disk\n+ max_eddington_fraction", r"$\dot M/\dot M_{\rm Edd}$",
  "Allow (slim disc) / cap super-Eddington growth -> sets whether BHs can exceed Eddington (relevant for high-z, LRD)."],
]
draw_table(p_frac, p_wrap, p_header, agn_rows, fs=12, hfs=13.5,
           header_bg="#34495e", halign=["left", "center", "left"], mono_cols=(0,))

note("Co-calibration:  SN feedback sets the LOW-mass end of the z=0 galaxy stellar mass function; AGN (coupling_efficiency, seed mass) sets the HIGH-mass end + BH masses.  Tune together.")
note("            Priority order (COLIBRE used ~4 emulated knobs):  (1) f_kin   (2) coupling_efficiency   (3) subgrid_seed_mass_Msun   (4) f_E scale (min/max, P_pivot).   On the eEOS, P is locked to density, so P_pivot behaves like a density pivot.")
advance(0.6)
finish(fig2, "/private/tmp/claude-501/-Volumes-Elements-scratch-SWIFT-LRD/01d829b5-b735-4b1e-a243-74b74350d281/scratchpad/table2_dualchannel.png")
