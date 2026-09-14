# model_interface.py
import numpy as np
from scipy.integrate import solve_ivp
import kinetics

# Define state ordering (must match everywhere)
# state = [Glc_ext, Glc, Glc6P, Fru6P, Fru16P2, Pyruvate,
#          AcCoA_Glc, Cit, LCFA_ext, LCFA_CoA_cyto, LCFA_CoA_mito,
#          AcCoA_Fat, Mal]
STATE_NAMES = [
    "Glc_ext","Glc","Glc6P","Fru6P","Fru16P2","Pyruvate",
    "AcCoA_Glc","Cit","LCFA_ext","LCFA_CoA_cyto","LCFA_CoA_mito",
    "AcCoA_Fat","Mal"
]

# initial state
DEFAULT_STATE = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

def derivatives(t, x, inputs):
    """
    Compose derivatives using functions from kinetics.py.
    inputs: dict with keys 'glu_in', 'fat_in', 'atp_draw'
    """
    # unpack state
    (Glc_ext, Glc, Glc6P, Fru6P, Fru16P2, Pyruvate,
     AcCoA_Glc, Cit, LCFA_ext, LCFA_CoA_cyto, LCFA_CoA_mito,
     AcCoA_Fat, Mal) = x

    # compute each derivative by calling kinetics functions
    dGlc_ext = kinetics.dGlc_ext_dt(Glc_ext, inputs["glu_in"])
    dGlc     = kinetics.dGlc_dt(Glc, Glc_ext)
    dGlc6P   = kinetics.dGlc6P_dt(Glc6P, Glc)
    dFru6P   = kinetics.dFru6P_dt(Fru6P, Glc6P, Cit)
    dFru16P2 = kinetics.dFru16P2_dt(Fru16P2, Fru6P, Cit)
    dPyruv   = kinetics.dPyruvate_dt(Pyruvate, Fru16P2, AcCoA_Glc)
    dAcCoA_G = kinetics.dAcetylCoA_Glc_dt(AcCoA_Glc, Pyruvate, inputs["atp_draw"])
    dCit     = kinetics.dCit_dt(Cit, AcCoA_Glc, AcCoA_Fat, Glc)
    dLCFA_ext = kinetics.dLCFA_ext_dt(LCFA_ext, inputs["fat_in"])
    dLCFA_cyto = kinetics.dLCFA_CoA_cyto_dt(LCFA_ext, LCFA_CoA_cyto, Mal, Glc)
    dLCFA_mito = kinetics.dLCFA_CoA_mito_dt(LCFA_CoA_cyto, LCFA_CoA_mito, Mal)
    dAcCoA_F = kinetics.dAcCoA_Fat_dt(LCFA_CoA_mito, AcCoA_Fat, inputs["atp_draw"])
    dMal     = kinetics.dMal_dt(Cit, Mal, Glc)


    return [
        dGlc_ext, dGlc, dGlc6P, dFru6P, dFru16P2, dPyruv,
        dAcCoA_G, dCit, dLCFA_ext, dLCFA_cyto, dLCFA_mito,
        dAcCoA_F, dMal
    ]



def compute_all_fluxes(x, inputs):
    """
    Compute instantaneous fluxes for visualization.

    Parameters
    - x: state vector in the same ordering as model_interface.STATE_NAMES
    - inputs: dict with keys 'glu_in', 'fat_in', 'atp_draw'

    Returns
    A dict with keys:
      'GLUT4', 'Glycolysis', 'PDH', 'Citrate_prod', 'ACL_ACC',
      'CD36', 'CPT1', 'BetaOx', 'FAS', 'Net_Cit'
    """
    # Unpack state (must match STATE_NAMES)
    (Glc_ext, Glc, Glc6P, Fru6P, Fru16P2, Pyruvate,
     AcCoA_Glc, Cit, LCFA_ext, LCFA_CoA_cyto, LCFA_CoA_mito,
     AcCoA_Fat, Mal) = x

    # Inputs
    glu_in = float(inputs.get("glu_in", 0.0))
    fat_in = float(inputs.get("fat_in", 0.0))
    atp_draw = float(inputs.get("atp_draw", 0.0))

    # --- GLUT4 uptake (kinetics.py constants) ---
    Vmax_GLUT4 = 1.0
    Km_GLUT4 = 5.0
    GLUT4 = Vmax_GLUT4 * Glc_ext / (Km_GLUT4 + Glc_ext + 1e-12)

    # --- Glycolysis proxy (lumped downstream flux v_LG) ---
    Vmax_LG = 3.0
    Km_LG = 0.4
    Glycolysis = Vmax_LG * Fru16P2 / (Km_LG + Fru16P2 + 1e-12)

    # --- PDH flux (inhibited by AcCoA_Glc) ---
    Vmax_PDH = 2.0
    Km_PDH = 0.3
    alpha_AcCoA = 0.6
    PDH = Vmax_PDH * Pyruvate / (Km_PDH + Pyruvate + 1e-12)
    PDH *= 1.0 / (1.0 + alpha_AcCoA * AcCoA_Glc)

    # --- Citrate production from AcCoA pools (same as kinetics.dCit_dt) ---
    k_Cit_prod_Glu = 0.4
    k_Cit_prod_Fat = 0.2
    Citrate_prod = k_Cit_prod_Glu * AcCoA_Glc + k_Cit_prod_Fat * AcCoA_Fat

    # --- ACL + ACC (Cit -> Mal) clearance (glucose-activated) ---
    Vmax_ACL_ACC = 1.5
    Km_ACL_ACC = 0.10
    K_ins = 5.0
    ACL_ACC = Vmax_ACL_ACC * Cit / (Km_ACL_ACC + Cit + 1e-12) * (Glc / (K_ins + Glc + 1e-12))

    # --- CD36 uptake for LCFA_ext ---
    Vmax_CD36 = 1.0
    Km_CD36 = 0.1
    CD36 = Vmax_CD36 * LCFA_ext / (Km_CD36 + LCFA_ext + 1e-12)

    # --- CPT1 transport (inhibited by Mal) ---
    Vmax_CPT1 = 1.0
    Km_CPT1 = 0.1
    alpha_Mal = 0.1
    CPT1 = Vmax_CPT1 * LCFA_CoA_cyto / (Km_CPT1 + LCFA_CoA_cyto + 1e-12)
    CPT1 *= 1.0 / (1.0 + alpha_Mal * max(Mal, 0.0))

    # --- Beta-oxidation in mitochondria ---
    Vmax_beta = 3.0
    Km_beta = 0.05
    BetaOx = Vmax_beta * LCFA_CoA_mito / (Km_beta + LCFA_CoA_mito + 1e-12)

    # --- FAS (Mal -> LCFA_CoA_cyto) for completeness ---
    Vmax_FAS = 1.0
    Km_FAS = 0.1
    K_FAS_Glc = 5.0
    Mal_pos = max(Mal, 0.0)
    FAS = Vmax_FAS * (Mal_pos / (Km_FAS + Mal_pos + 1e-12)) * (Glc / (K_FAS_Glc + Glc + 1e-12))

    # --- Net citrate (production minus ACL/ACC clearance) ---
    Net_Cit = Citrate_prod - ACL_ACC
  

    return {
        "GLUT4": float(max(0.0, GLUT4)),
        "Glycolysis": float(max(0.0, Glycolysis)),
        "PDH": float(max(0.0, PDH)),
        "Citrate_prod": float(max(0.0, Citrate_prod)),
        "ACL_ACC": float(max(0.0, ACL_ACC)),
        "CD36": float(max(0.0, CD36)),
        "CPT1": float(max(0.0, CPT1)),
        "BetaOx": float(max(0.0, BetaOx)),
        "FAS": float(max(0.0, FAS)),
        "Net_Cit": float(Net_Cit)
    }



def compute_instant_fluxes(state, inputs):
    """
    Backwards-compatible wrapper: return the single citrate-production flux.
    """
    fluxes = compute_all_fluxes(state, inputs)
    return float(fluxes.get("Citrate_prod", fluxes.get("Net_Cit", 0.0)))



def step_ivp_once(state, glu_in, fat_in, atp_draw, dt=0.05, method="RK45"):
    """
    Advance the full state by dt and return (next_state, flux).
    """
    inputs = {"glu_in": float(glu_in), "fat_in": float(fat_in), "atp_draw": float(atp_draw)}
    sol = solve_ivp(lambda t, y: derivatives(t, y, inputs),
                    t_span=(0.0, dt),
                    y0=state,
                    method=method,
                    atol=1e-6, rtol=1e-3)
    next_state = sol.y[:, -1]
    next_state = np.maximum(next_state, 0.0)  # clamp negatives
    flux = compute_instant_fluxes(next_state, inputs)
    return next_state, flux
