# model_interface.py

import numpy as np
from scipy.integrate import solve_ivp

import kinetics


# ------------------------------------------------------------
# Canonical model definitions
#
# Ownership now resides in kinetics.py.
# Re-exported here for backwards compatibility.
# ------------------------------------------------------------

STATE_NAMES = kinetics.STATE_NAMES
STATE_INDEX = kinetics.STATE_INDEX
DEFAULT_STATE = kinetics.DEFAULT_STATE


# ------------------------------------------------------------
# Derivative assembly
# ------------------------------------------------------------

def derivatives(t, x, inputs):
    """
    Compose state derivatives using kinetics.py.

    Parameters
    ----------
    t : float
        Solver time (unused but required by solve_ivp)

    x : ndarray
        State vector

    inputs : dict
        {
            "glu_in": float,
            "fat_in": float,
            "atp_draw": float
        }
    """

    (
        Glc_ext,
        Glc,
        Glc6P,
        Fru6P,
        Fru16P2,
        Pyruvate,
        AcCoA,
        Cit,
        LCFA_ext,
        LCFA_CoA_cyto,
        LCFA_CoA_mito,
        Mal,
        TAG
    ) = x

    dGlc_ext = kinetics.dGlc_ext_dt(
        Glc_ext,
        Cit,
        inputs["glu_in"]
    )

    dGlc = kinetics.dGlc_dt(
        Glc,
        Glc_ext,
        Cit
    )

    dGlc6P = kinetics.dGlc6P_dt(
        Glc6P,
        Glc
    )

    dFru6P = kinetics.dFru6P_dt(
        Fru6P,
        Glc6P,
        Cit
    )

    dFru16P2 = kinetics.dFru16P2_dt(
        Fru16P2,
        Fru6P,
        Cit
    )

    dPyruv = kinetics.dPyruvate_dt(
        Pyruvate,
        Fru16P2,
        AcCoA
    )

    dAcCoA = kinetics.dAcetylCoA_dt(
        AcCoA,
        Pyruvate,
        LCFA_CoA_mito,
        inputs["atp_draw"]
    )

    dCit = kinetics.dCit_dt(
        Cit,
        AcCoA,
        Glc,
        inputs["atp_draw"]
    )

    dLCFA_ext = kinetics.dLCFA_ext_dt(
        LCFA_ext,
        inputs["fat_in"]
    )

    dLCFA_cyto = kinetics.dLCFA_CoA_cyto_dt(
        LCFA_ext,
        LCFA_CoA_cyto,
        Mal,
        Glc
    )

    dLCFA_mito = kinetics.dLCFA_CoA_mito_dt(
        LCFA_CoA_cyto,
        LCFA_CoA_mito,
        Mal
    )

    dMal = kinetics.dMal_dt(
        Cit,
        Mal,
        Glc
    )

    dTAG = kinetics.dTAG_dt(
        LCFA_CoA_cyto,
        Glc
    )

    return [
        dGlc_ext,
        dGlc,
        dGlc6P,
        dFru6P,
        dFru16P2,
        dPyruv,
        dAcCoA,
        dCit,
        dLCFA_ext,
        dLCFA_cyto,
        dLCFA_mito,
        dMal,
        dTAG
    ]


# ------------------------------------------------------------
# Simulation step
# ------------------------------------------------------------

def step_ivp_once(
    state,
    glu_in,
    fat_in,
    atp_draw,
    dt=0.05,
    method="RK45"
):
    """
    Advance model state by one integration step.

    Returns
    -------
    next_state : ndarray
    fluxes : dict

    Fluxes are now computed by kinetics.py,
    which is the single source of truth.
    """

    inputs = {
        "glu_in": float(glu_in),
        "fat_in": float(fat_in),
        "atp_draw": float(atp_draw)
    }

    sol = solve_ivp(
        lambda t, y: derivatives(t, y, inputs),
        t_span=(0.0, dt),
        y0=state,
        method=method,
        atol=1e-6,
        rtol=1e-3
    )

    next_state = sol.y[:, -1]

    # Prevent tiny negative numerical values
    next_state = np.maximum(next_state, 0.0)

    fluxes = kinetics.compute_fluxes(
        next_state,
        inputs
    )

    return next_state, fluxes
