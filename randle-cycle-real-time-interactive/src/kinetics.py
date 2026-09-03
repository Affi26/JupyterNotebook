def dGlc_ext_dt(Glc_ext, Glc_in):
    """
    ODE for extracellular glucose.
    Kinetic parameters and v_GLUT4 are defined inline.
    """

    # --- kinetic parameters (explicit here) ---
    Vmax_GLUT4 = 1.0
    Km_GLUT4   = 5.0

    # --- kinetic expression ---
    v_GLUT4 = Vmax_GLUT4 * Glc_ext / (Km_GLUT4 + Glc_ext)

    # --- ODE ---
    return Glc_in - v_GLUT4



def dGlc_dt(Glc, Glc_ext):
    """
    ODE for cytosolic glucose.
    Kinetic parameters and fluxes are defined inline.
    """

    # --- kinetic parameters ---
    Vmax_GLUT4 = 1.0
    Km_GLUT4   = 5.0

    Vmax_HK = 2.0
    Km_HK   = 0.5

    # --- kinetic expressions ---
    v_GLUT4 = Vmax_GLUT4 * Glc_ext / (Km_GLUT4 + Glc_ext)
    v_HK    = Vmax_HK    * Glc     / (Km_HK    + Glc)

    # --- ODE ---
    return v_GLUT4 - v_HK



def dGlc6P_dt(Glc6P, Glc):
    """
    ODE for glucose-6-phosphate.
    All kinetic parameters and fluxes defined inline.
    """

    # --- kinetic parameters ---
    Vmax_HK = 2.0
    Km_HK   = 0.5

    Vmax_PGI = 1.5
    Km_PGI   = 0.3

    Vmax_GS = 0.8
    Km_GS   = 0.4

    # --- kinetic expressions ---
    v_HK  = Vmax_HK  * Glc    / (Km_HK  + Glc)
    v_PGI = Vmax_PGI * Glc6P  / (Km_PGI + Glc6P)
    v_GS  = Vmax_GS  * Glc6P  / (Km_GS  + Glc6P)

    # --- ODE ---
    return v_HK - v_PGI - v_GS



def dFru6P_dt(Fru6P, Glc6P, Cit):
    """
    ODE for fructose-6-phosphate.
    All kinetic parameters and fluxes defined inline.
    """

    # --- kinetic parameters ---
    Vmax_PGI = 1.5
    Km_PGI   = 0.3

    Vmax_PFK = 2.0
    Km_PFK   = 0.2
    alpha_Cit = 0.5   # citrate inhibition strength

    # --- kinetic expressions ---
    v_PGI = Vmax_PGI * Glc6P / (Km_PGI + Glc6P)

    # PFK inhibited by citrate
    v_PFK = (
        Vmax_PFK * Fru6P / (Km_PFK + Fru6P)
        * (1.0 / (1.0 + alpha_Cit * Cit))
    )

    # --- ODE ---
    return v_PGI - v_PFK



def dFru16P2_dt(Fru16P2, Fru6P, Cit):
    """
    ODE for fructose-1,6-bisphosphate.
    All kinetic parameters and fluxes defined inline.
    """

    # --- kinetic parameters ---
    Vmax_PFK = 2.0
    Km_PFK   = 0.2
    alpha_Cit = 0.5   # citrate inhibition strength

    Vmax_LG = 3.0
    Km_LG   = 0.4

    # --- kinetic expressions ---

    # PFK (inhibited by citrate)
    v_PFK = (
        Vmax_PFK * Fru6P / (Km_PFK + Fru6P)
        * (1.0 / (1.0 + alpha_Cit * Cit))
    )

    # Lumped glycolysis (ALD + downstream)
    v_LG = Vmax_LG * Fru16P2 / (Km_LG + Fru16P2)

    # --- ODE ---
    return v_PFK - v_LG



def dPyruvate_dt(Pyruvate, Fru16P2, AcCoA_Glc):
    """
    ODE for pyruvate.
    All kinetic parameters and fluxes defined inline.
    """

    # --- kinetic parameters ---
    Vmax_LG = 3.0
    Km_LG   = 0.4

    Vmax_PDH = 2.0
    Km_PDH   = 0.3
    alpha_AcCoA = 0.6   # PDH inhibition by Acetyl-CoA_Glc

    Vmax_LDH = 1.5
    Km_LDH   = 0.2

    # --- kinetic expressions ---

    # Lumped glycolysis (ALD + downstream)
    v_LG = Vmax_LG * Fru16P2 / (Km_LG + Fru16P2)

    # PDH (inhibited by Acetyl-CoA_Glc)
    v_PDH = (
        Vmax_PDH * Pyruvate / (Km_PDH + Pyruvate)
        * (1.0 / (1.0 + alpha_AcCoA * AcCoA_Glc))
    )

    # LDH (pyruvate → lactate)
    v_LDH = Vmax_LDH * Pyruvate / (Km_LDH + Pyruvate)

    # --- ODE ---
    return 2.0 * v_LG - v_PDH - v_LDH



def dAcetylCoA_Glc_dt(AcCoA_Glc, Pyruvate, ATP_draw):
    """
    ODE for glucose-derived acetyl-CoA.
    Includes ATP_draw_eff.
    """

    # --- enforce non-negative pool for kinetics ---
    AcCoA_pos = max(AcCoA_Glc, 0.0)

    # --- kinetic parameters ---
    Vmax_PDH    = 2.0
    Km_PDH      = 0.3
    alpha_AcCoA = 0.6   # PDH inhibition by Acetyl-CoA_Glc

    # --- effective ATP draw (saturating) ---
    K_draw = 0.5
    ATP_draw_eff = ATP_draw * (AcCoA_pos / (K_draw + AcCoA_Glc))

    # --- PDH flux (inhibited by AcCoA) ---
    v_PDH = (
        Vmax_PDH * Pyruvate / (Km_PDH + Pyruvate)
        * (1.0 / (1.0 + alpha_AcCoA * AcCoA_Glc))
    )

    # --- ODE ---
    return v_PDH - ATP_draw_eff



def dCit_dt(Cit, AcCoA_Glc):
    """
    ODE for citrate (single signaling pool).
    All kinetic parameters and fluxes defined inline.
    """

    # --- kinetic parameters ---
    k_Cit_prod  = 0.4
    k_Cit_clear = 0.3

    # --- ODE ---
    return k_Cit_prod * AcCoA_Glc - k_Cit_clear * Cit