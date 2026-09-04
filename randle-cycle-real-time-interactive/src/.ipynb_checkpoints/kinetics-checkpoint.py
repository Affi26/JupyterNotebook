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

    # --- kinetic parameters ---
    Vmax_PDH    = 2.0
    Km_PDH      = 0.3
    alpha_AcCoA = 0.6   # PDH inhibition by Acetyl-CoA_Glc

    # --- effective ATP draw (saturating) ---
    K_draw = 0.5
    ATP_draw_eff = ATP_draw * (AcCoA_Glc / (K_draw + AcCoA_Glc))

    # --- PDH flux (inhibited by AcCoA) ---
    v_PDH = (
        Vmax_PDH * Pyruvate / (Km_PDH + Pyruvate)
        * (1.0 / (1.0 + alpha_AcCoA * AcCoA_Glc))
    )

    # --- ODE ---
    return v_PDH - ATP_draw_eff



def dCit_dt(Cit, AcCoA_Glu, AcCoA_Fat, Glc):
    """
    ODE for citrate (single signaling pool).
    Production from glucose- and fat-derived AcCoA,
    clearance via lumped ACL+ACC (glucose-dependent).
    """

    # --- citrate production from AcCoA pools ---
    k_Cit_prod_Glu = 0.4   # from Acetyl-CoA_Glu
    k_Cit_prod_Fat = 0.2   # from Acetyl-CoA_Fat

    prod_term = (
        k_Cit_prod_Glu * AcCoA_Glu
        + k_Cit_prod_Fat * AcCoA_Fat
    )

    # --- ACL+ACC clearance (Cit -> Mal), glucose-activated ---
    Vmax_ACL_ACC = 1.5
    Km_ACL_ACC   = 0.10
    K_ins        = 5.0    # glucose/insulin sensitivity

    v_ACL_ACC = (
        Vmax_ACL_ACC
        * Cit / (Km_ACL_ACC + Cit)
        * Glc / (K_ins + Glc)
    )

    # --- ODE ---
    return prod_term - v_ACL_ACC




def dLCFA_ext_dt(LCFA_ext, LCFA_in):
    """
    ODE for extracellular LCFA.
    Kinetic parameters and v_CD36 are defined inline.
    """

    # --- kinetic parameters (explicit here) ---
    Vmax_CD36 = 1.0
    Km_CD36   = 0.1

    # --- kinetic expression ---
    v_CD36 = Vmax_CD36 * LCFA_ext / (Km_CD36 + LCFA_ext)

    # --- ODE ---
    return LCFA_in - v_CD36



def dLCFA_CoA_cyto_dt(LCFA_ext, LCFA_CoA_cyto, Mal, Glc):
    """
    ODE for cytosolic LCFA-CoA.
    All kinetic parameters and fluxes (v_CD36, v_FAS, v_CPT1, v_TAG)
    are defined explicitly inside this function.
    """

    # --- CD36 uptake kinetics ---
    Vmax_CD36 = 1.0
    Km_CD36   = 0.1
    v_CD36 = Vmax_CD36 * LCFA_ext / (Km_CD36 + LCFA_ext)

    # --- FAS kinetics (Mal -> LCFA-CoA_cyto) ---
    Vmax_FAS = 1.0
    Km_FAS   = 0.1
    K_FAS_Glc = 5.0
    v_FAS = Vmax_FAS * (Mal / (Km_FAS + Mal)) * (Glc / (K_FAS_Glc + Glc))

    # --- CPT1 transport (LCFA-CoA_cyto -> mito) ---
    Vmax_CPT1 = 1.0
    Km_CPT1   = 0.1
    alpha_Mal = 0.1   # malonyl-CoA inhibition coefficient (1/mM)

    v_CPT1 = Vmax_CPT1 * LCFA_CoA_cyto / (Km_CPT1 + LCFA_CoA_cyto)
    v_CPT1 *= 1.0 / (1.0 + alpha_Mal * Mal)

    # --- TAG synthase (LCFA-CoA_cyto -> TAG) ---
    Vmax_TAG = 0.8
    Km_TAG   = 0.1
    K_TAG_Glc = 5.0
    v_TAG = Vmax_TAG * LCFA_CoA_cyto / (Km_TAG + LCFA_CoA_cyto)
    v_TAG *= Glc / (K_TAG_Glc + Glc)

    # --- ODE ---
    return v_CD36 + v_FAS - v_CPT1 - v_TAG



def dLCFA_CoA_mito_dt(LCFA_CoA_cyto, LCFA_CoA_mito, Mal):
    """
    ODE for mitochondrial LCFA-CoA.
    Kinetic parameters and fluxes (v_CPT1, v_beta_ox)
    are defined explicitly inside this function.
    """

    # --- CPT1 transport (cyto -> mito) ---
    Vmax_CPT1 = 1.0
    Km_CPT1   = 0.1
    alpha_Mal = 0.1   # malonyl-CoA inhibition coefficient (1/mM)

    v_CPT1 = Vmax_CPT1 * LCFA_CoA_cyto / (Km_CPT1 + LCFA_CoA_cyto)
    v_CPT1 *= 1.0 / (1.0 + alpha_Mal * Mal)

    # --- beta-oxidation (mito LCFA-CoA -> Acetyl-CoA_Fat) ---
    Vmax_beta = 3.0
    Km_beta   = 0.05

    v_beta_ox = Vmax_beta * LCFA_CoA_mito / (Km_beta + LCFA_CoA_mito)

    # --- ODE ---
    return v_CPT1 - v_beta_ox




def dAcCoA_Fat_dt(LCFA_CoA_mito, AcCoA_Fat, ATP_draw):
    """
    ODE for fat-derived acetyl-CoA.
    Includes ATP_draw_eff (saturating sink), matching glucose AcCoA style.
    """

    # --- beta-oxidation (LCFA-CoA_mito -> Acetyl-CoA_Fat) ---
    Vmax_beta = 3.0
    Km_beta   = 0.05

    v_beta_ox = Vmax_beta * LCFA_CoA_mito / (Km_beta + LCFA_CoA_mito)

    # --- effective ATP draw (saturating) ---
    K_draw = 0.5
    ATP_draw_eff = ATP_draw * (AcCoA_Fat / (K_draw + AcCoA_Fat))

    # --- ODE ---
    return v_beta_ox - ATP_draw_eff



def dMal_dt(Cit, Mal, Glc):
    """
    ODE for malonyl-CoA.
    Sources: v_ACL_ACC (from citrate)
    Sinks:   v_FAS (to LCFA-CoA_cyto)
    """

    # --- enforce non-negative pool for kinetics ---
    Mal_pos = max(Mal, 0.0)

    # --- ACL+ACC kinetics (Cit -> Mal) ---
    Vmax_ACL_ACC = 1.5
    Km_ACL_ACC   = 0.10
    K_ins        = 5.0   # glucose/insulin sensitivity

    v_ACL_ACC = (
        Vmax_ACL_ACC
        * Cit / (Km_ACL_ACC + Cit)
        * Glc / (K_ins + Glc)
    )

    # --- FAS kinetics (Mal -> LCFA-CoA_cyto) ---
    Vmax_FAS   = 1.0
    Km_FAS     = 0.10
    K_FAS_Glc  = 5.0

    v_FAS = (
        Vmax_FAS
        * Mal_pos / (Km_FAS + Mal_pos)
        * Glc / (K_FAS_Glc + Glc)
    )

    # --- ODE ---
    return v_ACL_ACC - v_FAS

