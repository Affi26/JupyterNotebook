# Randle Cycle Model Documentation

## Overview

This document describes the current conceptual model behind the Randle Cycle Interactive Simulator. It is distilled from the original modelling notebook and serves as a human-readable reference.

## State Variables

### Glucose Path

- Glc_ext — extracellular glucose
- Glc — cytosolic glucose
- Glc6P — glucose-6-phosphate
- Fru6P — fructose-6-phosphate
- Fru16P2 — fructose-1,6-bisphosphate
- Fru26P2 — fructose-2,6-bisphosphate (conceptual, not currently implemented)
- Pyruvate — pyruvate pool
- Glycogen — glycogen store (conceptual, not currently implemented)

### Fat Path

- LCFA_ext — extracellular long-chain fatty acids
- LCFA_CoA_cyto — cytosolic LCFA-CoA
- LCFA_CoA_mito — mitochondrial LCFA-CoA
- TAG — triglyceride store (conceptual state, currently represented as a sink)

### Shared Pools

- AcCoA — acetyl-CoA
- Cit — citrate
- Mal — malonyl-CoA

### Inputs

- Glc_in
- LCFA_in
- ATP_draw

---

# ODE Structure

## Glucose Uptake

d(Glc_ext)/dt = Glc_in - v_GLUT4

d(Glc)/dt = v_GLUT4 - v_HK

d(Glc6P)/dt = v_HK - v_PGI - v_GS

## Glycolysis

d(Fru6P)/dt = v_PGI - v_PFK


d(Fru16P2)/dt = v_PFK - v_LG


d(Pyruvate)/dt = 2*v_LG - v_PDH - v_LDH

## Fat Uptake

d(LCFA_ext)/dt = LCFA_in - v_CD36


d(LCFA_CoA_cyto)/dt = v_CD36 + v_FAS - v_CPT1 - v_TAG


d(LCFA_CoA_mito)/dt = v_CPT1 - v_BetaOx

## Shared Substrates

### Acetyl-CoA

d(AcCoA)/dt = v_PDH + v_BetaOx - v_CitrateProd - ATP_draw

### Citrate

d(Cit)/dt = v_CitrateProd - v_ACL_ACC

### Malonyl-CoA

d(Mal)/dt = v_ACL_ACC - v_FAS

### TAG Storage

d(TAG)/dt = v_TAG

---

# Kinetic Terms

## v_GLUT4

Glucose uptake from extracellular space.

v_GLUT4 = Vmax_GLUT4 * Glc_ext / (Km_GLUT4 + Glc_ext)

## v_HK

Hexokinase: Glc -> Glc6P

## v_PGI

Phosphoglucose isomerase: Glc6P <-> Fru6P

## v_GS

Glycogen synthesis from Glc6P.

## v_PFK

PFK: Fru6P -> Fru16P2

Includes citrate inhibition.

## v_LG

Lumped lower glycolysis.

Fru16P2 -> Pyruvate

## v_PDH

Pyruvate -> AcCoA

Inhibited by AcCoA.

## v_LDH

Pyruvate -> Lactate sink.

## v_CD36

Fatty acid uptake.

## v_CPT1

LCFA_CoA_cyto -> LCFA_CoA_mito

Inhibited by malonyl-CoA.

## v_BetaOx

Beta oxidation.

LCFA_CoA_mito -> AcCoA

## v_CitrateProd

AcCoA -> Citrate

## v_ACL_ACC

Cit -> Mal

Activated by glucose/insulin proxy.

## v_FAS

Mal -> LCFA_CoA_cyto

Activated by glucose/insulin proxy.

## v_TAG

LCFA_CoA_cyto -> TAG storage

---

# Regulatory Logic

## Citrate

Citrate suppresses glycolysis via PFK inhibition.

## Acetyl-CoA

Acetyl-CoA suppresses PDH.

## Malonyl-CoA

Malonyl-CoA suppresses CPT1 and therefore fatty acid oxidation.

---

# Current Architecture

PARAMS
    ↓
Flux Functions (v_*)
    ↓
ODE Functions
    ↓
compute_fluxes()
    ↓
compute_inhibitions()
    ↓
WebSocket Transport
    ↓
PixiJS Frontend

---

# Current Simplifications

Implemented:

- Glucose uptake and glycolysis
- PDH regulation
- Fat uptake and beta oxidation
- Citrate pool
- Malonyl-CoA pool
- CPT1 inhibition
- ATP-demand sink terms

Not yet fully implemented:

- Explicit Fru-2,6-P2 dynamics
- Explicit glycogen store
- Explicit TAG state pool
- Hormonal signalling
- AMPK signalling
- Mitochondrial redox state
- TCA cycle intermediates beyond citrate
