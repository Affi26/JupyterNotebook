# Randle Cycle ODE Modelling

![[randle-1.png]]

![[randle-2.png]]

## 0. Notation and State Variables

We use reaction rates $v_{\text{X}}$ as kinetic terms, to be defined later.

### Glucose Path – Main State Variables

- $\text{Glc}_{\text{ext}}(t)$ – extracellular glucose  
- $\text{Glc}(t)$ – cytosolic glucose  
- $\text{Glc-6-P}(t)$ – glucose 6-phosphate  
- $\text{Fru-6-P}(t)$ – fructose 6-phosphate  
- $\text{Fru-1,6-P}_2(t)$ – fructose 1,6-bisphosphate  
- $\text{Fru-2,6-P}_2(t)$ – fructose 2,6-bisphosphate  
- $\text{Pyruvate}(t)$ – pyruvate (cytosol/mito entry pool)       
- $\text{Glycogen}(t)$ – glycogen store  

### Fat Path – Main State Variables

- $\text{LCFA}_{\text{ext}}(t)$ – extracellular long-chain fatty acids  
- $\text{LCFA-CoA}_{\text{cyto}}(t)$ – cytosolic LCFA-CoA  
- $\text{LCFA-CoA}_{\text{mito}}(t)$ – mitochondrial LCFA-CoA  
- $\text{TAG}(t)$ – triglyceride store  


### Shared State Variables
- $\text{Acetyl-CoA}(t)$ – acetyl-CoA  
- $\text{Cit}(t)$ – citrate pool from Acetyl-CoA
- $\text{Mal}(t)$ – malonyl-CoA from citrate pool

We also have $\text{Glc}_{\text{in}}(t)$, $\text{LCFA}_{\text{in}}(t)$, $\text{ATP}_{\text{draw}}(t)$ as user inputs.

---

## 1. Glucose Path ODEs

### 1.1 Glucose Uptake and Phosphorylation

Extracellular to cytosolic glucose via GLUT4:

$$
\frac{d\,\text{Glc}_{\text{ext}}}{dt}
= \text{Glc}_{\text{in}} - v_{\text{GLUT4}}
$$

where $\text{Glc}_{\text{in}}$ is the user-defined glucose input rate.

$$
\frac{d\,\text{Glc}}{dt}
= v_{\text{GLUT4}} - v_{\text{HK}}
$$

Hexokinase (HK) converts Glc to Glc-6-P:

$$
\frac{d\,\text{Glc-6-P}}{dt}
= v_{\text{HK}} - v_{\text{PGI}} - v_{\text{Glycogen synthase}}
$$

### 1.2 PFK and Glycolytic Intermediates

Isomerization to Fru-6-P:

$$
\frac{d\,\text{Fru-6-P}}{dt}
= v_{\text{PGI}} - v_{\text{PFK}}
$$

PFK produces Fru-1,6-P$_2$:

$$
\frac{d\,\text{Fru-1,6-P}_2}{dt}
= v_{\text{PFK}} - v_{\text{LG}}
$$

Downstream glycolytic flux (lumping ALD + downstream steps) to pyruvate:

$$
\frac{d\,\text{Pyruvate}}{dt}
= 2*v_{\text{LG}} - v_{\text{PDH}} - v_{\text{LDH}}
$$

LDH is step from Pyruvate >> Lactate

### 1.3 Fru-2,6-P$_2$ Regulation

PFK-2 produces Fru-2,6-P$_2$:

$$
\frac{d\,\text{Fru-2,6-P}_2}{dt}
= v_{\text{PFK2,prod}} - v_{\text{PFK2,deg}}
$$

Fru-2,6-P$_2$ will appear as an **activator** in $v_{\text{PFK}}$.

In this reduced Randle-cycle model, we do not track Fru-2,6-P$_2$ or PFK-2 explicitly.

Their effects on PFK activity are assumed to be absorbed into the effective
parameters of $v_{\text{PFK}}$ (e.g., $V_{\max,\text{PFK}}$, $K_{m,\text{PFK}}$).

The dominant regulatory signal we retain is **citrate inhibition of PFK**,
linking mitochondrial fat oxidation to suppression of glycolysis.

### 1.4 Glycogen Synthesis (optional)

From Glc-6-P:

$$
\frac{d\,\text{Glycogen}}{dt}
= v_{\text{Glycogen synthase}}
$$

---

## 2. Fat Path ODEs

### 2.1 LCFA Uptake and Activation

Extracellular LCFA via CD36:

$$
\frac{d\,\text{LCFA}_{\text{ext}}}{dt}
= \text{LCFA}_{\text{in}} - v_{\text{CD36}}
$$

Activation to LCFA-CoA in cytosol:

$$
\frac{d\,\text{LCFA-CoA}_{\text{cyto}}}{dt}
= v_{\text{CD36}} + v_{\text{FAS}} - v_{\text{CPT1}} - v_{\text{TAG synthase}}
$$

### 2.2 CPT1 Transport

Transport into mitochondria:

$$
\frac{d\,\text{LCFA-CoA}_{\text{mito}}}{dt}
= v_{\text{CPT1}} - v_{\beta\text{-ox}}
$$


## 3. Shared substrates

### 3.1 Acetyl-CoA

1. Pyruvate converts to acetyl-CoA via PDH.

2. Mitochondrial LCFA converts to acetyl-CoA via β-oxidation.

3. Acetyl-CoA is consumed to produce ATP (ATP draw)

4. Excess Acetyl-CoA beyond ATP requirements is converted to citrate

$$
\frac{d\,\text{Acetyl-CoA}}{dt}
= v_{\text{PDH}} + v_{\beta\text{-ox}} - v_{\text{cit}} - \text{ATP}_{\text{draw}}
$$


### 3.2 Citrate

Let $\text{Cit}(t)$ be the effective citrate pool (signal).

$$
\frac{d\,\text{Cit}}{dt}
= v_{\text{cit}} - v_{\text{ACL+ACC}}
$$

- $v_{\text{ACL+ACC}}$: lumped citrate clearance (ACL+ACC enzymes)



### 3.3 Malonyl-CoA

Malonyl-CoA from (ACL+ACC) (using citrate):

$$
\frac{d\,\text{Mal}}{dt}
= v_{\text{ACL+ACC}}
- v_{\text{FAS}}
$$

### TAG synthesis from cytosolic LCFA-CoA (optional):

$$
\frac{d\,\text{TAG}}{dt}
= v_{\text{TAG synthase}}
$$

---

We can later refine each $v_{\text{X}}$ into explicit kinetic rate laws (Michaelis–Menten, reversible, Hill, etc.) and start assigning parameter values.

---

# Kinetic Terms – $v_{GLUT4}$

We define the GLUT4 uptake flux with citrate inhibition:

## Full Expression

$$
v_{\text{GLUT4}}
= V_{\max,\text{GLUT4}}
\cdot \frac{\text{Glc}_{\text{ext}}}{K_{m,\text{GLUT4}} + \text{Glc}_{\text{ext}}}
\cdot \frac{1}{1 + \alpha_{\text{Cit,GLUT4}}\,\text{Cit}}
$$


<ul>
  <li><strong>$V_{\max,\text{GLUT4}}$</strong> — maximal GLUT4 transport capacity
    <ul>
      <li>Sets the upper limit of glucose uptake.</li>
      <li>Units: concentration/time.</li>
    </ul>
  </li>

  <li><strong>$K_{m,\text{GLUT4}}$</strong> — Michaelis constant
    <ul>
      <li>External glucose level at half-maximal uptake.</li>
      <li>Controls sensitivity to $\text{Glc}_{\text{ext}}$.</li>
    </ul>
  </li>

  <li><strong>$\alpha_{\text{Cit,GLUT4}}$</strong> — citrate inhibition coefficient
    <ul>
      <li>Strength of citrate suppression on GLUT4 flux.</li>
      <li>Units: inverse concentration of citrate.</li>
    </ul>
  </li>
</ul>

---

# Kinetic Terms – $v_{HK}$

We now define the **hexokinase flux** converting cytosolic glucose to Glc‑6‑P.

## Full Expression

$$
v_{\text{HK}}
= V_{\max,\text{HK}}
\cdot \frac{\text{Glc}}{K_{m,\text{HK}} + \text{Glc}}
\cdot \frac{1}{1 + \alpha_{\text{Glc6P,HK}} \,\text{Glc-6-P}}
$$


- $V_{\max,\text{HK}}$ – maximal HK rate (concentration/time)  
- $K_{m,\text{HK}}$ – HK Michaelis constant for glucose (concentration)  
- $\alpha_{\text{Glc6P,HK}}$ – strength of Glc‑6‑P product inhibition  
  - Units: inverse concentration of Glc‑6‑P

---

# Kinetic Terms – $v_{\text{PGI}}$

Phosphoglucose isomerase converts **Glc‑6‑P → F6P**.  
This step is fast and near‑equilibrium, so we use a simple reversible MM‑style term.

## Full Expression

$$
v_{\text{PGI}}
= V_{\max,\text{PGI}}
\cdot 
\left(
\frac{\text{Glc6P}}{K_{m,\text{PGI}} + \text{Glc6P}}
-
\frac{\text{F6P}}{K_{m,\text{PGI}} + \text{F6P}}
\right)
$$

- $V_{\max,\text{PGI}}$ – maximal PGI catalytic rate (concentration/time)  
- $K_{m,\text{PGI}}$ – Michaelis constant for both Glc‑6‑P and F6P  
- No inhibition term (PGI is not strongly regulated in Randle‑cycle dynamics)

---

# Kinetic Terms – $v_{\text{GS}}$

Glycogen Synthase converts **Glc‑6‑P → Glycogen**.  
We model it with activation by Glc‑6‑P (allosteric) and a simple MM dependence on substrate.

## Full Expression

$$
v_{\text{GS}}
= V_{\max,\text{GS}}
\cdot \frac{\text{Glc6P}}{K_{m,\text{GS}} + \text{Glc6P}}
\cdot \left( 1 + \beta_{\text{Glc6P,GS}} \cdot \text{Glc6P} \right)
$$

- $V_{\max,\text{GS}}$ — maximal glycogen synthase rate (concentration/time)  
- $K_{m,\text{GS}}$ — Michaelis constant for Glc‑6‑P  
- $\beta_{\text{Glc6P,GS}}$ — activation strength by Glc‑6‑P  
  - Units: inverse concentration of Glc‑6‑P  
  - Higher $\beta$ → stronger activation

This is a **phenomenological activation term**, capturing the fact that Glc‑6‑P strongly stimulates GS activity.

---

# Kinetic Terms – $v_{\text{PFK}}$

PFK‑1 converts **F6P → F1,6BP** and is inhibited by citrate.

## Full Expression

$$
v_{\text{PFK}}
= V_{\max,\text{PFK}}
\cdot \frac{\text{F6P}}{K_{m,\text{PFK}} + \text{F6P}}
\cdot \frac{1}{1 + \alpha_{\text{Cit,PFK}} \,\text{Cit}}
$$

- $V_{\max,\text{PFK}}$ – maximal PFK rate  
- $K_{m,\text{PFK}}$ – Michaelis constant for F6P  
- $\alpha_{\text{Cit,PFK}}$ – citrate inhibition strength (units: inverse citrate concentration)

---

# Kinetic Terms – Lumped Lower Glycolysis ($v_{\text{LG}}$)

We lump the steps from F1,6BP to pyruvate into a single effective flux.

## Full Expression

$$
v_{\text{LG}}
= V_{\max,\text{LG}}
\cdot \frac{\text{F1,6BP}}{K_{m,\text{LG}} + \text{F1,6BP}}
$$

- $V_{\max,\text{LG}}$ – maximal lower‑glycolysis rate (concentration/time)  
- $K_{m,\text{LG}}$ – effective Michaelis constant for F1,6BP  

---

# Kinetic Terms – $v_{\text{PDH}}$

Pyruvate dehydrogenase (PDH) converts **Pyruvate → Acetyl‑CoA** and is a key control point in the Randle cycle.  
It is inhibited by mitochondrial **acetyl‑CoA** (reflecting fat oxidation).

## Full Expression

$$
v_{\text{PDH}}
= V_{\max,\text{PDH}}
\cdot \frac{\text{Pyruvate}}{K_{m,\text{PDH}} + \text{Pyruvate}}
\cdot \frac{1}{1 + \alpha_{\text{Acetyl-CoA,PDH}} \,\text{Acetyl-CoA}}
$$

- $V_{\max,\text{PDH}}$ – maximal PDH rate (concentration/time)  
- $K_{m,\text{PDH}}$ – Michaelis constant for pyruvate  
- $\alpha_{\text{Acetyl-CoA,PDH}}$ – inhibition strength by acetyl‑CoA (inverse acetyl‑CoA concentration)

---

# Kinetic Terms – $v_{\text{LDH}}$

Lactate dehydrogenase (LDH) converts **Pyruvate ↔ Lactate**.  
Here we focus on the forward direction (pyruvate → lactate) as a sink for pyruvate.

## Full Expression

$$
v_{\text{LDH}}
= V_{\max,\text{LDH}}
\cdot \frac{\text{Pyruvate}}{K_{m,\text{LDH}} + \text{Pyruvate}}
$$

- $V_{\max,\text{LDH}}$ – maximal LDH rate (concentration/time)  
- $K_{m,\text{LDH}}$ – Michaelis constant for pyruvate  

(If you later add an explicit lactate pool, you can make this reversible.)

---

# Kinetic Terms – $v_{\text{CD36}}$

CD36 mediates **LCFA uptake** from the extracellular space and their activation to LCFA‑CoA in the cytosol.

## Full Expression

$$
v_{\text{CD36}}
= V_{\max,\text{CD36}}
\cdot \frac{\text{LCFA}_{\text{ext}}}{K_{m,\text{CD36}} + \text{LCFA}_{\text{ext}}}
$$

- $V_{\max,\text{CD36}}$ – maximal LCFA uptake rate via CD36 (concentration/time)  
- $K_{m,\text{CD36}}$ – Michaelis constant for extracellular LCFA  

(You can later multiply this by an insulin‑dependent factor if you want CD36 to be hormonally regulated.)

---

# Kinetic Terms – $v_{\text{CPT1}}$

CPT1 transports **LCFA‑CoA from cytosol into mitochondria**, and is inhibited by malonyl‑CoA.

## Full Expression

A simple form with malonyl‑CoA inhibition:

$$
v_{\text{CPT1}}
= V_{\max,\text{CPT1}}
\cdot \frac{\text{LCFA-CoA}_{\text{cyto}}}{K_{m,\text{CPT1}} + \text{LCFA-CoA}_{\text{cyto}}}
\cdot \frac{1}{1 + \alpha_{\text{Mal}} \cdot \text{Mal}}
$$

- $V_{\max,\text{CPT1}}$ – maximal LCFA‑CoA transport rate  
- $K_{m,\text{CPT1}}$ – Michaelis constant for $\text{LCFA-CoA}_{\text{cyto}}$  
- $\alpha_{\text{Mal}}$ – strength of malonyl‑CoA inhibition  

---

# Kinetic Terms – $v_{\beta\text{-ox}}$

β‑oxidation converts **LCFA‑CoA (mito)** into **acetyl‑CoA** inside the mitochondria.

## Full Expression

A simple saturable form:

$$
v_{\beta\text{-ox}}
= V_{\max,\beta}
\cdot
\frac{\text{LCFA-CoA}_{\text{mito}}}
     {K_{m,\beta} + \text{LCFA-CoA}_{\text{mito}}}
$$

- $V_{\max,\beta}$ — maximal β‑oxidation rate  
- $K_{m,\beta}$ — Michaelis constant for mitochondrial LCFA‑CoA  

(You can later multiply this by an AMPK factor or an NADH feedback term if desired.)

---

# Kinetic Terms – $v_{\text{cit}}$

Citrate production converts **mitochondrial acetyl‑CoA** into **citrate** as a signaling pool.

## Full Expression

We use a simple saturable dependence on acetyl‑CoA:

$$
v_{\text{cit}}
= V_{\max,\text{cit}}
\cdot
\frac{\text{Acetyl-CoA}}
     {K_{m,\text{cit}} + \text{Acetyl-CoA}}
$$

- $V_{\max,\text{cit}}$ — maximal citrate production rate  
- $K_{m,\text{cit}}$ — Michaelis constant for acetyl‑CoA in citrate synthesis  

(You can later modulate this by ATP status or NADH if desired.)

---

# Kinetic Terms – $v_{\text{ACL+ACC}}$

The lumped ACL+ACC step converts **citrate (Cit)** into **malonyl-CoA (Mal)** in a **glucose-dependent** manner.

## Full Expression

A simple saturable form with glucose activation:

$$
v_{\text{ACL+ACC}}
= V_{\max,\text{ACL+ACC}}
\cdot
\frac{\text{Cit}}
     {K_{m,\text{ACL+ACC}} + \text{Cit}}
\cdot
\frac{\text{Glc}}
     {K_{\text{ins}} + \text{Glc}}
$$

- $V_{\max,\text{ACL+ACC}}$ — maximal citrate → malonyl-CoA rate  
- $K_{m,\text{ACL+ACC}}$ — citrate Michaelis constant  
- $K_{\text{ins}}$ — effective glucose/insulin sensitivity constant  

(You can later replace the simple glucose term with a more detailed insulin/GLUT4 factor if desired.)

---

# Kinetic Terms – $v_{\text{FAS}}$

Fatty acid synthase (FAS) converts **malonyl-CoA (Mal)** into **LCFA-CoA (cyto)**, promoted by **glucose/insulin**.

## Full Expression

A simple saturable form with glucose activation:

$$
v_{\text{FAS}}
= V_{\max,\text{FAS}}
\cdot
\frac{\text{Mal}}
     {K_{m,\text{FAS}} + \text{Mal}}
\cdot
\frac{\text{Glc}}
     {K_{\text{FAS,Glc}} + \text{Glc}}
$$

- $V_{\max,\text{FAS}}$ — maximal Mal → LCFA-CoA$_{\text{cyto}}$ rate  
- $K_{m,\text{FAS}}$ — malonyl-CoA Michaelis constant  
- $K_{\text{FAS,Glc}}$ — glucose/insulin sensitivity constant for FAS  

(You can later couple this glucose term more tightly to your insulin/GLUT4 proxy if needed.)

---

# Kinetic Terms – $v_{\text{TAG synthase}}$

TAG synthase converts **LCFA-CoA (cyto)** into **triacylglycerol (TAG)**, optionally modulated by **glucose/insulin**.

## Full Expression

A simple saturable form (with optional glucose activation):

$$
v_{\text{TAG synthase}}
= V_{\max,\text{TAG}}
\cdot
\frac{\text{LCFA-CoA}_{\text{cyto}}}
     {K_{m,\text{TAG}} + \text{LCFA-CoA}_{\text{cyto}}}
\cdot
\frac{\text{Glc}}
     {K_{\text{TAG,Glc}} + \text{Glc}}
$$

- $V_{\max,\text{TAG}}$ — maximal LCFA-CoA → TAG rate  
- $K_{m,\text{TAG}}$ — LCFA-CoA Michaelis constant  
- $K_{\text{TAG,Glc}}$ — glucose/insulin sensitivity constant (set large or remove if you want TAG to be mostly passive)

(If you prefer TAG synthesis to be glucose‑independent, you can drop the glucose factor entirely.)

---

