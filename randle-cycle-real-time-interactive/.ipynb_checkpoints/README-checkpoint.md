# Randle Cycle Real-Time Interactive Simulator

Interactive visualization of a simplified Randle Cycle model.

The project combines:

- Python ODE simulation
- Real-time WebSocket communication
- PixiJS visualization
- Dynamic flux, state and regulation display

---

# Project Structure

randle-cycle-real-time-interactive/

│

├── doc/
│   ├── randle-cycle-ODE-modelling.ipynb
│   │       Full ODE derivations and modelling notes
│   │
│   ├── randle-cycle-simulation-TEST.ipynb
│   │       Early simulation experiments
│   │
│   ├── randle-1.png
│   ├── randle-2.png
│
├── backend/
│   │
│   ├── kinetics.py
│   │       Canonical model definition
│   │       - state metadata
│   │       - ODE equations
│   │       - flux calculations
│   │       - inhibition calculations
│   │
│   ├── model_interface.py
│   │       Numerical integration layer
│   │       - solve_ivp wrapper
│   │       - derivative assembly
│   │
│   ├── server.py
│   │       WebSocket server
│   │       Streams simulation state to frontend
│   │
│   └── test-kinetics-sim.py
│           Standalone backend testing
│
├── frontend/
│   │
│   ├── index.html
│   │       UI controls and canvas
│   │
│   ├── app.js
│   │       Real-time visualization
│   │       - node rendering
│   │       - flux rendering
│   │       - inhibition rendering
│   │       - WebSocket client
│   │
│   ├── layout.json
│   │       Network topology
│   │       - node positions
│   │       - flux edges
│   │       - inhibition edges
│   │
│   └── pixi.min.js
│
└── README.md

---

# Backend Architecture

## kinetics.py

Single source of truth for biological behaviour.

Owns:

- STATE_NAMES
- STATE_INDEX
- DEFAULT_STATE
- ODE equations
- compute_fluxes()
- compute_inhibitions()
- state_to_dict()

## model_interface.py

Responsible for numerical integration.

Uses:

- scipy.solve_ivp()

Returns:

- updated state vector
- current flux dictionary

## server.py

Responsible for communication.

Streams simulation data through WebSockets.

---

# Frontend Architecture

The frontend receives:

```json
{
  "state": {},
  "dstate": {},
  "flux": {},
  "inhibition": {}
}
```

## state

Current metabolite concentrations.

Examples:

```json
{
  "Cit": 1.24,
  "AcCoA": 0.03
}
```

## dstate

Current derivatives (dx/dt).

Examples:

```json
{
  "Cit": 0.21,
  "AcCoA": -0.01
}
```

Positive values indicate accumulation.

Negative values indicate depletion.

## flux

Reaction rates used for pathway visualization.

Examples:

```json
{
  "GLUT4": 0.43,
  "PDH": 0.21
}
```

## inhibition

Normalized inhibition strengths.

Range:

```text
0 = no inhibition
1 = maximal inhibition
```

Examples:

```json
{
  "Glycolysis": 0.75,
  "PDH": 0.12
}
```

---

# Visualization

## Nodes

Display:

- metabolite name
- concentration
- derivative (dx/dt)

Node size scales with concentration.

## Flux Edges

Represent material flow.

Color:

```text
white -> low flux
green -> high flux
```

Thickness:

```text
proportional to flux magnitude
```

## Inhibition Edges

Represent regulation.

Color:

```text
white -> low inhibition
red -> high inhibition
```

Thickness:

```text
proportional to inhibition strength
```

---

# layout.json

Network topology is defined entirely in layout.json.

Example:

```json
{
  "flux_edges": [
    ["Pyruvate", "AcCoA"]
  ],

  "inhibition_edges": [
    ["Cit", "Fru6P"]
  ]
}
```

This allows biological structure and visualization structure to remain separated from rendering logic.

---

# Current Biological Assumptions

Current simplified model includes:

- Glucose uptake via GLUT4
- Glycolysis
- PDH conversion
- Citrate production
- ACL/ACC lipogenesis branch
- Fatty acid uptake via CD36
- CPT1 transport
- Beta oxidation
- ATP-demand driven AcCoA consumption
- ATP-demand driven citrate consumption

The model is intentionally simplified and intended primarily for interactive visualization and systems-level exploration.

