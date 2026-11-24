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




# Running the Project

The project consists of two separate processes:

1. Backend simulation server
2. Frontend web server

Both must be running simultaneously.

---

## Start Backend Simulation Server

Open a terminal:

```bash
cd backend
python server.py
```

Expected output:

```text
WebSocket server running on ws://localhost:8765
```

The backend is responsible for:

- ODE integration
- Flux calculations
- Inhibition calculations
- WebSocket data streaming

---

## Start Frontend Web Server

Open a second terminal:

```bash
cd frontend
python -m http.server 8000
```

Expected output:

```text
Serving HTTP on 0.0.0.0 port 8000
```

Open browser:

```text
http://localhost:8000
```

---

## Browser Notes

Because layout.json is loaded using fetch():

```javascript
fetch("layout.json")
```

the frontend must be served through a web server.

Opening:

```text
index.html
```

directly using:

```text
file://
```

will trigger browser CORS restrictions and will not work.

---

## Updating Frontend Files

Browsers frequently cache:

- app.js
- layout.json

After making frontend changes:

```text
Ctrl + Shift + R
```

or

```text
Ctrl + F5
```

to force a reload.

If debugging:

```text
F12
→ Network
→ Disable Cache
```

can be useful.

---

## Typical Development Workflow

Terminal 1:

```bash
cd backend
python server.py
```

Terminal 2:

```bash
cd frontend
python -m http.server 8000
```

Browser:

```text
http://localhost:8000
```

Edit:

- kinetics.py
- model_interface.py
- server.py
- app.js
- layout.json

Then:

```text
Restart backend if Python changed
Refresh browser if frontend changed
```