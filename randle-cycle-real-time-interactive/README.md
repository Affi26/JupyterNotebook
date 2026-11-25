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
│   ├── randle-1.png
│   └── randle-2.png
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
│   └── test-kinetics-sim.ipynb
│           Standalone backend testing with jupyter notebook
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
│   └── layout.json
│           Network topology
│           - node positions
│           - flux edges
│           - inhibition edges
│
└── README.md

---

## Start the browser visualization

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