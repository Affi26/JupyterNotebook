contents \& structure:



randle-cycle-real-time-interactive/

│

├── doc/

│   ├── randle-cycle-ODE-modelling.ipynb        # Full ODE derivations + modelling notes

│   ├── randle-cycle-simulation-TEST.ipynb      # Early simulation tests (Python-only)

│   ├── randle-1.png                            # Reference diagram

│   ├── randle-2.png                            # Reference diagram

│

├── backend/

│   ├── kinetics.py                             # All dX/dt equations (core ODE system)

│   ├── test-kinetics-sim.py                    # numerical verification of kinetic ODEs

│   ├── model\_interface.py                      # Converts ODE state → flux dictionary

│   └── server.py                               # WebSocket server streaming fluxes to frontend

│

├── frontend/

│   ├── index.html                              # PixiJS canvas + script loader

│   ├── app.js                                  # Real-time metabolic animation (WebGL)

│   └── pixi.min.js                             # Local PixiJS build (avoids Edge blocking)

│

└── README.md

