# server.py

import asyncio
import json

import websockets

import kinetics

from model_interface import (
    DEFAULT_STATE,
    step_ivp_once,
    derivatives
)

# ------------------------------------------------------------
# Simulation parameters (modifiable from frontend)
# ------------------------------------------------------------

inputs = {
    "glu_in": 0.0,
    "fat_in": 0.0,
    "atp_draw": 0.0
}

# ------------------------------------------------------------
# Current simulation state
# ------------------------------------------------------------

state = DEFAULT_STATE.copy()

# ------------------------------------------------------------
# WebSocket handler
# ------------------------------------------------------------

async def handler(websocket):

    global state
    global inputs

    print("Client connected.")

    try:

        while True:

            # ------------------------------------------------
            # Receive frontend input updates
            # ------------------------------------------------

            try:

                msg = await asyncio.wait_for(
                    websocket.recv(),
                    timeout=0.001
                )

                data = json.loads(msg)

                for key in [
                    "glu_in",
                    "fat_in",
                    "atp_draw"
                ]:

                    if key in data:
                        inputs[key] = float(data[key])

            except asyncio.TimeoutError:
                pass

            # ------------------------------------------------
            # Advance simulation
            # ------------------------------------------------

            state, flux = step_ivp_once(
                state,
                inputs["glu_in"],
                inputs["fat_in"],
                inputs["atp_draw"],
                dt=0.05
            )

            # ------------------------------------------------
            # Compute inhibition metrics
            # ------------------------------------------------

            inhibition = kinetics.compute_inhibitions(
                state
            )

            # ------------------------------------------------
            # Build payload
            # ------------------------------------------------

            dxdt = derivatives(
                0,
                state,
                inputs
            )

            payload = {

                # Actual concentrations
                "state": kinetics.state_to_dict(
                    state
                ),

                # State derivatives (dx/dt)
                "dstate": kinetics.state_to_dict(
                    dxdt
                ),

                # Reaction rates
                "flux": flux,

                # Normalized inhibition strengths
                "inhibition": inhibition
            }

            await websocket.send(
                json.dumps(payload)
            )

            # ~20 FPS

            await asyncio.sleep(0.05)

    except websockets.ConnectionClosed:

        print("Client disconnected.")


# ------------------------------------------------------------
# Server entry point
# ------------------------------------------------------------

async def main():

    async with websockets.serve(
        handler,
        "localhost",
        8765
    ):

        print(
            "WebSocket server running on "
            "ws://localhost:8765"
        )

        await asyncio.Future()


if __name__ == "__main__":
    asyncio.run(main())