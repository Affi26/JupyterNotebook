import asyncio
import json
import websockets

from model_interface import (
    DEFAULT_STATE,
    step_ivp_once,
    compute_all_fluxes
)

# ------------------------------------------------------------
# Simulation parameters (modifiable from frontend)
# ------------------------------------------------------------
inputs = {
    "glu_in": 0.1,
    "fat_in": 0.1,
    "atp_draw": 0.1
}

# Current ODE state
state = DEFAULT_STATE.copy()


# ------------------------------------------------------------
# WebSocket handler
# ------------------------------------------------------------
async def handler(websocket):
    global state, inputs

    print("Client connected.")

    try:
        while True:
            # ------------------------------------------------
            # Receive input updates from frontend (optional)
            # ------------------------------------------------
            try:
                msg = await asyncio.wait_for(websocket.recv(), timeout=0.001)
                data = json.loads(msg)

                # Update inputs if provided
                for key in ["glu_in", "fat_in", "atp_draw"]:
                    if key in data:
                        inputs[key] = float(data[key])

            except asyncio.TimeoutError:
                pass  # no message received this frame

            # ------------------------------------------------
            # Step ODE system
            # ------------------------------------------------
            state, _ = step_ivp_once(
                state,
                inputs["glu_in"],
                inputs["fat_in"],
                inputs["atp_draw"],
                dt=0.05
            )

            # ------------------------------------------------
            # Compute flux dictionary
            # ------------------------------------------------
            flux = compute_all_fluxes(state, inputs)

            # ------------------------------------------------
            # Send fluxes to frontend
            # ------------------------------------------------
            await websocket.send(json.dumps(flux))

            # 20 FPS update rate
            await asyncio.sleep(0.05)

    except websockets.ConnectionClosed:
        print("Client disconnected.")


# ------------------------------------------------------------
# WebSocket server entry point
# ------------------------------------------------------------
async def main():
    async with websockets.serve(handler, "localhost", 8765):
        print("WebSocket server running on ws://localhost:8765")
        await asyncio.Future()  # run forever


if __name__ == "__main__":
    asyncio.run(main())
