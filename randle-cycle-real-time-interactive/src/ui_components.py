# ui_components.py

import plotly.graph_objects as go



# helper function
def norm(v, vmin=0.0, vmax=1.0):
    return max(0, min(1, (v - vmin) / (vmax - vmin)))


# ---------------------------------
# Color gradient for concentration
# ---------------------------------
def conc_color(value):
    v = max(0, min(1, value))
    
    # green → yellow → red
    r = int(255 * v)
    g = int(255 * (1 - v))
    b = 0

    return f"rgba({r},{g},{b},1.0)"


# -------------------------
# Color gradient for flux
# -------------------------
def flux_color(flux):
    f = max(0, min(1, flux))
    
    # green → yellow → red
    r = int(255 * f)
    g = int(255 * (1 - f))
    b = 0

    return f"rgba({r},{g},{b},1.0)"


# -------------------------
# Node constructor
# -------------------------
def make_node(name, x, y, value):
    size = 5 + 40 * value
    value_norm = norm(value)
    color =  conc_color(value)
    return {
        "type": "node",
        "name": name,
        "x": x,
        "y": y,
        "size": size,
        "color": color
    }


# -------------------------
# Arrow constructor
# -------------------------
def make_arrow(x0, y0, x1, y1, flux):
    width = 1 + 10 * flux
    flux_norm = norm(flux)
    color = flux_color(flux_norm)
    return {
        "type": "arrow",
        "x0": x0,
        "y0": y0,
        "x1": x1,
        "y1": y1,
        "width": width,
        "color": color
    }


# -------------------------
# Node renderer
# -------------------------
def render_node(node):
    return go.Scatter(
        x=[node["x"]],
        y=[node["y"]],
        mode="markers+text",
        marker=dict(size=node["size"], color=node["color"]),
        text=node["name"],
        textposition="top center"
    )


# -------------------------
# Arrow renderer
# -------------------------
def render_arrow_shape(arrow):
    return dict(
        type="line",
        x0=arrow["x0"],
        y0=arrow["y0"],
        x1=arrow["x1"],
        y1=arrow["y1"],
        line=dict(width=arrow["width"], color=arrow["color"]),
        layer="below"
    )
