// Flux color gradient helper
function fluxColor(f) {
    const clamped = Math.min(f, 1.0);

    const r = Math.floor(255 * clamped);
    const g = Math.floor(255 * (1 - clamped));
    const b = 50;

    return (r << 16) + (g << 8) + b;
}


// ------------------------------------------------------------
// Load metabolic layout (node positions + edges)
// ------------------------------------------------------------
async function loadLayout() {
    const res = await fetch("layout.json");
    return await res.json();
}

// ------------------------------------------------------------
// Create PixiJS application
// ------------------------------------------------------------
const app = new PIXI.Application({
    width: 1300,
    height: 600,
    backgroundColor: 0x1e1e1e,
});
document.body.appendChild(app.view);

// Containers
let nodes = {};
let edges = [];

// ------------------------------------------------------------
// Create nodes from layout.json
// ------------------------------------------------------------
function createNodes(nodePositions) {
    nodes = {};

    for (const [name, [x, y]] of Object.entries(nodePositions)) {
        const g = new PIXI.Graphics();
        g.beginFill(0x00aaff);
        g.drawCircle(0, 0, 20);
        g.endFill();

        g.x = x;
        g.y = y;

        // --- Add label ---
        const label = new PIXI.Text(name, {
            fontFamily: "Arial",
            fontSize: 14,
            fill: 0xffffff
        });
        label.x = x + 25;
        label.y = y - 10;

        nodes[name] = { gfx: g, label: label };
        app.stage.addChild(g);
        app.stage.addChild(label);
    }
}

// ------------------------------------------------------------
// Create edges (arrows) from layout.json
// ------------------------------------------------------------
function createEdges(edgeList) {
    edges = [];

    for (const [src, dst] of edgeList) {
        const line = new PIXI.Graphics();
        line.lineStyle(4, 0xffffff, 1.0);

        const s = nodes[src];
        const d = nodes[dst];

        line.moveTo(s.x, s.y);
        line.lineTo(d.x, d.y);

        edges.push({ src, dst, gfx: line });
        app.stage.addChild(line);
    }
}




const fluxMap = {
    "Glc_ext→Glc": "GLUT4",
    "Glc→Glc6P": "Glycolysis",
    "Fru16P2→Pyruvate": "Glycolysis",
    "Pyruvate→AcCoA_Glc": "PDH",
    "AcCoA_Glc→Cit": "Citrate_prod",
    "Cit→Mal": "ACL_ACC",

    "LCFA_ext→LCFA_CoA_cyto": "CD36",
    "LCFA_CoA_cyto→LCFA_CoA_mito": "CPT1",
    "LCFA_CoA_mito→AcCoA_Fat": "BetaOx"
};


// ------------------------------------------------------------
// WebSocket connection to backend
// ------------------------------------------------------------
const ws = new WebSocket("ws://localhost:8765");

const gluSlider = document.getElementById("gluSlider");
const fatSlider = document.getElementById("fatSlider");
const atpSlider = document.getElementById("atpSlider");

function sendInputs() {
    const msg = {
        glu_in: parseFloat(gluSlider.value),
        fat_in: parseFloat(fatSlider.value),
        atp_draw: parseFloat(atpSlider.value)
    };
    ws.send(JSON.stringify(msg));
}

gluSlider.oninput = sendInputs;
fatSlider.oninput = sendInputs;
atpSlider.oninput = sendInputs;


ws.onmessage = (event) => {
    const flux = JSON.parse(event.data);

    // --------------------------------------------------------
    // Node pulsing based on incoming flux
    // --------------------------------------------------------
    for (const [name, nodeObj] of Object.entries(nodes)) {

        // Find edges that end at this node
        const incoming = edges.filter(e => e.dst === name);

        let f_total = 0.0;

        for (const edge of incoming) {
            const key = `${edge.src}→${edge.dst}`;
            const fluxKey = fluxMap[key];
            if (fluxKey) {
                f_total += flux[fluxKey] || 0.0;
            }
        }

        const amplified = Math.min(f_total * 2.5, 1.0);

        nodeObj.gfx.scale.set(1 + 0.3 * amplified);
    }

    // --------------------------------------------------------
    // Edge thickness + color based on flux
    // --------------------------------------------------------
    for (const edge of edges) {
        const key = `${edge.src}→${edge.dst}`;
        const fluxKey = fluxMap[key];
        const f = fluxKey ? flux[fluxKey] || 0.0 : 0.0;

        const amplified = Math.min(f * 3.0, 1.0);

        const color = fluxColor(amplified);

        edge.gfx.clear();
        edge.gfx.lineStyle(2 + 10 * amplified, color, 0.9);

        const s = nodes[edge.src].gfx;
        const d = nodes[edge.dst].gfx;

        edge.gfx.moveTo(s.x, s.y);
        edge.gfx.lineTo(d.x, d.y);
    }
};


// ------------------------------------------------------------
// Initialize everything
// ------------------------------------------------------------
loadLayout().then(layout => {
    createNodes(layout.nodes);
    createEdges(layout.edges);
});
