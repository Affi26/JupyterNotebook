// ------------------------------------------------------------
// Flux color gradient helper
// ------------------------------------------------------------

function fluxColor(f) {

    const clamped =
        Math.min(f, 1.0);

    const r =
        Math.floor(
            255 * (1 - clamped)
        );

    const g = 255;

    const b =
        Math.floor(
            255 * (1 - clamped)
        );

    return (r << 16)
         + (g << 8)
         + b;
}

// ------------------------------------------------------------
// Inhibition color helper
// ------------------------------------------------------------

function inhibitionColor(x) {

    const c = Math.min(x, 1.0);

    const r = 255;
    const g = Math.floor(100 * (1 - c));
    const b = Math.floor(100 * (1 - c));

    return (r << 16) + (g << 8) + b;
}

// ------------------------------------------------------------
// Concentration color helper
// ------------------------------------------------------------

function concentrationColor(c) {

    const x = Math.min(
        Math.log10(1 + 10 * c),
        1.0
    );

    const r = Math.floor(255 * x);
    const g = Math.floor(180 * (1 - x));
    const b = 80;

    return (r << 16) + (g << 8) + b;
}

// ------------------------------------------------------------
// Load metabolic layout
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

// ------------------------------------------------------------
// Containers
// ------------------------------------------------------------

let nodes = {};
let fluxEdges = [];
let inhibitionEdges = [];

// ------------------------------------------------------------
// Create nodes
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

        const label = new PIXI.Text(name, {

            fontFamily: "Arial",
            fontSize: 14,
            fill: 0xffffff

        });

        label.x = x + 25;
        label.y = y - 10;

        nodes[name] = {

            gfx: g,
            label: label

        };

        app.stage.addChild(g);
        app.stage.addChild(label);
    }
}

// ------------------------------------------------------------
// Create edges (flux and inhibition are separate)
// ------------------------------------------------------------

function createFluxEdges(edgeList) {

    fluxEdges = [];

    for (const [src, dst] of edgeList) {

        const line = new PIXI.Graphics();

        const s = nodes[src];
        const d = nodes[dst];

        line.lineStyle(2, 0xffffff, 0.5);

        line.moveTo(
            s.gfx.x,
            s.gfx.y
        );

        line.lineTo(
            d.gfx.x,
            d.gfx.y
        );

        fluxEdges.push({
            src,
            dst,
            gfx: line
        });

        app.stage.addChild(line);
    }
}

function createInhibitionEdges(edgeList) {

    inhibitionEdges = [];

    for (const [src, dst] of edgeList) {

        const line = new PIXI.Graphics();

        const s = nodes[src];
        const d = nodes[dst];

        line.lineStyle(2, 0x550000, 0.5);

        line.moveTo(
            s.gfx.x,
            s.gfx.y
        );

        line.lineTo(
            d.gfx.x,
            d.gfx.y
        );

        inhibitionEdges.push({
            src,
            dst,
            gfx: line
        });

        app.stage.addChild(line);
    }
}

// ------------------------------------------------------------
// Flux mapping
// ------------------------------------------------------------

const fluxMap = {

    "Glc_ext→Glc": "GLUT4",
    "Glc→Glc6P": "Glycolysis",
    "Fru6P→Fru16P2": "Glycolysis",
    "Fru16P2→Pyruvate": "Glycolysis",
    "Pyruvate→AcCoA": "PDH",

    "LCFA_CoA_mito→AcCoA": "BetaOx",

    "AcCoA→Cit": "Citrate_prod",

    "Cit→Mal": "ACL_ACC",

    "LCFA_ext→LCFA_CoA_cyto": "CD36",
    "LCFA_CoA_cyto→LCFA_CoA_mito": "CPT1",

    "Cit→Pyruvate": "Citrate_prod",
    "Cit→Fru6P": "Citrate_prod",
    "Cit→Glc_ext": "Citrate_prod",

    "Mal→LCFA_CoA_mito": "FAS",
    "Mal→LCFA_CoA_cyto": "FAS"
};

const inhibitionMap = {

    "Cit→Pyruvate": "PDH",

    "Cit→Fru6P": "Glycolysis",

    "Cit→Glc_ext": "Glycolysis",

    "Mal→LCFA_CoA_mito": "CPT1",

    "Mal→LCFA_CoA_cyto": "BetaOx"
};

// ------------------------------------------------------------
// WebSocket
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

// ------------------------------------------------------------
// Receive updates
// ------------------------------------------------------------

ws.onmessage = (event) => {

    const data = JSON.parse(event.data);

    const state = data.state || {};
    const dstate = data.dstate || {};
    const flux = data.flux || {};
    const inh = data.inhibition || {};

    // --------------------------------------------------------
    // Update nodes
    // --------------------------------------------------------

    for (const [name, nodeObj] of Object.entries(nodes)) {

        const conc = state[name] || 0.0;
        
        const rate = dstate[name] || 0.0;

        const incoming =
            fluxEdges.filter(
                e => e.dst === name
            );

        let f_total = 0.0;

        for (const edge of incoming) {

            const key =
                `${edge.src}→${edge.dst}`;

            const fluxKey =
                fluxMap[key];

            if (fluxKey) {

                f_total +=
                    flux[fluxKey] || 0.0;
            }
        }

        const amplified =
            Math.min(
                f_total * 2.5,
                1.0
            );

        const inhVal =
            inh[name] || 0.0;

        const shrink =
            1 -
            0.3 *
            Math.min(inhVal, 1.0);

        const concScale =
            1 +
            0.5 *
            Math.log10(
                1 + 10 * conc
            );

        nodeObj.gfx.tint =
            concentrationColor(conc);

        nodeObj.gfx.scale.set(
            concScale * shrink
        );

        nodeObj.label.text = `${name}, ${conc.toFixed(3)}, ${rate.toExponential(2)}`;
    }

    // --------------------------------------------------------
    // Update edges
    // --------------------------------------------------------

    for (const edge of fluxEdges) {
    
        const key =
            `${edge.src}→${edge.dst}`;
    
        const fluxKey =
            fluxMap[key];
    
        const f =
            fluxKey
                ? flux[fluxKey] || 0.0
                : 0.0;
    
        const amplified =
            Math.min(f * 3.0, 1.0);
    
        edge.gfx.clear();
    
        edge.gfx.lineStyle(
            2 + 10 * amplified,
            fluxColor(amplified),
            0.9
        );
    
        edge.gfx.moveTo(
            nodes[edge.src].gfx.x,
            nodes[edge.src].gfx.y
        );
    
        edge.gfx.lineTo(
            nodes[edge.dst].gfx.x,
            nodes[edge.dst].gfx.y
        );
    }


    for (const edge of inhibitionEdges) {
    
        const key =
            `${edge.src}→${edge.dst}`;
    
        const inhKey =
            inhibitionMap[key];
    
        const inhVal =
            inhKey
                ? inh[inhKey] || 0.0
                : 0.0;
    
        edge.gfx.clear();
    
        edge.gfx.lineStyle(
            2 + 10 * inhVal,
            inhibitionColor(inhVal),
            0.9
        );
    
        edge.gfx.moveTo(
            nodes[edge.src].gfx.x,
            nodes[edge.src].gfx.y
        );
    
        edge.gfx.lineTo(
            nodes[edge.dst].gfx.x,
            nodes[edge.dst].gfx.y
        );
    }
    
};

// ------------------------------------------------------------
// Initialize
// ------------------------------------------------------------

loadLayout().then(layout => {

    createNodes(layout.nodes);

    console.log(layout);
    console.log(layout.flux_edges);
    console.log(layout.inhibition_edges);

    createFluxEdges(
        layout.flux_edges
    );

    createInhibitionEdges(
        layout.inhibition_edges
    );

});