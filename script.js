let fem;
let slider;

function setup() {
    createCanvas(window.innerWidth,window.innerHeight);
    background(25);

    fem = new FEM();
    fem.viewMargins();
    fem.setGeometry(5,1);
    fem.setNumberOfElements(10,2);
    fem.createMesh();
    fem.setMaterial(200e9,0.3,"stress");
    fem.createMatrices();

    let constraints = {
        points: [{x: 0, y: 1, indicator: "xy"}],
        lines: [{xy: 0, orientation: "v", indicator: "xy"}],
    };

    let loads = {
        points: [{x: 5, y: 1, indicator: "y", val: 0}],
        lines: [{xy: 5, orientation: "v", indicator: "y", val: -10e4 }],
    };

    fem.applyBoundaryConditions(constraints,loads);
    fem.solve();

    slider = createSlider(-1, 5, 0);
    slider.position(10,10);
}

function draw() {
    background(25);
    let fact = 10**slider.value();
    
    fem.viewMargins();
    fem.viewMesh(fem.mesh.nelx,fem.mesh.nely);
    fem.viewDisplacements(fact);
}