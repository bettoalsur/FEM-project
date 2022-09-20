let mesh = createMesh(5,1,10,2);

let material = {
    E0: 200e9,
    v : 0.3,
    caseOf: "stress"
}

let stiffnessMatrix = globalStiffnessMatrix();
let Kg = stiffnessMatrix.Kg;
let Kgm = stiffnessMatrix.Kgm;
let F = new Matrix(mesh.ndof,1);

let boundaryConditions = new Object;

boundaryConditions.DOFconstrained = {
    dofy: [],
    dofx: [],
    dofxy: [1,11,22]
};

boundaryConditions.loads = {
    Fx: [
        []
    ],
    Fy: [
        [mesh.nodx*mesh.nody-1 , -10e4]
    ]
}

applyBoundaryConditions();

let U = Kgm.solve(F);

let slider;
function setup(){
    slider = createSlider(0, 5, 0);
    slider.position(10,10);
}

let fact;
function draw() {
    fact = 10**slider.value();
    createCanvas(window.innerWidth,window.innerHeight);
    background(25);

    viewMargins();
    viewMesh();
    viewDisplacements();
}