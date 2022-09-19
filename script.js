
// let mesh = createMesh(1,5,4,20);
let mesh = createMesh(5,1,20,4);

let E = elasticMaterialMatrix({
    E0: 200e9,
    v : 0.3
},"stress");

let Ke = elementStiffnessMatrix();
let Kg = globalStiffnessMatrix();

/* let Kgm = applyBoundaryConditions({
    dofx: [0,1,2,3,4],
    dofy: [0,1,2,3,4]
}); */
let Kgm = applyBoundaryConditions({
    dofy: [0,21,42,63,84].concat([20,41,62,83,104]),
    dofx: [42,62]
}); 
/* let Kgm = applyBoundaryConditions({
    dofy: [0,21,42,63,84],
    dofx: [0,21,42,63,84]
});  */

let Fx = [];
let Fy = [];
let F = new Matrix(mesh.ndof,1);
/* F.set(mesh.ndof,1,-2e4);
F.set(mesh.ndof-2,1,-2e4);
F.set(mesh.ndof-4,1,-2e4);
F.set(mesh.ndof-6,1,-2e4);
F.set(mesh.ndof-8,1,-2e4); */

F.set(mesh.ndof-9*2,1,-10e4);
// F.set(mesh.ndof,1,-10e4);

let U = Kgm.solve(F);

let slider;
function setup(){
    slider = createSlider(1, 1e4, 1);
    slider.position(10,window.innerHeight-50);
}

let fact;
function draw() {
    fact = slider.value();
    createCanvas(window.innerWidth,window.innerHeight);
    background(25);
    viewMesh();
    viewDisplacements();
}
