
let mesh = createMesh(4,2,30,10);



let E = elasticMaterialMatrix({
    E0: 200e9,
    v : 0.3
},"stress");

let Ke = elementStiffnessMatrix();
let Kg = globalStiffnessMatrix();

let Kgm = applyBoundaryConditions({
    dofx: [0,mesh.nodx-1],
    dofy: [0,mesh.nodx-1]
});

let Fx = [];
let Fy = [];
let F = new Matrix(mesh.ndof,1);
// F.set(11*2+2,1,-13e6);
F.set(mesh.ndof-2*mesh.nelx/2,1,-13e6);
// F.set(13*2+2,1,-13e6);

let U = Kgm.solve(F);

function setup() {
    let tol = 50;
    let fact = 5e1;
    createCanvas(window.innerWidth, window.innerHeight);
    background(25);
    mesh.x.forEach((x,index) => {
        strokeWeight(2);

        let xp = map(x,0,mesh.L,tol,width-tol);
        let yp = map(mesh.y[index],0,mesh.H,height-tol,tol);
        stroke(255);
        point(xp,yp);
        
        let xpf = map(x+fact*U.vals[index*2],0,mesh.L,tol,width-tol);
        let ypf = map(mesh.y[index]+fact*U.vals[index*2+1],0,mesh.H,height-tol,tol);
        stroke(255,100,100);
        point(xpf,ypf);
        
        fill(255);
        noStroke();
        //text(index,xp+5,yp);
    });
}