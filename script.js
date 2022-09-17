
let mesh = createMesh(10,1,20,2)
let coor = new Matrix(4,2);
coor.vals = [0,0, mesh.l,0, mesh.l,mesh.h, 0,mesh.h];

let E0 = 200e9;
let v = 0.3;
let E = ElasticMaterialMatrix(E0,v);
let Ke = ElementStiffnessMatrix(coor);

let F = new Matrix(8,1);
F.set(6,1,-13e6);

[1,2,7,8].forEach(dof => {
    for(let i = 0; i < 8 ; i++) {
        let row = (dof-1)*8 + i;
        let col = i*8 + (dof-1);
        Ke.vals[row] = 0;
        Ke.vals[col] = 0;
    }
    Ke.vals[(dof-1)*(8+1)]=1;
})

let U = Ke.solve(F);
U.disp();

