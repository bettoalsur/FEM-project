function createMesh(L,H,nelx,nely) {
    let l = L/nelx;
    let h = H/nely;

    let nodx = nelx+1;
    let nody = nely+1;
    let ndof = nodx*nody*2;

    let x = Array(nodx*nody);
    let y = Array(nodx*nody);
    let index = 0;
    for (let j = 0; j < nody; j++){
        for (let i = 0; i < nodx; i++){
            x[index] = i*l;
            y[index] = j*h;
            index++;
        }
    }

    let con = new Matrix(nelx*nely,4);
    let aux = [];
    for(let j = 0 ; j < nely ; j++) {
        for(let i = 0 ; i < nelx ; i++) {
            aux = aux.concat( [j*nodx+i,j*nodx+i+1,(j+1)*nodx+i+1,(j+1)*nodx+i] );
        }
    }
    con.vals = aux;
    return {L,H,l,h,nelx,nely,nodx,nody,ndof,con,x,y};
}

function elasticMaterialMatrix() {
    if (material.caseOf != "stress" && material.caseOf != "strain") {
        console.log("unexpected material case");
        return;
    }
    let E0 = material.E0;
    let v  = material.v;
    let E = new Matrix(3,3);
    if (material.caseOf == "stress") {
        E.set(1,1,1);
        E.set(2,2,1);
        E.set(2,1,v);
        E.set(1,2,v);
        E.set(3,3,0.5*(1-v));
        E = E.mult(E0/(1-v*v));
    } else {
        E.set(1,1,1-v);
        E.set(2,2,1-v);
        E.set(2,1,v);
        E.set(1,2,v);
        E.set(3,3,0.5*(1-2*v));
        E = E.mult(E0/(1+v)/(1-2*v));
    }
    return E;
}

function elementStiffnessMatrix() {
    let coor = new Matrix(4,2);
    coor.vals = [0,0, mesh.l,0, mesh.l,mesh.h, 0,mesh.h];
    let xi = [ -Math.sqrt(3/5) , 0 , Math.sqrt(3/5) ];
    let wi = [5/9 , 8/9 , 5/9];
    let Ke = new Matrix(8,8);
    let E = elasticMaterialMatrix();
    xi.forEach((s,indexS) => {
        xi.forEach((t,indexT) => {
            let dNdst = new Matrix(2,4);
            dNdst.vals = [-(1-t), (1-t),(1+t),-(1+t),
                          -(1-s),-(1+s),(1+s), (1-s)];
            dNdst = dNdst.mult(1/4);
            
            let J = dNdst.mult(coor);
            let dNdxy = J.solve(dNdst);
            
            let B = new Matrix(3,8);
            for (let col = 0; col < 4 ; col++) {
                let indexX = col;
                let indexY = 4 + col;
                B.set(1, col*2+1, dNdxy.vals[indexX]);
                B.set(2, col*2+2, dNdxy.vals[indexY]);
                B.set(3, col*2+1, dNdxy.vals[indexY]);
                B.set(3, col*2+2, dNdxy.vals[indexX]);
            }
            let escalar = J.determinant()*wi[indexS]*wi[indexT];
            Ke = Ke.add( B.transpose().mult(E).mult(B).mult( escalar ) );
        })
    });
    return Ke;
}

function globalStiffnessMatrix() {
    let Ke = elementStiffnessMatrix();
    let Kg = new Matrix(mesh.ndof,mesh.ndof);
    let Kgm = new Matrix(mesh.ndof,mesh.ndof);
    for (let i = 0 ; i < mesh.con.rows ; i++) {
        let nods = mesh.con.vals.slice(i*4,(i+1)*4);
        let dofs = Array(8);
        nods.forEach((nod,index) => {
            dofs[index*2] = nod*2;
            dofs[index*2+1] = nod*2+1;
        });
        dofs.forEach((row,indexRow) => {
            dofs.forEach((col,indexCol) => {
                let globalIndex = row*mesh.ndof + col;
                let localIndex = indexRow*8 + indexCol;
                Kg.vals[globalIndex] += Ke.vals[localIndex];
                Kgm.vals[globalIndex] += Ke.vals[localIndex];
            });
        });
    }
    return {Kg,Kgm};
}

function applyConstraints() {
    let DOFconstrained = boundaryConditions.DOFconstrained;
    let dofx = DOFconstrained.dofx.concat(DOFconstrained.dofxy).map(dof => dof*2);
    let dofy = DOFconstrained.dofy.concat(DOFconstrained.dofxy).map(dof => dof*2+1);
    dofx.concat(dofy).forEach(dof => {
        for(let i = 0; i < mesh.ndof ; i++) {
            let row = dof*mesh.ndof + i;
            let col = i*mesh.ndof + dof;
            Kgm.vals[row] = 0;
            Kgm.vals[col] = 0;
        }
        Kgm.vals[dof*(mesh.ndof+1)]=1;
    });
}

function applyLoads() {
    let loads = boundaryConditions.loads;
    loads.Fx.forEach(row => {
        F.vals[row[0]*2] = row[1];
    });
    loads.Fy.forEach(row => {
        F.vals[row[0]*2+1] = row[1];
    });
    let DOFconstrained = boundaryConditions.DOFconstrained;
    let dofx = DOFconstrained.dofx.concat(DOFconstrained.dofxy).map(dof => dof*2);
    let dofy = DOFconstrained.dofy.concat(DOFconstrained.dofxy).map(dof => dof*2+1);
    dofx.concat(dofy).forEach(dof => {
        F.vals[dof] = 0;
    });
}

function applyBoundaryConditions() {
    applyConstraints();
    applyLoads();
}