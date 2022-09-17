function ElasticMaterialMatrix(E0,v) {
    let E = new Matrix(3,3);
    E.set(1,1,1);
    E.set(2,2,1);
    E.set(2,1,v);
    E.set(1,2,v);
    E.set(3,3,0.5*(1-v));
    E = E.mult(E0/(1-v*v));
    return E;
}

function ElementStiffnessMatrix(coor) {
    let xi = [ -Math.sqrt(3/5) , 0 , Math.sqrt(3/5) ];
    let wi = [5/9 , 8/9 , 5/9];
    let Ke = new Matrix(8,8);
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

function createMesh(L,H,nelx,nely) {
    let l = L/nelx;
    let h = H/nely;

    let nodx = nelx+1;
    let nody = nely+1;

    let con = new Matrix(nelx*nely,4);
    let aux = [];
    for(let j = 0 ; j < nely ; j++) {
        for(let i = 0 ; i < nelx ; i++) {
            aux = aux.concat( [j*nodx+i,j*nodx+i+1,(j+1)*nodx+i+1,(j+1)*nodx+i] );
        }
    }
    con.vals = aux;
    return {L,H,l,h,nelx,nely,nodx,nody,con};
}