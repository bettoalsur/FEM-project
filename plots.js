
let minMargin = 50;
let margins = Array(4);

function viewMesh() {
    noFill();
    stroke(255);
    strokeWeight(1);

    createMargins();
    let l = (margins[2]-margins[0])/mesh.nelx;
    let h = (margins[1]-margins[3])/mesh.nely;
    mesh.con.vals.forEach((nod,index)=>{
        if (index%4==3) {
            let xp = map(mesh.x[nod],0,mesh.L,margins[0],margins[2]);
            let yp = map(mesh.y[nod],0,mesh.H,margins[1],margins[3]);
            rect(xp,yp,l,h);
        }
    });
}

function viewMargins() {
    let margin = minMargin - 5;
    noFill();
    stroke(255);
    strokeWeight(1);
    rect(margin,margin,width-2*margin,height-2*margin);
}

function createMargins() {
    let marginY = height/2 - (mesh.H/mesh.L)*(width -2*minMargin)/2;
    let marginX = width/2  - (mesh.L/mesh.H)*(height-2*minMargin)/2
    if (marginY > minMargin || marginX < minMargin) {
        margins[0] = minMargin;
        margins[2] = width - minMargin;
        margins[1] = height - marginY;
        margins[3] = marginY;
    } else {
        margins[0] = marginX;
        margins[2] = width - marginX;
        margins[1] = height - minMargin;
        margins[3] = minMargin;
    }
}

function viewDisplacements() {
    noFill();
    stroke(255,100,100);
    strokeWeight(1);

    createMargins();

    for(let row = 0 ; row < mesh.con.rows ; row++) {
        let element = mesh.con.vals.slice(row*4,(row+1)*4);
        element.forEach((nod,index) => {

            let u1 = mesh.x[nod] + fact*U.vals[nod*2];
            let u2 = mesh.x[element[(index+1)%4]] + fact*U.vals[element[(index+1)%4]*2]

            let v1 = mesh.y[nod] + fact*U.vals[nod*2+1];
            let v2 = mesh.y[element[(index+1)%4]] + fact*U.vals[element[(index+1)%4]*2+1]

            let x1 = map(u1,0,mesh.L,margins[0],margins[2]);
            let x2 = map(u2,0,mesh.L,margins[0],margins[2]);

            let y1 = map(v1,0,mesh.H,margins[1],margins[3]);
            let y2 = map(v2,0,mesh.H,margins[1],margins[3]);

            line(x1,y1,x2,y2);
        })
    }
}