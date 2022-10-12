class FEM {
    
    //CANVAS
    constructor() {
        this.canvas = {
            minMargin: 50,
            margins: Array(4),
        } 
    }

    viewMargins() {
        let margin = this.canvas.minMargin - 10;
        noFill();
        stroke(255);
        strokeWeight(1);
        rect(margin,margin,width-2*margin,height-2*margin);
    }

    // GEOMETRY
    setGeometry(L,H) {
        this.geometry = {L,H};
        this.scaleGeometry();
        this.viewGeometry();
    }

    scaleGeometry() {
        let minMargin = this.canvas.minMargin;
        let marginY = height/2 - (this.geometry.H/this.geometry.L)*(width -2*minMargin)/2;
        let marginX = width/2  - (this.geometry.L/this.geometry.H)*(height-2*minMargin)/2
        if (marginY > minMargin || marginX < minMargin) {
            this.canvas.margins[0] = minMargin;
            this.canvas.margins[2] = width - minMargin;
            this.canvas.margins[1] = height - marginY;
            this.canvas.margins[3] = marginY;
        } else {
            this.canvas.margins[0] = marginX;
            this.canvas.margins[2] = width - marginX;
            this.canvas.margins[1] = height - minMargin;
            this.canvas.margins[3] = minMargin;
        }
    }

    viewGeometry() {
        noFill();
        stroke(255);
        strokeWeight(1);
        let margins = this.canvas.margins;
        rect(margins[0],margins[3],margins[2]-margins[0],margins[1]-margins[3]);
    }

    // MESH
    setNumberOfElements(nelx,nely) {
        this.mesh = {nelx,nely}
        this.viewMesh();
    }

    viewMesh() {
        noFill();
        stroke(255);
        strokeWeight(1);
        let margins = this.canvas.margins;
        let dx = (margins[2]-margins[0])/this.mesh.nelx;
        let dy = (margins[1]-margins[3])/this.mesh.nely;
        for (let i = 0 ; i < this.mesh.nelx + 1 ; i++) {
            line(margins[0]+i*dx,margins[3],margins[0]+i*dx,margins[1]);
        }
        for (let j = 0 ; j < this.mesh.nely + 1 ; j++) {
            line(margins[0],margins[3]+j*dy,margins[2],margins[3]+j*dy);
        }
    }

    createMesh() {
        let nelx = this.mesh.nelx;
        let nely = this.mesh.nely;

        let l = this.geometry.L/nelx;
        let h = this.geometry.H/nely;
    
        let nodx = nelx+1;
        let nody = nely+1;
        let nDOF = nodx*nody*2;
    
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
        index = 0;
        for(let j = 0 ; j < nely ; j++) {
            for(let i = 0 ; i < nelx ; i++) {
                con.vals[index  ] =     j*nodx+i;
                con.vals[index+1] =     j*nodx+i+1;
                con.vals[index+2] = (j+1)*nodx+i+1;
                con.vals[index+3] = (j+1)*nodx+i;
                index+=4;
            }
        }
        this.mesh.l = l;
        this.mesh.h = h;
        this.mesh.nelx = nelx;
        this.mesh.nely = nely;
        this.mesh.nodx = nodx;
        this.mesh.nody = nody;
        this.mesh.nDOF = nDOF;
        this.mesh.con = con;
        this.mesh.x = x;
        this.mesh.y = y;
    }

    // MATERIAL
    setMaterial(E0,v,caseOf) {
        this.material = {E0,v,caseOf};
    }

    // MATRICES
    createMatrices() {
        this.matrices = {};
        this.matrices.E = this.elasticMaterialMatrix();
        this.matrices.Ke = this.elementStiffnessMatrix();
        let globalStiffness = this.globalStiffnessMatrix();
        this.matrices.Kg = globalStiffness.Kg;
        this.matrices.Kgm = globalStiffness.Kgm;
        this.matrices.F = new Matrix(this.mesh.nDOF,1);
    }

    elasticMaterialMatrix() {
        if (this.material.caseOf != "stress" && this.material.caseOf != "strain") {
            console.log("unexpected material case");
            return;
        }
        let v  = this.material.v;
        let E = new Matrix(3,3);
        if (this.material.caseOf == "stress") {
            E.set(1,1,1/(1-v*v));
            E.set(2,2,1/(1-v*v));
            E.set(2,1,v/(1-v*v));
            E.set(1,2,v/(1-v*v));
            E.set(3,3,0.5*(1-v)/(1-v*v));
            return E;
        }
        E.set(1,1,(1-v)/(1+v)/(1-2*v));
        E.set(2,2,(1-v)/(1+v)/(1-2*v));
        E.set(2,1,v/(1+v)/(1-2*v));
        E.set(1,2,v/(1+v)/(1-2*v));
        E.set(3,3,0.5*(1-2*v)/(1+v)/(1-2*v));
        return E;
    }

    elementStiffnessMatrix() {
        let coor = new Matrix(4,2);
        coor.vals = [0,0, this.mesh.l,0, this.mesh.l,this.mesh.h, 0,this.mesh.h];
        let xi = [ -Math.sqrt(3/5) , 0 , Math.sqrt(3/5) ];
        let wi = [5/9 , 8/9 , 5/9];
        let dNdst = new Matrix(2,4);
        let B = new Matrix(3,8);
        let Ke = new Matrix(8,8);
        xi.forEach((s,indexS) => {
            xi.forEach((t,indexT) => {
                dNdst.vals = [-(1-t)/4, (1-t)/4,(1+t)/4,-(1+t)/4,
                              -(1-s)/4,-(1+s)/4,(1+s)/4, (1-s)/4];
                let J = dNdst.mult(coor);
                let dNdxy = J.solve(dNdst);
                for (let col = 0; col < 4 ; col++) {                
                    B.set(1, col*2+1, dNdxy.vals[col]);
                    B.set(2, col*2+2, dNdxy.vals[col+4]);
                    B.set(3, col*2+1, dNdxy.vals[col+4]);
                    B.set(3, col*2+2, dNdxy.vals[col]);
                }
                Ke = Ke.add( B.transpose().mult(this.matrices.E).mult(B).mult(
                    J.determinant()*wi[indexS]*wi[indexT]
                ));
            })
        });
        return Ke;
    }

    globalStiffnessMatrix() {
        let Kg = new Matrix(this.mesh.nDOF,this.mesh.nDOF);
        let Kgm = new Matrix(this.mesh.nDOF,this.mesh.nDOF);
        for (let i = 0 ; i < this.mesh.con.rows ; i++) {
            let nods = this.mesh.con.vals.slice(i*4,(i+1)*4);
            let dofs = Array(8);
            nods.forEach((nod,index) => {
                dofs[index*2] = nod*2;
                dofs[index*2+1] = nod*2+1;
            });
            dofs.forEach((row,indexRow) => {
                dofs.forEach((col,indexCol) => {
                    let globalIndex = row*this.mesh.nDOF + col;
                    let localIndex = indexRow*8 + indexCol;
                    Kg.vals[globalIndex] += this.material.E0*this.matrices.Ke.vals[localIndex];
                    Kgm.vals[globalIndex] += this.material.E0*this.matrices.Ke.vals[localIndex];
                });
            });
        }
        return {Kg,Kgm};
    }

    // BOUNDARY CONDITIONS
    applyBoundaryConditions(constraints,loads) {
        let loadedDOF = this.getDOF(loads,"load");
        if (loadedDOF.length <= 0) {
            console.log("No-loads were applied");
            return;
        };
        let constrainedDOF = this.getDOF(constraints,"constraint");
        if (constrainedDOF.length <= 0) {
            console.log("No-constraints were applied");
            return;
        };
        this.applyLoads(loadedDOF,constrainedDOF);
        this.applyConstraints(constrainedDOF);
        this.viewConstraints(constrainedDOF);
        this.viewLoads();
    }

    getDOF(condition,typeOfCondition) {
        let output = [];
        condition.points.forEach(point=>{
            let nod = this.point2NOD(point);
            if(point.indicator.includes("x")) {
                (typeOfCondition=="load") ? output.push({dof: nod*2  , val: point.val}) : output.push({dof: nod*2})
            }
            if(point.indicator.includes("y")) {
                (typeOfCondition=="load") ? output.push({dof: nod*2+1, val: point.val}) : output.push({dof: nod*2+1})
            }
        })
        condition.lines.forEach(line=>{
            let nods = this.line2NOD(line);
            nods.forEach(nod=>{
                if(line.indicator.includes("x")) {
                    (typeOfCondition=="load") ? output.push({dof: nod*2  , val: line.val/nods.length}) : output.push({dof: nod*2})
                }
                if(line.indicator.includes("y")) {
                    (typeOfCondition=="load") ? output.push({dof: nod*2+1, val: line.val/nods.length}) : output.push({dof: nod*2+1})
                }
            });
        });
        return output
    }

    point2NOD(point) {
        return this.mesh.x.reduce((acc,currX,index)=>{
            let dis1 = (currX - point.x)**2 + (this.mesh.y[index] - point.y)**2;
            let dis2 = (this.mesh.x[acc] - point.x)**2 + (this.mesh.y[acc] - point.y)**2;
            return dis1 < dis2 ? index : acc
        },0);
    }

    line2NOD(line) {
        let tol = 0.03;
        let dxy  = line.orientation == "v" ? this.mesh.l/2 : this.mesh.h/2;
        let arr = line.orientation == "v" ? this.mesh.x   : this.mesh.y;
        return arr.reduce((acc,curr,index)=>{
            if (curr >= line.xy-dxy*(1+tol) && curr <= line.xy+dxy*(1+tol)) acc.push(index);
            return acc;
        },[]);
    }

    applyConstraints(constrainedDOF) {
        constrainedDOF.forEach(condition => {
            let dof = condition.dof;
            for(let i = 0; i < this.mesh.nDOF ; i++) {
                let row = dof*this.mesh.nDOF + i;
                let col = i*this.mesh.nDOF + dof;
                this.matrices.Kgm.vals[row] = 0;
                this.matrices.Kgm.vals[col] = 0;
            }
            this.matrices.Kgm.vals[dof*(this.mesh.nDOF+1)]=1;
        });
    }

    applyLoads(loadedDOF,constrainedDOF) {
        constrainedDOF = constrainedDOF.map(condition=>condition.dof);
        loadedDOF.forEach(condition=>{
            if (!constrainedDOF.includes(condition.dof)) {
                this.matrices.F.vals[condition.dof] += condition.val;
            }
        });
    }

    viewConstraints(constrainedDOF) {
        noFill();
        stroke(255,100,100);
        strokeWeight(2);
        let scale = (this.canvas.margins[2]-this.canvas.margins[0])/this.geometry.L;
        let dl = Math.min(this.mesh.l,this.mesh.h)/10*scale;
        constrainedDOF.forEach(condition=>{
            let index;
            if(condition.dof%2==0) {
                index = condition.dof/2;
                let x = map(this.mesh.x[index],0,this.geometry.L,this.canvas.margins[0],this.canvas.margins[2]);
                let y = map(this.mesh.y[index],0,this.geometry.H,this.canvas.margins[1],this.canvas.margins[3]);
                triangle(x,y,x-dl*0.87,y+dl*0.5,x-dl*0.87,y-dl*0.5);
            } else {
                index = (condition.dof-1)/2;
                let x = map(this.mesh.x[index],0,this.geometry.L,this.canvas.margins[0],this.canvas.margins[2]);
                let y = map(this.mesh.y[index],0,this.geometry.H,this.canvas.margins[1],this.canvas.margins[3]);
                triangle(x,y,x-dl*0.5,y+dl*0.87,x+dl*0.5,y+dl*0.87);
            }
        });
    }
    
    viewLoads() {
        noFill();
        stroke(100,255,100);
        let dx = (this.canvas.margins[2]-this.canvas.margins[0])/this.mesh.nelx*0.85;
        let dy = (this.canvas.margins[1]-this.canvas.margins[3])/this.mesh.nely*0.85;
        let dl = Math.min(dx,dy)*0.10;
        let maxForce = this.matrices.F.vals.filter(force => force != 0).map(force => Math.abs(force)).sort((a,b)=>b-a)[0];
        let xScale = dx/maxForce;
        let yScale = dy/maxForce;
        this.matrices.F.vals.forEach((force,dof)=>{
            if (force != 0) {
                let index = (dof%2==0) ? dof/2 : (dof-1)/2;
                let dir = dl*Math.sign(force);
                let x = map(this.mesh.x[index],0,this.geometry.L,this.canvas.margins[0],this.canvas.margins[2]);
                let y = map(this.mesh.y[index],0,this.geometry.H,this.canvas.margins[1],this.canvas.margins[3]);
                strokeWeight(2);
                ellipse(x,y,dl);
                strokeWeight(2.25);
                if(dof%2==0) {
                    line(x-dir/2,y,x-dir/2-force*xScale,y);
                    line(x-dir/2,y,x-dir,y-dir/2);
                    line(x-dir/2,y,x-dir,y+dir/2);
                } else {
                    line(x,y+dir/2,x,y+dir/2+force*yScale);
                    line(x,y+dir/2,x-dir/2,y+dir);
                    line(x,y+dir/2,x+dir/2,y+dir);
                }
            }
        });
    }

    // SOLVE
    solve() {
        let U = this.matrices.Kgm.solve(this.matrices.F);
        this.solved = {U};
    }

    viewDisplacements(fact) {
        noFill();
        stroke(255,100,100);
        strokeWeight(1);
        for(let row = 0 ; row < this.mesh.con.rows ; row++) {
            let element = this.mesh.con.vals.slice(row*4,(row+1)*4);
            element.forEach((nod,index) => {
    
                let u1 = this.mesh.x[nod] + fact*this.solved.U.vals[nod*2];
                let u2 = this.mesh.x[element[(index+1)%4]] + fact*this.solved.U.vals[element[(index+1)%4]*2]
    
                let v1 = this.mesh.y[nod] + fact*this.solved.U.vals[nod*2+1];
                let v2 = this.mesh.y[element[(index+1)%4]] + fact*this.solved.U.vals[element[(index+1)%4]*2+1]
    
                let x1 = map(u1,0,this.geometry.L,this.canvas.margins[0],this.canvas.margins[2]);
                let x2 = map(u2,0,this.geometry.L,this.canvas.margins[0],this.canvas.margins[2]);
    
                let y1 = map(v1,0,this.geometry.H,this.canvas.margins[1],this.canvas.margins[3]);
                let y2 = map(v2,0,this.geometry.H,this.canvas.margins[1],this.canvas.margins[3]);
    
                line(x1,y1,x2,y2);
            })
        }
    }
}