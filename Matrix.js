class Matrix {
    constructor(rows,cols) {
        this.rows = rows;
        this.cols = cols;
        this.vals = new Array(rows*cols);
        for(let i = 0; i < this.vals.length; i++) this.vals[i] = 0;
    }

    disp() {
        console.log("[ ]:");
        for(let i = 0; i < this.rows; i++) console.log(this.vals.slice(i*this.cols,(i+1)*this.cols));
    }

    set(row,col,val) {
        if (row <= 0 || row > this.rows || col <= 0 || col > this.cols) {
            console.log("position out of range");
            return;
        }
        let index = (row-1)*(this.cols) + (col-1);
        this.vals[index] = val;
    }

    transpose() {
        let C = new Matrix(this.cols,this.rows);
        let indexC = 0;
        for (let i = 0 ; i < this.cols ; i++) {
            for (let j = 0; j < this.rows ; j++) {
                let indexA = j*this.cols + i;
                C.vals[indexC] = this.vals[indexA];
                indexC++;
            }
        }
        return C;
    }

    determinant() {
        if (this.rows != this.cols ) {
            console.log("Matrix must be square");
            return;
        }
        let size = this.rows;
        if (size == 1) return this.vals[0];
        if (size == 2) return this.vals[0]*this.vals[3]-this.vals[1]*this.vals[2];
        let sign = 1;
        let result = 0;
        for (let k = 0 ; k < size ; k++){
            let cofactors = new Matrix(size-1,size-1);
            cofactors.vals = this.vals.slice(size,size*size).filter((_,index)=>index%size!=k);
            result += sign * this.vals[k] * cofactors.determinant();
            sign*=(-1);
        }
        return result;
    }

    add(B) {
        if ( !(B instanceof Matrix) && typeof(B) != "number") {
            console.log("invalid input");
            return;
        }
        if ( typeof(B) == "number" ) {
            let C = new Matrix(this.rows,this.cols);
            C.vals = this.vals.map(x => x+B);
            return C;
        } 
        if(this.rows != B.rows || this.cols != B.cols) {
            console.log("sizes are inconsistent");
            return;
        }
        let C = new Matrix(this.rows,this.cols);
        this.vals.forEach((val,index) => {
            C.vals[index] = val + B.vals[index];
        });
        return C;
    }

    mult(B) {
        if ( !(B instanceof Matrix) && typeof(B) != "number") {
            console.log("invalid input");
            return;
        }
        if ( typeof(B) == "number" ) {
            let C = new Matrix(this.rows,this.cols);
            C.vals = this.vals.map(x => x*B);
            return C;
        } 
        if (this.cols != B.rows) {
            console.log("sizes are inconsistent");
            return;
        }
        let C = new Matrix(this.rows,B.cols);
        let size = this.cols;
        for (let i = 0 ; i < C.rows ; i++) {
            for (let j = 0 ; j < C.cols ; j++) {
                let indexC = i*C.cols + j;
                for (let k = 0; k < size ; k++) {
                    let indexA = i*this.cols + k;
                    let indexB = k*B.cols + j;
                    C.vals[indexC] += this.vals[indexA]*B.vals[indexB];
                }
            }
        }
        return C;
    }

    fillrandi(N) {
        this.vals = this.vals.map(x => Math.ceil(Math.random()*N) );
    }

    filldiag(N) {
        if (this.rows != this.cols ) {
            console.log("Matrix must be square");
            return;
        }
        let size = this.rows;
        for (let i = 0 ; i < size ; i++){
            let index = i*this.cols + i;
            this.vals[index] = N;
        }
    }

    LU() {
        if (this.rows != this.cols ) {
            console.log("Matrix must be square");
            return;
        }
        let size = this.rows;
        let L = new Matrix(size,size);
        let U = new Matrix(size,size);
        let Acopy = new Matrix(size,size);
        Acopy.vals = this.vals.map(x=>x);
        L.filldiag(1);
        for (let k = 0 ; k < size ; k++) {
            U.vals[k*size + k] = Acopy.vals[k*size + k];
            for (let i = k+1 ; i < size ; i++) {
                L.vals[i*size + k] = Acopy.vals[i*size + k] / U.vals[k*size + k];
                U.vals[k*size + i] = Acopy.vals[k*size + i];
            }
            for (let i = k+1 ; i < size ; i++) {
                for (let j = k+1 ; j < size ; j++) {
                    Acopy.vals[i*size + j] -= L.vals[i*size + k]*U.vals[k*size + j];
                }
            }
        }
        return {L,U};
    }

    solve(b) {
        if (this.rows != this.cols ) {
            console.log("Matrix must be square");
            return;
        }
        if ( !(b instanceof Matrix) ) {
            console.log("invalid input");
            return;
        }
        if (this.cols != b.rows) {
            console.log("sizes are inconsistent");
            return;
        }
        let res = this.LU();
        let L = res.L;
        let U = res.U;
        let result = new Matrix(b.rows,b.cols);
        if (b.cols == 1) {
            result.vals = solveLinear(L,U,b.vals);
            return result;
        }
        for (let i = 0 ; i < b.cols ; i++) {
            let vec = b.vals.filter((_,index)=>index%b.cols==i);
            solveLinear(L,U,vec).forEach((val,index)=>{
                result.vals[index*b.cols + i] = val;
            });
        }
        return result;
    }
}

function solveLinear(L,U,vec) {
    let size = vec.length;
    let y = [];
    let x = [];
    for (let k = 0 ; k < size ; k++ ) {
        y[k] = vec[k];
        for (let i = 0 ; i < k ; i++ ) {
            y[k] -= L.vals[k*size + i]*y[i];
        }
    }
    for (let k = size - 1 ; k >= 0 ; k--) {
        for (let i = k+1 ; i < size ; i++) {
            y[k] -= x[i]*U.vals[k*size + i];
        }
        x[k] = y[k] / U.vals[k*size + k];
    }
    return x;
}