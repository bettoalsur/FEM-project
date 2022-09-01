class Matrix {
    constructor(rows,cols) {
        this.rows = rows;
        this.cols = cols;
        this.vals = [];
        for (let i = 0 ; i < this.rows*this.cols ; i++) {
            this.vals.push(0);
        }
    }

    disp() {
        console.log("[ ]:");
        for (let i = 0 ; i < this.rows ; i++) {
            let row = [];
            for (let j = 0 ; j < this.cols ; j++) {
                let index = i*this.cols + j;
                row.push(this.vals[index]);
            }
            console.log(row);
        }
    }

    set(row,col,val) {
        if (row >= 0 && row <= this.rows && col >= 0 && col <= this.cols) {
            let index = (row-1)*(this.cols) + (col-1);
            this.vals[index] = val;
        } else {
            console.log("position out of range");
        }
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
        let size = this.rows*this.cols;
        for (let i = 0 ; i < size ; i++) {
            C.vals[i] = this.vals[i] + B.vals[i];
        }
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
        if (b.cols != 1) {
            console.log("input is not a vector");
            return;
        } // must be removed in future
        let res = this.LU();
        let L = res.L;
        let U = res.U;
        let size = b.rows;
        let y = [];
        let x = [];
        for (let k = 0 ; k < size ; k++ ) {
            y[k] = b.vals[k];
            for (let i = 0 ; i < k ; i++ ) {
                y[k] -= L.vals[k*size + i]*y[i];
            }
        }
        for (let k = size - 1 ; k>= 0 ; k--) {
            for (let i = k+1 ; i < size ; i++) {
                y[k] -= x[i]*U.vals[k*size + i];
            }
            x[k] = y[k] / U.vals[k*size + k];
        }
        let C = new Matrix(size,1);
        C.vals = x;
        return C;
    }
}