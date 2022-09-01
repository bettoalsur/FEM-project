
let A = new Matrix(4,4);
A.vals = [-1,-1,6,9,-5,5,-3,6,7,-3,5,-6,3,-3,-3,3];
A.disp();

let b = new Matrix(4,1);
b.vals = [-29,-54,38,41];
b.disp()

let C = A.solve(b);
C.disp();

A.mult(C).disp();