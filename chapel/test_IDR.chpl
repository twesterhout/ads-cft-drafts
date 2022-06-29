use IDR;
use Random;
use LinearAlgebra;

config const n : int = 100;

proc testDense(type eltType) {
  var A : [0 ..# n, 0 ..# n] eltType = noinit;
  fillRandom(A);
  forall i in 0 ..# n do
    A[i, i] += n:eltType;

  var b : [0 ..# n] eltType = noinit;
  fillRandom(b);

  idrsMethod(A, b, relTol=sqrt(epsilon(eltType)));
}



proc main() {

  var A : [0 ..# 4, 0 ..# 4] real;
  A[0, ..] = [0.95261881, 0.03465117, 0.87306725, 0.61793314];
  A[1, ..] = [0.2121299 , 0.12656886, 0.94639018, 0.78809681];
  A[2, ..] = [0.82831116, 0.66228212, 0.48960417, 0.4083788];
  A[3, ..] = [0.08738027, 0.63494534, 0.72405378, 0.50726822];
  var C : [0 ..# 4] real = [0.62064746, 0.62386813, 0.73991291, 0.35641101];
  var X : [0 ..# 4] real;

  idrsMethod(X, A, C, 4, 1e-8, 1e-6, 8);
  writeln("X = ", X);

  testDense(real(32));
  testDense(real(64));
  testDense(complex(64));
  testDense(complex(128));
}
