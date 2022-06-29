module IDR {

  private use LinearAlgebra;
  private use Random;

  config param idrVerbose : bool = true;

  proc computeOmega(t, s) {
      const angle = sqrt(2.0) / 2;
      const ns = norm(s);
      const nt = norm(t);
      const ts = inner(t, s);
      const rho = abs(ts / (nt * ns));
      var om = ts / (nt * nt);
      if rho < angle then
          om = om * angle:om.type / rho;
      return om;
  }

  /* Machine epsilon for real(64) */
  proc epsilon(type t: real(64)) : real {
    extern const DBL_EPSILON: real;
    return DBL_EPSILON;
  }

  /* Machine epsilon for real(32) */
  proc epsilon(type t: real(32)) : real {
    extern const FLT_EPSILON: real;
    return FLT_EPSILON;
  }

  proc epsilon(type t: complex(64)) : real {
    return epsilon(real(32));
  }

  proc epsilon(type t: complex(128)) : real {
    return epsilon(real(64));
  }

  proc idrsMethod(A, const ref C : [] ?eltType,
                  s : int = 8,
                  absTol : real = 0,
                  relTol : real = sqrt(epsilon(eltType)),
                  maxIter : int = C.size) {
    var X : [C.domain] eltType = 0;
    idrsMethod(X, A, C, s, absTol, relTol, maxIter);
    return X;
  }
  proc idrsMethod(ref X : [] ?eltType, A, const ref C : [] eltType,
                  s : int = 8,
                  absTol : real = 0,
                  relTol : real = sqrt(epsilon(eltType)),
                  maxIter : int = C.size) {
    if idrVerbose {
      writeln("=== IDR(", s, ") ===");
      writeln("iter\tresnorm");
    }

    var R = C - dot(A, X);
    var normR = norm(R);
    var iteration = 1;
    const tol = max(relTol * normR, absTol);

    // writeln(R);
    // writeln(normR);
    // writeln(tol);

    // Initial guess is a good enough solution
    if normR <= tol then return;

    const dim = C.size;
    type eltType = C.eltType;

    var P : [0 ..# s, 0 ..# dim] eltType = noinit;
    fillRandom(P);
    // P[0, ..] = [0.5331830160438613, 0.4540291355871424, 0.017686826714964354, 0.17293302893695128];
    // P[1, ..] = [0.9589258763297348, 0.9735659798036858, 0.30386961108678245, 0.17690859963285543];
    // writeln("P = ", P);

    var U : [0 ..# s, 0 ..# dim] eltType;
    var G : [0 ..# s, 0 ..# dim] eltType;
    var Q : [0 ..# dim] eltType;
    var V : [0 ..# dim] eltType;
    var M = eye(s, s, eltType=eltType);
    var f : [0 ..# s] eltType;

    var omega : eltType = 1;
    while normR > tol && iteration <= maxIter {
      forall i in 0 ..# s do
        f[i] = inner(P[i, ..], R);
      // writeln("f = ", f);

      for k in 0 ..# s {
        // Solve small system and make v orthogonal to P
        const c = solve_tril(M[k .., k ..].reindex(0 ..# s - k, 0 ..# s - k),
                             f[k ..].reindex(0 ..# s - k),
                             unit_diag=false);
        // writeln("c = ", c);

        // Constructing linear combinations
        V = c[0] * G[k, ..];
        Q = c[0] * U[k, ..];
        for i in k + 1 ..< s {
            V += c[i - k] * G[i, ..];
            Q += c[i - k] * U[i, ..];
        }
        // writeln("V = ", V);
        // writeln("Q = ", Q);

        // Compute new U[:,k] and G[:,k], G[:,k] is in space G_j
        V = R - V;

        // Preconditioning
        // ldiv!(Pl, V)

        U[k, ..] = Q + omega * V;
        G[k, ..] = dot(A, U[k, ..]); // mul!(G[k], A, U[k])
        // writeln("U = ", U);
        // writeln("G = ", G);

        // Bi-orthogonalise the new basis vectors
        for i in 0 ..< k {
          const alpha = inner(P[i, ..], G[k, ..]) / M[i, i];
          // writeln("alpha = ", alpha);
          G[k, ..] -= alpha * G[i, ..];
          U[k, ..] -= alpha * U[i, ..];
        }
        // writeln("U = ", U);
        // writeln("G = ", G);


        // New column of M = P'*G  (first k-1 entries are zero)
        forall i in k ..< s do
          M[i, k] = inner(P[i, ..], G[k, ..]);

        // Make r orthogonal to q_i, i = 1..k
        const beta = f[k] / M[k, k];
        R -= beta * G[k, ..];
        X += beta * U[k, ..];

        normR = norm(R);
        if idrVerbose then writeln(iteration, "\t", normR);
        if normR < tol || iteration == maxIter then return;

        if k < s then
          f[k + 1 ..< s] -= beta * M[k + 1 ..< s, k];

        // writeln("f = ", f);
        // writeln("M = ", M);
        iteration += 1;
      }

      // Now we have sufficient vectors in G_j to compute residual in G_j+1
      // Note: r is already perpendicular to P so v = r
      V = R;

      // Preconditioning
      // ldiv!(Pl, V)

      Q = dot(A, V);
      omega = computeOmega(Q, R);
      R -= omega * Q;
      X += omega * V;

      normR = norm(R);
      iteration += 1;
    }

    // if idrVerbose then
    //   writeln();
  }

} // end module IDR
