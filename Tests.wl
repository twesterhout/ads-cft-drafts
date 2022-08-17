(* ::Package:: *)

SetOptions[EvaluationNotebook[], CellContext -> Notebook];
Get["EDCRGTCcode.m"] (* Riemann Geometry Tensor Calculus *);
Needs["Utilities`"];
On[Assert];


xGrid = {0.4986655005698084, 0.9973310011396168, 1.4959965017094252, 1.9946620022792336, 2.4933275028490423, 2.9919930034188504};
yGrid = {0.0};
zGrid = {0.9045084971874737, 0.6545084971874737, 0.34549150281252633, 0.09549150281252633};
zGridFull = Catenate[{{1.0}, zGrid, {0.0}}];


Clear[combineGrids];
combineGrids[xGrid_, yGrid_, zGrid_] :=
    Module[{x, y, z},
        x = Table[xx, {zz, zGrid}, {yy, yGrid}, {xx, xGrid}] // Flatten;
        y = Table[yy, {zz, zGrid}, {yy, yGrid}, {xx, xGrid}] // Flatten;
        z = Table[zz, {zz, zGrid}, {yy, yGrid}, {xx, xGrid}] // Flatten;
        {x, y, z}
    ];


{bulkX, bulkY, bulkZ} = combineGrids[xGrid, yGrid, zGrid];
{conformalX, conformalY, conformalZ} = combineGrids[xGrid, yGrid, {0.0}];
{horizonX, horizonY, horizonZ} = combineGrids[xGrid, yGrid, {1.0}];
fullGridLength = Length[xGrid] * Length[yGrid] * Length[zGridFull];
bulkGridLength = Length[bulkX];
conformalGridLength = Length[conformalX];
horizonGridLength = Length[horizonX];


SeedRandom[42];
fullQs = Table[RandomReal[], {zz, zGridFull}, {yy, yGrid}, {xx, xGrid}, {p, Parameters}];
fullDQs = Table[RandomReal[], {zz, zGridFull}, {yy, yGrid}, {xx, xGrid}, {direction, 1, 4}, {p, Parameters}];
fullDDQs = Table[RandomReal[], {zz, zGridFull}, {yy, yGrid}, {xx, xGrid}, {direction1, 1, 4}, {direction2, 1, 4}, {p, Parameters}];
bulkQs = ArrayReshape[fullQs[[2 ;; -2]], {bulkGridLength, Length[Parameters]}];
bulkDQs = ArrayReshape[fullDQs[[2 ;; -2]], {bulkGridLength, 4, Length[Parameters]}];
bulkDDQs = ArrayReshape[fullDDQs[[2 ;; -2]], {bulkGridLength, 4, 4, Length[Parameters]}];
conformalQs = ArrayReshape[fullQs[[-1]],{conformalGridLength, Length[Parameters]}];
conformalDQs = ArrayReshape[fullDQs[[-1]],{conformalGridLength, 4, Length[Parameters]}];
conformalDDQs = ArrayReshape[fullDDQs[[-1]],{conformalGridLength, 4, 4, Length[Parameters]}];
horizonQs = ArrayReshape[fullQs[[1]],{horizonGridLength, Length[Parameters]}];
horizonDQs = ArrayReshape[fullDQs[[1]],{horizonGridLength, 4, Length[Parameters]}];
horizonDDQs = ArrayReshape[fullDDQs[[1]],{horizonGridLength, 4, 4, Length[Parameters]}];


ClearAll[parameterReplacement, parameterIndex, derivativeReplacement];

parameterReplacement[Q_] :=
    MapThread[#2[x, z_] -> Q[[#1]]&, {Range[Length[Parameters]], Parameters}];

parameterIndex[s_Symbol] :=
    FirstPosition[Parameters, s][[1]];

derivativeReplacement[DQ_, DDQ_][x_] :=
    x;
derivativeReplacement[DQ_, DDQ_][Derivative[1, 0][s_Symbol][x_Symbol, z_]] :=
    DQ[[2, parameterIndex[s]]];
derivativeReplacement[DQ_, DDQ_][Derivative[0, 1][s_Symbol][x_Symbol, z_]] :=
    DQ[[4, parameterIndex[s]]];
derivativeReplacement[DQ_, DDQ_][Derivative[2, 0][s_Symbol][x_Symbol, z_]] :=
    DDQ[[2, 2, parameterIndex[s]]];
derivativeReplacement[DQ_, DDQ_][Derivative[1, 1][s_Symbol][x_Symbol, z_]] :=
    DDQ[[2, 4, parameterIndex[s]]];
derivativeReplacement[DQ_, DDQ_][Derivative[0, 2][s_Symbol][x_Symbol, z_]] :=
    DDQ[[4, 4, parameterIndex[s]]];

ClearAll[repeatedRewrite, computeQuantity, evalOnBulkGrid, evalOnConformalGrid, evalOnHorizonGrid];
repeatedRewrite[expr_, f_] :=
    Map[Replace[#, x_ :> f[x]]&, expr /. {x_ :> f[x]}, Infinity];

computeQuantity[expr_, {xx_, yy_, zz_, Q_, DQ_, DDQ_}] :=
    Module[{e = expr},
        e = repeatedRewrite[e, derivativeReplacement[DQ, DDQ]];
        e = e //. parameterReplacement[Q];
        e = e /. {L -> 1, \[Mu] -> 2.3, k0 -> 2.1, \[Theta] -> 0.785, V0 -> 5.0, x -> xx, y -> yy, z -> zz};
        e
    ];

evalOnBulkGrid[expr_] :=
    N[computeQuantity[expr, #]& /@ Transpose[{bulkX, bulkY, bulkZ, bulkQs, bulkDQs, bulkDDQs}]];
evalOnConformalGrid[expr_] :=
    N[computeQuantity[expr, #]& /@ Transpose[{conformalX, conformalY, conformalZ, conformalQs, conformalDQs, conformalDDQs}]];
evalOnHorizonGrid[expr_] :=
    N[computeQuantity[expr, #]& /@ Transpose[{horizonX, horizonY, horizonZ, horizonQs, horizonDQs, horizonDDQs}]];


dCoord = {dt, dx, dy, dz};
P[z_] :=
    1 + z + z^2 - (\[Mu]^2 z^3) / 2;
ds = (L^2) / z^2 (-(1 - z) P[z] Qtt[x, z] dt^2 + (Qzz[x, z] dz^2) / (P[z] (1 - z)) + Qxx[x, z] (dx + z^2 Qxz[x, z] dz)^2 + Qyy[x, z] dy^2);
gdd = Map[Coefficient[ds, #]&, TensorProduct[dCoord, dCoord], {2}];
gdd = UpperTriangularize[gdd];
gdd = (gdd + Transpose[gdd]) / 2 // Simplify;

RGtensors[gdd, Coord, {1, 0, 0}];

Ad = {(1 - z) \[Psi][x, z], 0, 0, 0};
Fdd = (Transpose[#] - #&) @ covD[Ad] // Simplify;
FUd = Raise[Fdd, {1}];
FUU = Raise[FUd, {2}];
divFd = covDiv[FUd, {1}];
EOM\[Psi] = divFd[[1]] // Simplify;

\[Phi] = z * \[Phi]1[x, z];
\[Chi] = z * \[Chi]1[x, z];
V = -2 (\[Phi]^2 + \[Chi]^2) / L^2;
DV\[Phi] = -4 \[Phi] / L^2;
DV\[Chi] = -4 \[Chi] / L^2;
covD\[Phi] = covD[\[Phi]];
covD\[Chi] = covD[\[Chi]];
EOM\[Phi] = covDiv[Raise[covD\[Phi], {1}], {1}] - DV\[Phi] // Simplify;
EOM\[Chi] = covDiv[Raise[covD\[Chi], {1}], {1}] - DV\[Chi] // Simplify;

GUddRef = GUdd /. {Qtt -> (1&), Qzz -> (1&), Qxx -> (1&), Qyy -> (1&), Qxz -> (0&)} // Simplify;
\[Xi]U = multiDot[gUU, (GUdd - GUddRef), {1, 2}, {2, 3}];
\[Xi]d = Lower[\[Xi]U, 1];
div\[Xi]dd = (Transpose[#] + #) / 2& @ covD[\[Xi]d] // Simplify;
Gdd = Rdd + (3 / L^2 - V) * gdd -
    2 (Outer[Times, covD\[Phi], covD\[Phi]] + Outer[Times, covD\[Chi], covD\[Chi]]) -
    (multiDot[Fdd, Raise[Fdd, {2}], {2, 2}] -
        gdd / 4 * multiDot[Fdd, FUU, {1, 1}, {2, 2}]
    );
EOMQ = Gdd - div\[Xi]dd // Simplify;

bulkEquations = {
    EOMQ[[1, 1]],
    EOMQ[[2, 2]],
    EOMQ[[3, 3]],
    EOMQ[[4, 4]],
    EOMQ[[2, 4]],
    EOM\[Psi],
    EOM\[Phi],
    EOM\[Chi]
};
bulkEquations = evalOnBulkGrid[bulkEquations];
bulkEquations = ArrayReshape[bulkEquations, {Length[zGrid], Length[yGrid], Length[xGrid], Length[Parameters]}];

conformalEquations = constructConformalBoundaryConditions[];
conformalEquations = evalOnConformalGrid[conformalEquations];
conformalEquations = ArrayReshape[conformalEquations, {1, Length[yGrid], Length[xGrid], Length[Parameters]}];

horizonEquations = constructHorizonBoundaryConditions[minimalOrder[gdd]];
horizonEquations = evalOnHorizonGrid[horizonEquations];
horizonEquations = ArrayReshape[horizonEquations, {1, Length[yGrid], Length[xGrid], Length[Parameters]}];

combinedEquations =
    ArrayReshape[
        Catenate[{horizonEquations, bulkEquations, conformalEquations}],
        {Length[xGrid] * Length[yGrid] * (Length[zGrid] + 2), Length[Parameters]}
    ];


(*Export["test_data.h5",
    {"x" \[Rule] xGrid,
        "y" \[Rule] yGrid,
        "z" \[Rule] zGridFull,
        "Qs" \[Rule] ArrayReshape[fullQs, {fullGridLength, Length[Parameters]}],
        "DQs" \[Rule] ArrayReshape[fullDQs, {fullGridLength, 4, Length[Parameters]}],
        "DDQs" \[Rule] ArrayReshape[fullDDQs, {fullGridLength, 4, 4, Length[Parameters]}],
        "equations" \[Rule] combinedEquations
    },
    "Datasets"
];*)
