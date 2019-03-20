(****************************************************************)
(* SpinWeightedSpheroidalEigenvalue                             *)
(****************************************************************)
VerificationTest[
    SpinWeightedSpheroidalEigenvalue[2, 3, 2, 0]
    ,
    6
    ,
    TestID->"SpinWeightedSpheroidalEigenvalue with zero spheroidicity."
]

VerificationTest[
        SpinWeightedSpheroidalEigenvalue[2, 10, 10, 1.2345],
        79.76477183487565,
        TestID -> "SpinWeightedSpheroidalEigenvalue with positive spheroidicity"
]

Module[{l=2,m=2,c=SetPrecision[8.58583+1.98352 I,50]}, VerificationTest[
        N[SpheroidalEigenvalue[l, m, c],10],
        N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*c] + 2*m*I*c,10],
        TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 against built-in Mathematica function"
]]

Module[{l=2,m=2,c=SetPrecision[2.28831+5.351950 I,50]}, VerificationTest[
        N[SpheroidalEigenvalue[l, m, c],10],
        N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*c] + 2*m*I*c,10],
        TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 against built-in Mathematica function"
]]

