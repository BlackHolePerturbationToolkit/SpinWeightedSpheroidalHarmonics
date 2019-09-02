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

(****************************************************************)
(* SpinWeightedSpheroidalHarmonicS s=phi=0                      *)
(* Test that for s=phi=0, function reduces to                   *)
(* Sqrt[(2 l + 1)/(4 \[Pi])] Sqrt[(l - m)!/(l + m)!]            *)
(*   SpheroidalPS[l, m, I \[Gamma], Cos[\[Theta]]]              *)
(****************************************************************)
Module[
    {s = 0, l = 1, m = 0, \[Gamma] = 0.5, \[Theta] = 0.5, \[Phi] = 0},
    VerificationTest[
        N[Sqrt[(2 l + 1)/(4 \[Pi])] Sqrt[(l - m)!/(l + m)!] SpheroidalPS[l, m, I \[Gamma], Cos[\[Theta]]]],
        N[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]]],
        TestID -> "Compare SpinWeightedSpheroidalEigenvalue with s=phi=0 to equivalent formula in documentation with l=1"
    ]
]

Module[
    {s = 0, l = 2, m = -1, \[Gamma] = 0.8, \[Theta] = 0.2, \[Phi] = 0},
    VerificationTest[
        N[Sqrt[(2 l + 1)/(4 \[Pi])] Sqrt[(l - m)!/(l + m)!] SpheroidalPS[l, m, I \[Gamma], Cos[\[Theta]]]],
        N[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]]],
        TestID -> "Compare SpinWeightedSpheroidalEigenvalue with s=phi=0 to equivalent formula in documentation with l=2"
    ]
]

(*complex spheroidicity tests near expected QNM values for s=0*)
(*based on QNMs as computed with qnm package made available by Leo Stein*)

Module[{testc=SetPrecision[{-0.0287516-0.161094 I,-0.149234-0.148211 I,-0.426033-0.115162 I,-0.0246915-0.993013 I},30],l=2,m=2},
VerificationTest[N[SpheroidalEigenvalue[l, m, testc], 10], 
 N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*testc] + 2*m*I*testc, 10], 
 TestID -> 
  "SpinWeightedSpheroidalEigenvalue with s=0 against (Mathematica) SpheroidalEigenvalue for complex (approximate QNM) spheroidicity (l=2,m=2)"]]

Module[{testc=SetPrecision[{-0.028664-0.849146 I,-0.0860387-0.848094 I,-0.426033-0.115162 I,-1.15863-1.82117 I},30],l=14,m=1},
VerificationTest[N[SpheroidalEigenvalue[l, m, testc], 10], 
 N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*testc] + 2*m*I*testc, 10], 
 TestID -> 
  "SpinWeightedSpheroidalEigenvalue with s=0 against (Mathematica) SpheroidalEigenvalue for complex (approximate QNM) spheroidicity (l=14,m=1)"]]

Module[{testc=SetPrecision[{-0.0286038-1.1459 I,-0.0860387-0.848094 I,-0.637869-3.03133 I,-1.10607-3.0076 I},30],l=17,m=17},
VerificationTest[N[SpheroidalEigenvalue[l, m, testc], 10], 
 N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*testc] + 2*m*I*testc, 10], 
 TestID -> 
  "SpinWeightedSpheroidalEigenvalue with s=0 against (Mathematica) SpheroidalEigenvalue for complex (approximate QNM) spheroidicity (l=17,m=17)"]]


