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
        TestID -> "SpinWeightedSpheroidalEigenvalue with positive spheroidicity",
]

(*complex spheroidicity tests near expected QNM values for s=0*)

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

(*general complex spheroidicity comparison branch choice SpheroidalEigenvalue for s=0*)

Module[{l = 2, m = 2, c = 2.288310998002696550202017533592879772186279296875`49.74505997240401 + 5.3519509384954062625183723866939544677734375`50.11405704809964 I}, 
 VerificationTest[N[SpheroidalEigenvalue[l, m, c], 10], 
  N[SpinWeightedSpheroidalEigenvalue[0, l, m, I*c] + 2*m*I*c, 10], 
  TestID -> 
   "SpinWeightedSpheroidalEigenvalue with s=0 branch choice with respect to (Mathematica) SpheroidalEigenvalue"]]


