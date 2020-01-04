BeginTestSection["SpinWeightedSpheroidalEigenvalue"]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 3, 2, 0]
  ,
  6
  ,
  TestID->"SpinWeightedSpheroidalEigenvalue with zero spheroidicity."
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 10, 10, 1.2345]
  ,
  79.76477183487565
  ,
  TestID->"SpinWeightedSpheroidalEigenvalue with positive machine-precision spheroidicity"
]

idStringEigenvalue[s_, l_, m_, \[Gamma]_, method_:""] :=
  "SpinWeightedSpheroidalEigenvalue["<>ToString[s]<>", "<>ToString[l]<>", "<>ToString[m]<>", "<>ToString[\[Gamma], InputForm]<>"]";

Get[PacletResource["SpinWeightedSpheroidalHarmonics", "EigenvalueLeaverData.m"]];

Table[
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]]
    ,
    EigenvalueLeaverData[s, l, m, \[Gamma]]
    ,
    TestID -> idStringEigenvalue[s, l, m, \[Gamma], "Leaver"]
  ],
  {\[Gamma], {0.000000001`32}}, {s, -3, 3}, {l, Abs[s], 4}, {m, -l, l}
]

VerificationTest[
    SpinWeightedSpheroidalEigenvalue[2, 2, 1, 0.000000001],
    -3.333333442653828*^-9,
    SameTest -> withinRoundoff,
    TestID -> "SpinWeightedSpheroidalEigenvalue with machine precision near-zero spheroidicity."
    
]

Module[{s = -2, l = 2, m = 1, \[Gamma] = 0.1 I}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] + 2 m \[Gamma] - \[Gamma]^2,
    4.005766385035429 - 0.13328628274049953*I,
    SameTest -> withinRoundoff,
    TestID -> idStringEigenvalue[s, l, m, \[Gamma]]
  ]
]

Module[{s = 0, l = 2, m = 1, \[Gamma] = 0.1 + 0.2 I}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] + 2 m \[Gamma] - \[Gamma]^2,
    6.012859862254978 - 0.01713352832954182 I,
    SameTest -> withinRoundoff,
    TestID -> idStringEigenvalue[s, l, m, \[Gamma]]
  ]
]

EndTestSection[]

BeginTestSection["SpinWeightedSpheroidalHarmonicS"]

Get[PacletResource["SpinWeightedSpheroidalHarmonics", "SlmLeaverData.m"]];

idString[s_, l_, m_, \[Gamma]_, method_:""] :=
  "SpinWeightedSpheroidalHarmonicS["<>ToString[s]<>", "<>ToString[l]<>", "<>ToString[m]<>", "<>ToString[\[Gamma], InputForm]<>"]";

Table[
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma],  Method -> "Leaver"][0.3`20, 0.4`20]
    ,
    SlmLeaverData[s, l, m, \[Gamma]]
    ,
    TestID -> idString[s, l, m, \[Gamma], "Leaver"]
  ],
  {\[Gamma], {0.000001`20}}, {s, -3, 3}, {l, Abs[s], 4}, {m, -l, l}
]

Module[{s = 0, l = 1, m = 0, \[Gamma] = 0.5, \[Theta] = 0.5, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]]
    ,
    0.4305973913818948
    ,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

Module[{s = 0, l = 2, m = -1, \[Gamma] = 0.8, \[Theta] = 0.2, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]]
    ,
    0.15410463040391498
    ,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

Module[{s = -2, l = 10, m = 5, \[Gamma] = 1.0, \[Theta] = 0.1, \[Phi] = 4}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi], Method -> "Leaver"],
    -0.010529405169055157 - 0.02355597401012832*I,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

Module[{s = 1, l = 2, m = 2, \[Gamma] = 2.12`20, \[Theta] = \[Pi]/3., \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    0.053503866650272916,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

Module[{s = -2, l = 7, m = 2, \[Gamma] = 0.01`40, \[Theta] = \[Pi]/2, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi], Method -> "Leaver"],
    -0.17001322692047803109041563873499247769`26.19971750849647,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]


EndTestSection[]

BeginTestSection["Complex spheroidicity near expected QNM values"]

(* Based on QNMs as computed with qnm package made available by Leo Stein *)

Module[{testc = {-0.0287516-0.161094 I,-0.149234-0.148211 I,-0.426033-0.115162 I,-0.0246915-0.993013 I}, l=2, m=2},
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc
    ,
    {6.02153433010513 - 0.007939169387951778*I, 5.99974298637217 - 0.03791683322948014*I, 
     5.855756412693891 - 0.0841715869120514*I, 6.842766741987868 - 0.041839426991681715*I}
    ,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=2,m=2)"
  ]
]

Module[{testc = {-0.028664-0.849146 I,-0.0860387-0.848094 I,-0.426033-0.115162 I,-1.15863-1.82117 I}, l=14, m=1}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc
    ,
    {210.36147982153605 - 0.024437405606186034*I,  210.35727704692184 - 0.07326079871853541*I, 209.91558030447734 - 0.049233841074118545*I, 210.9886233104435 - 2.1200593243717076*I}
    ,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=14,m=1)"
  ]
]

Module[{testc = {-0.0286038 - 1.1459*I, -0.0860387 - 0.848094*I, -0.637869 - 3.03133*I, -1.10607 - 3.0076*I}, l = 17, m = 17},
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc
    ,
    {307.2767705582399 - 0.06377931094717815*I, 306.6926124407844 - 0.14198975395701074*I, 314.5435916995019 - 3.7614072964319583*I, 313.61055899767786 - 6.471497169350965*I}
    ,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=17,m=17)"
  ]
]

EndTestSection[]
