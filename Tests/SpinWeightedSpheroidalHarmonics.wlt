BeginTestSection["SpinWeightedSpheroidalEigenvalue"]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 3, 2, 0]
  ,
  6
  ,
  TestID->"SpinWeightedSpheroidalEigenvalue with 0 spheroidicity."
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 3, 2, 0.0]
  ,
  6
  ,
  TestID->"SpinWeightedSpheroidalEigenvalue with 0.0 spheroidicity."
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 3, 2, 0.0`32]
  ,
  6
  ,
  TestID->"SpinWeightedSpheroidalEigenvalue with 0.0``32 spheroidicity."
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

Get[FileNameJoin[{$testDir, "EigenvalueLeaverData.m"}]];

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

With[{s = -2, l = 2, m = 1, \[Gamma] = 0.1 I}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] + 2 m \[Gamma] - \[Gamma]^2,
    4.005766385035429 - 0.13328628274049953*I,
    SameTest -> withinRoundoff,
    TestID -> idStringEigenvalue[s, l, m, \[Gamma]]
  ]
]

With[{s = 0, l = 2, m = 1, \[Gamma] = 0.1 + 0.2 I}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] + 2 m \[Gamma] - \[Gamma]^2,
    6.012859862254978 - 0.01713352832954182 I,
    SameTest -> withinRoundoff,
    TestID -> idStringEigenvalue[s, l, m, \[Gamma]]
  ]
]

(* Suboptions *)
VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32],
  -0.66407962199376178419950018744266,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[...]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> "Leaver"],
  -0.66407962199376178419950018744266,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> \"Leaver\"]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"Leaver"}],
  -0.66407962199376178419950018744266,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"Leaver\"}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"Leaver", "InitialGuess" -> 0.1}],
  -0.66407962199376178419950018744266,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"Leaver\", <<valid suboption>>}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"Leaver", "NumInversions" -> 3}],
  -0.66407962199376178419950018744266,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"Leaver\", \"NumInversions\" -> 3}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"Leaver", "NumTerms" -> 1}],
  -0.66407962199376178419950018744266,
  {SpinWeightedSpheroidalEigenvalue::optx},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"Leaver\", <<invalid suboption>>}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> "SphericalExpansion"],
  -0.664079621993761784199500178632,
  {SpinWeightedSpheroidalEigenvalue::numterms},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> \"SphericalExpansion\"]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"SphericalExpansion"}],
  -0.664079621993761784199500178632,
  {SpinWeightedSpheroidalEigenvalue::numterms},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"SphericalExpansion\"}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"SphericalExpansion", "InitialGuess" -> 0.1}],
  -0.664079621993761784199500178632,
  {SpinWeightedSpheroidalEigenvalue::optx, SpinWeightedSpheroidalEigenvalue::numterms},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., Method -> {\"SphericalExpansion\", <<invalid subobtion>>}]"
]

VerificationTest[
  SpinWeightedSpheroidalEigenvalue[2, 2, 2, 0.1`32, Method -> {"SphericalExpansion", "NumTerms" -> 4}],
  -0.664079621993761781254003787957,
  {},
  TestID->"SpinWeightedSpheroidalEigenvalue[..., {\"SphericalExpansion\", <<valid suboption>>}]"
]

EndTestSection[]

BeginTestSection["SpinWeightedSpheroidalHarmonicS"]

Get[FileNameJoin[{$testDir, "SlmLeaverData.m"}]];

idString[s_, l_, m_, \[Gamma]_, method_:""] :=
  "SpinWeightedSpheroidalHarmonicS["<>ToString[s]<>", "<>ToString[l]<>", "<>ToString[m]<>", "<>ToString[\[Gamma], InputForm]<>"]";

idStringThetaDeriv[s_, l_, m_, \[Gamma]_, method_:""] :=
  "D[SpinWeightedSpheroidalHarmonicS["<>ToString[s]<>", "<>ToString[l]<>", "<>ToString[m]<>", "<>ToString[\[Gamma], InputForm]<>method<>"], \[Theta]]";

idStringPhiDeriv[s_, l_, m_, \[Gamma]_, method_:""] :=
  "D[SpinWeightedSpheroidalHarmonicS["<>ToString[s]<>", "<>ToString[l]<>", "<>ToString[m]<>", "<>ToString[\[Gamma], InputForm]<>method<>"], \[Phi]]";

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

With[{s = 0, l = 1, m = 0, \[Gamma] = 0.5, \[Theta] = 0.5, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]][\[Theta], \[Phi]],
    0.4305973913818948,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = 0, l = 2, m = -1, \[Gamma] = 0.8, \[Theta] = 0.2, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]][\[Theta], \[Phi]],
    0.15410463040391498,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = -2, l = 10, m = 5, \[Gamma] = 1.0, \[Theta] = 0.1, \[Phi] = 4}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    -0.010529282862103993 - 0.023555700389793967*I,
    SameTest -> withinRoundoff,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = 1, l = 2, m = 2, \[Gamma] = 2.12`20, \[Theta] = \[Pi]/3., \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    0.053503866650272916,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = -2, l = 7, m = 2, \[Gamma] = 0.01`40, \[Theta] = \[Pi]/2, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    -0.17001322692047803109041563873499247769`26.19971750849647,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = 0, l = 1, m = 0, \[Gamma] = 0.01`30, \[Theta] = 0, \[Phi] = 0}, 
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    0.48860446631177094306032408791530492821`27.827288196470665,
    TestID -> idString[s, l, m, \[Gamma]]
  ]
]

With[{s = 0, l = 2, m = 2, \[Gamma] = 0.00000001, \[Theta] = 0.3, \[Phi] = 0.2},
  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]]][\[Theta], \[Phi]],
    0.20088968520084047 + 0.0849347966173594*I,
    SameTest -> withinRoundoff,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma]]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> Automatic]][\[Theta], \[Phi]],
    0.20088968520084047 + 0.0849347966173594*I,
    SameTest -> withinRoundoff,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> Automatic"]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"]][\[Theta], \[Phi]],
    0.20088968520084047 + 0.0849347966173594*I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    SameTest -> withinRoundoff,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> \"SphericalExpansion\""]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"]][\[Theta], \[Phi]],
    0.20088968520084058 + 0.08493479661735945*I,
    SameTest -> withinRoundoff,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> \"Leaver\""]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]]][\[Theta], \[Phi]],
    -0.026273411446970076 + 0.062142461805285834*I,
    SameTest -> withinRoundoff,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma]]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> Automatic]][\[Theta], \[Phi]],
    -0.026273411446970076 + 0.062142461805285834*I,
    SameTest -> withinRoundoff,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> Automatic"]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"]][\[Theta], \[Phi]],
    -0.026273411446970076 + 0.062142461805285834*I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    SameTest -> withinRoundoff,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> \"SphericalExpansion\""]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"]][\[Theta], \[Phi]],
    -0.026273411446970086 + 0.062142461805285855*I,
    SameTest -> withinRoundoff,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> \"Leaver\""]
  ]
]


With[{s = 0, l = 2, m = 2, \[Gamma] = 0.00000001`32, \[Theta] = 0.3`32, \[Phi] = 0.2`32},
  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]]][\[Theta], \[Phi]],
    0.2008896852008405311614794231142 + 0.08493479661735942846637968004680 I,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma]]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> Automatic]][\[Theta], \[Phi]],
    0.2008896852008405311614794231142 + 0.08493479661735942846637968004680 I,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> Automatic"]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"]][\[Theta], \[Phi]],
    0.2008896852008405311614794231142 + 0.08493479661735942846637968004680 I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> \"SphericalExpansion\""]
  ]

  VerificationTest[
    Derivative[1, 0][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"]][\[Theta], \[Phi]],
    0.2008896852008405311614794231142 + 0.08493479661735942846637968004680 I,
    TestID -> idStringThetaDeriv[s, l, m, \[Gamma], ", Method -> \"Leaver\""]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]]][\[Theta], \[Phi]],
    -0.02627341144697007921403373259026 + 0.06214246180528584097100501520202 I,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> Automatic"]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> Automatic]][\[Theta], \[Phi]],
    -0.02627341144697007921403373259026 + 0.06214246180528584097100501520202 I,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma]]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"]][\[Theta], \[Phi]],
    -0.02627341144697007921403373259026 + 0.06214246180528584097100501520202 I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> \"SphericalExpansion\""]
  ]

  VerificationTest[
    Derivative[0, 1][SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"]][\[Theta], \[Phi]],
    -0.02627341144697007921403373259026 + 0.06214246180528584097100501520202 I,
    TestID -> idStringPhiDeriv[s, l, m, \[Gamma], ", Method -> \"Leaver\""]
  ]
]

(* Suboptions *)
With[{s = 2, l = 2, m = 2, \[Gamma] = 0.1`32, \[Theta] = 0.3`32, \[Phi] = 0.4`32},
  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]][\[Theta], \[Phi]],
    0.00019649578982284742659052547689 + 0.00020231964149966822640907627571 I,
    {},
    TestID->"SpinWeightedSpheroidalHarmonicS[...]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "Leaver"][\[Theta], \[Phi]],
    0.00019649578982284742659052547689 + 0.00020231964149966822640907627571 I,
    {},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> \"Leaver\"]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"Leaver"}][\[Theta], \[Phi]],
    0.00019649578982284742659052547689 + 0.00020231964149966822640907627571 I,
    {},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> {\"Leaver\"}]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"Leaver", "NumTerms" -> 2}][\[Theta], \[Phi]],
    0.0001841212861596120040037305869324 + 0.0001895783754036400622852984961101 I,
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> {\"Leaver\", NumTerms -> 2}]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"Leaver", "Invalid" -> 1}][\[Theta], \[Phi]],
    0.00019649578982284742659052547689 + 0.00020231964149966822640907627571 I,
    {SpinWeightedSpheroidalHarmonicS::optx},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> {\"Leaver\", <<invalid suboption>>}]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"][\[Theta], \[Phi]],
    0.0001964957898238118597970926718157 + 0.0002023196415006612440242575779433 I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> \"SphericalExpansion\"]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"SphericalExpansion"}][\[Theta], \[Phi]],
    0.0001964957898238118597970926718157 + 0.0002023196415006612440242575779433 I,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> {\"SphericalExpansion\"}]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"SphericalExpansion", "InitialGuess" -> 0.1}][\[Theta], \[Phi]],
    0.0001964957898238118597970926718157 + 0.0002023196415006612440242575779433 I,
    {SpinWeightedSpheroidalHarmonicS::optx, SpinWeightedSpheroidalHarmonicS::numterms},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., Method -> {\"SphericalExpansion\", <<invalid subobtion>>}]"
  ]

  VerificationTest[
    SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> {"SphericalExpansion", "NumTerms" -> 4}][\[Theta], \[Phi]],
    0.0001964957988921799016469279216708 + 0.0002023196508378026294361569589402 I,
    {},
    TestID->"SpinWeightedSpheroidalHarmonicS[..., {\"SphericalExpansion\", <<valid suboption>>}]"
  ]

]

(* Zero spheroidicity *)
VerificationTest[
  SpinWeightedSpheroidalHarmonicS[2, 2, 2, 0.0][\[Theta], \[Phi]],
  1/2 E^(2 I \[Phi]) Sqrt[5/\[Pi]] Sin[\[Theta]/2]^4,
  TestID->"SpinWeightedSpheroidalHarmonicS[...] with 0.0 spheroidicity"
]

VerificationTest[
  SpinWeightedSpheroidalHarmonicS[2, 2, 2, 0][\[Theta], \[Phi]],
  1/2 E^(2 I \[Phi]) Sqrt[5/\[Pi]] Sin[\[Theta]/2]^4,
  TestID->"SpinWeightedSpheroidalHarmonicS[...] with 0 spheroidicity"
]

VerificationTest[
  SpinWeightedSpheroidalHarmonicS[2, 2, 2, 0.0``32][\[Theta], \[Phi]],
  1/2 E^(2 I \[Phi]) Sqrt[5/\[Pi]] Sin[\[Theta]/2]^4,
  TestID->"SpinWeightedSpheroidalHarmonicS[...] with 0.0``32 spheroidicity"
]

With[{s = -2, l = 5, m = 2, \[Gamma] = 0.1}, 
  VerificationTest[
    SetAccuracy[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]]["ExpansionCoefficients"][[1;;5]], 6],
    <|0 -> -8.818050086267075, 1 -> 21.65089643682929, 2 -> -17.171313724018557, 
      3 -> 4.847739931702922, 4 -> -0.2937583721632349|>,
    TestID -> "ExpansionCoefficients (Leaver)"
  ]

  VerificationTest[
    SetAccuracy[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], Method -> "SphericalExpansion"]["ExpansionCoefficients"][[2;;6]], 6],
    <|-2 -> 0.000055654127758144056, -1 -> -0.0169795873354771, 0 -> 0.9997443511214629, 
      1 -> 0.014929457524949523, 2 -> 0.00018496617816892035|>,
    {SpinWeightedSpheroidalHarmonicS::numterms},
    TestID -> "ExpansionCoefficients (SphericalExpansion)"
  ]
]

EndTestSection[]

BeginTestSection["Complex spheroidicity near expected QNM values"]

(* Based on QNMs as computed with qnm package made available by Leo Stein *)

With[{testc = {-0.0287516-0.161094 I,-0.149234-0.148211 I,-0.426033-0.115162 I,-0.0246915-0.993013 I}, l=2, m=2},
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc,
    {6.02153433010513 - 0.007939169387951778*I, 5.99974298637217 - 0.03791683322948014*I, 
     5.855756412693891 - 0.0841715869120514*I, 6.842766741987868 - 0.041839426991681715*I},
    SameTest -> withinRoundoff,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=2,m=2)"
  ]
]

With[{testc = {-0.028664-0.849146 I,-0.0860387-0.848094 I,-0.426033-0.115162 I,-1.15863-1.82117 I}, l=14, m=1}, 
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc,
    {210.36147982153605 - 0.024437405606186034*I,  210.35727704692184 - 0.07326079871853541*I, 209.91558030447734 - 0.049233841074118545*I, 210.9886233104435 - 2.1200593243717076*I},
    SameTest -> withinRoundoff,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=14,m=1)"
  ]
]

With[{testc = {-0.0286038 - 1.1459*I, -0.0860387 - 0.848094*I, -0.637869 - 3.03133*I, -1.10607 - 3.0076*I}, l = 17, m = 17},
  VerificationTest[
    SpinWeightedSpheroidalEigenvalue[0, l, m, I testc] + 2 m I testc,
    {307.2767705582399 - 0.06377931094717815*I, 306.6926124407844 - 0.14198975395701074*I, 314.5435916995019 - 3.7614072964319583*I, 313.61055899767786 - 6.471497169350965*I},
    SameTest -> withinRoundoff,
    TestID -> "SpinWeightedSpheroidalEigenvalue with s=0 for complex (approximate QNM) spheroidicity (l=17,m=17)"
  ]
]

EndTestSection[]
