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
