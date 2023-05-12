{
  "More Information" -> {
    "A SpinSpinWeightedSpheroidalHarmonicSFunction is created when using SpinSpinWeightedSpheroidalHarmonicS without specifying values for \[Theta] and \[Phi].",
    "A SpinSpinWeightedSpheroidalHarmonicSFunction may be treated as a numeric function of two arguments, \[Theta] and \[Phi].",
    "The parameters s, l, m and \[Gamma], the method used, and and pre-computed expansion coefficients are stored inside the SpinSpinWeightedSpheroidalHarmonicSFunction.",
    "Numerical evaluation for specific \[Theta] and \[Phi] is done by summing the pre-computed expansion coefficients times the basis functions (the particular basis depends on the Method)."
    },
  "Numerical Evaluation" -> {
    "Slm = SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32]",
    "Slm[\[Pi]/4, 0]",
    "Slm[\"ExpansionCoefficients\"]"
    },
  "See Also" -> {"SpinWeightedSpheroidalHarmonicS"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}