{
  "More Information" -> {
    "A SpinSpinWeightedSpheroidalHarmonicSFunction may be treated as a numeric function of two arguments, \[Theta] and \[Phi].",
    "The parameters s, l, m and \[Gamma] and pre-computed expansion coefficients are stored inside the SpinSpinWeightedSpheroidalHarmonicSFunction.",
    "Numerical evaluation for specific \[Theta] and \[Phi] is done by summing the pre-computed expansion coefficients times the basis functions (the particular basis depends on the Method)."
    }
  "Numerical Evaluation" -> {
    "Slm = SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32]",
    "Slm[\[Pi]/4, 0]"
    },
  "Option Descriptions" -> {
    "Method" -> "The method to use to numerically evaluate the harmonic. Possible choices are \"Leaver\" and \"SphericalExpansion\"."
   },
  "See Also" -> {"SpinWeightedSpheroidalHarmonicS"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}