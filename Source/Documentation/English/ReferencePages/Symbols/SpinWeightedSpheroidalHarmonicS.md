{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalHarmonicS computes the oblate spin-weighted spheroidal harmonic with oblateness \[Gamma].",
    "For \[Gamma]=0, SpinWeightedSpheroidalHarmonicS reduces the SpinWeightedSphericalHarmonicY.",
    "For s=0 and \[Gamma]=0, SpinWeightedSpheroidalHarmonicS is equivalent to SphericalHarmonicY.",
    "For s=0, SpinWeightedSpheroidalHarmonicS[0, l, m, \[Gamma], \[Theta], \[Phi]] is equivalent to SpheroidalPS[l, m, I \[Gamma], \[Theta]] Exp[I m \[Phi]].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalHarmonicS produces a numerical value of the same precision.",
    "Series expansions in \[Gamma] can be computed to arbitrary order for generic s, l and m."
    },
  "Numerical Evaluation" -> {
    "SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1, \[Pi]/4, 0]",
    "SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32, \[Pi]/4, 0]",
    },
  "Series Expansion" -> {
    "Series[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]], {\[Gamma], 0, 2}]",
    "Series[SpinWeightedSpheroidalHarmonicS[-2, 2, 2, \[Gamma], \[Theta], \[Phi]], {\[Gamma], 0, 2}]",
    },
  "Option Descriptions" -> {
    "MaxIterations" -> "The maximum number of iterations to use when computing numerical values."
   },
  "See Also" -> {"SpinWeightedSpheroidalEigenvalue", "SpinWeightedSphericalHarmonicY"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}