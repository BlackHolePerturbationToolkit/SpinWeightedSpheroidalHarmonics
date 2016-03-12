{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalHarmonicS computes the oblate spin-weighted spheroidal harmonic with oblateness \[Gamma].",
    "For \[Gamma]=0, SpinWeightedSpheroidalHarmonicS reduces to SpinWeightedSphericalHarmonicY.",
    "For s=0 and \[Gamma]=0, SpinWeightedSpheroidalHarmonicS is equivalent to SphericalHarmonicY.",
    "For s=0, SpinWeightedSpheroidalHarmonicS[0, l, m, \[Gamma], \[Theta], \[Phi]] is equivalent to Sqrt[(2 l + 1)/(4 \[Pi])] Sqrt[(l - m)!/(l + m)!] SpheroidalPS[l, m, I \[Gamma], Cos[\[Theta]]].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalHarmonicS produces a numerical value of the same precision.",
    "Numerical values are computed using a series expansion about \[Gamma] = 0, with the number of terms in the expansion determined automatically to ensure the desired precision is reached.",
    "The MaxIterations option controls the maximum number of terms to include in a numerical approximation. For large \[Gamma]/l it may be necessary to set this to a large value in order to achieve the desired accuracy.",
    "For sufficiently large \[Gamma]/l the numerical evaluation of SpinWeightedSpheroidalHarmonicS may fail to converge for any value of MaxIterations.",
    "Series expansions about \[Gamma] = 0 can be computed to arbitrary order for generic s, l and m."
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