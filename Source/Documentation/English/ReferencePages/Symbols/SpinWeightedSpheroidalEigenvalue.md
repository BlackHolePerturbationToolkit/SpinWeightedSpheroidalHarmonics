{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalEigenvalue computes the eigenvalue for the oblate spheroidal equation with oblateness \[Gamma].",
    "For \[Gamma]=0, SpinWeightedSpheroidalEigenvalue reduces to l(l+1) - s(s+1).",
    "For s=0, SpinWeightedSpheroidalEigenvalue[0, l, m, \[Gamma]] is equivalent to SpheroidalEigenvalue[l, m, I \[Gamma]] - 2 m \[Gamma].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalEigenvalue produces a numerical value of the same precision.",
    "Series expansions in \[Gamma] can be computed to arbitrary order for generic s, l and m."
    },
  "Numerical Evaluation" -> {
    "SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 0.1]",
    "SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 0.1`32]",
    },
  "Series Expansion" -> {
    "Series[SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]], {\[Gamma], 0, 2}]",
    "Series[SpinWeightedSpheroidalEigenvalue[-2, 2, 2, \[Gamma]], {\[Gamma], 0, 2}]",
    },
  "Option Descriptions" -> {
    "MaxIterations" -> "The maximum number of iterations to use when computing numerical values."
   },
  "See Also" -> {"SpinWeightedSpheroidalHarmonicS", "SpinWeightedSphericalHarmonicY"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}