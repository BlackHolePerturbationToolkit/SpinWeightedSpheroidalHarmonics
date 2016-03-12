{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalEigenvalue computes the eigenvalue for the oblate spheroidal equation with oblateness \[Gamma].",
    "SpinWeightedSpheroidalEigenvalue is the eigenvalue, \[Lambda], appearing in the spin-weighted spheroidal equation, D[Sin[\[Theta]] D[ SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]], \[Theta]], \[Theta]]/ Sin[\[Theta]] + (-\[Gamma]^2 Sin[\[Theta]]^2 - (m + s Cos[\[Theta]])^2/Sin[\[Theta]]^2 - 2 s \[Gamma] Cos[\[Theta]] + s + 2 m \[Gamma] + \[Lambda]) SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] == 0.",
    "For \[Gamma]=0, SpinWeightedSpheroidalEigenvalue reduces to l(l+1) - s(s+1).",
    "For s=0, SpinWeightedSpheroidalEigenvalue[0, l, m, \[Gamma]] is equivalent to SpheroidalEigenvalue[l, m, I \[Gamma]] - 2 m \[Gamma].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalEigenvalue produces a numerical value of the same precision.",
    "Numerical values are computed using a series expansion about \[Gamma] = 0, with the number of terms in the expansion determined automatically to ensure the desired precision is reached.",
    "The MaxIterations option controls the maximum number of terms to include in a numerical approximation. For large \[Gamma]/l it may be necessary to set this to a large value in order to achieve the desired accuracy.",
    "For sufficiently large \[Gamma]/l the numerical evaluation of SpinWeightedSpheroidalEigenvalue may fail to converge for any value of MaxIterations.",
    "Series expansions about \[Gamma] = 0 can be computed to arbitrary order for generic s, l and m."
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