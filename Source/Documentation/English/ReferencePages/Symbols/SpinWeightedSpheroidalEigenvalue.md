{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalEigenvalue computes the eigenvalue for the oblate spheroidal equation with oblateness \[Gamma].",
    "SpinWeightedSpheroidalEigenvalue is the eigenvalue, \[Lambda], appearing in the spin-weighted spheroidal equation, D[Sin[\[Theta]] D[ SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]], \[Theta]], \[Theta]]/ Sin[\[Theta]] + (-\[Gamma]^2 Sin[\[Theta]]^2 - (m + s Cos[\[Theta]])^2/Sin[\[Theta]]^2 - 2 s \[Gamma] Cos[\[Theta]] + s + 2 m \[Gamma] + \[Lambda]) SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] == 0.",
    "The definition for \[Lambda] is consistent with that of [Teukolsky, S. A., Astrophys. J. 185, 635–647 (1973). doi:10.1086/152444] and that of [Sasaki, M. & Tagoshi, H. Living Rev. Relativ. (2003) 6: 6. doi:10.12942/lrr-2003-6]."
    "For \[Gamma]=0, SpinWeightedSpheroidalEigenvalue reduces to l(l+1) - s(s+1).",
    "For s=0, SpinWeightedSpheroidalEigenvalue[0, l, m, \[Gamma]] is equivalent to SpheroidalEigenvalue[l, m, I \[Gamma]] - 2 m \[Gamma].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalEigenvalue produces a numerical value of the same precision.",
    "Numerical values are computed using a spectral/eigenvalue approach for an initial guess followed by Leaver's method to obtain an accurate answer. The precision of the result will be approximately equal to the precision of \[Gamma].",
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
  "See Also" -> {"SpinWeightedSpheroidalHarmonicS", "SpinWeightedSphericalHarmonicY"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}