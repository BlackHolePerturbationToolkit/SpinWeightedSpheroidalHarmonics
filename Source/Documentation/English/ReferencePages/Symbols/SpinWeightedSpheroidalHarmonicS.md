{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSpheroidalHarmonicS computes the oblate spin-weighted spheroidal harmonic with oblateness \[Gamma].",
    "SpinWeightedSpheroidalHarmonicS satisfies the spin-weighted spheroidal equation, D[Sin[\[Theta]] D[ SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]], \[Theta]], \[Theta]]/ Sin[\[Theta]] + (-\[Gamma]^2 Sin[\[Theta]]^2 - (m + s Cos[\[Theta]])^2/Sin[\[Theta]]^2 - 2 s \[Gamma] Cos[\[Theta]] + s + 2 m \[Gamma] + \[Lambda]) SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] == 0.",
    "SpinWeightedSpheroidalHarmonicS is unit-normalised on the two-sphere, Integrate[SpinWeightedSpheroidalHarmonicS[s, l1, m1, \[Gamma], \[Theta], \[Phi]] Conjugate[SpinWeightedSpheroidalHarmonicS[s, l2, m2, \[Gamma], \[Theta], \[Phi]]] Sin[\[Theta]], {\[Theta], 0, \[Pi]}, {\[Phi], 0, 2 \[Pi]}] \[Equal] KroneckerDelta[l1, l2] KroneckerDelta[m1, m2]. This is consistent with the Meixner-Sch\[ADoubleDot]fke convention for the spheroidal Legendre functions in the same sense that the normalisation of SphericalHarmonicY is consistent with the Legendre functions.",
    "For \[Gamma]=0, SpinWeightedSpheroidalHarmonicS reduces to SpinWeightedSphericalHarmonicY.",
    "For s=0 and \[Gamma]=0, SpinWeightedSpheroidalHarmonicS is equivalent to SphericalHarmonicY.",
    "For s=0, SpinWeightedSpheroidalHarmonicS[0, l, m, \[Gamma], \[Theta], 0] is equivalent to Sqrt[(2 l + 1)/(4 \[Pi])] Sqrt[(l - m)!/(l + m)!] SpheroidalPS[l, m, I \[Gamma], Cos[\[Theta]]].",
    "For numerical values of \[Gamma], SpinWeightedSpheroidalHarmonicS produces a numerical value of the same precision.",
    "Numerical values are computed using either Leaver's method or an expansion in terms of spin-weighted spherical harmonics.",
    "The number of terms used for numerical evaluation is such that the normalisation is determined as accurately as possible given the precision of the input \[Gamma].",
    "Series expansions about \[Gamma] = 0 can be computed to arbitrary order for generic s, l and m.",
    "SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma]] produces a SpinWeightedSpheroidalHarmonicSFunction, which may be efficiently evaluated for multiple values of \[Theta] and \[Phi].",
    "For Methods that support suboptions, they may be passed as a list. For example Method \[Rule] {\"Leaver\", \"InitialGuess\" \[Rule] 0.1} sets the \"InitialGuess\" suboption to have value 0.1.",
    "The possible suboptions that the \"Lever\" method accepts are: \"NumTerms\" \[Rule] n.",
    "The possible suboptions that the \"SphericalExpansion\" method accepts are: \"NumTerms\" \[Rule] n."
    },
  "Numerical Evaluation" -> {
    "SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1, \[Pi]/4, 0]",
    "SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32, \[Pi]/4, 0]",
    "Slm = SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32]; Slm[\[Pi]/4, 0]"
    },
  "Series Expansion" -> {
    "Series[SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]], {\[Gamma], 0, 2}]",
    "Series[SpinWeightedSpheroidalHarmonicS[-2, 2, 2, \[Gamma], \[Theta], \[Phi]], {\[Gamma], 0, 2}]",
    },
  "Option Descriptions" -> {
    "Method" -> "The method to use to numerically evaluate the harmonic. Possible choices are \"Leaver\" and \"SphericalExpansion\"."
   },
  "See Also" -> {"SpinWeightedSpheroidalEigenvalue", "SpinWeightedSpheroidalHarmonicSFunction", "SpinWeightedSphericalHarmonicY"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}