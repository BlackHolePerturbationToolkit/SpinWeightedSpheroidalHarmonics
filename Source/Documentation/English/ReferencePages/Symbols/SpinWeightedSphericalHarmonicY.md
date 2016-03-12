{
  "More Information" -> {
    "For real \[Gamma], SpinWeightedSphericalHarmonicY computes the spin-weighted spherical harmonic.",
    "For s=0, SpinWeightedSphericalHarmonicY is equivalent to SphericalHarmonicY.",
    "SpinWeightedSphericalHarmonicY is unit-normalised on the two-sphere, Integrate[SpinWeightedSphericalHarmonicY[s, l1, m1, \[Theta], \[Phi]] Conjugate[SpinWeightedSphericalHarmonicY[s, l2, m2, \[Theta], \[Phi]]] Sin[\[Theta]], {\[Theta], 0, \[Pi]}, {\[Phi], 0, 2 \[Pi]}] \[Equal] KroneckerDelta[l1, l2] KroneckerDelta[m1, m2].",
  },
  "Numerical Evaluation" -> {
    "SpinWeightedSphericalHarmonicY[-2, 2, 2, 0.1, 0.3]"
    },
  "Analytic expressions" -> {
    "SpinWeightedSphericalHarmonicY[-2, 2, 2, \[Theta], \[Phi]]",
    },
  "See Also" -> {"SpinWeightedSpheroidalEigenvalue", "SpinWeightedSpheroidalHarmonicS"},
  "More About" -> {"SpinWeightedSpheroidalHarmonics"},
  "Tutorials" -> {"SpinWeightedSpheroidalHarmonics"}
}