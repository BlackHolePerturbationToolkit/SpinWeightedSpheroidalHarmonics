VerificationTest[
  SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 0.1`32];
  SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 0.1`32, Method -> "Leaver"];
  SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 0.1`32, Method -> {"SphericalExpansion", "NumTerms" -> 6}];
  SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32][0.1`32, 0.2`32];
  SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32, Method -> "Leaver"][0.1`32, 0.2`32];
  SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 0.1`32, Method -> {"SphericalExpansion", "NumTerms" -> 6}][0.1`32, 0.2`32];
  Series[SpinWeightedSpheroidalEigenvalue[2, 2, 2, \[Gamma]], {\[Gamma], 0, 1}];
  Series[SpinWeightedSpheroidalEigenvalue[2, 2, 2, \[Gamma]], {\[Gamma], \[Infinity], 1}];
  Names["SpinWeightedSpheroidalHarmonics`Private`*$" ~~ DigitCharacter ..],
  {},
  TestID -> "Memory leaks"
]
