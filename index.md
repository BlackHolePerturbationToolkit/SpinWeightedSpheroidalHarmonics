{% include head.html %}

<p>
 <h1 style="display:inline">SpinWeightedSpheroidalHarmonics</h1> <span style="float:right;"><a href="{{ site.github.repository_url }}" class = "code_btn">Get the code!</a></span>
</p>

The SpinWeightedSpheroidalHarmonics package for Mathematica provides functions for computing
spin-weighted spheroidal harmonics, spin-weighted spherical harmonics and their associated eigenvalues.
Support is included for both arbitrary-precision numerical evaluation, and for series expansions.

<p align="center"><img  src="swsh.png" alt="S(s=-2, l=2, gamma=1.9)"/></p>

Explicitly the package computes solutions to:

## Numerical evaluation


## Small frequency expansions

Small frequency expansions can be performed for the eigenfunctions and eigenvalues


```
Series[SpinWeightedSpheroidalEigenvalue[s, l, m, γ], {γ, 0, 6}]
```

## High frequency expansions

You can easily calculate series expansions of the eigenvalue about $\gamma = \infty$. For example, for the $s=2,l=2,m=2$ case, the command:
```
Series[SpinWeightedSpheroidalEigenvalue[2, 2, 2, γ], {γ, ∞, 6}]
```
returns  

$-6 \gamma  -1 + \frac{3}{4 \gamma } -\frac{15}{64 \gamma ^3} -\frac{3}{16 \gamma ^4} +\frac{3}{512 \gamma ^5} + \frac{27}{128 \gamma ^6} +O\left[\frac{1}{\gamma}\right]^7$
