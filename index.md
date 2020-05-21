{% include head.html %}

<p>
 <h1 style="display:inline">SpinWeightedSpheroidalHarmonics</h1> <span style="float:right;"><a href="https://bhptoolkit.org/mathematica-install.html" class = "code_btn">Install this package!</a></span>
</p>

The SpinWeightedSpheroidalHarmonics package for Mathematica provides functions for computing
spin-weighted spheroidal harmonics, spin-weighted spherical harmonics and their associated eigenvalues.
Support is included for both arbitrary-precision numerical evaluation, and for series expansions.

<p align="center"><img width="50%" src="swsh.png" alt="S(s=-2, l=2, gamma=1.9)"/></p>

Explicitly the package computes solutions to the equation

$$\left[\frac{1}{\sin\theta} \frac{d}{d\theta} \left( \sin\theta \frac{d}{d\theta}\right) - \gamma^2 \sin^2\theta - \frac{(m+s\cos\theta)^2}{\sin^2\theta} - 2\gamma s \cos\theta + s+ 2m\gamma + {}_s\lambda_{lm} \right] { }_{s}S_{lm} = 0 \nonumber $$

where

$S$ is the spin-weighted spheroidal harmonic  
$s$ is the spin-weight  
$l,m$ are the multipolar indices  
$\gamma$ is the spheroidicity   
$\lambda$ is the spheroidal eigenvalue

## Numerical evaluation

Evaluating the spheroidal eigenfunction for $s=-2,l=2,m=2,\gamma=1.5$ is done via

```
SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 1.5`30, π/2, 0]
0.066929509191673988268462820
```

Likewise the eigenvalue, $\lambda$, is computed by

```
SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 1.5`30]
-5.5776273646788261891028539326
```

Note the output tracks the precision of the input so it is very easy to compute high precision results.

## Small frequency expansions

Small frequency expansions can be performed for the eigenfunctions and eigenvalues

For given $s,l,m$ the eigenfunction can be expanded via, e.g.,

```
Series[SpinWeightedSpheroidalaHarmonicS[2, 3, 2, γ], {γ, 0, 1}]
```

This returns

$\frac{1}{2} \sqrt{\frac{7}{\pi }} e^{2 i \phi } \sin ^4\left(\frac{\theta }{2}\right) (3 \cos (\theta )+2)-\frac{1}{72} \gamma \left(\sqrt{\frac{7}{\pi }} e^{2 i \phi } \sin^4\left(\frac{\theta }{2}\right) (54 \cos (\theta )+27 \cos (2 \theta)+29)\right)+O\left(\gamma ^2\right)$
  
The series expansion above can also be performed for generic values of $s,l,m$ which returns a series in `SpinWeightedSphericalHarmonicsY` functions.

For generic $s,l,m$ the eigenvalues are expanded via:
```
Series[SpinWeightedSpheroidalEigenvalue[s, l, m, γ], {γ, 0, 1}]
```

returns

$\left(l^2+l-s (s+1)\right) + \gamma  \left(-\frac{2 m s^2}{l (l+1)}-2 m\right) + O\left(\gamma ^2\right).$

The code is fast but note that in general the series expansion will be quicker to compute if you provide explicit values of $s,l,m$.

## High frequency expansions

You can also calculate series expansions of the eigenvalue about $\gamma = \infty$. For example, for the $s=2,l=2,m=2$ case, the command:
```
Series[SpinWeightedSpheroidalEigenvalue[2, 2, 2, γ], {γ, ∞, 6}]
```
returns  

$-6 \gamma  -1 + \frac{3}{4 \gamma } -\frac{15}{64 \gamma ^3} -\frac{3}{16 \gamma ^4} +\frac{3}{512 \gamma ^5} + \frac{27}{128 \gamma ^6} +O\left[\frac{1}{\gamma}\right]^7$

Currently, for the high frequency expansions, you have to given explicit values for $s,l,m$.

## Further examples

See the Documentation Centre for a tutorial and documentation on individual commands. See also the [Mathematica Toolkit Examples](https://github.com/BlackHolePerturbationToolkit/MathematicaToolkitExamples) for further example notebooks.