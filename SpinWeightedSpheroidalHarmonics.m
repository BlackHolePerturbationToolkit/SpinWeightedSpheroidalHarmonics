(* ::Package:: *)

BeginPackage["SpinWeightedSpheroidalHarmonics`"];

SpinWeightedSphericalHarmonicY::usage = "SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]] gives the spin-weighted spherical harmonic with of spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicS::usage = "SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] gives the spin-weighted oblate spheroidal harmonic with spheroidicity \[Gamma], spin-weight s, degree l and order m.";
SpinWeightedSpheroidalEigenvalue::usage = "SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] gives the spin-weighted oblate spheroidal eigenvalue with spheroidicity \[Gamma], spin-weight s, degree l and order m.";

Begin["`Private`"];

(**********************************************************)
(* Internal functions                                     *)
(**********************************************************)

(* Series expansion coefficients *)
\[Alpha][s_,l_,m_]:=1/l (Sqrt[l^2-m^2] Sqrt[l^2-s^2])/(Sqrt[2l-1] Sqrt[2l+1]);
\[Alpha][s_,l_,m_]/;l<Abs[s]=0;
\[Alpha][0,l_,m_]:= Sqrt[l^2-m^2]/(Sqrt[2l-1] Sqrt[2l+1]);
\[Beta][s_,l_,m_]:=-((m s)/(l(l+1)));
\[Beta][s_,l_,m_]/;l<Abs[s]=0;
\[Beta][0,l_,m_]=0;
s\[Lambda]lm[s_,l_,m_][0]:=l(l+1)-s(s+1); 
s\[Lambda]lm[s_,l_,m_][i_Integer?Positive]:=-(1/d[s,l,m][0,0])(Sum[d[s,l,m][i-k,0] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+l,m] \[Alpha][s,l,m] d[s,l,m][-2+i,-2]+\[Alpha][s,l,m] (\[Beta][s,-1+l,m]+\[Beta][s,l,m]) d[s,l,m][-2+i,-1]+(-1+\[Alpha][s,l,m]^2+\[Alpha][s,1+l,m]^2+\[Beta][s,l,m]^2) d[s,l,m][-2+i,0]+\[Alpha][s,1+l,m] (\[Beta][s,l,m]+\[Beta][s,1+l,m]) d[s,l,m][-2+i,1]+\[Alpha][s,1+l,m] \[Alpha][s,2+l,m] d[s,l,m][-2+i,2]-2 s \[Alpha][s,l,m] d[s,l,m][-1+i,-1]+2 (m-s \[Beta][s,l,m]) d[s,l,m][-1+i,0]-2 s \[Alpha][s,1+l,m] d[s,l,m][-1+i,1](*+d[s,l,m][i,0] s\[Lambda]lm[s,l,m][0]-(l-s) (1+l+s) d[s,l,m][i,0]*))
d[s_,l_,m_][i_?Negative,j_]=0
d[s_,l_,m_][i_?NonNegative,j_]/;Abs[j]>i=0
d[s_,l_,m_][0,0]=1;
d[s_,l_,m_][1,1]:=-((s \[Alpha][s,1+l,m])/(1+l));
d[s_,l_,m_][1,0]=0;
d[s_,l_,m_][1,-1]:=(s \[Alpha][s,l,m])/l;
d[0,l_,m_][1,1]:=0;
d[0,l_,m_][1,-1]:=0;
d[s_,l_,m_][i_Integer?Positive,0]:=- 1/(2d[s,l,m][0,0]) Sum[d[s,l,m][i-k,j] d[s,l,m][k,j],{k,1,i-1},{j,-i,i}];
d[s_,l_,m_][i_Integer?Positive,j_Integer]:=1/(j (1+j+2 l)) (Sum[d[s,l,m][i-k,j] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+j+l,m] \[Alpha][s,j+l,m] d[s,l,m][-2+i,-2+j]+\[Alpha][s,j+l,m] (\[Beta][s,-1+j+l,m]+\[Beta][s,j+l,m]) d[s,l,m][-2+i,-1+j]+(-1+\[Alpha][s,j+l,m]^2+\[Alpha][s,1+j+l,m]^2+\[Beta][s,j+l,m]^2) d[s,l,m][-2+i,j]+\[Alpha][s,1+j+l,m] (\[Beta][s,j+l,m]+\[Beta][s,1+j+l,m]) d[s,l,m][-2+i,1+j]+\[Alpha][s,1+j+l,m] \[Alpha][s,2+j+l,m] d[s,l,m][-2+i,2+j]-2 s \[Alpha][s,j+l,m] d[s,l,m][-1+i,-1+j]+2 (m-s \[Beta][s,j+l,m]) d[s,l,m][-1+i,j]-2 s \[Alpha][s,1+j+l,m] d[s,l,m][-1+i,1+j]);
d[s_,l_,m_][i_Integer?Positive,j_Integer]/;l+j<Abs[s]:=0;

simplify[expr_] := Collect[expr, {HoldPattern[\[Alpha][__]], HoldPattern[\[Beta][__]]}, Simplify];

(* Spectral method for seed in numerical evaluation *)
kHat[s_, 0, m_, \[Gamma]_] := -l (1+l)+\[Gamma]^2/3;
kHat[s_, l_, m_, \[Gamma]_] := -l (1+l)+(2 m s^2 \[Gamma])/(l+l^2)+1/3 (1+(2 (l+l^2-3 m^2) (l+l^2-3 s^2))/(l (-3+l+8 l^2+4 l^3))) \[Gamma]^2;
k1[s_, l_, m_, \[Gamma]_] /; l<2 := 0;
k1[s_, l_, m_, \[Gamma]_] := (Sqrt[((-1+l-m) (l-m) (-1+l+m) (l+m) (-1+l-s) (l-s) (-1+l+s) (l+s))/((-3+2 l) (1+2 l))] \[Gamma]^2)/((-1+l) l (-1+2 l));
k2[s_, l_, m_, \[Gamma]_] := (Sqrt[((1+l-m) (2+l-m) (1+l+m) (2+l+m) (1+l-s) (2+l-s) (1+l+s) (2+l+s))/((1+2 l) (5+2 l))] \[Gamma]^2)/((1+l) (2+l) (3+2 l));
kTilde1[s_, 0, m_, \[Gamma]_] := 0;
kTilde1[s_, 1, m_, \[Gamma]_] := -((2 s Sqrt[((l^2-m^2) (l^2-s^2))/(-1+4 l^2)] \[Gamma])/l);

kTilde1[s_, l_, m_, \[Gamma]_] := -((2 s Sqrt[((l^2-m^2) (l^2-s^2))/(-1+4 l^2)] \[Gamma] (-1+l^2+m \[Gamma]))/(l (-1+l^2)));
kTilde2[s_, 0, m_, \[Gamma]_] := -((2 s Sqrt[((1+2 l+l^2-m^2) (1+2 l+l^2-s^2))/(3+8 l+4 l^2)] \[Gamma])/(1+l));
kTilde2[s_, l_, m_, \[Gamma]_] := -((2 s Sqrt[((1+2 l+l^2-m^2) (1+2 l+l^2-s^2))/(3+8 l+4 l^2)] \[Gamma] (2 l+l^2+m \[Gamma]))/(l (2+3 l+l^2)));

SWSHEigenvalueSpectral[s_, l_, m_, \[Gamma]_] :=
 Module[{n, Matrix, Eigens, lmin},
  n = l - m + 15; (* FIXME: This should be dependent on \[Gamma] *)
  lmin = Max[Abs[s], Abs[m]];
  Matrix = SparseArray[{
    {i_, i_} :> kHat[s, lmin+i-1, m, \[Gamma]],
    {i_, j_} /; j-i==-2 :> k2[s, lmin+i-3, m, \[Gamma]],
    {i_, j_} /; j-i==-1 :> kTilde1[s, lmin+i-1, m, \[Gamma]],
    {i_, j_} /; j-i==1 :> kTilde2[s, lmin+i-1, m, \[Gamma]],
    {i_, j_} /; j-i==2 :> k1[s, lmin+i+1, m, \[Gamma]]
    }, {n, n}];

  Eigens = Sort[-Eigenvalues[Matrix], Greater];
  Eigens[[-(l-lmin)-1]]-s(s+1)
];

(* Leaver's Method *)
SWSHEigenvalueLeaver[s_, l_, m_, \[Gamma]_, Aini_] :=
 Module[{Myprec, Nmax, nInv, \[Alpha], \[Beta], \[Alpha]n, \[Beta]n, \[Gamma]n, n, LHS, RHS, Eq, A, Aval, Avar},
  Myprec = Floor[Precision[\[Gamma]]];
  Nmax = 100;(* FIXME: This should be dependent on \[Gamma] *)
  nInv = l-Abs[m];
  \[Alpha] = Abs[m+s];
  \[Beta] = Abs[m-s];
  \[Alpha]n[n_] := (-4\[Gamma](n+\[Alpha]+1)(n+\[Beta]+1)(n+(\[Alpha]+\[Beta])/2+1+s ))/((2n+\[Alpha]+\[Beta]+2)(2n+\[Alpha]+\[Beta]+3));
  \[Beta]n[n_, A_] := A + s(s+1)+\[Gamma]^2 -(n +(\[Alpha]+\[Beta])/2)(n +(\[Alpha]+\[Beta])/2+1) + If[s!=0,(8m s^2 \[Gamma])/((2n+\[Alpha]+\[Beta])(2n+\[Alpha]+\[Beta]+2)),0];
  \[Gamma]n[n_] := (4\[Gamma] n(n+\[Alpha]+\[Beta])(n+(\[Alpha]+\[Beta])/2-s))/((2n+\[Alpha]+\[Beta]-1)(2n+\[Alpha]+\[Beta]));
  RHS[Ax_] := -ContinuedFractionK[-\[Alpha]n[n-1] \[Gamma]n[n], \[Beta]n[n,Ax], {n, nInv+1, Nmax}];
  LHS[Ax_] := \[Beta]n[nInv, Ax] + ContinuedFractionK[-\[Alpha]n[nInv-n] \[Gamma]n[nInv-n+1], \[Beta]n[nInv-n, Ax], {n, 1, nInv}];
  Eq[A_?NumericQ] := LHS[A] - RHS[A];
  Aval = Avar /. FindRoot[Eq[Avar]==0, {Avar, Aini}, AccuracyGoal -> Myprec-3, WorkingPrecision -> Myprec, Method -> "Secant"];
  Aval
]

(**********************************************************)
(* SpinWeightedSpheroidalEigenvalue                       *)
(**********************************************************)

SpinWeightedSpheroidalEigenvalue::ncvb = "Failed to converge after `1` iterations. SpinWeightedSpheroidalEigenvalue obtained `2` and `3` for the result and error estimates. Increasing the value of the MaxIterations option may yield a more accurate answer.";

SyntaxInformation[SpinWeightedSpheroidalEigenvalue] =
 {"ArgumentsPattern" -> {_, _, _, _}};
SetAttributes[SpinWeightedSpheroidalEigenvalue, {NumericFunction, Listable}];

SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_] /; l < Abs[s] := 0;

SpinWeightedSpheroidalEigenvalue[s_, l_, m_, (0|0.)] :=
  l(l+1) - s(s+1);

SpinWeightedSpheroidalEigenvalue[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ] :=
 Module[{Aini},
  Aini = SWSHEigenvalueSpectral[s, l, m, N[\[Gamma]]];
  SWSHEigenvalueLeaver[s, l, m, \[Gamma], SetPrecision[Aini, Precision[\[Gamma]]]] - 2 m \[Gamma] + \[Gamma]^2
];

SpinWeightedSpheroidalEigenvalue /: N[SpinWeightedSpheroidalEigenvalue[s_Integer, l_Integer, m_Integer, \[Gamma]_?NumericQ], Nopts___] :=
  SpinWeightedSpheroidalEigenvalue[s, l, m, N[\[Gamma], Nopts]];

SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_], {\[Gamma]_, 0, order_}] :=
 Module[{i, j, coeffs}, Internal`InheritedBlock[{d, s\[Lambda]lm}, Block[{\[Alpha], \[Beta]},
  Do[
    Do[
      d[s, l, m][i, j] = simplify[d[s, l, m][i, j]], {j, -i, i}]; 
    s\[Lambda]lm[s, l, m][i] = simplify[s\[Lambda]lm[s, l, m][i]];
  , {i, 0, order}];
  coeffs = Table[s\[Lambda]lm[s, l, m][i], {i, 0, order}];
  SeriesData[\[Gamma], 0, coeffs, 0, order + 1, 1]
]]];

(**********************************************************)
(* SpinWeightedSpheroidalHarmonicS                        *)
(**********************************************************)

SpinWeightedSpheroidalHarmonicS::ncvb = "Failed to converge after `1` iterations. SpinWeightedSpheroidalHarmonicS obtained `2` and `3` for the result and error estimates. Increasing the value of the MaxIterations option may yield a more accurate answer.";

SyntaxInformation[SpinWeightedSpheroidalHarmonicS] =
 {"ArgumentsPattern" -> {_, _, _, _, _, _, ___}};
Options[SpinWeightedSpheroidalHarmonicS] = {MaxIterations -> Automatic};
SetAttributes[SpinWeightedSpheroidalHarmonicS, {NumericFunction, Listable}];

SpinWeightedSpheroidalHarmonicS[s_, l_, m_, (0|0.), \[Theta]_, \[Phi]_, opts___] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]];

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ, OptionsPattern[]] :=
 Module[{Si, S, maxOrder, i, j}, Internal`InheritedBlock[{d, s\[Lambda]lm, SpinWeightedSphericalHarmonicY},
  maxOrder = OptionValue["MaxIterations"];
  If[maxOrder == Automatic, maxOrder = Max[32, Precision[\[Gamma]]]];
  S = 0;
  Do[
    Do[
      d[s, l, m][i, j] = Simplify[d[s, l, m][i, j]];
      SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]] = SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]];
    , {j, -i, i}];
    s\[Lambda]lm[s, l, m][i] = Simplify[s\[Lambda]lm[s, l, m][i]];
    Si = Sum[d[s, l, m][i, j] SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]] \[Gamma]^i, {j, -i, i}];
    S += Si;
    Which[
      Si != 0 && S == (S-Si),
        Break[],
      i == maxOrder,
        Message[SpinWeightedSpheroidalHarmonicS::ncvb, maxOrder, S, Si];
    ];
    ,{i, 0, maxOrder}];
  S
]];

SpinWeightedSpheroidalHarmonicS /: N[SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?NumericQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ, opts___:OptionsPattern[]], Nopts___] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, N[\[Gamma], Nopts], \[Theta], \[Phi], opts];

SpinWeightedSpheroidalHarmonicS /: 
  Series[SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_, \[Theta]_, \[Phi]_], {\[Gamma]_, 0, order_}] :=
 Module[{i, j, coeffs}, Internal`InheritedBlock[{d, s\[Lambda]lm}, Block[{\[Alpha], \[Beta]},
  Do[
    Do[
      d[s, l, m][i, j] = simplify[d[s, l, m][i, j]], {j, -i, i}]; 
    s\[Lambda]lm[s, l, m][i] = simplify[s\[Lambda]lm[s, l, m][i]];
  , {i, 0, order}];
  coeffs = Table[Sum[d[s, l, m][i, j] SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]], {j, -i, i}], {i, 0, order}];
  SeriesData[\[Gamma], 0, coeffs, 0, order + 1, 1]
]]];


(**********************************************************)
(* SpinWeightedSphericalHarmonicY                         *)
(**********************************************************)

SyntaxInformation[SpinWeightedSphericalHarmonicY] =
 {"ArgumentsPattern" -> {_, _, _, _, _, ___}};
SetAttributes[SpinWeightedSphericalHarmonicY, {NumericFunction, Listable}];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] :=
  (-1)^m Sqrt[((l+m)!(l-m)!(2l+1))/(4\[Pi] (l+s)!(l-s)!)] Sum[Binomial[l-s,r] Binomial[l+s,r+s-m] (-1)^(l-r-s) If[(2 l - 2 r - s + m)==0, 1, Sin[\[Theta]/2]^(2 l - 2 r - s + m)] Cos[\[Theta]/2]^(2 r + s - m),{r,Max[m-s,0],Min[l-s,l+m]}]Exp[I m \[Phi]];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] /; Abs[m]>l = 0;

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, 0.] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0];

SpinWeightedSphericalHarmonicY[0, l_, m_, \[Theta]_, \[Phi]_] :=
  SphericalHarmonicY[l, m, \[Theta], \[Phi]];

End[];

EndPackage[];
