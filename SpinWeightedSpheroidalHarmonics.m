(* ::Package:: *)

BeginPackage["SpinWeightedSpheroidalHarmonics`"];

SpinWeightedSphericalHarmonicY::usage = "SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]] gives the spin-weighted spherical harmonic with of spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicS::usage = "SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] gives the spin-weighted oblate spheroidal harmonic with spheroidicity \[Gamma], spin-weight s, degree l and order m.";
SpinWeightedSpheroidalEigenvalue::usage = "SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] gives the spin-weighted oblate spheroidal eigenvalue with spheroidicity \[Gamma], spin-weight s, degree l and order m.";

Begin["`Private`"];

(**********************************************************)
(* Internal functions                                     *)
(**********************************************************)

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

(**********************************************************)
(* SpinWeightedSpheroidalEigenvalue                       *)
(**********************************************************)

SpinWeightedSpheroidalEigenvalue::ncvb = "Failed to converge after `1` iterations. SpinWeightedSpheroidalEigenvalue obtained `2` and `3` for the result and error estimates. Increasing the value of the MaxIterations option may yield a more accurate answer.";

Options[SpinWeightedSpheroidalEigenvalue] = {MaxIterations -> Automatic};
SetAttributes[SpinWeightedSpheroidalEigenvalue, {NumericFunction, Listable}];

SpinWeightedSpheroidalEigenvalue[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, OptionsPattern[]] :=
 Module[{\[Lambda]i, \[Lambda], maxOrder, i, j}, Internal`InheritedBlock[{d, s\[Lambda]lm},
  maxOrder = OptionValue["MaxIterations"];
  If[maxOrder == Automatic, maxOrder = Max[32, Precision[\[Gamma]]]];
  \[Lambda] = 0;
  Do[
    Do[
      d[s, l, m][i, j] = Simplify[d[s, l, m][i, j]], {j, -i, i}];
    s\[Lambda]lm[s, l, m][i] = Simplify[s\[Lambda]lm[s, l, m][i]];
    \[Lambda]i = s\[Lambda]lm[s, l, m][i] \[Gamma]^i;
    \[Lambda] += \[Lambda]i;
    Which[
      \[Lambda]i != 0 && \[Lambda] == (\[Lambda]-\[Lambda]i),
      Break[],
      i == maxOrder,
      Message[SpinWeightedSpheroidalEigenvalue::ncvb, maxOrder, \[Lambda], \[Lambda]i];
    ];
  ,{i, 0, maxOrder}];
  \[Lambda]
]];

SpinWeightedSpheroidalEigenvalue /: N[SpinWeightedSpheroidalEigenvalue[s_Integer, l_Integer, m_Integer, \[Gamma]_?NumericQ, opts___:OptionsPattern[]], Nopts___] :=
  SpinWeightedSpheroidalEigenvalue[s, l, m, N[\[Gamma], Nopts], opts];

SpinWeightedSpheroidalEigenvalue[s_, l_, m_, (0|0.), opts___:OptionsPattern[]] :=
  l(l+1) - s(s+1);

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

Options[SpinWeightedSpheroidalHarmonicS] = {MaxIterations -> Automatic};
SetAttributes[SpinWeightedSpheroidalHarmonicS, {NumericFunction, Listable}];

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

SpinWeightedSpheroidalHarmonicS[s_, l_, m_, (0|0.), \[Theta]_, \[Phi]_, opts___] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]];

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

SetAttributes[SpinWeightedSphericalHarmonicY, {NumericFunction, Listable}];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] :=
  (-1)^m Sqrt[((l+m)!(l-m)!(2l+1))/(4\[Pi] (l+s)!(l-s)!)] Sin[\[Theta]/2]^(2l) Sum[Binomial[l-s,r] Binomial[l+s,r+s-m] (-1)^(l-r-s)  Cot[\[Theta]/2]^(2 r+s-m),{r,Max[m-s,0],Min[l-s,l+m]}]Exp[I m \[Phi]];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] /; Abs[m]>l = 0;

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, 0.] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0];

End[];

EndPackage[];
