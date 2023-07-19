(* ::Package:: *)

(* ::Title:: *)
(*SpinWeightedSpheoridalHarmonics package*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["SpinWeightedSpheroidalHarmonics`"];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{SpinWeightedSphericalHarmonicY, SpinWeightedSpheroidalHarmonicS, SpinWeightedSpheroidalHarmonicSFunction, SpinWeightedSpheroidalEigenvalue}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


SpinWeightedSphericalHarmonicY::usage = "\!\(\*RowBox[{\"SpinWeightedSphericalHarmonicY\", \"[\", RowBox[{StyleBox[\"s\", \"TI\"], \",\", StyleBox[\"l\", \"TI\"], \",\", StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"\[Theta]\", \"TR\"], \",\", StyleBox[\"\[Phi]\", \"TR\"]}], \"]\"}]\) gives the spin-weighted spherical harmonic \!\(\*RowBox[{SubscriptBox[\"\[InvisiblePrefixScriptBase]\", \"s\"], SubscriptBox[\"Y\", RowBox[{\"l\", \"m\"}]],\"(\", RowBox[{\"\[Theta]\", \",\", \"\[Phi]\"}], \")\"}]\).";
SpinWeightedSpheroidalHarmonicSFunction::usage = "\!\(\*RowBox[{\"SpinWeightedSpheroidalHarmonicSFunction\", \"[\", StyleBox[\"data\", \"TI\"], \"]\"}]\) represents a solution to the spin-weighted spheroidal equation.";
SpinWeightedSpheroidalHarmonicS::usage = "\!\(\*RowBox[{\"SpinWeightedSpheroidalHarmonicS\", \"[\", RowBox[{StyleBox[\"s\", \"TI\"], \",\", StyleBox[\"l\", \"TI\"], \",\", StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"\[Gamma]\", \"TR\"]}], \"]\", \"[\", RowBox[{StyleBox[\"\[Theta]\", \"TR\"], \",\", StyleBox[\"\[Phi]\", \"TR\"]}], \"]\"}]\) gives the spin-weighted spheroidal harmonic \!\(\*RowBox[{SubscriptBox[\"\[InvisiblePrefixScriptBase]\", \"s\"], SubscriptBox[\"S\", RowBox[{\"l\", \"m\"}]],\"(\", RowBox[{\"\[Theta]\", \",\", \"\[Phi]\", \";\", \"\[Gamma]\"}], \")\"}]\).
\!\(\*RowBox[{\"SpinWeightedSpheroidalHarmonicS\", \"[\", RowBox[{StyleBox[\"s\", \"TI\"], \",\", StyleBox[\"l\", \"TI\"], \",\", StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"\[Gamma]\", \"TR\"]}], \"]\"}]\) constructs a SpinWeightedSpheroidalHarmonicSFunction that can be evaluated for specific values of \!\(\*StyleBox[\"\[Theta]\", \"TR\"]\) and \!\(\*StyleBox[\"\[Phi]\", \"TR\"]\).";
SpinWeightedSpheroidalEigenvalue::usage = "\!\(\*RowBox[{\"SpinWeightedSpheroidalEigenvalue\", \"[\", RowBox[{StyleBox[\"s\", \"TI\"], \",\", StyleBox[\"l\", \"TI\"], \",\", StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"\[Gamma]\", \"TR\"]}], \"]\"}]\) gives the spin-weighted spheroidal eigenvalue with spin-weight s, degree l and order m.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


SpinWeightedSpheroidalEigenvalue::numterms = "Automatic determination of the number of terms to use in the SphericalExpansion method may be unreliable in certain cases. Currently using `1` terms. It is recommended to verify the result is unchanged with a different number of terms.";
SpinWeightedSpheroidalEigenvalue::optx = "Unknown options in `1`";
SpinWeightedSpheroidalEigenvalue::params = "Invalid parameters s=`1`, l=`2`, m=`3`";
SpinWeightedSpheroidalEigenvalue::infseriesparams = "Series expansion about \[Gamma]=\[Infinity] only supported for integer or half-integer parameters, but specified parameters were s=`1`, l=`2`, m=`3`";
SpinWeightedSpheroidalHarmonicS::numterms = "Automatic determination of the number of terms to use in the SphericalExpansion method may be unreliable in certain cases. Currently using `1` terms. It is recommended to verify the result is unchanged with a different number of terms.";
SpinWeightedSpheroidalHarmonicS::optx = "Unknown options in `1`";
SpinWeightedSpheroidalHarmonicS::normprec = "Normalisation cannot be determined for \"Leaver\" method so it will be set to 1. To obtain an accurate result either provide a higher precision input spheroidicity or use the \"SphericalExpansion\" method instead.";
SpinWeightedSpheroidalHarmonicS::prec = "Spin-weighted spheroidal harmonic cannot be computed using \"Leaver\" method with the given working precision. To obtain an accurate result either provide a higher precision input spheroidicity or use the \"SphericalExpansion\" method instead.";
SpinWeightedSpheroidalHarmonicS::params = "Invalid parameters s=`1`, l=`2`, m=`3`";
SpinWeightedSphericalHarmonicY::params = "Invalid parameters s=`1`, l=`2`, m=`3`";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];



(* ::Section::Closed:: *)
(*Useful functions*)


(* ::Subsection::Closed:: *)
(*Series expansion coefficients*)


\[Alpha][s_,l_,m_]:=1/l (Sqrt[l^2-m^2] Sqrt[l^2-s^2])/(Sqrt[2l-1] Sqrt[2l+1]);
\[Alpha][0,l_,m_]:= Sqrt[l^2-m^2]/(Sqrt[2l-1] Sqrt[2l+1]);
\[Alpha][s_,l_,0]:= Sqrt[l^2-s^2]/(Sqrt[2l-1] Sqrt[2l+1]);
\[Alpha][s_,0,m_]:= 0;
\[Alpha][s_,1/2,m_]:= 0;
\[Alpha][s_,-1/2,m_]:= 0;
\[Beta][s_,l_,m_]:=-((m s)/(l(l+1)));
\[Beta][s_,0,m_]=0;
\[Beta][s_,-1,m_]=0;
s\[Lambda]lm[s_,l_,m_][0]:=l(l+1)-s(s+1); 
s\[Lambda]lm[s_,l_,m_][i_Integer?Positive]:=-(1/d[s,l,m][0,0])(Sum[d[s,l,m][i-k,0] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+l,m] \[Alpha][s,l,m] d[s,l,m][-2+i,-2]+\[Alpha][s,l,m] (\[Beta][s,-1+l,m]+\[Beta][s,l,m]) d[s,l,m][-2+i,-1]+(-1+\[Alpha][s,l,m]^2+\[Alpha][s,1+l,m]^2+\[Beta][s,l,m]^2) d[s,l,m][-2+i,0]+\[Alpha][s,1+l,m] (\[Beta][s,l,m]+\[Beta][s,1+l,m]) d[s,l,m][-2+i,1]+\[Alpha][s,1+l,m] \[Alpha][s,2+l,m] d[s,l,m][-2+i,2]-2 s \[Alpha][s,l,m] d[s,l,m][-1+i,-1]+2 (m-s \[Beta][s,l,m]) d[s,l,m][-1+i,0]-2 s \[Alpha][s,1+l,m] d[s,l,m][-1+i,1](*+d[s,l,m][i,0] s\[Lambda]lm[s,l,m][0]-(l-s) (1+l+s) d[s,l,m][i,0]*))
d[s_,l_,m_][i_?Negative,j_]=0
d[s_,l_,m_][i_?NonNegative,j_]/;Abs[j]>i=0
d[s_,l_,m_][0,0]=1;
d[s_,l_,m_][1,1]:=-((s \[Alpha][s,1+l,m])/(1+l));
d[s_,l_,m_][1,0]=0;
d[s_,l_,m_][1,-1]:=(s \[Alpha][s,l,m])/l;
d[s_,0,m_][1,-1]:=0;
d[s_,l_,m_][i_Integer?Positive,0]:=- 1/(2d[s,l,m][0,0]) Sum[d[s,l,m][i-k,j] d[s,l,m][k,j],{k,1,i-1},{j,-i,i}];
d[s_,l_,m_][i_Integer?Positive,j_Integer]:=1/(j (1+j+2 l)) (Sum[d[s,l,m][i-k,j] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+j+l,m] \[Alpha][s,j+l,m] d[s,l,m][-2+i,-2+j]+\[Alpha][s,j+l,m] (\[Beta][s,-1+j+l,m]+\[Beta][s,j+l,m]) d[s,l,m][-2+i,-1+j]+(-1+\[Alpha][s,j+l,m]^2+\[Alpha][s,1+j+l,m]^2+\[Beta][s,j+l,m]^2) d[s,l,m][-2+i,j]+\[Alpha][s,1+j+l,m] (\[Beta][s,j+l,m]+\[Beta][s,1+j+l,m]) d[s,l,m][-2+i,1+j]+\[Alpha][s,1+j+l,m] \[Alpha][s,2+j+l,m] d[s,l,m][-2+i,2+j]-2 s \[Alpha][s,j+l,m] d[s,l,m][-1+i,-1+j]+2 (m-s \[Beta][s,j+l,m]) d[s,l,m][-1+i,j]-2 s \[Alpha][s,1+j+l,m] d[s,l,m][-1+i,1+j]);
d[s_,l_,m_][i_Integer?Positive,j_Integer]/; 1+j+2 l == 0:= 0;

simplify[expr_] := Collect[expr, {HoldPattern[\[Alpha][__]], HoldPattern[\[Beta][__]]}, Simplify];


(* ::Subsection::Closed:: *)
(*Functions for spectral method for seed in numerical evaluation*)


kHat[s_, 0, 0, \[Gamma]_] := \[Gamma]^2/3;
kHat[s_, l_, m_, \[Gamma]_] := -l (1+l)+(2 m s^2 \[Gamma])/(l+l^2)+1/3 (1+(2 (l+l^2-3 m^2) (l+l^2-3 s^2))/(l (-3+l+8 l^2+4 l^3))) \[Gamma]^2;
kHat[s_,1/2,m_,\[Gamma]_]/;Abs[s]==1/2:=1/864 (-648+576 m \[Gamma]+216 \[Gamma]^2+288 m^2 \[Gamma]^2);
k2[s_, l_, m_, \[Gamma]_] := (Sqrt[((1+l-m) (2+l-m) (1+l+m) (2+l+m) (1+l-s) (2+l-s) (1+l+s) (2+l+s))/((1+2 l) (5+2 l))] \[Gamma]^2)/((1+l) (2+l) (3+2 l));
kTilde2[s_, 0, 0, \[Gamma]_] := -2 s Sqrt[(1-s^2)/3] \[Gamma];
kTilde2[s_, l_, m_, \[Gamma]_] := -((2 s Sqrt[((1+2 l+l^2-m^2) (1+2 l+l^2-s^2))/(3+8 l+4 l^2)] \[Gamma] (2 l+l^2+m \[Gamma]))/(l (2+3 l+l^2)));


(* ::Subsection::Closed:: *)
(*Continued fraction*)


(* ::Text:: *)
(*This is preferred over Mathematica's ContinedFractionK function as we can get an error estimate on the result using this function.*)


CF[a_, b_, {n_, n0_}] := 
  Module[{A, B, ak, bk, res = Indeterminate, j = n0},
   A[n0 - 2] = 1;
   B[n0 - 2] = 0;
   ak[k_] := ak[k] = (a /. n -> k);
   bk[k_] := bk[k] = (b /. n -> k);
   A[n0 - 1] = 0(*bk[n0-1]*);
   B[n0 - 1] = 1;
   A[k_] := A[k] = bk[k] A[k - 1] + ak[k] A[k - 2];
   B[k_] := B[k] = bk[k] B[k - 1] + ak[k] B[k - 2];
   While[res =!= (res = A[j]/B[j]), j++];
   Clear[A, B, ak, bk];
   res
   ];




(* ::Section::Closed:: *)
(*SpinWeightedSpheroidalEigenvalue*)


(* ::Subsection::Closed:: *)
(*Spherical expansion method*)


Options[SWSHEigenvalueSpectral] = {"NumTerms" -> Automatic};


SWSHEigenvalueSpectral[s_, l_, m_, \[Gamma]_, OptionsPattern[]]:=
 Module[{nDown, nUp, Matrix, Eigens, lmin},
  (* FIXME: Improve the estimate of nmax. It should depend on the accuarcy sought. *)
  Switch[OptionValue["NumTerms"],
   Automatic,
    nUp = Ceiling[Abs[3/2\[Gamma]]]+5;
    Message[SpinWeightedSpheroidalEigenvalue::numterms, nUp];,
   _Integer,
    nUp = OptionValue["NumTerms"];,
   _,
    Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> {"SphericalExpansion", opts}];
    Return[$Failed];
  ];
  lmin=Max[Abs[s],Abs[m]];
  nDown = Min[nUp, l-lmin];

  Matrix=SparseArray[
	{{i_,i_}:>kHat[s,l-nDown-1+i,m,\[Gamma]],
	{i_,j_}/;j-i==-2:>k2[s,l-nDown-3+i,m,\[Gamma]],
	{i_,j_}/;j-i==-1:>kTilde2[s,l-nDown+i-2,m,\[Gamma]],
	{i_,j_}/;j-i==1:>kTilde2[s,l-nDown+i-1,m,\[Gamma]],
	{i_,j_}/;j-i==2:>k2[s,l-nDown+i+-1,m,\[Gamma]]}
  ,{nUp+nDown+1,nUp+nDown+1}];

  (* To choose the eigenvalue corrsponding to the desired l, we assume that the real
     part of eigenvalue is a monotonic function of l. *)
  Eigens=-Sort[Quiet[Eigenvalues[Matrix],{Eigenvalues::arh,Eigenvalues::arhm}]];

  Eigens[[-(nDown+1)]]-s(s+1)
];


(* ::Subsection::Closed:: *)
(*Leaver's method*)


Options[SWSHEigenvalueLeaver] = {"InitialGuess" -> "SphericalExpansion", "NumInversions" -> Automatic};


SWSHEigenvalueLeaver[s_, l_, m_, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{Myprec, Nmax, nInv, \[Alpha], \[Beta], \[Alpha]n, \[Beta]n, \[Gamma]n, n, LHS, RHS, Eq, A, Aval, Avar, \[Lambda]ini, Aini},
  \[Lambda]ini = OptionValue["InitialGuess"];
  If[\[Lambda]ini === "SphericalExpansion",
    \[Lambda]ini = Quiet[SetPrecision[SWSHEigenvalueSpectral[s, l, m, N[\[Gamma]]] - 2 m \[Gamma] + \[Gamma]^2, Precision[\[Gamma]]], SpinWeightedSpheroidalEigenvalue::numterms];
  ];
  Aini = \[Lambda]ini + 2 m \[Gamma] - \[Gamma]^2;
  Switch[OptionValue["NumInversions"],
   Automatic,
    nInv = l-Max[Abs[m],Abs[s]];,
   _Integer,
    nInv = OptionValue["NumInversions"];,
   _,
    Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> {"Leaver", opts}];
    Return[$Failed];
  ];
  Myprec = Max[Precision[\[Gamma]], 3];
  \[Alpha] = Abs[m+s];
  \[Beta] = Abs[m-s];
  \[Alpha]n[n_] := (-4\[Gamma](n+\[Alpha]+1)(n+\[Beta]+1)(n+(\[Alpha]+\[Beta])/2+1+s ))/((2n+\[Alpha]+\[Beta]+2)(2n+\[Alpha]+\[Beta]+3));
  \[Beta]n[n_, A_] := A + s(s+1)+\[Gamma]^2 -(n +(\[Alpha]+\[Beta])/2)(n +(\[Alpha]+\[Beta])/2+1) + If[s!=0,(8m s^2 \[Gamma])/((2n+\[Alpha]+\[Beta])(2n+\[Alpha]+\[Beta]+2)),0];
  \[Gamma]n[n_] := (4\[Gamma] n(n+\[Alpha]+\[Beta])(n+(\[Alpha]+\[Beta])/2-s))/((2n+\[Alpha]+\[Beta]-1)(2n+\[Alpha]+\[Beta]));
  RHS[Ax_] := -CF[-\[Alpha]n[n-1] \[Gamma]n[n], \[Beta]n[n,Ax], {n, nInv+1}];
  LHS[Ax_] := \[Beta]n[nInv, Ax] + ContinuedFractionK[-\[Alpha]n[nInv-n] \[Gamma]n[nInv-n+1], \[Beta]n[nInv-n, Ax], {n, 1, nInv}];
  Eq[A_?NumericQ] := LHS[A] - RHS[A];
  Aval = Avar /. Quiet[Check[FindRoot[Eq[Avar]==0, {Avar, Aini}, AccuracyGoal -> Myprec-3, WorkingPrecision -> Myprec, Method -> "Secant"], Avar -> $Failed, {Power::infy, FindRoot::nlnum}], {Power::infy, FindRoot::nlnum}];
  Clear[\[Alpha]n, \[Beta]n, \[Gamma]n, LHS, RHS, Eq];
  Aval
];


(* ::Subsection::Closed:: *)
(*SpinWeightedSpheroidalEigenvalue*)


SyntaxInformation[SpinWeightedSpheroidalEigenvalue] =
 {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};


Options[SpinWeightedSpheroidalEigenvalue] = {Method -> Automatic};


SetAttributes[SpinWeightedSpheroidalEigenvalue, {NumericFunction, Listable, NHoldAll}];


SpinWeightedSpheroidalEigenvalue[s_?NumericQ, l_?NumericQ, m_?NumericQ, \[Gamma]_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] :=
 (Message[SpinWeightedSpheroidalEigenvalue::params, s, l, m]; $Failed);


SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; \[Gamma] == 0 :=
  l(l+1) - s(s+1);


SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_?InexactNumberQ, OptionsPattern[]] := 
 Module[{opts, Aini, \[Lambda]},
  Switch[OptionValue[Method],
    "SphericalExpansion",
      \[Lambda] = SWSHEigenvalueSpectral[s, l, m, \[Gamma]] - 2 m \[Gamma] + \[Gamma]^2,
    {"SphericalExpansion", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHEigenvalueSpectral]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      ];
      \[Lambda] = SWSHEigenvalueSpectral[s, l, m, \[Gamma], opts] - 2 m \[Gamma] + \[Gamma]^2,
    Automatic | "Leaver",
      \[Lambda] = SWSHEigenvalueLeaver[s, l, m, \[Gamma]] - 2 m \[Gamma] + \[Gamma]^2;,
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHEigenvalueLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      ];
      \[Lambda] = SWSHEigenvalueLeaver[s, l, m, \[Gamma], opts] - 2 m \[Gamma] + \[Gamma]^2,
      _,
      Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      \[Lambda] = $Failed;
  ];
  \[Lambda]
];


SpinWeightedSpheroidalEigenvalue /: N[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_?NumericQ, opts:OptionsPattern[]], Nopts___] :=
  SpinWeightedSpheroidalEigenvalue[s, l, m, N[\[Gamma], Nopts], opts];


(* ::Subsection::Closed:: *)
(*Small-\[Gamma] expansion*)


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


(* ::Subsection::Closed:: *)
(*Large-\[Gamma] expansion*)


(* ::Text:: *)
(*For details see "High-order asymptotics for the Spin-Weighted Spheroidal Equation at large real frequency", M. Casals, A. C. Ottewill, N. Warburton, arXiv:1810.00432*)


SWSHq[s_, l_, m_] := Module[{slm,z0},
	slm = Abs[m+Abs[s]] + Abs[s];
	z0  = If[EvenQ[l+m], 0, 1];

	If[l  >= slm, l+1-z0,  2l+1 - slm ]
]


SpinWeightedSpheroidalEigenvalue /:
 Series[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_], {\[Gamma]_, \[Infinity], order_}] /; !AllTrue[{2s, 2l, 2m}, IntegerQ] :=
 (Message[SpinWeightedSpheroidalEigenvalue::infseriesparams, s, l, m]; $Failed);


SpinWeightedSpheroidalEigenvalue /:
 Series[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, -\[Gamma]_], {\[Gamma]_, \[Infinity], order_}] :=
   Series[SpinWeightedSpheroidalEigenvalue[s, l, -m, \[Gamma]], {\[Gamma], \[Infinity], order}];


SpinWeightedSpheroidalEigenvalue /:
 Series[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_], {\[Gamma]_, \[Infinity], order_}] :=
Module[{slm,z0,q,aFgen,AFgen,Asgen,\[Delta]gen,\[Nu]gen,RecRelgen,n,c,p,Serngen,Asum,sElm,s\[Lambda]lm,res},
  aFgen[0,_]=0;
  aFgen[0,0]=1;

  AFgen[0]=1;
  AFgen[r_]=AFgen[0]Sum[aFgen[r,k]c^-k,{k,Abs[r], order+1}];
  
  Asgen=Sum[\[Delta]gen[k]c^-k, {k, 1, order+1}];
  \[Nu]gen=1/2 (-1+SWSHq[s,l,m]-s-Abs[m+s]);
  RecRelgen=1/2  (1+n+\[Nu]gen+Abs[m+s]) (2(1+n+\[Nu]gen+s)-Abs[m-s]+Abs[m+s])AFgen[n+1]+(-1-m-m^2-2 (n+ \[Nu]gen)+4 c n-2 (n+\[Nu]gen)^2-s-m s-2 (n+\[Nu]gen) s-4 c Asgen-(1+2 n+2 \[Nu]gen+s) Abs[m+s])AFgen[n] + 1/2  (n+\[Nu]gen) (2 (n+\[Nu]gen+s)+Abs[m-s]+Abs[m+s])AFgen[n-1];

  Serngen[nn_]:=Serngen[nn]=Collect[CoefficientList[Normal[Series[RecRelgen/.n->nn,{c, \[Infinity], order+1}]],c^-1],aFgen[_,_]];

  Do[
    \[Delta]gen[p]=\[Delta]gen[p]/.Solve[Serngen[0][[p]]==0,\[Delta]gen[p]][[1]];
    Do[If[k!=0,aFgen[k,p]=aFgen[k,p]/.Solve[Serngen[k][[p]]==0,aFgen[k,p]][[1]]],{k, -p, p}]
  ,{p, 1, order+1}];
  
  Asum = 4\[Gamma] Sum[\[Delta]gen[k]\[Gamma]^-k,{k, 2, order+1}];

  sElm = -\[Gamma]^2+ 2 SWSHq[s,l,m] \[Gamma]-1/2 (SWSHq[s,l,m]^2-m^2-2s^2+1)+Asum;

  s\[Lambda]lm = sElm-s(s+1)+\[Gamma]^2-2m \[Gamma];
  
  res = SeriesData[\[Gamma], \[Infinity], Reverse[CoefficientList[Expand[s\[Lambda]lm \[Gamma]^order],\[Gamma]]], -1, order+1, 1];

  Remove[aFgen, AFgen, \[Delta]gen, Serngen, c, n];

  res
]


(* ::Section::Closed:: *)
(*SpinWeightedSpheroidalHarmonicS*)


(* ::Subsection::Closed:: *)
(*Spherical expansion method*)


Options[SWSHSSpectral] = {"NumTerms" -> Automatic};


SWSHSSpectral[s_, l_, m_, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{lmin, nUp, nDown, A, esys,evec,eval,sign,pos},
  (* FIXME: Improve the estimate of nmax. It should depend on the accuarcy sought. *)
  Switch[OptionValue["NumTerms"],
   Automatic,
    nUp = Ceiling[Abs[3/2\[Gamma]]]+5;
    Message[SpinWeightedSpheroidalHarmonicS::numterms, nUp];,
   _Integer,
    nUp = OptionValue["NumTerms"];,
   _,
    Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> {"SphericalExpansion", opts}];
    Return[$Failed];
  ];
  lmin = Max[Abs[s],Abs[m]];
  nDown = Min[l-lmin,nUp];

  A = SparseArray[
        {{i_,i_} :> kHat[s, l-nDown-1+i, m, \[Gamma]],
         {i_,j_} /; j-i==-2 :> k2[s, l-nDown-3+i, m, \[Gamma]],
         {i_,j_} /; j-i==-1 :> kTilde2[s, l-nDown+i-2, m, \[Gamma]],
         {i_,j_} /; j-i==1 :> kTilde2[s, l-nDown+i-1, m, \[Gamma]],
         {i_,j_} /; j-i==2 :> k2[s, l-nDown+i+-1, m, \[Gamma]]},
        {nUp+nDown+1, nUp+nDown+1}];
  esys = Quiet[Eigensystem[A], {Eigenvalues::arh,Eigenvalues::arhm}];
  eval = -Sort[esys[[1]]][[-(nDown+1)]];
  pos  = Position[esys[[1]], -eval][[1]];
  evec = First[esys[[2,pos]]];

  sign=Sign[evec[[Min[l-lmin+1,(nUp+nDown)/2+1]]]];
  SpinWeightedSpheroidalHarmonicSFunction[<|"s" -> s, "l" -> l, "m" -> m, "\[Gamma]"-> \[Gamma], "Eigenvalue" -> eval-s(s+1) - 2 m \[Gamma] + \[Gamma]^2, "Method" -> {"SphericalExpansion", "Coefficients" -> sign*evec, "\!\(\*SubscriptBox[\(l\), \(min\)]\)" -> l-nDown, "\!\(\*SubscriptBox[\(l\), \(max\)]\)" -> l+nUp}|>]
];


(* ::Subsection::Closed:: *)
(*Leaver's method*)


Options[SWSHSLeaver] = {"NumTerms" -> Automatic};


SWSHSLeaver[s_, l_, m_, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], k1, k2, \[Alpha]n, \[Beta]n, \[Gamma]n, an, n, sign, norm, anTab, nmin, nmax, normterm, prec=Precision[\[Gamma]]},
  nmin = 0;
  nmax = 1;
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]];

  (* Leaver's recurrence formula, Eq. 20 of Leaver 1985 *)
  k1 = Abs[m-s]/2;
  k2 = Abs[m+s]/2;
  \[Alpha]n[n_] := -2(n+1)(n+2k1+1);
  \[Beta]n[n_] := n(n-1)+2n(k1+k2+1-2\[Gamma])-(2\[Gamma](2k1+s+1)-(k1+k2)(k1+k2+1))-(s(s+1)+\[Lambda]+2m \[Gamma]);
  \[Gamma]n[n_] := 2\[Gamma](n+k1+k2+s);
  an[nmin] = 1;
  an[nmin+1] = -\[Beta]n[nmin]/\[Alpha]n[nmin];
  an[n_] := an[n] = (-\[Beta]n[n-1]an[n-1]-\[Gamma]n[n-1]an[n-2])/\[Alpha]n[n-1];

  (* Compute coefficients until we reach the desired tolerance for computing the norm *)
  normterm[i_] := 2^i Pochhammer[i+2 (1+k1+k2), -2k2-1] Hypergeometric1F1[1+i+2 k1, i+2 (1+k1+k2), 4 \[Gamma]] Sum[an[j]an[i-nmin-j], {j, 0, i-nmin}];
  Switch[OptionValue["NumTerms"],
   Automatic,
    norm = normterm[nmin];
    While[norm != (norm += normterm[nmax]), nmax++];,
   _Integer,
    nmax = OptionValue["NumTerms"] - 1;
    norm = Sum[normterm[i], {i, nmin, nmax}];,
   _,
    Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> {"Leaver", opts}];
    Return[$Failed];
  ];
  
  anTab = Table[an[n], {n, nmin, nmax}];

  (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)
  norm = Sqrt[2\[Pi]] (2^(1+2 k1+2 k2) E^(-2 \[Gamma]) Gamma[1+2 k2] norm)^(1/2);
  If[Precision[norm] == 0.,
    Message[SpinWeightedSpheroidalHarmonicS::normprec];
    norm = 1;
  ];

  (* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
  sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1];
  
  (* Clear Symbols with DownValues to prevent memory leaks *)
  Clear[\[Alpha]n, \[Beta]n, \[Gamma]n, an, normterm];

  (* Return a SpinWeightedSpheroidalHarmonicSFunction which can be evaluated for arbitratry \[Theta], \[Phi] *)
  SpinWeightedSpheroidalHarmonicSFunction[<|"s" -> s, "l" -> l, "m" -> m, "\[Gamma]"-> \[Gamma], "Eigenvalue" -> \[Lambda], "Method" -> {"Leaver", "Coefficients" -> sign anTab/norm, "\!\(\*SubscriptBox[\(n\), \(max\)]\)" -> nmax}|>]
];


(* ::Subsection::Closed:: *)
(*SpinWeightedSpheroidalHarmonicS*)


SyntaxInformation[SpinWeightedSpheroidalHarmonicS] =
 {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};


Options[SpinWeightedSpheroidalHarmonicS] = {Method -> Automatic};


SetAttributes[SpinWeightedSpheroidalHarmonicS, {Listable, NHoldAll}];


SpinWeightedSpheroidalHarmonicS[s_?NumericQ, l_?NumericQ, m_?NumericQ, \[Gamma]_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] := 
 (Message[SpinWeightedSpheroidalHarmonicS::params, s, l, m]; $Failed);


SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; \[Gamma] == 0 :=
  SpinWeightedSpheroidalHarmonicSFunction[<|"s" -> s, "l" -> l, "m" -> m, "\[Gamma]"-> \[Gamma], "Eigenvalue" -> SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]], "Method" -> {"SphericalExact"}|>] /; OptionValue[Method] == Automatic;


SpinWeightedSpheroidalHarmonicS[s_?NumericQ, l_?NumericQ, m_?NumericQ, \[Gamma]:(_?InexactNumberQ|0), OptionsPattern[]] :=
 Module[{opts, Slm},
  Switch[OptionValue[Method],
    "Eigenvalue"|"SphericalExpansion",
      Slm = SWSHSSpectral[s, l, m, \[Gamma]],
    {"Eigenvalue"|"SphericalExpansion", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHSSpectral]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
      ];
      Slm = SWSHSSpectral[s, l, m, \[Gamma], opts],
    Automatic | "Leaver",
      Slm = SWSHSLeaver[s, l, m, \[Gamma]];,
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHSLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
      ];
      Slm = SWSHSLeaver[s, l, m, \[Gamma], opts],
    _,
     Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
     Slm = $Failed;
  ];
  Slm
];


(* ::Subsection::Closed:: *)
(*Small-\[Gamma] expansion*)


SpinWeightedSpheroidalHarmonicS /: 
  Series[SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_][\[Theta]_, \[Phi]_], {\[Gamma]_, 0, order_}] :=
 Module[{i, j, coeffs}, Internal`InheritedBlock[{d, s\[Lambda]lm}, Block[{\[Alpha], \[Beta]},
  Do[
    Do[
      d[s, l, m][i, j] = simplify[d[s, l, m][i, j]], {j, -i, i}]; 
    s\[Lambda]lm[s, l, m][i] = simplify[s\[Lambda]lm[s, l, m][i]];
  , {i, 0, order}];
  coeffs = Table[Sum[d[s, l, m][i, j] If[TrueQ[l+j < Abs[s] || l+j < Abs[m]], 0, SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]]], {j, -i, i}], {i, 0, order}];
  SeriesData[\[Gamma], 0, coeffs, 0, order + 1, 1]
]]];


(* ::Section::Closed:: *)
(*SpinWeightedSpheroidalHarmonicSFunction*)


SyntaxInformation[SpinWeightedSpheroidalHarmonicSFunction] =
 {"ArgumentsPattern" -> {_}};


SetAttributes[SpinWeightedSpheroidalHarmonicSFunction, {NHoldAll}];


(* ::Subsection::Closed:: *)
(*Output format*)


SpinWeightedSpheroidalHarmonicSFunction /:
 MakeBoxes[swshf:SpinWeightedSpheroidalHarmonicSFunction[a_Association], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended, s = a["s"], l = a["l"], m = a["m"], \[Gamma] = a["\[Gamma]"], \[Lambda] = a["Eigenvalue"], method = First[a["Method"]], opts = Association[Rest[a["Method"]]]},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", s}], "  ",
                  BoxForm`SummaryItem[{"l: ", l}], "  ",
                  BoxForm`SummaryItem[{"m: ", m}], "  ",
                  BoxForm`SummaryItem[{"\[Gamma]: ", \[Gamma]}]}],
                  BoxForm`SummaryItem[{"Eigenvalue: ", \[Lambda]}]
             };
  extended = Join[{BoxForm`SummaryItem[{"Method: ", method}]},
                   Switch[method,
                    "SphericalExpansion",
                     {BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(l\), \(min\)]\): ", opts["\!\(\*SubscriptBox[\(l\), \(min\)]\)"]}],
                      BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(l\), \(max\)]\): ", opts["\!\(\*SubscriptBox[\(l\), \(max\)]\)"]}],
                      BoxForm`SummaryItem[{"Expansion coefficients: ", Column[opts["Coefficients"]]}]},
                    "Leaver",
                     {BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(n\), \(max\)]\): ", opts["\!\(\*SubscriptBox[\(n\), \(max\)]\)"]}],
                      BoxForm`SummaryItem[{"Expansion coefficients: ", Column[opts["Coefficients"]]}]},
                    _,
                     {}
                    ]];
  BoxForm`ArrangeSummaryBox[
    SpinWeightedSpheroidalHarmonicSFunction,
    swshf,
    PolarPlot[{swshf[\[Theta],0],swshf[\[Theta],\[Pi]]}, {\[Theta],0,\[Pi]}, Axes->False, PlotRange->All, ImagePadding->All, PlotStyle->ColorData[97,1],
      ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
    summary,
    extended,
    form]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


SpinWeightedSpheroidalHarmonicSFunction[assoc_][key_String] /; KeyExistsQ[assoc, key] := assoc[key];


Keys[SpinWeightedSpheroidalHarmonicSFunction[assoc_]] ^:= Join[Keys[assoc], {"ExpansionCoefficients"}];


(* ::Subsection::Closed:: *)
(*Expansion coefficients*)


SpinWeightedSpheroidalHarmonicSFunction[assoc_]["ExpansionCoefficients"] /; First[assoc["Method"]] == "Leaver" :=
 Module[{opts = Association[Rest[assoc["Method"]]], coeffs, nmax},
  {nmax, coeffs} = {opts["\!\(\*SubscriptBox[\(n\), \(max\)]\)"], opts["Coefficients"]};
  Association[Thread[Range[0, nmax] -> coeffs]]
 ];


SpinWeightedSpheroidalHarmonicSFunction[assoc_]["ExpansionCoefficients"] /; First[assoc["Method"]] == "SphericalExpansion" :=
 Module[{opts = Association[Rest[assoc["Method"]]], lmin, lmax, coeffs},
  {lmin, lmax, coeffs} = {opts["\!\(\*SubscriptBox[\(l\), \(min\)]\)"], opts["\!\(\*SubscriptBox[\(l\), \(max\)]\)"], opts["Coefficients"]};
  Association[Thread[Range[lmin, lmax] -> coeffs]]
 ];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


(* ::Subsubsection::Closed:: *)
(*SphericalExpansion method*)


SpinWeightedSpheroidalHarmonicSFunction[assoc_][\[Theta]_?NumericQ, \[Phi]_?NumericQ] /; First[assoc["Method"]] == "SphericalExpansion" :=
 Module[{s = assoc["s"], m = assoc["m"], opts = Association[Rest[assoc["Method"]]], lmin, lmax, coeffs},
  {lmin, lmax, coeffs} = {opts["\!\(\*SubscriptBox[\(l\), \(min\)]\)"], opts["\!\(\*SubscriptBox[\(l\), \(max\)]\)"], opts["Coefficients"]};
  Sum[coeffs[[l-lmin+1]]SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0], {l, lmin, lmax}]Exp[I m \[Phi]]
];


(* ::Subsubsection::Closed:: *)
(*Leaver's method*)


SpinWeightedSpheroidalHarmonicSFunction[assoc_][\[Theta]_?NumericQ, \[Phi]_?NumericQ] /; First[assoc["Method"]] == "Leaver" :=
 Module[{u = Cos[\[Theta]], k1, k2, oneplusu, oneminusu, res, s = assoc["s"], l = assoc["l"], m = assoc["m"], \[Gamma] = assoc["\[Gamma]"], opts = Association[Rest[assoc["Method"]]], coeffs, nUp},
  {nUp, coeffs} = {opts["\!\(\*SubscriptBox[\(n\), \(max\)]\)"], opts["Coefficients"]};
  k1 = Abs[m-s]/2;
  k2 = Abs[m+s]/2;
  oneplusu = 2 Cos[\[Theta]/2]^2;
  oneminusu = 2 Sin[\[Theta]/2]^2;
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *)
  res = E^(\[Gamma] u) If[k1==0, 1, oneplusu^k1] If[k2==0, 1, oneminusu^k2] coeffs . Join[{1},oneplusu^Range[1,nUp]] Exp[I m \[Phi]];

  (* Print a warning message if the harmonic cannot be determined *)
  If[Precision[res] == 0. && Accuracy[res] < Accuracy[{\[Gamma], \[Theta], \[Phi]}], Message[SpinWeightedSpheroidalHarmonicS::prec]];

  res
];


(* ::Subsubsection::Closed:: *)
(*SphericalExact method*)


SpinWeightedSpheroidalHarmonicSFunction[assoc_][\[Theta]_, \[Phi]_] /; First[assoc["Method"]] == "SphericalExact" :=
  SpinWeightedSphericalHarmonicY[assoc["s"], assoc["l"], assoc["m"], \[Theta], \[Phi]];


(* ::Subsubsection::Closed:: *)
(*Derivatives*)


Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[assoc_]][\[Theta]_?NumericQ, \[Phi]_?NumericQ] /; First[assoc["Method"]] == "SphericalExpansion" :=
 Module[{\[Theta]1,\[Phi]1, res, s = assoc["s"], m = assoc["m"],  opts = Association[Rest[assoc["Method"]]], lmin, lmax, coeffs},
  {lmin, lmax, coeffs} = {opts["\!\(\*SubscriptBox[\(l\), \(min\)]\)"], opts["\!\(\*SubscriptBox[\(l\), \(max\)]\)"], opts["Coefficients"]};
  res = Sum[coeffs[[l-lmin+1]]D[SpinWeightedSphericalHarmonicY[s, l, m, \[Theta]1, 0],{\[Theta]1,d1}], {l, lmin, lmax}]D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]};
  Remove[\[Theta]1,\[Phi]1];
  res
];


Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[assoc_]][\[Theta]_?NumericQ, \[Phi]_?NumericQ] /; First[assoc["Method"]] == "Leaver":=
 Module[{\[Theta]1, \[Phi]1, u, k1, k2, oneplusu, oneminusu, res, s = assoc["s"], m = assoc["m"], \[Gamma] = assoc["\[Gamma]"], opts = Association[Rest[assoc["Method"]]], coeffs, nUp},
  {nUp, coeffs} = {opts["\!\(\*SubscriptBox[\(n\), \(max\)]\)"], opts["Coefficients"]};
  k1 = Abs[m-s]/2;
  k2 = Abs[m+s]/2;
  u=Cos[\[Theta]1];
  oneplusu = 2 Cos[\[Theta]1/2]^2;
  oneminusu = 2 Sin[\[Theta]1/2]^2;
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *)
  res = D[E^(\[Gamma] u) If[k1==0, 1, oneplusu^k1] If[k2==0, 1, oneminusu^k2] coeffs . oneplusu^Range[0,nUp],{\[Theta]1,d1}] D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]};
  Remove[\[Theta]1,\[Phi]1];
  res
];
  


(* ::Section::Closed:: *)
(*SpinWeightedSphericalHarmonicY*)


SyntaxInformation[SpinWeightedSphericalHarmonicY] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


SetAttributes[SpinWeightedSphericalHarmonicY, {NumericFunction, Listable}];


SpinWeightedSphericalHarmonicY[s_?NumericQ, l_?NumericQ, m_?NumericQ, \[Theta]_, \[Phi]_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] := 
 (Message[SpinWeightedSphericalHarmonicY::params, s, l, m]; $Failed);


SpinWeightedSphericalHarmonicY[s_?NumericQ, l_?NumericQ, m_?NumericQ, \[Theta]_, \[Phi]_] :=
  (-1)^m Sqrt[((l+m)!(l-m)!(2l+1))/(4\[Pi] (l+s)!(l-s)!)] Sum[Binomial[l-s,r] Binomial[l+s,r+s-m] (-1)^(l-r-s) If[(2 l - 2 r - s + m)===0, 1, Sin[\[Theta]/2]^(2 l - 2 r - s + m)] If[(2 r + s - m)===0, 1, Cos[\[Theta]/2]^(2 r + s - m)],{r,Max[m-s,0],Min[l-s,l+m]}]Exp[I m \[Phi]];


FunctionExpand[SpinWeightedSphericalHarmonicY[s_, l_, m_, \[Theta]_, \[Phi]_]] ^:=
  ((-1)^m*E^(I*m*\[Phi])*Sqrt[((1 + 2*l)*(l - m)!*(l + m)!)/((l - s)!*(l + s)!)]*
  Piecewise[{{(-1)^(l - s)*Binomial[l + s, -m + s]*Cos[\[Theta]/2]^(-m + s)*
      Hypergeometric2F1[-l - m, -l + s, 1 - m + s, -Cot[\[Theta]/2]^2]*
      Sin[\[Theta]/2]^(2*l + m - s), m <= s}}, (-1)^(l - m)*Binomial[l - s, m - s]*
    Cos[\[Theta]/2]^(m - s)*Hypergeometric2F1[-l + m, -l - s, 1 + m - s, -Cot[\[Theta]/2]^2]*
    Sin[\[Theta]/2]^(2*l - m + s)])/(2*Sqrt[Pi]);


SpinWeightedSphericalHarmonicY[s_, l_, m_, \[Theta]_, 0.] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0];


SpinWeightedSphericalHarmonicY[0, l_, m_, \[Theta]_, \[Phi]_] :=
  SphericalHarmonicY[l, m, \[Theta], \[Phi]];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{SpinWeightedSphericalHarmonicY, SpinWeightedSpheroidalHarmonicS, SpinWeightedSpheroidalHarmonicSFunction, SpinWeightedSpheroidalEigenvalue}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[]
EndPackage[];
