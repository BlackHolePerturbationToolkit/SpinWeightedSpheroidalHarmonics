(* ::Package:: *)

(* ::Title:: *)
(*SpinWeightedSpheoridalHarmonics package*)


(* ::Section::Closed:: *)
(*Public function definitions and preamble*)


BeginPackage["SpinWeightedSpheroidalHarmonics`"];

SpinWeightedSphericalHarmonicY::usage = "SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]] gives the spin-weighted spherical harmonic with of spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicS::usage = "SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] gives the spin-weighted oblate spheroidal harmonic with spheroidicity \[Gamma], spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicSFunction::usage = "SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], data, method] represents a function for computing numerical values of spin-weighted spheroidal harmonics.";
SpinWeightedSpheroidalEigenvalue::usage = "SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] gives the spin-weighted oblate spheroidal eigenvalue with spheroidicity \[Gamma], spin-weight s, degree l and order m.";

$SpinWeightedSpheroidalHarmonicsInformation::usage = "$SpinWeightedSpheroidalHarmonicsInformation is a list of rules that gives information about the version of the SpinWeightedSpheroidalHarmonics package you are running.";
$SpinWeightedSpheroidalHarmonicsInstallationDirectory::usage = "$SpinWeightedSpheroidalHarmonicsInstallationDirectory gives the top-level directory in which the SpinWeightedSpheroidalHarmonics package is installed.";

$SpinWeightedSpheroidalHarmonicsVersionNumber::usage = "$SpinWeightedSpheroidalHarmonicsVersionNumber is a real number which gives the current version number for the SpinWeightedSpheroidalHarmonics package.";
$SpinWeightedSpheroidalHarmonicsReleaseNumber::usage = "$SpinWeightedSpheroidalHarmonicsReleaseNumber is an integer which gives the current release number for the SpinWeightedSpheroidalHarmonics package.";
$SpinWeightedSpheroidalHarmonicsVersion::usage = "$SpinWeightedSpheroidalHarmonicsVersionNumber is a string that gives the version of the SpinWeightedSpheroidalHarmonics package you are running.";

Begin["`Private`"];

(**********************************************************)
(* Package version information                            *)
(**********************************************************)

$SpinWeightedSpheroidalHarmonicsInstallationDirectory = FileNameDrop[FindFile["SpinWeightedSpheroidalHarmonics`"], -2];

$SpinWeightedSpheroidalHarmonicsVersionNumber        = 1.0;
$SpinWeightedSpheroidalHarmonicsReleaseNumber        = 0;

$SpinWeightedSpheroidalHarmonicsVersion :=
 Module[{path, version, release, buildid, gitrev, gitdir},
  path = $SpinWeightedSpheroidalHarmonicsInstallationDirectory;
  version = ToString[NumberForm[$SpinWeightedSpheroidalHarmonicsVersionNumber, {Infinity, 1}]];
  release = ToString[$SpinWeightedSpheroidalHarmonicsReleaseNumber];

  buildid = Quiet@ReadList[FileNameJoin[{path, "BUILD_ID"}], "String"];
  If[SameQ[buildid, $Failed],
    buildid = "";
  ,
    buildid = " (" <> First[buildid] <> ")";
  ];

  (* First, check for a GIT_REVISION file. If it exists, use its contents as the revision. *)
  gitrev = Quiet@ReadList[FileNameJoin[{path, "GIT_REVISION"}],"String"];

  (* Otherwise, try to determine the git revision directly *)
  If[SameQ[gitrev, $Failed],
    gitdir = FileNameJoin[{path, ".git"}];
    If[FileType[gitdir] === Directory,
      gitrev = Quiet@ReadList["!git --git-dir "<>gitdir<>" rev-parse HEAD", String];
      If[gitrev === {}, gitrev = $Failed];
    ];
  ];

  (* If it worked, ReadList returns a list but we just want the first element (line) *)
  If[Head[gitrev] === List, gitrev = First[gitrev]];

  (* Check we have a git revision and otherwise give up trying *)
  If[Head[gitrev] === String && StringMatchQ[gitrev, RegularExpression["[0-9a-f]{5,40}"]], gitrev = " (" <> gitrev <> ")", gitrev = ""];

  version <> "." <> release <> buildid <> gitrev
]

$SpinWeightedSpheroidalHarmonicsInformation :=
  {"InstallationDirectory" -> $SpinWeightedSpheroidalHarmonicsInstallationDirectory,
   "Version" -> $SpinWeightedSpheroidalHarmonicsVersion,
   "VersionNumber" -> $SpinWeightedSpheroidalHarmonicsVersionNumber,
   "ReleaseNumber" -> $SpinWeightedSpheroidalHarmonicsReleaseNumber}



(* ::Section::Closed:: *)
(*Useful functions*)


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

(* Functions for spectral method for seed in numerical evaluation *)
kHat[s_, 0, 0, \[Gamma]_] := \[Gamma]^2/3;
kHat[s_, l_, m_, \[Gamma]_] := -l (1+l)+(2 m s^2 \[Gamma])/(l+l^2)+1/3 (1+(2 (l+l^2-3 m^2) (l+l^2-3 s^2))/(l (-3+l+8 l^2+4 l^3))) \[Gamma]^2;
k2[s_, l_, m_, \[Gamma]_] := (Sqrt[((1+l-m) (2+l-m) (1+l+m) (2+l+m) (1+l-s) (2+l-s) (1+l+s) (2+l+s))/((1+2 l) (5+2 l))] \[Gamma]^2)/((1+l) (2+l) (3+2 l));
kTilde2[s_, 0, 0, \[Gamma]_] := -2 s Sqrt[(1-s^2)/3] \[Gamma];
kTilde2[s_, l_, m_, \[Gamma]_] := -((2 s Sqrt[((1+2 l+l^2-m^2) (1+2 l+l^2-s^2))/(3+8 l+4 l^2)] \[Gamma] (2 l+l^2+m \[Gamma]))/(l (2+3 l+l^2)));


(*Module for computing a contined fraction. This is preferred over Mathematica's ContinedFractionK function as we can get an error estimate on the result using this function*)
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
   res
   ];




(* ::Section::Closed:: *)
(*SpinWeightedSpheroidalEigenvalue*)


(* ::Subsection::Closed:: *)
(*Spectral and Leaver's method*)


SWSHEigenvalueSpectral[s_,l_,m_,\[Gamma]_]:=
 Module[{nmin,nmax,Matrix,Eigens,lmin},
  nmax=Ceiling[Abs[1.5\[Gamma]-\[Gamma]^2/250]]+5;(*FIXME: Improve the estimate of nmax*)
  lmin=Max[Abs[s],Abs[m]];
  nmin=Min[nmax,l-lmin];

  Matrix=SparseArray[
	{{i_,i_}:>kHat[s,l-nmin-1+i,m,\[Gamma]],
	{i_,j_}/;j-i==-2:>k2[s,l-nmin-3+i,m,\[Gamma]],
	{i_,j_}/;j-i==-1:>kTilde2[s,l-nmin+i-2,m,\[Gamma]],
	{i_,j_}/;j-i==1:>kTilde2[s,l-nmin+i-1,m,\[Gamma]],
	{i_,j_}/;j-i==2:>k2[s,l-nmin+i+-1,m,\[Gamma]]}
  ,{nmax+nmin+1,nmax+nmin+1}];
  Eigens=-Sort[Eigenvalues[Matrix]];

  Eigens[[-(nmin+1)]]-s(s+1)
];

SWSHEigenvalueLeaver[s_, l_, m_, \[Gamma]_, Aini_] :=
 Module[{Myprec, Nmax, nInv, \[Alpha], \[Beta], \[Alpha]n, \[Beta]n, \[Gamma]n, n, LHS, RHS, Eq, A, Aval, Avar},
  Myprec = Max[Precision[\[Gamma]], MachinePrecision];
  nInv = l-Abs[m];
  \[Alpha] = Abs[m+s];
  \[Beta] = Abs[m-s];
  \[Alpha]n[n_] := (-4\[Gamma](n+\[Alpha]+1)(n+\[Beta]+1)(n+(\[Alpha]+\[Beta])/2+1+s ))/((2n+\[Alpha]+\[Beta]+2)(2n+\[Alpha]+\[Beta]+3));
  \[Beta]n[n_, A_] := A + s(s+1)+\[Gamma]^2 -(n +(\[Alpha]+\[Beta])/2)(n +(\[Alpha]+\[Beta])/2+1) + If[s!=0,(8m s^2 \[Gamma])/((2n+\[Alpha]+\[Beta])(2n+\[Alpha]+\[Beta]+2)),0];
  \[Gamma]n[n_] := (4\[Gamma] n(n+\[Alpha]+\[Beta])(n+(\[Alpha]+\[Beta])/2-s))/((2n+\[Alpha]+\[Beta]-1)(2n+\[Alpha]+\[Beta]));
  RHS[Ax_] := -CF[-\[Alpha]n[n-1] \[Gamma]n[n], \[Beta]n[n,Ax], {n, nInv+1}];
  LHS[Ax_] := \[Beta]n[nInv, Ax] + ContinuedFractionK[-\[Alpha]n[nInv-n] \[Gamma]n[nInv-n+1], \[Beta]n[nInv-n, Ax], {n, 1, nInv}];
  Eq[A_?NumericQ] := LHS[A] - RHS[A];
  Aval = Avar /. FindRoot[Eq[Avar]==0, {Avar, Aini}, AccuracyGoal -> Myprec-3, WorkingPrecision -> Myprec, Method -> "Secant"];
  Aval
]



(* ::Subsection::Closed:: *)
(*SpinWeightedSpheroidalEigenvalue (uses a combination of the spectral and Leaver's method)*)


(**********************************************************)
(* SpinWeightedSpheroidalEigenvalue                       *)
(**********************************************************)

(*We only compute the eigenvalues/eigenfunctions when (s,l,m) are either all integers or all half-integers*)
ValidParamQ[s_Integer,l_Integer,m_Integer]:=True
ValidParamQ[s_,l_,m_]:=Mod[{s,l,m},1]=={1/2,1/2,1/2}

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






(* ::Subsection::Closed:: *)
(*Small-\[Gamma] expansion*)


(*Small \[Gamma] expansions*)
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


(* ::Subsection:: *)
(*Large-\[Gamma] expansion*)


(*For details see "High-order asymptotics for the Spin-Weighted Spheroidal Equation at large real frequency", M. Casals, A. C. Ottewill, N. Warburton, arXiv:1810.00432*)

(*Large \[Gamma] expansions*) 
SWSHq[s_?NumberQ, l_?NumberQ, m_?NumberQ] := Module[{slm,z0},
	slm = Abs[m+Abs[s]] + Abs[s];
	z0  = If[EvenQ[l+m], 0, 1];

	If[l  >= slm, l+1-z0,  2l+1 - slm ]
]

(*\!\(\(Teukolsky - Starobinsky\ identities\ 
\*SubscriptBox[\(\(imply\)\(\ \)\), \(s\)]
\*SubscriptBox[\(\[Lambda]\), \(lm\((\(-c\))\)\)]\)
\*SubscriptBox[\(=\), \(s\)]
\*SubscriptBox[\(\[Lambda]\), \(l\((\(-m\))\)c\)]\)*)

SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_?NumberQ, l_?NumberQ, m_?NumberQ, -\[Gamma]_], {\[Gamma]_, \[Infinity], order_}]:= Series[SpinWeightedSpheroidalEigenvalue[s, l, -m, \[Gamma]], {\[Gamma], \[Infinity], order}]

SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_?NumberQ, l_?NumberQ, m_?NumberQ, \[Gamma]_], {\[Gamma]_, \[Infinity], order_}] :=
Module[{slm,z0,q,aFgen,AFgen,Asgen,\[Delta]gen,\[Nu]gen,RecRelgen,n,c,p,Serngen,Asum,sElm,s\[Lambda]lm},
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
  
  SeriesData[\[Gamma], \[Infinity], Reverse[CoefficientList[Expand[s\[Lambda]lm \[Gamma]^order],\[Gamma]]], -1, order+1, 1]
  
]



(* ::Section:: *)
(*SpinWeightedSpheroidalHarmonicS*)


(**********************************************************)
(* SpinWeightedSpheroidalHarmonicS                        *)
(**********************************************************)

SyntaxInformation[SpinWeightedSpheroidalHarmonicS] =
 {"ArgumentsPattern" -> {_, _, _, _, ___}};

Options[SpinWeightedSpheroidalHarmonicS] = {Method -> Automatic};

SetAttributes[SpinWeightedSpheroidalHarmonicS, {NumericFunction, Listable}];

SpinWeightedSpheroidalHarmonicS[s_, l_, m_, (0|0.), \[Theta]_, \[Phi]_, opts___] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]];
  
SpinWeightedSpheroidalHarmonicS[s_, l_, m_, 0] :=
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, 0, Method->"SphericalExact"]

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?MachineNumberQ] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma],
    Method -> (OptionValue[SpinWeightedSpheroidalHarmonicS, "Method"] /. Automatic -> "Eigenvalue")];

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma],
    Method -> (OptionValue[SpinWeightedSpheroidalHarmonicS, "Method"] /. Automatic -> "Leaver")];

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ, opts___:OptionsPattern[]] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], opts][\[Theta], \[Phi]];

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, Method->"Eigenvalue"] :=
 Module[{lmin, nmax, nmin, A, esys,evec,eval,sign,pos},
  lmin = Max[Abs[s],Abs[m]];
  nmax = Ceiling[Abs[3/2\[Gamma]-\[Gamma]^2/250]]+50;(*FIXME: Improve the estimate of nmax*)
  If[EvenQ[nmax],nmax+=1];
  nmin = Min[l-lmin,nmax];

  A = -SparseArray[
        {{i_,i_} :> kHat[s, l-nmin-1+i, m, \[Gamma]],
         {i_,j_} /; j-i==-2 :> k2[s, l-nmin-3+i, m, \[Gamma]],
         {i_,j_} /; j-i==-1 :> kTilde2[s, l-nmin+i-2, m, \[Gamma]],
         {i_,j_} /; j-i==1 :> kTilde2[s, l-nmin+i-1, m, \[Gamma]],
         {i_,j_} /; j-i==2 :> k2[s, l-nmin+i+-1, m, \[Gamma]]},
        {nmax+nmin+1, nmax+nmin+1}];
  esys = Eigensystem[A,Method->If[Precision[\[Gamma]]==MachinePrecision,"Banded","Automatic"]];
  eval = Sort[esys[[1]],Greater][[-(nmin+1)]];
  pos  = Position[esys[[1]], eval][[1]];
  evec = First[esys[[2,pos]]];

  sign=Sign[evec[[Min[l-lmin+1,(nmax+nmin)/2+1]]]];
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], {sign*evec, nmin, nmax}, Method -> "Eigenvalue"]
];

SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, Method -> "Leaver"] :=
 Module[{\[Lambda], k1, k2, \[Alpha]n, \[Beta]n, \[Gamma]n, an, n, norm, anTab, nmin, nmax, maxcoeff, prec=Precision[\[Gamma]]},
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

  (* Compute coefficients until we reach the desired tolerance *)
  maxcoeff = Max[Abs[an[nmin]], Abs[an[nmin+1]]];
  While[Abs[an[nmax]]/maxcoeff > 10^-prec,
    maxcoeff = Max[maxcoeff, Abs[an[nmax++]]];
  ];
  anTab = Table[an[n], {n, nmin, nmax}];

  (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)
  norm = Sqrt[2\[Pi]] (2^(1+2 k1+2 k2) E^(-2 \[Gamma]) Gamma[1+2 k2] Sum[anTab[[1;;i-nmin+1]].anTab[[i-nmin+1;;1;;-1]] 2^i Pochhammer[i+2 (1+k1+k2),-2k2-1] Hypergeometric1F1[1+i+2 k1,i+2 (1+k1+k2),4 \[Gamma]], {i, nmin, nmax}])^(1/2);

  (* Return a SpinWeightedSpheroidalHarmonicSFunction which can be evaluated for arbitratry \[Theta], \[Phi] *)
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], {anTab/norm, nmin, nmax}, Method -> "Leaver"]
];

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


(* ::Section:: *)
(*SpinWeightedSpheroidalHarmonicSFunction*)


(**********************************************************)
(* SpinWeightedSpheroidalHarmonicSFunction                *)
(**********************************************************)

SyntaxInformation[SpinWeightedSpheroidalHarmonicSFunction] =
 {"ArgumentsPattern" -> {_, _, _, _, {{__}, _, _}, ___}};

Options[SpinWeightedSpheroidalHarmonicSFunction] = {Method -> "Leaver"};

SetAttributes[SpinWeightedSpheroidalHarmonicSFunction, {NumericFunction}];

Format[SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, coeffs_ /;(Head[coeffs]=!=SequenceForm), Method -> method_]] :=
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], SequenceForm@@{"<<", coeffs[[2]]+coeffs[[3]]+1, ">>"}, Method -> method];

SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {dn_List, nmin_Integer, nmax_Integer}, Method -> "Eigenvalue"][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
 Sum[dn[[k+nmin+1]]SpinWeightedSphericalHarmonicY[s, l+k, m, \[Theta], 0], {k, -nmin, nmax}]Exp[I m \[Phi]];

SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {an_List, nmin_Integer, nmax_Integer}, Method -> "Leaver"][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
 Module[{u = Cos[\[Theta]], k1 = Abs[m-s]/2, k2 = Abs[m+s]/2, oneplusu, oneminusu},
  oneplusu = 2 Cos[\[Theta]/2]^2;
  oneminusu = 2 Sin[\[Theta]/2]^2;
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *)
  E^(\[Gamma] u) oneplusu^k1 oneminusu^k2 an.oneplusu^Range[0,nmax] Exp[I m \[Phi]]
];

(*Derivatives*)
Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {dn_List, nmin_Integer, nmax_Integer}, Method -> "Eigenvalue"]][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=Module[{\[Theta]1,\[Phi]1},
  Sum[dn[[k+nmin+1]]D[SpinWeightedSphericalHarmonicY[s, l+k, m, \[Theta]1, 0],{\[Theta]1,d1}], {k, -nmin, nmax}]D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]}
];

Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {an_List, nmin_Integer, nmax_Integer}, Method -> "Leaver"]][\[Theta]_?NumericQ, \[Phi]_?NumericQ]:=
 Module[{\[Theta]1, \[Phi]1, u, k1 = Abs[m-s]/2, k2 = Abs[m+s]/2, oneplusu, oneminusu},
  u=Cos[\[Theta]1];
  oneplusu = 2 Cos[\[Theta]1/2]^2;
  oneminusu = 2 Sin[\[Theta]1/2]^2;
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *)
  D[E^(\[Gamma] u) oneplusu^k1 oneminusu^k2 an.oneplusu^Range[0,nmax],{\[Theta]1,d1}] D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]}
];
  
SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, 0, Method -> "SphericalExact"][\[Theta]_,\[Phi]_] := SpinWeightedSphericalHarmonicY[s,l,m,\[Theta],\[Phi]]


(* ::Section::Closed:: *)
(*SpinWeightedSphericalHarmonicY*)


(**********************************************************)
(* SpinWeightedSphericalHarmonicY                         *)
(**********************************************************)

SyntaxInformation[SpinWeightedSphericalHarmonicY] =
 {"ArgumentsPattern" -> {_, _, _, _, _, ___}};
SetAttributes[SpinWeightedSphericalHarmonicY, {NumericFunction, Listable}];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] :=
  (-1)^m Sqrt[((l+m)!(l-m)!(2l+1))/(4\[Pi] (l+s)!(l-s)!)] Sum[Binomial[l-s,r] Binomial[l+s,r+s-m] (-1)^(l-r-s) If[(2 l - 2 r - s + m)==0, 1, Sin[\[Theta]/2]^(2 l - 2 r - s + m)] If[(2 r + s - m)==0, 1, Cos[\[Theta]/2]^(2 r + s - m)],{r,Max[m-s,0],Min[l-s,l+m]}]Exp[I m \[Phi]];

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] /; Abs[m]>l = 0;

SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, 0.] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0];

SpinWeightedSphericalHarmonicY[0, l_, m_, \[Theta]_, \[Phi]_] :=
  SphericalHarmonicY[l, m, \[Theta], \[Phi]];

End[];

EndPackage[];
