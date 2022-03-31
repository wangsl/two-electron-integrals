
BeginPackage["RysChebyshevFitting`"]

(* Reference: M. M. Shepherd and J. G. Laframboise Mathematics of Computation, v36, p249 (1981) *)

RysChebyshevFitting::usage = "Chebyshev polynomails fitting for Rys polynomials roots and weights"

Begin["`private`"]

<<RysRootsWeights`

$MinPrecision = 512;
$MaxExtraPrecision = 1024;
$WorkingPrecision = 512;
SetOptions[$Output, PageWidth->100];

RysChebyshevFitting[np_, a_, b_, m_] :=
  Module[{t, k, x, T, W, TW, i, idx, cT, cW, j, FitT, FitW},
    t = Table[Cos[(2k+1)/(m+1) Pi/2], {k, 0, m}]; 
    x = 1/2(a+b) + 1/2(b-a) t; 

    T = Table[0, {i, 1, m+1}];
    W = Table[0, {i, 1, m+1}];
    For[i=1, i<=Length[x], i++,
      TW = RsyRootsWeights[x[[i]], np];
      T[[i]] = TW[[1]];
      W[[i]] = TW[[2]];
    ]

    For[idx=1, idx<=np, idx++,
      cT = Table[Sum[T[[k+1, idx]] ChebyshevT[j, t[[k+1]]], {k, 0, m}]/(m+1) 2, {j, 0, m}];
      cT[[1]] = cT[[1]]/2;
      FitT[x_] := Sum[cT[[i+1]] ChebyshevT[i, (2x-(a+b))/(b-a)], {i, 0, m}];

      cW = Table[Sum[W[[k+1, idx]] ChebyshevT[j, t[[k+1]]], {k, 0, m}]/(m+1) 2, {j, 0, m}];
      cW[[1]] = cW[[1]]/2;
      FitW[x_] := Sum[cW[[i+1]] ChebyshevT[i, (2x-(a+b))/(b-a)], {i, 0, m}];

      Print[]
      For[i=0, i<=m, i++,
        Print[ "Coeffs: ", np, " ", idx, " ", i, " ", 
              NumberForm[cT[[i+1]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
              NumberForm[cW[[i+1]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)]
        ]
      ]
    ]
  ];

End[]
EndPackage[]
