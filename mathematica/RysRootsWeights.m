
BeginPackage["RysRootWeights`"]

RsyRootsWeights::usage = " "

Begin["`private`"]

(* J. Comp. Phys., 21, 144 (1976) *)

$MinPrecision = 256;
$MaxExtraPrecision = 1024;
$WorkingPrecision = 512;

L[k_, r_, x_, n_] := 
  Module[{s, l},
    s = Product[(x-r[[l]])/(r[[k]]-r[[l]]), {l, Drop[Table[l, {l, 1, n}], {k}]}];
    Return[s];
  ]

X[n_, x_] := 
  Module[{s},
    s = 1/2 x^(-(n+1)/2) (Gamma[(n+1)/2] - Gamma[(n+1)/2, x]);
    Return[s];
  ]

RsyRootsWeights[x1_, ntotal_, getU_:0] :=   
  Module[{x, t, P, i, n, c, p, s, R, k1, alpha1, l, s1, k2, alpha2, s2, R0, S0},
    x = SetPrecision[x1, 512];
    P = Table[0, {i, 1, 1000}]; 
    For[n = 0, n <= ntotal, n++,
      c = Table[0, {i, 1, 1000}]; 
      p = t^(2 n);
      For[i = 0, i < n, i++,
        c[[i+1]] = Integrate[Exp[-x t^2] p P[[i+1]], {t, 0, 1}];
      ];

      For[i = 0, i < n, i++,
        p = p - c[[i+1]] P[[i+1]];
      ];
      s = Integrate[Exp[-x t^2] p^2, {t, 0, 1}];
      p = p/Sqrt[s];
      P[[n+1]] = p;
    
      R = t /. NSolve[p == 0, t];

      R0 = Table[0, {i, 1, n}];
      S0 = Table[0, {i, 1, n}];

      For[i = 1, i <= n, i++,
        k1 = i;
        alpha1 = CoefficientList[L[k1, R, t, 2*n], t];
        s1 = Sum[alpha1[[l]] X[l-1, x], {l, 1, 2*n-1}];
        k2 = 2*n + 1 - k1;
        alpha2 = CoefficientList[L[k2, R, t, 2*n], t];
        s2 = Sum[alpha2[[l]] X[l-1, x], {l, 1, 2*n - 1}];
        R0[[n-i+1]] = R[[k2]];
        S0[[n-i+1]] = s1+s2;
      ];
    ];
    Return[ If[getU == 0, {R0, S0}, {R0^2/(1-R0^2), S0}] ];
  ];

End[]
EndPackage[]
