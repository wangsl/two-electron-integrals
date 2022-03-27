
<<RysRootsWeights`

n = 5;

TW = RsyRootsWeights[5.0, n];
For[i=1, i<=n, i++,
  If[i == 1, Print[]];
  Print[ "RW ", 
        NumberForm[TW[[1]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
        NumberForm[TW[[2]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)]
        ];
]

TW = RsyRootsWeights[5.0, n, 1];
For[i=1, i<=n, i++,
  If[i == 1, Print[]];
  Print[ "UW ", 
        NumberForm[TW[[1]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
        NumberForm[TW[[2]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)]
        ];
]
