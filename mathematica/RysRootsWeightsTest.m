
<<RysRootsWeights`

SetOptions[$Output, PageWidth->120];

n = 12;

(* H. F. King and M. Dupuis, J. Comput. Phys. 21, 144 (1976) 
  Table I: Rys Roots and Weights for n = 5
*)
x = 159;

TW = RsyRootsWeights[x, n];

U = TW[[1]]^2/(1-TW[[1]]^2);

For[i=1, i<=n, i++,
  If[i == 1, Print[]];
  Print[ "RW ", 
        NumberForm[TW[[1]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
        NumberForm[U[[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
        NumberForm[TW[[2]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)]
        ];
]

Quit[];

TW = RsyRootsWeights[x, n, 1];
For[i=1, i<=n, i++,
  If[i == 1, Print[]];
  Print[ "UW ", 
        NumberForm[TW[[1]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)], " "
        NumberForm[TW[[2]][[i]], {30, 26}, ScientificNotationThreshold->{0, 0}, NumberFormat->(Row[{#1, "E", #3}] &)]
        ];
]

(*
RW 1.20616479067479802744154227E-1  2.24047067536327285316724159E-1
RW 3.59993608978937694021013000E-1  1.23946126190236707856987793E-1
RW 5.91831842524405761143815416E-1  3.89863670917377183527900239E-2
RW 8.02341831901655824441309038E-1  7.63621205330385457423530003E-3
RW 9.56727269752466327509937093E-1  1.09653673890797594965577650E-3

UW 1.47631137474104001020588612E-2  2.24047067536327285316724159E-1
UW 1.48890985046712248883497665E-1  1.23946126190236707856987793E-1
UW 5.39088846984881919146386635E-1  3.89863670917377183527900239E-2
UW 1.80703657434083929347248146E    7.63621205330385457423530003E-3
UW 1.08101497669347635096310397E1   1.09653673890797594965577650E-3
*)