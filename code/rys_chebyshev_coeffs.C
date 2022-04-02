#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "gauss_hermite.h"
#include "rys_chebyshev_coeffs.h"
#include "fort.h"

#include "rys_chebyshev_coeffs_t.h"
#include "rys_chebyshev_coeffs_u.h"

#define RYS_CHEBYSHEV(X) \
  { X.rys_order, X.chebyshev_order, \
    X.x_min, X.x_max, \
    (const double *) X.roots_coefficients, \
    (const double *) X.weights_coefficients \
  }

void RysChebyshev::setup_parameters()
{
  if(has_setup_parameters) return;
  setup_parameters_u();
  has_setup_parameters = 1;
}

static void RysChebyshev::setup_parameters_t()
{
  int i = -1;
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_0);      //   0 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_1);      //   1 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_2);      //   2 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_3);      //   3 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_4);      //   4 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_5);      //   5 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_6);      //   6 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_7);      //   7 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_8);      //   8 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_9);      //   9 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_10);     //  10 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_11);     //  11 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_1_12);     //  12 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_0);      //  13 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_1);      //  14 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_2);      //  15 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_3);      //  16 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_4);      //  17 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_5);      //  18 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_6);      //  19 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_7);      //  20 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_8);      //  21 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_9);      //  22 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_10);     //  23 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_11);     //  24 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_2_12);     //  25 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_0);      //  26 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_1);      //  27 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_2);      //  28 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_3);      //  29 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_4);      //  30 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_5);      //  31 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_6);      //  32 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_7);      //  33 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_8);      //  34 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_9);      //  35 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_10);     //  36 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_11);     //  37 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_3_12);     //  38 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_0);      //  39 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_1);      //  40 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_2);      //  41 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_3);      //  42 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_4);      //  43 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_5);      //  44 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_6);      //  45 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_7);      //  46 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_8);      //  47 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_9);      //  48 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_10);     //  49 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_11);     //  50 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_4_12);     //  51 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_0);      //  52 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_1);      //  53 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_2);      //  54 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_3);      //  55 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_4);      //  56 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_5);      //  57 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_6);      //  58 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_7);      //  59 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_8);      //  60 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_9);      //  61 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_10);     //  62 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_11);     //  63 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_5_12);     //  64 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_0);      //  65 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_1);      //  66 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_2);      //  67 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_3);      //  68 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_4);      //  69 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_5);      //  70 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_6);      //  71 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_7);      //  72 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_8);      //  73 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_9);      //  74 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_10);     //  75 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_11);     //  76 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_6_12);     //  77 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_0);      //  78 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_1);      //  79 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_2);      //  80 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_3);      //  81 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_4);      //  82 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_5);      //  83 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_6);      //  84 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_7);      //  85 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_8);      //  86 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_9);      //  87 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_10);     //  88 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_11);     //  89 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_7_12);     //  90 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_0);      //  91 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_1);      //  92 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_2);      //  93 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_3);      //  94 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_4);      //  95 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_5);      //  96 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_6);      //  97 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_7);      //  98 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_8);      //  99 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_9);      // 100 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_10);     // 101 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_11);     // 102 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_8_12);     // 103 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_0);      // 104 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_1);      // 105 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_2);      // 106 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_3);      // 107 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_4);      // 108 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_5);      // 109 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_6);      // 110 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_7);      // 111 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_8);      // 112 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_9);      // 113 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_10);     // 114 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_11);     // 115 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_9_12);     // 116 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_0);     // 117 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_1);     // 118 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_2);     // 119 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_3);     // 120 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_4);     // 121 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_5);     // 122 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_6);     // 123 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_7);     // 124 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_8);     // 125 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_9);     // 126 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_10);    // 127 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_11);    // 128 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_10_12);    // 129 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_0);     // 130 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_1);     // 131 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_2);     // 132 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_3);     // 133 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_4);     // 134 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_5);     // 135 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_6);     // 136 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_7);     // 137 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_8);     // 138 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_9);     // 139 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_10);    // 140 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_11);    // 141 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_11_12);    // 142 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_0);     // 143 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_1);     // 144 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_2);     // 145 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_3);     // 146 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_4);     // 147 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_5);     // 148 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_6);     // 149 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_7);     // 150 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_8);     // 151 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_9);     // 152 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_10);    // 153 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_11);    // 154 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_t_12_12);    // 155 
}

static void RysChebyshev::setup_parameters_u()
{
  int i = -1;
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_0);      //   0 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_1);      //   1 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_2);      //   2 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_3);      //   3 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_4);      //   4 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_5);      //   5 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_6);      //   6 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_7);      //   7 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_8);      //   8 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_9);      //   9 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_10);     //  10 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_11);     //  11 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_1_12);     //  12 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_0);      //  13 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_1);      //  14 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_2);      //  15 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_3);      //  16 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_4);      //  17 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_5);      //  18 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_6);      //  19 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_7);      //  20 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_8);      //  21 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_9);      //  22 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_10);     //  23 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_11);     //  24 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_2_12);     //  25 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_0);      //  26 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_1);      //  27 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_2);      //  28 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_3);      //  29 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_4);      //  30 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_5);      //  31 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_6);      //  32 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_7);      //  33 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_8);      //  34 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_9);      //  35 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_10);     //  36 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_11);     //  37 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_3_12);     //  38 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_0);      //  39 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_1);      //  40 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_2);      //  41 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_3);      //  42 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_4);      //  43 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_5);      //  44 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_6);      //  45 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_7);      //  46 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_8);      //  47 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_9);      //  48 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_10);     //  49 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_11);     //  50 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_4_12);     //  51 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_0);      //  52 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_1);      //  53 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_2);      //  54 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_3);      //  55 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_4);      //  56 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_5);      //  57 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_6);      //  58 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_7);      //  59 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_8);      //  60 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_9);      //  61 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_10);     //  62 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_11);     //  63 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_5_12);     //  64 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_0);      //  65 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_1);      //  66 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_2);      //  67 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_3);      //  68 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_4);      //  69 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_5);      //  70 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_6);      //  71 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_7);      //  72 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_8);      //  73 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_9);      //  74 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_10);     //  75 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_11);     //  76 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_6_12);     //  77 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_0);      //  78 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_1);      //  79 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_2);      //  80 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_3);      //  81 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_4);      //  82 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_5);      //  83 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_6);      //  84 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_7);      //  85 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_8);      //  86 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_9);      //  87 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_10);     //  88 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_11);     //  89 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_7_12);     //  90 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_0);      //  91 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_1);      //  92 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_2);      //  93 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_3);      //  94 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_4);      //  95 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_5);      //  96 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_6);      //  97 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_7);      //  98 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_8);      //  99 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_9);      // 100 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_10);     // 101 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_11);     // 102 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_8_12);     // 103 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_0);      // 104 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_1);      // 105 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_2);      // 106 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_3);      // 107 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_4);      // 108 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_5);      // 109 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_6);      // 110 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_7);      // 111 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_8);      // 112 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_9);      // 113 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_10);     // 114 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_11);     // 115 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_9_12);     // 116 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_0);     // 117 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_1);     // 118 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_2);     // 119 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_3);     // 120 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_4);     // 121 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_5);     // 122 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_6);     // 123 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_7);     // 124 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_8);     // 125 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_9);     // 126 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_10);    // 127 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_11);    // 128 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_10_12);    // 129 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_0);     // 130 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_1);     // 131 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_2);     // 132 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_3);     // 133 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_4);     // 134 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_5);     // 135 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_6);     // 136 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_7);     // 137 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_8);     // 138 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_9);     // 139 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_10);    // 140 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_11);    // 141 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_11_12);    // 142 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_0);     // 143 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_1);     // 144 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_2);     // 145 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_3);     // 146 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_4);     // 147 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_5);     // 148 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_6);     // 149 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_7);     // 150 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_8);     // 151 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_9);     // 152 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_10);    // 153 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_11);    // 154 
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_12_12);    // 155 
}

#undef RYS_CHEBYSHEV

const RysChebyshevCoeffs *RysChebyshev::parameter(const int n_rys, const double x)
{
  for(int i = 0; i < parameters_length; i++) {
    const RysChebyshevCoeffs &p = parameters[i];
    if(n_rys == p.rys_order && p.x_min <= x && x <= p.x_max) 
      return &p;
  }
  return 0;
}

void RysChebyshev::initialize()
{
  RysChebyshev::setup_parameters();
  GaussHermite::setup_parameters();
}

inline void chebyshev_t(const int n, const double x, double *t)
{
  assert(-1.0 <= x && x <= 1.0);

  if(n >= 1) t[0] = 1.0;
  if(n >= 2) t[1] = x;
  
  if(n <= 2) return;

  const double two_x = x + x; 
  for(int i = 2; i <= n; i++)   
    t[i] = two_x*t[i-1] - t[i-2];
}

// Numerical Methods for Special Functions by Amparo Gil, Javier Segura, and Nico Temme
// Page 75 Algorithm 3.1. Clenshaw’s method for a Chebyshev sum

void RysChebyshev::calculate_rys_roots_and_weights(
  const int rys_order, const double x, 
  double *roots, double *weights,
  const int need_t)
{
  initialize();

  const GaussHermiteRootAndWeights *gh_parameger = GaussHermite::parameter(rys_order);
  if(x > gh_parameger->x_min) {
    const double one_over_sqrt_x = pow(x, -0.5);
    for(int k = 0; k < rys_order; k++) {
      const double t = one_over_sqrt_x*gh_parameger->roots[k];
      const double t2 = t*t;
      roots[k] = t2/(1.0-t2);
      weights[k] = one_over_sqrt_x*gh_parameger->weights[k];
    }
  } else {
    const RysChebyshevCoeffs *coeffs = RysChebyshev::parameter(rys_order, x);
    assert(coeffs && coeffs->rys_order == rys_order);

    const int &chebyshev_order = coeffs->chebyshev_order;
    const double &x_min = coeffs->x_min;
    const double &x_max = coeffs->x_max;
    
    double *T = new double [chebyshev_order+1];
    assert(T);
    chebyshev_t(chebyshev_order, (2*x-x_min-x_max)/(x_max-x_min), T);

    for(int k = 0; k < rys_order; k++) {
      const double *r_coeffs = coeffs->roots_coefficients + k*(chebyshev_order+1);
      const double *w_coeffs = coeffs->weights_coefficients + k*(chebyshev_order+1);
      double &r_sum = roots[k];
      double &w_sum = weights[k];
      r_sum = 0.0;
      w_sum = 0.0;
      for(int i = 0; i < chebyshev_order+1; i++) {
        r_sum += r_coeffs[i]*T[i];
        w_sum += w_coeffs[i]*T[i];
      }
    }
    if(T) { delete [] T; T = 0; }
  }

  if(need_t) {
    for(int k = 0; k < rys_order; k++) {
      const double &u = roots[k];
      roots[k] = sqrt(u/(1+u));
    }
  }
}                                                    

// Reference: J. Comput. Phys. 21, 144 (1976) Table 1 for n=5
/*
    0.00   1  0.148874338981631216  0.295524224714752870E+00
    0.00   2  0.433395394129247213  0.269266719309996461E+00
    0.00   3  0.679409568299024436  0.219086362515982069E+00
    0.00   4  0.865063366688984536  0.149451349150580587E+00
    0.00   5  0.973906528517172076  0.666713443086881657E-01

    5.00   1  0.120616479067479812  0.224047067536327305E+00
    5.00   2  0.359993608978937729  0.123946126190236719E+00
    5.00   3  0.591831842524405793  0.389863670917377156E-01
    5.00   4  0.802341831901655866  0.763621205330385338E-02
    5.00   5  0.956727269752466358  0.109653673890797621E-02

   10.00   1  0.101239395075509955  0.182931707896803408E+00
   10.00   2  0.304887573873516171  0.809149721111935422E-01
   10.00   3  0.511824380911261478  0.152238928314807464E-01
   10.00   4  0.722257402679194249  0.114136525649064825E-02
   10.00   5  0.921039112895102874  0.354524106744085883E-04
*/

void RysChebyshev::test()
{
  const int rys_order = 5;

  double *roots = new double [rys_order];
  assert(roots);
  double *weights = new double [rys_order];
  assert(weights);

  std::cout << std::endl;
  for(double x = 0.0; x < 120.1; x += 5.0) {
    RysChebyshev::calculate_rys_roots_and_weights(rys_order, x, roots, weights, 1);
    for(int k = 0; k < rys_order; k++) {
      std::cout << " " 
                << k << " "   
                << std::fixed << std::setprecision(2) << x << " " 
                << std::scientific << std::setprecision(20) << roots[k] << " " << weights[k]
                << std::endl;
    }
    std::cout << std::endl;
  }

  if(roots) { delete [] roots; roots = 0; }
  if(weights) { delete [] weights; weights = 0; }
}

// Fortran versions
extern "C" {
  // CalculateRysRootsAndWeights
  void FORT(calculaterysrootsandweights)(
    const int &rys_order, const double &x, 
    double *roots, double *weights, 
    const int &need_u)
  {
    RysChebyshev::calculate_rys_roots_and_weights(rys_order, x, roots, weights, need_u);
    std::cout.flush();
  }

  // RysChebyshevInitialize
  void FORT(ryschebyshevinitialize)() { RysChebyshev::initialize(); }
}
