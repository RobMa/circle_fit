#!/usr/bin/env python3

import sympy as sp
from codegenerator import codegen

Mzz, Mxz, Myz, Mxx, Mxy, Myy, n, eta = sp.symbols('Mzz Mxz Myz Mxx Mxy Myy n eta', positive=True)
Mz, Mx, My = sp.symbols('Mz Mx My', real=True)

M = sp.Matrix([
    [Mzz, Mxz, Myz, Mz],
    [Mxz, Mxx, Mxy, Mx],
    [Myz, Mxy, Myy, My],
    [ Mz,  Mx,  My,  n]
])
Msym = sp.MatrixSymbol('M_out', 4, 4)

C = sp.Matrix([
        [4*Mz, 2*Mx, 2*My, 0],
        [2*Mx, n,    0,    0],
        [2*My, 0,    n,    0],
        [0,    0,    0,    0]
])
Csym = sp.MatrixSymbol('C_out', 4, 4)

Q3 = sp.simplify((M-eta*C).det())
Q3 = sp.collect(Q3, eta)
Q3_monic = sp.Matrix.zeros(4, 1)
Q3_monic[0] = Q3.coeff(eta, 0)
Q3_monic[1] = Q3.coeff(eta, 1)
Q3_monic[2] = Q3.coeff(eta, 2)
Q3_monic[3] = Q3.coeff(eta, 3)
Q3_monic_sym = sp.MatrixSymbol('Q3_out', *Q3_monic.shape)

codegen([ sp.Eq(Csym, C), sp.Eq(Msym, M) ], 'src/generated/M_C.hpp.template', 'src/generated/M_C.hpp')

codegen([sp.Eq(Q3_monic_sym, Q3_monic)], 'src/generated/Q3.hpp.template', 'src/generated/Q3.hpp')

