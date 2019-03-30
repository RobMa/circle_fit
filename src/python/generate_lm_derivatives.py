#!/usr/bin/python3
import sympy as sp
print('Deriving formulas...')
x, y, A, D, Theta = sp.symbols('x y A D Theta', real=True)

P = A * (x**2 + y**2) + sp.sqrt(1.0+4*A*D) * (x*sp.cos(Theta) + y*sp.sin(Theta)) + D

d = (2.0*P)/(1.0+sp.sqrt(1.0+4.0*A*P))

d_dA = sp.diff(d, A)
d_dD = sp.diff(d, D)
d_dTheta = sp.diff(d, Theta)

d_dA_dA = sp.diff(d_dA, A)
d_dA_dD = sp.diff(d_dA, D)
d_dA_dTheta = sp.diff(d_dA, Theta)

d_dD_dD = sp.diff(d_dD, D)
d_dD_dTheta = sp.diff(d_dD, Theta)

d_dTheta_dTheta = sp.diff(d_dTheta, Theta)

from codegenerator import codegen

print('Generating code...')

codegen([
    sp.Eq(sp.Symbol('d'), d),

    sp.Eq(sp.Symbol('d_dA'), d_dA),
    sp.Eq(sp.Symbol('d_dD'), d_dD),
    sp.Eq(sp.Symbol('d_dTheta'), d_dTheta),

    sp.Eq(sp.Symbol('d_dA_dA'), d_dA_dA),
    sp.Eq(sp.Symbol('d_dA_dD'), d_dA_dD),
    sp.Eq(sp.Symbol('d_dA_dTheta'), d_dA_dTheta),

    sp.Eq(sp.Symbol('d_dD_dD'), d_dD_dD),
    sp.Eq(sp.Symbol('d_dD_dTheta'), d_dD_dTheta),

    sp.Eq(sp.Symbol('d_dTheta_dTheta'), d_dTheta_dTheta),
], 'src/generated/lm_derivatives.hpp.template', 'src/generated/lm_derivatives.hpp')

