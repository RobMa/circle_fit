#!/usr/bin/env python3

# Three point exact fitting Algorithm described in "Direct Least-Squares Fitting of Algebraic Surfaces" from Vaughan Pratt (1987).
# Works well only in the absence of noise.

import sympy as sp
from math import pi, sin, cos, sqrt, fabs
from random import uniform

""" 
Derive the formula for the three point exact circle fit
"""
def derive_q():
    x1, y1, x2, y2, x3, y3 = sp.symbols('x1 y1 x2 y2 x3 y3', real=True)
    Aplus = sp.Matrix([
    [1,x1,y1,x1**2+y1**2],
    [1,x2,y2,x2**2+y2**2],
    [1,x3,y3,x3**2+y3**2],
    [1,x1,y1,x1**2+y1**2]])
    q = Aplus.nullspace()[0]
    return sp.simplify(q), (x1,y1,x2,y2,x3,y3)

q,q_syms = derive_q()

def circle_fit_three_points(x1,y1,x2,y2,x3,y3):
    substitutions = { q_syms[0]: x1, q_syms[1]: y1, q_syms[2]:x2, q_syms[3]:y2, q_syms[4]:x3, q_syms[5]:y3 }
    q_sol = q.subs(substitutions)
    x0 = -q_sol[1]/2/q_sol[3]
    y0 = -q_sol[2]/2/q_sol[3]
    R = sqrt(fabs( x0**2 + y0**2 - q_sol[0]/q_sol[3] ))
    return R, x0, y0

class Circle():
    def __init__(self, R=0,x=0,y=0):
        self.R = R
        self.x=x
        self.y=y

    def __str__(self):
        return 'Circle(R={0},x={1},y={2}'.format(self.R, self.x, self.y)

test_circle = Circle(R=uniform(1,10), x=uniform(-10,10), y=uniform(-10,10))
print('Test Circle:', test_circle)

random_angles = [uniform(0,pi*2) for it in range(3)]
random_points = [[test_circle.x+test_circle.R*cos(angle), test_circle.y+test_circle.R*sin(angle)] for angle in random_angles]

fit_R, fit_x, fit_y = circle_fit_three_points(random_points[0][0], random_points[0][1], random_points[1][0], random_points[1][1], random_points[2][0], random_points[2][1])

fit_circle = Circle(R=fit_R, x=fit_x, y=fit_y)
print('Fit Circle:', fit_circle)



