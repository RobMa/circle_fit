#!/usr/bin/env python3
import os
from ctypes import cdll, c_int, c_double, pointer
libpath = os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'libcircle_fit.so')
circle_fit_lib = cdll.LoadLibrary(libpath)

def estimate_circle_taubin(x, y):
    number_of_points = c_int(len(x))
    DoubleArray = c_double * len(x)
    x_arr = DoubleArray(*x)
    y_arr = DoubleArray(*y)
    a = c_double(0)
    b = c_double(0)
    r = c_double(0)
    circle_fit_lib.estimate_circle_taubin(x_arr, y_arr, number_of_points, pointer(a), pointer(b), pointer(r))
    return (r.value,a.value,b.value)


def estimate_circle_lm(x,y):
    number_of_points = c_int(len(x))
    DoubleArray = c_double * len(x)
    x_arr = DoubleArray(*x)
    y_arr = DoubleArray(*y)
    a = c_double(0)
    b = c_double(0)
    r = c_double(0)
    circle_fit_lib.estimate_circle_lm(x_arr, y_arr, number_of_points, pointer(a), pointer(b), pointer(r))
    return (r.value,a.value,b.value)


def estimate_circle(x,y):
    number_of_points = c_int(len(x))
    DoubleArray = c_double * len(x)
    x_arr = DoubleArray(*x)
    y_arr = DoubleArray(*y)
    a = c_double(0)
    b = c_double(0)
    r = c_double(0)
    circle_fit_lib.estimate_circle(x_arr, y_arr, number_of_points, pointer(a), pointer(b), pointer(r))
    return (r.value,a.value,b.value)

