#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib
from math import pi, sin, cos, sqrt
from random import uniform, gauss, seed
import os
from ctypes import cdll, c_int, c_double, pointer, POINTER
libpath = os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'libcircle_fit_dynamic.so')
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


def estimate_circle_lm(x, y, rxy_init):
    number_of_points = c_int(len(x))
    DoubleArray = c_double * len(x)
    x_arr = DoubleArray(*x)
    y_arr = DoubleArray(*y)
    rxy_init_arr = (c_double*3)(*rxy_init)
    rxy_out_arr = (c_double*3)(0, 0, 0)
    circle_fit_lib.estimate_circle_lm(x_arr, y_arr, number_of_points, rxy_init_arr, rxy_out_arr)
    return (rxy_out_arr[0],rxy_out_arr[1],rxy_out_arr[2])


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


class LmTrace:
    def __init__(self, numiters):
        self.x = [0.0]*numiters
        self.y = [0.0]*numiters
        self.r = [0.0]*numiters
        self.mse = [0.0]*numiters
        self.grad_norm = [0.0]*numiters
        self.numiters = numiters

def estimate_circle_lm_trace(x, y, rxy_init):
    number_of_points = c_int(len(x))
    DoubleArray = c_double * len(x)
    x_arr = DoubleArray(*x)
    y_arr = DoubleArray(*y)
    rxy_init_arr = (c_double*3)(*rxy_init)
    out_numiters = c_int(0)
    DoublePointer = POINTER(c_double)
    out_x = DoublePointer()
    out_y = DoublePointer()
    out_r = DoublePointer()
    out_grad = DoublePointer()
    out_mse = DoublePointer()
    circle_fit_lib.estimate_circle_lm_trace(x_arr, y_arr, number_of_points, rxy_init_arr,
        pointer(out_numiters), pointer(out_x), pointer(out_y), pointer(out_r), pointer(out_grad), pointer(out_mse))
    trace = LmTrace(out_numiters.value)
    for i in range(out_numiters.value):
        trace.x[i] = out_x[i]
        trace.y[i] = out_y[i]
        trace.r[i] = out_r[i]
        trace.mse[i] = out_mse[i]
        trace.grad_norm[i] = out_grad[i]
    return trace


def circle_to_points(thetas, a=0.5, b=0.5, r=0.5):
    points = [(cos(theta)*r+a, sin(theta)*r+b) for theta in thetas]
    x, y = [point[0] for point in points], [point[1] for point in points]
    return x, y


def sample_circle(n=3, theta_range=(-pi,pi), noise_sigma=0, a=0.5, b=0.5, r=0.5, noise_uniform=0):
    thetas = [uniform(theta_range[0], theta_range[1]) for i in range(n)]
    points = [(lambda theta: [
        cos(theta) * r + a + gauss(0, noise_sigma) + uniform(-noise_uniform, noise_uniform),
        sin(theta) * r + b + gauss(0, noise_sigma) + uniform(-noise_uniform, noise_uniform)
        ])(theta) for theta in thetas]
    return [point[0] for point in points], [point[1] for point in points]


def plot_circle(ax,label,x=0.5,y=0.5,r=0.5,color='black',linestyle='-',linewidth=1):
    circ_patch = plt.Circle([x, y], r)
    circ_patch.set_fill(False)
    circ_patch.set_linestyle(linestyle)
    circ_patch.set_color(color)
    circ_patch.set_linewidth(linewidth)
    circ_patch.set_label(label)
    ax.add_patch(circ_patch)
    harrow, = ax.plot([x,x+r], [y,y], linestyle=linestyle, linewidth=linewidth, color=color)
    varrow, = ax.plot([x,x], [y,y+r], linestyle=linestyle, linewidth=linewidth, color=color)
    # harrow = plt.Arrow(x,y,r,0, width=0)
    # harrow.set_linestyle(linestyle)
    # harrow.set_color(color)
    # harrow.set_linewidth(linewidth)
    # varrow = plt.Arrow(x,y,0,r, width=0)
    # varrow.set_linestyle(linestyle)
    # varrow.set_color(color)
    # varrow.set_linewidth(linewidth)
    # ax.add_patch(harrow)
    # ax.add_patch(varrow)
    return (circ_patch, harrow, varrow)


def plot_circle_update(rxy, patches):
    circ_patch, harrow, varrow = patches
    circ_patch.center = rxy[1:3]
    circ_patch.set_radius(rxy[0])
    harrow.set_xdata([rxy[1]+rxy[0], rxy[1]])
    harrow.set_ydata([rxy[2], rxy[2]])
    varrow.set_xdata([rxy[1], rxy[1]])
    varrow.set_ydata([rxy[2], rxy[2]+rxy[0]])



def mean_squared_error(x,y,rxy):
    E = 0
    for i, xi in enumerate(x):
        yi = y[i]
        E = E + (sqrt((xi-rxy[1])**2+(yi-rxy[2])**2)-rxy[0])**2
    E = E / len(x)
    return E

