import circle_fit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from math import pi, sin, cos, sqrt
from random import uniform, gauss, seed
import yaml
from circle_fit import circle_to_points, sample_circle, plot_circle, mean_squared_error


# # Three points
# x, y = circle_to_points([0, pi/2, pi])
# est = circle_fit.estimate_circle(x,y)
# fig, ax = plt.subplots()
# circ_patch = plt.Circle([0.5,0.5], 0.5)
# circ_patch.set_fill(False)
# ax.add_patch(circ_patch)
# ax.scatter(x, y)
# plt.show()


# 100 points with noise
seed(a=1)
orig = (2, 100, 0)
x, y = sample_circle(n=100, noise_sigma=0.1, noise_uniform=0.0, theta_range=(0,pi/3), a=orig[1], b=orig[2], r=orig[0])
with open('xy.yaml', 'w') as f:
    yaml.dump((x,y), f)

taubin = circle_fit.estimate_circle_taubin(x,y)
est = circle_fit.estimate_circle_lm(x,y,taubin)
print('Original circle MSE=', mean_squared_error(x, y, orig))
print('Taubin circle MSE=', mean_squared_error(x, y, taubin))
print('LM circle MSE=', mean_squared_error(x, y, est))

fig, ax = plt.subplots()
plot_circle(ax, 'original', color='black', linestyle='-',x=orig[1],y=orig[2],r=orig[0])
plot_circle(ax, 'taubin', x=taubin[1], y=taubin[2], r=taubin[0], color='green', linestyle='--')
plot_circle(ax, 'lm', x=est[1], y=est[2], r=est[0], color='red', linestyle='-.')
ax.scatter(x, y, label='input', marker=matplotlib.markers.MarkerStyle(marker='.'))
plt.legend()
plt.show()


