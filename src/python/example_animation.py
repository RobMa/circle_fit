import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib
import yaml
import circle_fit
from circle_fit import circle_to_points, sample_circle, plot_circle, mean_squared_error, plot_circle_update

from matplotlib import rc
rc('font',**{'family':'serif','sans-serif':['Computer Modern Roman']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

orig = (2, 100, 0)
with open("xy.yaml", 'r') as f:
    x, y =yaml.load(f)

orig_mse = mean_squared_error(x, y, orig)
print('orig_mse=', orig_mse)
taubin = circle_fit.estimate_circle_taubin(x, y)
taubin_mse = mean_squared_error(x, y, taubin)
print('taubin_mse=', taubin_mse)
trace = circle_fit.estimate_circle_lm_trace(x, y, taubin)

fig, (ax1, ax2) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3, 1]})
fig.set_size_inches([5*3/4, 5], forward=True)
fig.set_tight_layout(True)

ax1.set_ylim([-2, 2.1])
ax1.set_xlim([98, 102.1])
ax2.set_xlim([0, 20])
ax2.set_ylim([trace.mse[-1], trace.mse[0] + 5e-5])
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax2.set_xlabel('Iteration')
ax2.set_ylabel('MSE')

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))

ax1.scatter(x, y, label='Samples', marker=matplotlib.markers.MarkerStyle(marker='.'))
orig_circle = plot_circle(ax1, 'Original', color='black', linestyle='-',x=orig[1],y=orig[2],r=orig[0])
ax2.plot([0,20], [orig_mse, orig_mse], color='black', linestyle='-')
lm_mse, = ax2.plot([0,0], [0,0], color='red', linestyle='-.')

taubin_circle = plot_circle(ax1, 'Taubin', x=-100, y=-100, r=0.01, color='green', linestyle='--')
lm_circle = plot_circle(ax1, 'LM', x=-100, y=-100, r=0.01, color='red', linestyle='-.')
ax1.legend()
ax1.set_title('Input points')

def update(i):
    label = 'timestep {0}'.format(i)
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    if i == 2:
        plot_circle_update(taubin, taubin_circle)
        # plot_circle(ax1, 'taubin', x=taubin[1], y=taubin[2], r=taubin[0], color='green', linestyle='--')
        ax2.plot([0,20], [taubin_mse,taubin_mse], color='green', linestyle='--')
        ax1.set_title('Algebraic fit: GWAF Taubin')
        return ax1, ax2, taubin_circle
    if i >= 4 and i-4 < trace.numiters:
        ax1.set_title('LM-Algorithm: Iteration {}'.format(i-4))
        plot_circle_update([trace.r[i-4], trace.x[i-4], trace.y[i-4]], lm_circle)
        lm_mse.set_xdata(np.arange(0,i-3))
        lm_mse.set_ydata(trace.mse[0:i-3])
        return ax1, ax2, lm_circle
    # line.set_ydata(x - 5 + i)
    # ax.set_xlabel(label)
    return

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, trace.numiters+5), interval=2000)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('example_animation.gif', dpi=80, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()

