import numpy as np
from pylab import *

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator

# define how your plots look:

def setup_axes(fig, rect, theta, radius):

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()

    # Find grid values appropriate for the coordinate (degree).
    # The argument is an approximate number of grids.
    grid_locator1 = angle_helper.LocatorD(2)

    # And also use an appropriate formatter:
    tick_formatter1 = angle_helper.FormatterDMS()

    # set up number of ticks for the r-axis
    grid_locator2 = MaxNLocator(5)
    grid_locator1 = MaxNLocator(6)

    # the extremes are passed to the function
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(theta[0], theta[1], radius[0],
                                          radius[1]),
                                grid_locator1=grid_locator1,
                                grid_locator2=grid_locator2,
                                tick_formatter1=tick_formatter1,
                                tick_formatter2=None,
                                )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    # the axis artist lets you call axis with
    # "bottom", "top", "left", "right"
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")

    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")

    ax1.axis["left"].label.set_text("$D$ (kpc)")
    ax1.axis["top"].label.set_text(ur"Azimuthal Angle, $\Phi$")

    ax1.grid(True)
    
    # create a parasite axes
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                         # drawn twice, and possibly over some other
                         # artists. So, we decrease the zorder a bit to
                         # prevent this.

    return ax1, aux_ax


"""

The code above is from a google search.

To plot:

"""

fig = figure(1,figsize=(9,4))

ax1, aux_ax1 = setup_axes(fig,111,theta=[0.0,90.0],radius=[0,210])

aux_ax1.scatter(pa,D,s=vrange,marker='o',zorder=30,c='#602C6E',
                linewidth=0,alpha=0.9)

"""

Simple! Hah.

I read in the pa and D values, then set vrange equal to some parameter
(velocity spread of the absorber in my case).

To do 2 subplots, change the line to:

ax1, aux_ax1 = setup_axes(fig,121,theta=[0.0,90.0],radius=[0,210])

and add:

ax2, aux_ax2 = setup_axes(fig,122,theta=[0.0,90.0],radius=[0,210])

and so on.

"""
