import numpy as np
from pylab import *

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator

# define how your plots look:

def setup_axes(fig, rect, theta, radius, quad):


    # quad controls the quadrant of the plot and controls the orientation 
    # of the plot where quad=1 is upper left, 2 is upper right, 3 is 
    # lower left, 4 is lower right   
    if quad==1:
        tr_rotate = Affine2D().translate(np.pi/2.0, 0)
    elif quad==2:
        tr_rotate = Affine2D().translate(0, 0)
    elif quad==3:
        tr_rotate = Affine2D().translate(np.pi, 0)
    else:
        tr_rotate = Affine2D().translate(3.0*np.pi/2.0, 0)
        
    tr_scale = Affine2D().scale(np.pi/180., 1.) 

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = tr_scale + tr_rotate + PolarAxes.PolarTransform()

    # Find grid values appropriate for the coordinate (degree).
    # The argument is an approximate number of grids.
    grid_locator1 = angle_helper.LocatorD(2)

    # And also use an appropriate formatter:
    tick_formatter1 = angle_helper.FormatterDMS()

    # set up number of ticks for the r-axis
    grid_locator2 = MaxNLocator(5)
    grid_locator1 = MaxNLocator(6)

    # the extremes are passed to the function
    thetaMin = 0
    thetaMax = 90
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(thetaMin, thetaMax, radius[0],
                                          radius[1]),
                                grid_locator1=grid_locator1,
                                grid_locator2=grid_locator2,
                                tick_formatter1=tick_formatter1,
                                tick_formatter2=None,
                                )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    print ax1.get_xlim()
    print ax1.get_ylim()
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
    ax1.axis["top"].label.set_text(ur"$\Phi$")

    # Set the axis labels based on the quadrant
    if quad==1:
    
        # Set visibilities
        ax1.axis["right"].set_visible(True)
        ax1.axis["left"].set_visible(True)
        ax1.axis["right"].toggle(ticklabels=True, label=True)
        ax1.axis["left"].toggle(ticklabels=False, label=False)

        # Tick Labels
        #ax1.axis["right"].major_ticks.set_axis_direction('right')
        ax1.axis["right"].major_ticklabels.set_axis_direction('bottom')
        #ax1.axis["right"].major_ticklabels.set_pad(-20)
        #ax1.axis["right"].set_axis_direction('top')
        #ax1.axis["right"].set_ticklabel_direction('+')
        
        # Axis labels
        ax1.axis["right"].label.set_rotation(0)
        ax1.axis["right"].label.set_text("$D$ (kpc)")
        ax1.axis["right"].label.set_pad(20)

    elif quad==2: 
        ax1.axis["left"].set_visible(True)
        ax1.axis["left"].toggle(ticklabels=True, label=True)
        ax1.axis["left"].major_ticklabels.set_axis_direction("bottom")
        ax1.axis["left"].set_axis_direction("bottom")
        ax1.axis["left"].label.set_text("$D$ (kpc)")

    elif quad==3:
        # Set visibilities
        ax1.axis["right"].set_visible(True)
        ax1.axis["left"].set_visible(True)
        ax1.axis["right"].toggle(ticklabels=False, label=False)
        ax1.axis["left"].toggle(ticklabels=True, label=True)

        # Tick Labels
        #ax1.axis["right"].major_ticks.set_axis_direction('right')
        ax1.axis["left"].major_ticklabels.set_axis_direction('top')
        #ax1.axis["right"].major_ticklabels.set_pad(-20)
        #ax1.axis["right"].set_axis_direction('top')
        #ax1.axis["right"].set_ticklabel_direction('+')
        
        # Axis labels
        ax1.axis["left"].label.set_rotation(180)
        ax1.axis["left"].label.set_text("$D$ (kpc)")
        ax1.axis["left"].label.set_pad(20)

    else: 
        # Set visibilities
        ax1.axis["right"].set_visible(True)
        ax1.axis["left"].set_visible(True)
        ax1.axis["right"].toggle(ticklabels=True, label=True)
        ax1.axis["left"].toggle(ticklabels=False, label=False)

        # Tick Labels
        #ax1.axis["right"].major_ticks.set_axis_direction('right')
        ax1.axis["right"].major_ticklabels.set_axis_direction('top')
        #ax1.axis["right"].major_ticklabels.set_pad(-20)
        #ax1.axis["right"].set_axis_direction('top')
        #ax1.axis["right"].set_ticklabel_direction('+')
        
        # Axis labels
        ax1.axis["right"].label.set_rotation(180)
        ax1.axis["right"].label.set_text("$D$ (kpc)")
        ax1.axis["right"].label.set_pad(0)



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

pa, d, vrange = [], [], []
for i in range(20):
    pa.append(np.random.random()*90.0)
    d.append(np.random.random()*200.0)
    vrange.append(np.random.random()*25.0)

fig = figure(1,figsize=(10,10))
#fig, ((a1, a2),(a3,a4)) = subplots(2,2,figsize=(10,10))

#ax1, aux_ax1 = setup_axes(fig,a1,theta=[0.0,90.0],radius=[0,210], quad=1)
#ax2, aux_ax2 = setup_axes(fig,a2,theta=[0.0,90.0],radius=[0,210], quad=2)
#ax3, aux_ax3 = setup_axes(fig,a3,theta=[0.0,90.0],radius=[0,210], quad=3)
#ax4, aux_ax4 = setup_axes(fig,a4,theta=[0.0,90.0],radius=[0,210], quad=4)

ax1, aux_ax1 = setup_axes(fig,221,theta=[0.0,90.0],radius=[0,210], quad=1)
ax2, aux_ax2 = setup_axes(fig,222,theta=[0.0,90.0],radius=[0,210], quad=2)
ax3, aux_ax3 = setup_axes(fig,223,theta=[0.0,90.0],radius=[0,210], quad=3)
ax4, aux_ax4 = setup_axes(fig,224,theta=[0.0,90.0],radius=[0,210], quad=4)

aux_ax1.scatter(pa,d,s=vrange,marker='o',zorder=30,c='#602C6E',
                linewidth=0,alpha=0.9)
aux_ax2.scatter(pa,d,s=vrange,marker='o',zorder=30,c='#602C6E',
                linewidth=0,alpha=0.9)
aux_ax3.scatter(pa,d,s=vrange,marker='o',zorder=30,c='#602C6E',
                linewidth=0,alpha=0.9)
aux_ax4.scatter(pa,d,s=vrange,marker='o',zorder=30,c='#602C6E',
                linewidth=0,alpha=0.9)
#fig.subplots_adjust(wspace=-.5)
fig.savefig('polartest.pdf', bbox_inches='tight')
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
