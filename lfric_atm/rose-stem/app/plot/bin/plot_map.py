#!/usr/bin/env python
''' Quick plot of for lfric_atm global output '''

# Need to set a non-interactive backend for suites
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Index to the fields
varname=0
colbar_min=1
colbar_max=2

# Fields which are available to plot
theta           = ['theta',           265,  300]
m_v             = ['m_v',             1e-3, 15e-3]
m_cl            = ['m_cl',            0,    5e-5]
m_ci            = ['m_ci',            0,    1e-4]
u1              = ['u1',              -15,  15]
sw_heating_rate = ['sw_heating_rate', 0,    7e-5]


def load_cube_by_varname(filename, var):
   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

def gen_markersize(n_points):
    """
    The "points" based plot of the LFRic grid needs the size of the
    dots to be set sensibly or plotting artifacts are introduced.
    This method uses a rough exponential fit to pick a size.
    """
    c_value = int((n_points/6)**0.5)

    def exponential(x, a, b, c):
        return a*np.exp(-b*x)+c

    # C Value sample points
    x = [48, 108, 168, 192]
    # Markersize values
    y = [40, 12, 6, 1]

    # Create the fit
    pcoeffs, _ = curve_fit(exponential, x, y, p0=(1, 1e-6, 1))

    # Now derive the value (don't return a value less than 1)
    return max(round(exponential(c_value, *pcoeffs)), 1)

def do_plot(datapath, plotfield, plotpath='.'):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    lfric = load_cube_by_varname(datapath, plotfield[varname])
    lfric = lfric[-1, 0]

    # Get the x and y co-ordinates
    x_coord = np.around(lfric.coord('longitude').points, decimals=5)
    y_coord = np.around(lfric.coord('latitude').points,  decimals=5)

    # Save the min and max of the data
    field_min = np.around(np.min(lfric.data), decimals=5)
    field_max = np.around(np.max(lfric.data), decimals=5)

    # Set up the colourbar
    plt.set_cmap(plt.cm.RdYlBu_r)

    plt.figure(figsize=(8, 5))
    markersize = gen_markersize(lfric.data.shape[0])
    plt.scatter(x_coord, y_coord, c = lfric.data,
                edgecolor = "none", s = markersize,
                vmin = plotfield[colbar_min],
                vmax = plotfield[colbar_max])

    plt.title(plotfield[varname]+', min = '+str(field_min)
                                +', max = '+str(field_max) )
    plt.xlim([np.min(x_coord), np.max(x_coord)])
    plt.ylim([-90, 90])

    plt.savefig(plotpath+'/'+plotfield[varname]+'.png', bbox_inches='tight')


if __name__ == "__main__":

    import sys
    try:
        datapath, plotpath = sys.argv[1:3]
    except ValueError:
        print("Usage: {0} <datapath> <plotpath>".format(sys.argv[0]))
        exit(1)
    do_plot(datapath, theta,           plotpath)
    do_plot(datapath, m_v,             plotpath)
    do_plot(datapath, m_cl,            plotpath)
    do_plot(datapath, m_ci,            plotpath)
    do_plot(datapath, u1,              plotpath)
    do_plot(datapath, sw_heating_rate, plotpath)
