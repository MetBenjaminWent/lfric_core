#!/usr/bin/env python
##############################################################################
# (C) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

''' Plot the contents of IO_Dev diagnostic files - give the diagnostic file
    path as the argument and edit the list of plot_fields below to choose
    which fields will be plotted '''

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris

matplotlib.rc('font', **{'family': "Arial",
                         'size': 12})


def load_cube_by_varname(filename, var):
    '''Loads an iris cube from a file based on the varname'''
    variable_constraint = iris.Constraint(
        cube_func=(lambda c: c.var_name == var))
    return iris.load_cube(filename, constraint=variable_constraint)


def plot_last_domain(filename, varname, ax):
    '''Plots the last domain of a given field as a 2D scatter'''

    cube_out = load_cube_by_varname(filename, varname)

    # Strip out all but the last domain
    data = cube_out.data
    if len(np.shape(data)) == 2:
        data = data[-1]
    elif len(np.shape(data)) == 3:
        data = data[-1][-1]

    long = cube_out.coord('longitude').points
    lat = cube_out.coord('latitude').points

    ax.text(0.9, 0.9, varname, transform=ax.transAxes,
        verticalalignment ='top', horizontalalignment ='right',)

    ax.scatter(long, lat, c=data, cmap='gist_earth')
    ax.set_xlabel(r"X [$^\circ$]")
    ax.set_ylabel(r"Y [$^\circ$]", labelpad=-4)
    print("{0} domain plotted...".format(varname))


def plot_time_mean(filename, varname, ax):
    ''' Plots the time-series of the field mean'''

    cube_out = load_cube_by_varname(filename, varname)

    # Compute average of each time entry
    if len(np.shape(cube_out.data)) == 1:
        data = cube_out.data
    elif len(np.shape(cube_out.data)) == 2:
        data = np.nanmean(np.where(
            cube_out.data > -1e300, cube_out.data, np.nan), axis=1)
    elif len(np.shape(cube_out.data)) == 3:
        data = np.nanmean(np.where(
            cube_out.data > -1e300, cube_out.data, np.nan), axis=(1, 2))

    time = cube_out.coord('time').points
    ax.text(0.05, 0.85, varname, transform=ax.transAxes)

    ax.scatter(time, data)
    ax.set_ylim(0.8*min(data), 1.2*max(data))
    ax.set_xlabel(r"Time [s]")
    ax.set_ylabel(r"Field mean", labelpad=0)
    print("{0} time-series plotted...".format(varname))


def create_plots(config, datapath, plot_fields):
    '''Creates a series of plots based on an input config string and a list of
       fields to plot'''

    if config == "domains":

        xdim = int(np.ceil(np.sqrt(len(plot_fields))))
        ydim = int(np.ceil(len(plot_fields)/xdim))
        fig, axs = plt.subplots(ydim, xdim, figsize=(4*xdim, 2.5*ydim),
                                sharex=True, sharey=True)

        if ydim*xdim > 1:
            axs = axs.ravel()

        for i in range(len(plot_fields)):
            field = plot_fields[i]
            plot_last_domain(datapath, field, axs[i])

        remove_axes = np.arange(len(plot_fields)-len(axs), 0)
        for a in remove_axes:
            axs[a].set_axis_off()

        fig.suptitle("/".join(datapath.split("/")[-5::]))

    elif config == "time_var":

        ydim = 2
        xdim = len(plot_fields)

        fig, axs = plt.subplots(ydim, xdim, figsize=(4*xdim, 2.5*ydim))

        axs = axs.ravel()

        for i in range(len(plot_fields)):
            field = plot_fields[i]
            plot_time_mean(datapath, field, axs[i])
            plot_last_domain(datapath, field, axs[i+len(plot_fields)])

        fig.suptitle("/".join(datapath.split("/")[-5::]))

    else:
        print("Invalid config key '{0}' for IO_Dev plotting "
              "script.".format(config))
        sys.exit(1)

    plt.tight_layout(pad=2.0, h_pad=0.9, w_pad=0.9)

    output_path = datapath.replace(".nc", "_{0}_plot.png".format(config))
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print("Resulting plots saved in {0}".format(output_path))


if __name__ == "__main__":

    try:
        DATAPATH = sys.argv[1]
        CONFIGS = sys.argv[2].split(':')
        PLOT_FIELDS = sys.argv[3].split(':')
    except ValueError:
        print("Usage: {0} <datapath> <list:of:configs> "
              "<list:of:fields>".format(sys.argv[0]))
        sys.exit(1)

    for config in CONFIGS:
        create_plots(config, DATAPATH, PLOT_FIELDS)
