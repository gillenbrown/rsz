import model

import matplotlib.pyplot as plt
import matplotlib.colors as mplcol
import matplotlib.cm as cmx
import numpy as np

def cmd(cluster):
    """
    Plot a color-magnitude diagram for the cluster. Red sequence galaxies
    will be highlighted in red.

    :param cluster: cluster containing the galaxies to be plotted in the CMD.
    :return: figure AND axis objects for the plot
    """

    # throw out a few sources that don't need to be plotted.
    valid_sources = [source for source in cluster.sources_list if (
            source.near_center and (source.ch1 - source.ch2).error < 5.0)]

    # set up the plot
    fig, ax = plt.subplots(figsize=(9,6), tight_layout=True)

    # plot points one by one, so I don't have to make long lists for each
    # thing I want to plot. This is simpler
    for source in valid_sources:
        mag = source.ch2.value
        color = source.ch1 - source.ch2

        # red sequence members will be colored red, while non RS galaxies
        # will be colored black.
        if source.RS_member:
            point_color = "r"
        else:
            point_color = "k"

        ax.errorbar(x=mag, y=color.value, yerr=color.error, c=point_color)

    # label and clean up the axes
    ax.set_xlim(18, 23)
    ax.set_ylim(-1, 0.5)
    ax.set_xlabel("ch2")
    ax.set_ylabel("ch1 - ch2")

    # return both the figure and the axis so that other functions can add
    # more cool stuff to the plot.
    return fig, ax

def add_models(fig, ax):
    """
    Adds the RS models to the given axis, adding a colorbar to code redshift.

    The fig and ax should be obtained by calling cmd(). The models will be
    meaningless on any other plot.

    :param fig: Figure containing the plots.
    :param ax: CMD axis the models will be drawn onto. The colorbar will
               also steal room from this axis.
    :return: none, but fig and ax are modified by adding the models and
            colorbar.
    """

    # get the model predictions
    models = model.model_dict(0.05)

    # set the colormap, so we can color code lines by redshift
    spectral = plt.get_cmap("spectral")
    #normalize the colormap
    c_norm = mplcol.Normalize(vmin=np.float64(min(models)),
                              vmax=np.float64(max(models)))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=spectral)

    for z in models:
        # use the scalar map to get the color for this redshift.
        color_val = scalar_map.to_rgba(z)

        # make two points in the line for this RS at this redshfit
        ch2s = [10, 30]
        colors = [models[z].rs_color(ch2) for ch2 in ch2s]
        # plot that line
        ax.plot(ch2s, colors, color=color_val, linewidth=0.5, zorder=0)

        # also plot the points that corresponding to L*
        ax.scatter(models[z].mag_point, models[z].color_point, c=color_val,
                   s=30, linewidth=0, zorder=0)

    # add the colorbar
    scalar_map.set_array([]) # Don't know what this does
    cb = fig.colorbar(mappable=scalar_map, ax=ax)
    cb.set_label("Redshift")
    cb.set_ticks(np.arange(0, 3, 0.1))
