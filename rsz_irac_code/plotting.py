import matplotlib.pyplot as plt
import matplotlib.colors as mplcol
import matplotlib.cm as cmx
import numpy as np

import model
import data


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
    fig, ax = plt.subplots(figsize=(9,6))

    # plot points one by one, so I don't have to make long lists for each
    # thing I want to plot. This is simpler
    for source in valid_sources:
        mag = source.ch2.value
        color = source.ch1_m_ch2
        # red sequence members will be colored red, while non RS galaxies
        # will be colored black.
        if source.RS_member:
            point_color = "r"
        else:
            point_color = "k"

        ax.errorbar(x=mag, y=color.value, yerr=color.error, c=point_color,
                    fmt=".", elinewidth=0.35, capsize=0, markersize=5)

    # label and clean up the axes
    ax.set_xlim(18, 23)
    ax.set_ylim(-1, 0.5)
    ax.set_xlabel("ch2  [AB]")
    ax.set_ylabel("ch1 - ch2  [AB]")
    ax.text(0.02, 0.96, cluster.name.replace("_", " "), transform=ax.transAxes,
            horizontalalignment="left", verticalalignment="center",
            bbox=dict(facecolor="w", linewidth=0.0))

    # return both the figure and the axis so that other functions can add
    # more cool stuff to the plot.
    return fig, ax




def add_vega_labels(ax):
    """

    :param ax:
    :return:
    """
     # move ticks left and bottom, since we want Vega ticks on right
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    # put Vega mags/colors on the top and right
    # ch1: Vega - AB = -2.787
    ab_v_1 = -2.787
    # ch2: Vega - AB = -3.260
    ab_v_2 = -3.260
    # to do this, we need to make a second axis
    vega_color_ax = ax.twinx()
    vega_mag_ax = ax.twiny()
    vega_mag_ax.set_xlim(18 + ab_v_2, 23 + ab_v_2)
    vega_color_ax.set_ylim(-1 + (ab_v_1 - ab_v_2), 0.5 + (ab_v_1 - ab_v_2))
    vega_mag_ax.set_xlabel("ch2  [Vega]")
    vega_mag_ax.xaxis.set_label_position("top")
    vega_mag_ax.xaxis.tick_top()
    vega_color_ax.set_ylabel("ch1 - ch2  [Vega]")
    vega_color_ax.yaxis.set_label_position("right")
    vega_color_ax.yaxis.tick_right()

    return vega_color_ax, vega_mag_ax


def add_all_models(fig, ax, steal_axs):
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
    spectral = plt.get_cmap("RdYlBu_r")
    # some other decent colormaps: coolwarm, bwr, jet, RdBu_r, RdYlBu_r,
    # Spectral_r, rainbow
    #normalize the colormap
    c_norm = mplcol.Normalize(vmin=np.float64(min(models)),
                              vmax=np.float64(max(models)))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=spectral)

    for z in models:
        # use the scalar map to get the color for this redshift.
        color_val = scalar_map.to_rgba(z)
        # plot the line for this model
        add_one_model(ax, models[z], color_val)

    # add the colorbar
    scalar_map.set_array([]) # Don't know what this does
    cb = fig.colorbar(mappable=scalar_map, ax=steal_axs, pad=0.15,
                      fraction=0.05, anchor=(1.0, 1.0))
    cb.set_label("Redshift")
    cb.set_ticks(np.arange(0, 3, 0.1))


def add_one_model(ax, model, color):
    """
    Adds one model to the axis.

    Ax should be generated with the cmd function.

    :param ax: axis on which to plot the model
    :param model: model to plot on the axis
    :return: none, but ax will be modified
    """
    # make two points in the line for this RS
    ch2s = [10, 30]
    colors = [model.rs_color(ch2) for ch2 in ch2s]
    # plot that line
    ax.plot(ch2s, colors, color=color, linewidth=0.5, zorder=0)

    # also plot the points that corresponding to L*
    ax.scatter(model.mag_point, model.color_point, c=color,
               s=30, linewidth=0, zorder=0)


def add_redshift(ax, redshift):
    """Label the plot with a redshift in the upper right corner.

    Can be either just a value or one with errors, too.

    :param ax: axis to put the redshift of.
    :param redshift: redshift to be put in the corner. Can be either a
                    single float of a data object.
    :return: None, but the redshift will be added in a box in the upper
    right corner of ax.
    """
    if type(redshift) is data.AsymmetricData:
        # we need to enclose the value in braces so LaTeX can recognize them
        #  properly
        value = "{" + str(round(redshift.value, 2)) + "}"
        upper = "{+" + str(round(redshift.upper_error, 2)) + "}"
        lower = "{\,-" + str(round(redshift.lower_error, 2)) + "}"
        text = r"z=$\mathregular{{{value}^{upper}_{lower}}}$".format(
            value=value, upper=upper, lower=lower)
    else:
        text = "z=" + str(round(redshift, 2))


    ax.text(0.97, 0.96, text, transform=ax.transAxes,
            horizontalalignment="right", verticalalignment="center",
            bbox=dict(facecolor="w", linewidth=0.0))


def location(cluster):
    """Plot the location of all the sources in the cluster, with RS members
    highlighted.

    All sources are plotted. Ones that passed the location cut will be
    darker than those that didn't. Red sequence galaxies will be red
    regardless of whether or not they passed the location cut.

    :param cluster: cluster to plot the data for
    :return: fig, ax of the plot
    """

    fig, ax = plt.subplots(figsize=(9,8), tight_layout=True)
    rs_member_ra, rs_member_dec = [], []
    center_ra, center_dec = [], []
    rest_ra, rest_dec = [], []
    for source in cluster.sources_list:
        if source.RS_member:
            rs_member_ra.append(source.ra)
            rs_member_dec.append(source.dec)
        elif source.near_center:
            center_ra.append(source.ra)
            center_dec.append(source.dec)
        else:
            rest_ra.append(source.ra)
            rest_dec.append(source.dec)
    ax.scatter(rest_ra, rest_dec, c="0.7", linewidth=0)
    ax.scatter(center_ra, center_dec, c="k", linewidth=0,
               label="Location Cut")
    ax.scatter(rs_member_ra, rs_member_dec, c="r", linewidth=0,
               label="Red Sequence")

    # add labels and clean up the plot
    ax.set_xlabel("ra")
    ax.set_ylabel("dec")
    legend = ax.legend(loc=3)
    legend.get_frame().set_linewidth(0.5)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    # label the cluster name
    ax.text(0.02, 0.96, cluster.name.replace("_", " "), transform=ax.transAxes,
            horizontalalignment="left", verticalalignment="center",
            bbox=dict(facecolor="w", linewidth=0.0))
    ax.invert_xaxis()  # ra is backwards
    ax.set_aspect("equal", adjustable="box")  # we want ra and dec to be
    # scaled the same
    return fig, ax
