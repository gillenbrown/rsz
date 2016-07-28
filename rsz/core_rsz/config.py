import decimal

# ------------------------------- READ THIS -----------------------------------
#
# This file contains the information the code uses to fit the red sequence
# in different band combinations. To add a new band, create a new dictionary
# that holds all the same keys as the ones that already exist. Then add it to
# the `config_matches` dictionary. The key is the string with the name of the
# bands that are used, and the value is the congifuration dictionary containing
# all the parameters the code needs.
#
# This file also stores AB to Vega conversions, so if your filters aren't
# already in the ab_to_vega dictionary at the bottom, add them.
#
# --------------------------- Param documentation -----------------------------
#
# This section describes the parameters the code uses to fit the red sequence.
#
# color: Color being fitted. This needs to be of the format
#        "first_band-second_band" with no spaces. The name of the bands
#        need to match the names of the filters in the ezgal/data/filters
#        directory, otherwise the code won't be able to fine the filters.
#        Read http://www.baryons.org/ezgal/filters.php to see all the filters
#        ezgal has. You can also add your own filters if need be. Read
#        http://www.baryons.org/ezgal/manual/file_formats.html#filter-file-format
#        for more info.
# blue_band: Bluer of the two bands used above.
# red_band: Redder of the two bands used above.
#
# z_min: minimum possible redshift to be fit. This needs to be a decimal
#        object, since the code uses that format under the hood to avoid
#        floating point errors. Enter the redshift into the function as a
#        string, to avoid floating point errors when it's being created.
# z_max: maximum possible redshift to be fit. Sane format as z_min.
#    Note: Choose these two to span a region where color is monotonic as a
#    function of redshift.
# correction:  This is the fit, which comes from
#              numpy.polynomial.polynomial.polyfit
#              It is the coefficients of the polynomial, starting with the
#              lowest power of z. It takes uncorrected redshifts and turns
#              them into calibrated redshifts. To find this parameter, plot
#              rsz redshifts against well-measured redshifts, and find the
#              function that fits. This is what goes here.
# slope_fit: This is the fit (same as `correction`) that describes the slope
#            of the red sequence as a function of redshift. As an initial
#            guess, zero slope and intercept will normally be fine. To do this
#            properly, you need to figure out what the slope of the red
#            sequence is as a function of redshift, then plug that info in
#            here.
#
# plot_lims: In the plots we create, we need to know the limits. This parameter
#            is a list of 4 items: The minimum red magnitude, the maximum
#            red magnitude, the minimum color, and the maximum color. These
#            are all set in AB mags. The best way to choose these it to set
#            them to something very large, see where the points are, and hone
#            in from there.
#
# ------------  Redshift Fitting Parameters  -----------
#
# These parameters are very important, as they describe how the red sequence
# is determined. To do that properly, we need some background on how the code
# works. First, it counts up the number of galaxies near each RS model, and
# picks the model with the most nearby galaxies as an initial guess.
# It then selects those nearby galaxies as potential RS members, and performs
# chi-squared fitting on those to pick the model that best fits them. This is
# done twice, as an iterative process, to ensure that things converge.
# After that is done, we select the final red sequence members as those near
# the model selected as the redshift for the cluster.
#
# initial_mag: Magnitude cuts used in the initial guess fitting process. This
#              is a list of the number of magnitudes brighter and fainter than
#              the characteristic magnitude. Galaxies will pass this cut if
#              they have a magnitude m such that:
#              model_m* - initial_mag[0] < m < model_m* + initial_mag[1]
#              Galaxies that pass this cut will be counted. You definitely want
#              to choose a large bright cut, since the bright galaxies are the
#              most important to determining the red sequence. The code
#              performs a error cut, so the faint end cut isn't as important,
#              but choosing a lower value (less faint) can often increase the
#              signal-to-noise.
# initial_color: How much bluer, then redder, than the RS model that galaxies
#                can be to be counted during the initial fitting process.
#                Galaxies will pass this cut if they have a color c such that:
#                model_c - initial_color[0] < c < model_c + initial_color[1]
#                Note that both entries in the list need to be positive due to
#                how they are defined. Small numbers (0.1, 0.2) are best,
#                since they provide the tightest restriction.
# bluer_color_cut: During the main fitting process, the code does color cuts
#                  on each side of the potential red sequence models. These
#                  are done from the characteristic color of the RS models.
#                  This parameter specifies the blue side of those cuts. This
#                  needs to be a list, where the cuts are for each successive
#                  iteration. Start with a wide value so that the true RS will
#                  be included from the very beginning, even if the initial
#                  guess isn't very good. Don't make it too big, though, since
#                  including lots of blue galaxies (that often have smaller
#                  error bars) can drag the fit away from the true RS. The
#                  final value should be roughly the size of the real RS. To
#                  tune these parameters, enable the `fitting_procedure`
#                  plot in the parameter file for your run, and see whether the
#                  color cuts are working properly.
# redder_color_cut: This functions the same as `bluer_color_cut`, just on the
#                   red side of the RS. The advice there applies here too.
#                   That said, the red cut can often be more forgiving than the
#                   blue cut, since there are less foreground objects on the
#                   red side to drag things away from the true RS.
# brighter_mag_cut: Same as the previous two parameters, except this is a plain
#                   magnitude cut (based on distance from the characteristic
#                   magnitude of the RS model) This needs to be large
#                   enough to include all the galaxies in the RS.
# dimmer_mag_cut: Same as brighter_mag_cut, just on the dim side. Often a small
#                 value for this is good, since throwing away the faint
#                 galaxies improves your signal to noise. The RS is most
#                 visible on the bright end, so going too faint here will make
#                 it harder for the code to distinguish the real red sequence.
#  note on all of these: For a galaxy to be selected as a red sequence member,
#    it must have a mag m such that
#        model_m* - brighter_mag_cut < m < model_m* + dimmer_mag cut
#    then for each successsive color iterations, the color c must be
#        model_c - bluer_color_cut[i] < c < model_c + redder_color_cut[i]
#    Only galaxies that pass both of these cuts will be called red squence
#    galaies.
#
# -----------------------------------------------------------------------------


# IRAC ch1 - ch2
ch1_m_ch2 = dict()
ch1_m_ch2["color"] = "ch1-ch2"
ch1_m_ch2["blue_band"] = "ch1"
ch1_m_ch2["red_band"] = "ch2"

ch1_m_ch2["z_min"] = decimal.Decimal("0.7")
ch1_m_ch2["z_max"] = decimal.Decimal("1.7")
ch1_m_ch2["correction"] = [-0.17985356,  1.1423761]
ch1_m_ch2["slope_fit"] = [0, 0]

ch1_m_ch2["plot_lims"] = [18, 22, -1, 0.5]

ch1_m_ch2["initial_mag"] = [2.0, 0.0]
ch1_m_ch2["initial_color"] = [0.1, 0.1]
ch1_m_ch2["bluer_color_cut"] = [0.2, 0.1]
ch1_m_ch2["redder_color_cut"] = [0.2, 0.1]
ch1_m_ch2["brighter_mag_cut"] = 2.5
ch1_m_ch2["dimmer_mag_cut"] = 0
ch1_m_ch2["final_rs_mag"] = [2.0, 0.6]
ch1_m_ch2["final_rs_color"] = [0.15, 0.15]

# SDSS r - z
sloan_r_m_sloan_z = dict()
sloan_r_m_sloan_z["color"] = "sloan_r-sloan_z"
sloan_r_m_sloan_z["blue_band"] = "sloan_r"
sloan_r_m_sloan_z["red_band"] = "sloan_z"

sloan_r_m_sloan_z["z_min"] = decimal.Decimal("0.5")
sloan_r_m_sloan_z["z_max"] = decimal.Decimal("1.5")
sloan_r_m_sloan_z["correction"] = [0.01705775352432836, 1.0834470213733527]
sloan_r_m_sloan_z["slope_fit"] = [-0.00343316, -0.14489063]

sloan_r_m_sloan_z["plot_lims"] = [20, 23.5, 0, 3.5]

sloan_r_m_sloan_z["initial_mag"] = [2.0, 0.6]
sloan_r_m_sloan_z["initial_color"] = [0.2, 0.2]
sloan_r_m_sloan_z["bluer_color_cut"] = [0.25, 0.225]
sloan_r_m_sloan_z["redder_color_cut"] = [0.4, 0.3]
sloan_r_m_sloan_z["brighter_mag_cut"] = 1.4
sloan_r_m_sloan_z["dimmer_mag_cut"] = 0.6
sloan_r_m_sloan_z["final_rs_mag"] = [2.0, 0.6]
sloan_r_m_sloan_z["final_rs_color"] = [0.35, 0.35]

# ----------------- ADD NEW COLOR COMBO DICTS HERE ----------------------------
cfg_matches = {"ch1-ch2": ch1_m_ch2,
               "sloan_r-sloan_z": sloan_r_m_sloan_z}
# -----------------------------------------------------------------------------

# store Vega to AB conversions. This stores factor, such that
# Vega_mag = AB_mag + factor, or AB_mag = Vega_mag - factor
# These were obtained from http://www.baryons.org/ezgal/filters.php
ab_to_vega = {"ch1": -2.787,
              "ch2": -3.260,
              "sloan_u": -0.904,
              "sloan_g": 0.098,
              "sloan_r": -0.146,
              "sloan_i": -0.357,
              "sloan_z": -0.521}

# ------------------ Validating this file  ------------------------------------
#
# Don't add anything down here, the code here just checks that you filled
# things up properly.

all_keys = ["color", "blue_band", "red_band", "z_min", "z_max", "correction",
            "slope_fit", "plot_lims", "initial_mag", "initial_color",
            "bluer_color_cut", "redder_color_cut", "brighter_mag_cut",
            "dimmer_mag_cut", "final_rs_mag", "final_rs_color"]

for color, cfg in cfg_matches.items():
    # check that all the keys are there
    for key in all_keys:
        if key not in cfg:
            raise ValueError("Please add the {} parameter to the {} \n"
                             "\tconfiguration dictionary in config.py.\n"
                             "\n".format(key, color))

    # then check that there aren't any keys that shouldn't be there
    for key in cfg:
        if key not in all_keys:
            raise ValueError("The parameter {} that you added to the {}\n"
                             "\tconfiguration dictionary is not needed.\n"
                             "Please remove it.\n".format(key, color))

    # then validate the other parameters.

    # the color should match
    if color != cfg["color"]:
        raise ValueError("The color you specified in the {} dictionary\n"
                         "\tdoes not match the name you gave it in the \n"
                         "\t`cfg_matches` dictionary, which is {}. \n"
                         "\tPlease fix this.\n".format(cfg["color"], color))

    # the bands should match the name of the color
    if color != "-".join([cfg["blue_band"], cfg["red_band"]]):
        raise ValueError("The name of the bands you specified in the\n"
                         "\t`red_band` and `blue_band` parameters of the {}\n"
                         "\tdictionary doesn't match the color you gave: {}.\n"
                         "\tPlease fix this. The color should have the form:\n"
                         "\tblue_band-red_band.".format(color, color))

    # these bands should have AB-Vega conversions
    for band in [cfg["red_band"], cfg["blue_band"]]:
        if band not in ab_to_vega:
            raise ValueError("Please add {} to the `ab_to_vega` dictionary\n"
                             "\tin config.py.".format(band))

    # the z_max and min need to be in decimal format
    for z in ["z_max", "z_min"]:
        if type(cfg[z]) != decimal.Decimal:
            raise ValueError("The {} paramter in the {} config dictionary\n"
                             "\tneeds to be of type decimal. See the other\n"
                             "\timplemented ones for an example."
                             "\n".format(z, color))

    # the correction and slope need to be a list with nonzero entries
    for key in ["correction", "slope_fit",
                "bluer_color_cut", "redder_color_cut"]:
        if type(cfg[key]) != list or len(cfg[key]) < 1:
            raise ValueError("The {} parameter in the {} configuration\n"
                             "\tdictionary needs to be a list with at\n"
                             "\tleast one entry.".format(key, color))

    # the plot_lims needs to be a list of 4 values.
    if type(cfg["plot_lims"]) != list or len(cfg["plot_lims"]) != 4:
        raise ValueError("The parameter plot_lims in the {} configuration\n"
                         "\tdictionary needs to be a list with 4 items."
                         "".format(color))

    # the initial mag, initial color, final mag, and final color are all
    # lists with two values.
    for key in ["initial_mag", "initial_color",
                "final_rs_mag", "final_rs_color"]:
        if type(cfg[key]) != list or len(cfg[key]) != 2:
            raise ValueError("The parameter {} in the {} configuration\n"
                             "\tdictionary needs to be a list with 2 items,\n"
                             "\tas described in config.py.".format(key, color))

    for key in ["brighter_mag_cut", "dimmer_mag_cut"]:
        if type(cfg[key]) not in [float, int]:
            raise ValueError("The parameter {} in the configuration dict\n"
                             "\tneeds to be a single value (ie not\n"
                             "\ta list, even if it only has one item)."
                             "".format(key))

    # the bluer and redder color cuts need to have the same length
    if len(cfg["bluer_color_cut"]) != len(cfg["redder_color_cut"]):
        raise ValueError("The length of the `bluer_color_cut` and\n"
                         "\t`redder_color_cut` lists need to be the same\n"
                         "\tin the {} dictionary in config.py".format(color))




