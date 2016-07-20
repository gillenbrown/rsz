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

ch1_m_ch2 = dict()
ch1_m_ch2["color"] = "ch1-ch2"
ch1_m_ch2["blue_band"] = "ch1"
ch1_m_ch2["red_band"] = "ch2"

ch1_m_ch2["z_min"] = decimal.Decimal("0.7")
ch1_m_ch2["z_max"] = decimal.Decimal("1.7")
ch1_m_ch2["correction"] = [-0.17985356,  1.1423761]
ch1_m_ch2["slope_fit"] = [0, 0]

# plotting info
# mag_min, mag_max, color_min, color_max
ch1_m_ch2["plot_lims"] = [18, 22, -1, 0.5]

# redshift fitting
# ----------------
# initial redshift fitting mag and color cuts
# mag: how many mags brighter, then fainter than the M* point they can be
# color: how much bluer, then redder it can be
ch1_m_ch2["initial_mag"] = [2.0, 0.0]
ch1_m_ch2["initial_color"] = [0.1, 0.1]
ch1_m_ch2["bluer_color_cut"] = [0.2, 0.1]
ch1_m_ch2["redder_color_cut"] = [0.2, 0.1]
ch1_m_ch2["brighter_mag_cut"] = 2.5
ch1_m_ch2["dimmer_mag_cut"] = 0
ch1_m_ch2["final_rs_mag"] = [2.0, 0.6]  # brighter, fainter
ch1_m_ch2["final_rs_color"] = [0.15, 0.15]  # bluer, redder



sloan_r_m_sloan_z = dict()
sloan_r_m_sloan_z["color"] = "sloan_r-sloan_z"
sloan_r_m_sloan_z["blue_band"] = "sloan_r"
sloan_r_m_sloan_z["red_band"] = "sloan_z"

sloan_r_m_sloan_z["z_min"] = decimal.Decimal("0.5")
sloan_r_m_sloan_z["z_max"] = decimal.Decimal("1.5")
sloan_r_m_sloan_z["correction"] = [0.01705775352432836, 1.0834470213733527]
sloan_r_m_sloan_z["slope_fit"] = [-0.00343316, -0.14489063]

# plotting info
# mag_min, mag_max, color_min, color_max
sloan_r_m_sloan_z["plot_lims"] = [20, 23.5, 0, 3.5]

# redshift fitting
# ----------------
# initial redshift fitting mag and color cuts
# mag: how many mags brighter, then fainter than the M* point they can be
# color: how much bluer, then redder it can be
sloan_r_m_sloan_z["initial_mag"] = [2.0, 0.6]
sloan_r_m_sloan_z["initial_color"] = [0.2, 0.2]
sloan_r_m_sloan_z["bluer_color_cut"] = [0.25, 0.225]
sloan_r_m_sloan_z["redder_color_cut"] = [0.4, 0.3]
sloan_r_m_sloan_z["brighter_mag_cut"] = 1.4
sloan_r_m_sloan_z["dimmer_mag_cut"] = 0.6
sloan_r_m_sloan_z["final_rs_mag"] = [2.0, 0.6]  # brighter, fainter
sloan_r_m_sloan_z["final_rs_color"] = [0.35, 0.35]  # bluer, redder

cfg_matches = {"ch1-ch2": ch1_m_ch2,
               "sloan_r-sloan_z": sloan_r_m_sloan_z}

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
