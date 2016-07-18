fitting_combos = ["ch1-ch2", "sloan_r-sloan_z"]

import decimal

# Param documentation:
# z_min: minimum possible redshift to be fit.
# z_max: maximum possible redshift to be fit.
#        Note: Choose these two to span a region where color is monotonic
# correction:  This is the fit, which comes from
#              numpy.polynomial.polynomial.polyfit
#              It is the coefficients of the polynomial, starting with the
#              lowest power of z. It takes uncorrected redshifts and turns
#              them into calibrated redshifts. To find this parameter, plot
#              rsz redshifts against well-measured redshifts, and find the
#              function that fits. This is what goes here.

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

cfg_matches = {"ch1-ch2": ch1_m_ch2,
               "sloan_r-sloan_z": sloan_r_m_sloan_z}

# store Vega to AB conversions. This stores factor, such that
# AB_mag = Vega_mag + factor
vega_to_ab = {"ch1": 2.788,  # http://irsa.ipac.caltech.edu/data/COSMOS/tables/scosmos/scosmos_irac_200706_colDescriptions.html
              "ch2": 3.255,  # same
              "sloan_u": 0.91, # http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
              "sloan_g": -0.08, # same
              "sloan_r": 0.16, # same
              "sloan_i": 0.37, # same
              "sloan_z": 0.54} # same
