fitting_combos = ["ch1-ch2", "sloan_r-sloan_z"]

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
ch1_m_ch2["z_min"] = 0.7
ch1_m_ch2["z_max"] = 1.7
ch1_m_ch2["correction"] = [-0.17985356,  1.1423761]
ch1_m_ch2["slope_fit"] = [0, 0]


sloan_r_m_sloan_z = dict()
sloan_r_m_sloan_z["z_min"] = 0.5
sloan_r_m_sloan_z["z_max"] = 1.5
sloan_r_m_sloan_z["correction"] = [0.01705775352432836, 1.0834470213733527]
sloan_r_m_sloan_z["slope_fit"] = [-0.00343316, -0.14489063]

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
