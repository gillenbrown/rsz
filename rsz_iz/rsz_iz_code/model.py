import os
import decimal

import numpy as np
import ezgal


class _Slope(object):
    """
    Callable class used to calculate slopes of the red sequence.

    There is probably an easier way to make all this work, but I couldn't
    figure out the best way. Suggestions are definitely welcome.

    It currently returns zero as the slope no matter what, but this can
    easily be changed. To do this, have the slope(z) function return the
    slope given an actual redshift. Everything will work from there.
    """

    def __init__(self):
        """
        Find the slope of the red sequence as a function of redshfit.

        Uses the slopes of the red sequence as published in Eisenhardt 2007.

        :return: a function that returns the slope at a given redshift
        """

        # first we need to see at what redshift does r-z see the various
        # Eisenhardt colors.

        # Then do some interpolating here to find the slope at any redshift

        def slope(z):

            return 0  # just for now
            # # fit is a polynomial (highest power first) that describes the
            # # slope of the red sequence as a function of redshift
            # fit = [-0.14489063, -0.00343316]
            # # The work to get that is in the iPython notebook. It needs
            # # work to make it prettier and more satisfying
            #
            # # turn to decimal , so it can play nice with redshifts, which are
            # # of type decimal
            # fit = [decimal.Decimal(i) for i in fit]
            # # plug the redshift into the polynomial. This looks ugly, but
            # # that's all it's doing.
            # return float(sum([coeff * (z**i) for i, coeff in
            #                   enumerate(reversed(fit))]))

        self.slope_function = slope

    # I'll use a call method so I can call the object like a function to get
    #  the slope of the red sequence at the desired redshift.
    def __call__(self, redshift):
        return self.slope_function(redshift)


class RSModel(object):
    """
    Class storing data from the EzGal models.
    """

    # get an object that can return the slope at any redshift
    slope = _Slope()

    @staticmethod
    def correction(redshift):
        """
        Does the redshift correction.

        After comparing clusters with known
        redshifts to those produced by this program, this correction was
        devised to make this program return results more closely matching
        the actual redshift. It's just an empirical correction to make the
        redshifts better.

        This correction was found by running the code on SPT clusters that
        have spec or photo zs. The results were then compared, and a function
        was found that makes the redshifts line up the best.

        :param redshift: uncorrected redshift
        :returr: corrected redshift
        """
        return redshift  # just for now
        # # This is the fit, which is described by a linear relationship.
        # #  The first item is the slope, the second is the intercept.
        # #  If there were more items, they would go in descending power order.
        # fit = [1.0834470213733527, 0.01705775352432836]
        # # turn to decimal , so it can play nice with redshifts, which are
        # # of type decimal
        # fit = [decimal.Decimal(i) for i in fit]
        # # Applying the correction is just plugging in the redshift to the
        # # correction polynomial. This is a complicated way to do that.
        # corrected = sum([coeff * (redshift**i) for i, coeff in
        #                  enumerate(reversed(fit))])
        # # round and turn into type decimal, which is what redshifts need to be.
        # return decimal.Decimal(str(round(corrected, 3)))

    def __init__(self, redshift, i_mag, z_mag):
        """
        Initializes the useful parameters that will be used to calculate
        where the red sequence is.

        :param redshift: redshift at which these points are calculated
        :param r_mag: r_mag at that redshift
        :param z_mag: z_mag at that redshift
        :return: Model object.
        """

        self.mag_point = z_mag
        self.color_point = i_mag - z_mag
        self.z = RSModel.correction(redshift)
        self._slope = RSModel.slope(redshift)

    def rs_color(self, z_mag):
        """
        Calculate the color (r-z) of the red sequence at the given z
        magnitude.

        :param z_mag: z magnitude at which we want the color of the
                      red sequence
        :return: float value with the color of the red sequence

        Algorithm: Start with point slope form of a line
        y - y1 = m(x - x1)
        y = y1 + m(x - x1)
        Where x is z magnitude, y is r-z color, m is the slope of the
        red sequence, and x1 and y1 are the zeropoint that was returned by
        EzGal.
        """

        return self.color_point + self._slope * (z_mag - self.mag_point)


def model_dict(spacing):
    """
    Create a dictionary of model objects, that represent the red sequence at
    different redshifts.

    :param spacing: spacing of the redshifts at which the RS will be
    calculated.
    :return: dictionary where the keys are the redshifts of the models, and
             the values are the models themselves.
    """

    # get the EzGal model object
    model = _make_model()

    # decide the formation redshift and observed redshift
    zf = 3.0
    zs = np.arange(0.3, 1.900001, spacing)

    # normalize to Coma
    model.set_normalization(filter='ks', mag=10.9, apparent=True, vega=True,
                            z=0.023)

    # Calculate the observables we want in AB mags
    mags = model.get_apparent_mags(zf, filters=["sloan_z", "ch1"], zs=zs,
                                   vega=True)

    # initialize an empty dictionary that will be filled with RSModel objects
    rs_models = dict()

    # turn redshifts into decimal objects, to avoid floating point errors.
    # turn to strings first, then add zeros when necessary, then turn to
    # decimal objects. We add zeroes to get 1.10, rather than 1.1, for example.
    # Decimal objects are used because they are not subject to floating
    # point errors, and still work in simple algebra. Having floating point
    # errors in dictionary keys is a bad thing.

    # We don't do the correction here, so we only need to round to 2 digits.
    # the correction is done later, and after it we will wants 3 digits.
    zs = [str(round(z, 2)) for z in zs]
    z_2_digits = []
    for z in zs:
        if z[-2] == ".":
            z += "0"
        z_2_digits.append(z)
    decimal_zs = [decimal.Decimal(z) for z in z_2_digits]

    for z, m in zip(decimal_zs, mags):
        i_mag, z_mag = m  # split the magnitudes into i and z
        # The correction to the redshift is done when the RSModel
        # object is initialized
        this_model = RSModel(z, i_mag, z_mag)
        rs_models[this_model.z] = this_model

    return rs_models


def _make_model():
    """
    Make an EzGal object.

    This also takes care of making sure that the model exists and is
    pre-evolved in the way we like.

    :return: EzGal object.
    """

    evolved_model_name = "bc03_exp_0.1_z_0.02_chab_evolved_zf_3.0.model"
    default_model_name = "bc03_exp_0.1_z_0.02_chab.model"

    try:  # to open the evolved model
        model = ezgal.ezgal(evolved_model_name)
    except ValueError:  # the default model wasn't found
        try:  # to open the default model
            model = ezgal.model(default_model_name)
            _evolve_model(model, savename=evolved_model_name)

        except ValueError:  # the default model doesn't exist
            raise ValueError("Please download the default model, which is "
                             "bc03_exp_0.1_z_0.02_chab.model. This can be "
                             "found on the EzGal website at "
                             "http://www.baryons.org/ezgal/download.php")
    return model


def _evolve_model(model, savename):
    """
    Will do model evolution for a model, and save the resulting model.

    :param model: model to be evolved.
    :return:None, but model will be saved.
    """
    # find the place ezgal stores all the filters
    filters_dir = model.data_dir + "filters" + os.sep

    # make a list of all filters in that directory to use in the model.
    all_filters = [f for f in os.listdir(filters_dir) if f != "README"]
    for f in all_filters:
        model.add_filter(f)

    # add a bunch of formation redshifts, too
    model.set_zfs(np.arange(1.0, 4.5, 0.5))
    # This will calculate the evolution

    # Save the model.
    location = model.data_dir + "models" + os.sep + savename
    model.save_model(location)
