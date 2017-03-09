import os
import decimal

import numpy as np
import ezgal

import config


class _Slope(object):
    """
    Callable class used to calculate slopes of the red sequence.

    There is probably an easier way to make all this work, but I couldn't
    figure out the best way. Suggestions are definitely welcome.
    """

    def __init__(self, cfg):
        """
        Find the slope of the red sequence as a function of redshfit.

        Uses the slopes of the red sequence as published in Eisenhardt 2007.

        :return: a function that returns the slope at a given redshift
        """

        # Then do some interpolating here to find the slope at any redshift

        def slope(z):
            # fit is a polynomial (lowest power first) that describes the
            # slope of the red sequence as a function of redshift, and is
            # different for each band combo.
            fit = cfg["slope_fit"]

            # turn to decimal, so it can play nice with redshifts, which are
            # of type decimal
            fit = [decimal.Decimal(i) for i in fit]
            # plug the redshift into the polynomial. This looks ugly, but
            # that's all it's doing.
            return float(sum([coeff * (z ** i) for i, coeff in
                              enumerate(fit)]))
        # we set the object's slope_function object to be the function we just
        # defined, so we can call it later.
        self.slope_function = slope

    # I'll use a call method so I can call the object like a function to get
    #  the slope of the red sequence at the desired redshift.
    def __call__(self, redshift):
        return self.slope_function(redshift)


class RSModel(object):
    """
    Class storing data from the EzGal models.
    """

    @staticmethod
    def correction(redshift, cfg):
        """
        Does the redshift correction.

        After comparing clusters with known
        redshifts to those produced by this program, this correction was
        devised to make this program return results more closely matching
        the actual redshift. It's just an empirical correction to make the
        redshifts better.

        :param redshift: uncorrected redshift
        :param cfg: Configuration dictionary.
        :return: corrected redshift
        """
        # This is the fit, which comes from numpy.polynomial.polynomial.polyfit
        # It is the coefficients of the polynomial, starting with the
        # lowest power of z.
        fit = cfg["correction"]
        # turn to decimal , so it can play nice with redshifts, which are
        # of type decimal
        fit = [decimal.Decimal(i) for i in fit]
        # Applying the correction is just plugging in the redshift to the
        # correction polynomial. This is a complicated way to do that.
        corrected = sum([coeff * (redshift**i) for i, coeff in
                         enumerate(fit)])
        # round and turn into type decimal, which is what redshifts need to be.
        return decimal.Decimal(str(round(corrected, 3)))

    def __init__(self, redshift, blue_mag, red_mag, cfg):
        """
        Initializes the useful parameters that will be used to calculate
        where the red sequence is.

        These parameters should all come straight out of ezgal's
        get_apparent_mags() function.

        :param redshift: redshift at which these points are calculated
        :param blue_mag: blue_mag at that redshift
        :param red_mag: red magnitude at that redshift
        :return: Model object.
        """

        self.mag_point = red_mag
        self.color_point = blue_mag - red_mag
        self.z = RSModel.correction(redshift, cfg)  # use corrected redshift

        # we also put a slope object in the object, so we can easily access
        # the slope of the red sequence in these bands.
        self._slope_obj = _Slope(cfg)
        self._slope = self._slope_obj(redshift)

    def rs_color(self, red_mag):
        """
        Calculate the color of the red sequence at the given red_mag
        magnitude.

        :param red_mag: red magnitude at which we want the color of the
                        red sequence
        :return: float value with the color of the red sequence

        Algorithm: Start with point slope form of a line
        y - y1 = m(x - x1)
        y = y1 + m(x - x1)
        Where x is red magnitude, y is color, m is the slope of the
        red sequence, and x1 and y1 are the zero point that was returned by
        EzGal.
        """

        return self.color_point + self._slope * (red_mag - self.mag_point)


def model_dict(spacing):
    """
    Create a dictionary of model objects, that represent the red sequence at
    different redshifts.

    :param spacing: spacing of the redshifts at which the RS will be
                    modeled.
    :return: dictionary that holds dictionaries. The first set of keys are the
             different color combinations used (ex: ch1-ch2). The first values
             are dictionaries themselves, which have keys of redshift, and
             values of RSModel objects that represent the RS at that redshift
             in that color combo.
    """

    # We need to figure out which filters we need info for. We'll use a set
    # to avoid duplicates
    filters = set()
    for color in config.cfg_matches:
        band1, band2 = color.split("-")
        filters.add(band1)
        filters.add(band2)
    filters = list(filters)  # we go back to a list so EzGal will like it.

    # get the EzGal model object
    model = _make_model(filters)

    # decide the formation redshift and observed redshift
    zf = 3.0
    zs = np.arange(0.1, 2.500001, spacing)

    # normalize to Coma
    model.set_normalization(filter='ks', mag=10.9, apparent=True, vega=True,
                            z=0.023)

    # then get AB mags in those filters
    mags = model.get_apparent_mags(zf, filters=filters, zs=zs, ab=True)
    # this ^ is a 2d array. The first dimension is redshift, second is filters.

    # initialize an empty dictionary that will be filled with RSModel objects.
    # This will be a nested dictionary. The first set of keys will be the
    # different colors, the second set will be redshifts, and the values will
    # be the model objects at that color and redshift
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

    # we can then put things into the dictionary
    for color in config.cfg_matches:
        rs_models[color] = dict()
        band_1, band_2 = color.split("-")
        band_1_idx = filters.index(band_1)
        band_2_idx = filters.index(band_2)

        this_config = config.cfg_matches[color]

        for z, m in zip(decimal_zs, mags):
            # only do things if the redshift is in the right range
            if this_config["z_min"] <= z <= this_config["z_max"]:
                mag_1 = m[band_1_idx]  # split the magnitudes into bands
                mag_2 = m[band_2_idx]
                # The correction to the redshift is done when the RSModel
                # object is initialized
                this_model = RSModel(z, mag_1, mag_2, this_config)
                rs_models[color][this_model.z] = this_model

    return rs_models


def _make_model(filters):
    """
    Make an EzGal object.

    This also takes care of making sure that the model exists and is
    pre-evolved in the way we like.
    :param filters: List of filters that the model needs to contain
                    calculations for.
    :return: EzGal object.
    """
    # we have an evolved model in the local directory pre-evolved that we
    # want to grab if we can.
    code_dir = os.path.dirname(os.path.realpath(__file__))
    evolved_model = code_dir + os.sep +"bc03_exp_0.1_z_0.02_chab_evolved.model"
    # if it's not there, then get the un-evolved model.
    default_model = "bc03_exp_0.1_z_0.02_chab.model"

    try:  # to open the evolved model
        model = ezgal.ezgal(evolved_model)
        # see if all the filters we want are already in the model.
        if not all([filt in model.filters for filt in filters]):
            # if not, add them
            _evolve_model(model, filters, evolved_model)
    except ValueError:  # the default model wasn't found
        try:  # to open the default model
            model = ezgal.model(default_model)
            # add the filters we need
            _evolve_model(model, filters, save_name=evolved_model)

        except ValueError:  # the default model doesn't exist
            raise ValueError("Please download the default model, which is "
                             "bc03_exp_0.1_z_0.02_chab.model. This can be "
                             "found on the EzGal website at "
                             "http://www.baryons.org/ezgal/download.php")
    return model


def _evolve_model(model, filters, save_name):
    """
    Will do model evolution for a model, and save the resulting model.

    :param model: model to be evolved.
    :param filters: List of filters that need to be added to the model, and
                    have evolution calculated for.
    :param save_name: location where the evolved model will be saved.
    :return: None, but model will be saved.
    """
    # tell the user what they need to know
    print("Calculating model evolution, will take a while...")
    print("Only needs to be done once, unless you add more filters later.")
    print("Please ignore the warnings that will follow. The divide by zero\n" \
          "errors are from EzGal's internals and are fine.\n" \
          "The deprecation warnings are from EzGal using \n" \
          "Pyfits, which is old.\n")

    # add a bunch of formation redshifts, so we shouldn't have to do that later
    zfs = np.arange(1.0, 6.001, 0.05)
    model.set_zfs(zfs)

    # normalize to Coma
    model.set_normalization(filter='ks', mag=10.9, apparent=True, vega=True,
                            z=0.023)
    # then add all the filters to our model before calculating its evolution.
    try:
        for filt in filters:
            if filt not in model.filters:
                model.add_filter(filt, grid=True)
    except ValueError:
        raise IOError("The filter file {} was not found.\n"
                      "\tMake sure the name of the filter you gave in the\n"
                      "\tconfig file is in the `ezgal/data/filters` \n"
                      "\tdirectory. If it is a custom filter file, please \n"
                      "\tadd it there.".format(filt))

    # calculate mags to ensure things get calcuated. We don't do anything with
    # this information, we just need to make sure the model is fully evolved
    # and ready for the user to use without any further calculations.
    for zf in zfs:
        zs = np.arange(0.1, zf, 0.01)
        model.get_apparent_mags(zf, filters=filters, zs=zs, ab=True)

    # Save the model.
    model.save_model(save_name)
