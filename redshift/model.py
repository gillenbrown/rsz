import ezgal
import os
import numpy as np

class _Slope(object):
    """
    Callable class used to calculate slopes of the red sequence.
    """

    def __init__(self):
        """
        Find the slope of the red sequence as a function of redshfit.

        Uses the slopes of the red sequence as published in Eisenhardt 2007.

        :return: a function that returns the slope at a given redshift
        """

        # first we need to see at what redshift does ch1-ch2 see the various
        # Eisenhardt colors.
        # I did this in an iPython notebook that is in the
        # calculating_slopes folder, if you want to see the work. I'll just
        # manually add the results here, to avoid unecessary computation
        redshifts = [1, 2, 3]
        slopes = [-0.1, -0.04, -0.08]
        #TODO: These ^ are  made up! Find the actual ones.

        # TODO: do some interpolating here to find the slope at any redshift
        def f(z):
            return -0.05

        self.slope_function = f

    # I'll use a call method so I can call the object like a function to get
    #  the slope of the red sequence at the desired redshift.
    def __call__(self, redshift):
        return self.slope_function(redshift)


class RSModel(object):
    """
    Class storing data from the EzGal models.
    """


    # get slopes for all redshifts
    slope = _Slope()

    @staticmethod
    def correction(redshift):
        """
        Does the redshift correction.

        After comparing clusters with known
        redshifts to those produced by this program, this correction was
        devised to make this program return results more closely matching
        the actual redshift. It's just an impirical correction to make the
        redshifts better.

        :param z: uncorrected redshift
        :return: corrected redshift
        """
        # TODO: actually make this correction.
        # TODO: make an iPython notebook documenting the way I found this
        # correction, once I actually do that.
        return redshift


    def __init__(self, redshift, ch1_mag, ch2_mag):
        """
        Initializes the useful parameters that will be used to calculate
        where the red sequence is.

        :param redshift: redshift at which these points are calculated
        :param ch1_mag: ch1_mag at that redshift
        :param ch2_mag: ch2_mag at that redshift
        :return: Model object.
        """

        self.mag_point = ch2_mag
        self.color_point = ch1_mag - ch2_mag
        self.z = RSModel.correction(redshift)
        self._slope = RSModel.slope(redshift)

    def rs_color(self, ch2_mag):
        """
        Calculate the color (ch1-ch2) of the red sequence at the given ch2
        magnitude.

        :param ch1_mag: ch1 magnitude at which we want the color of the
                        red sequence
        :return: float value with the color of the red sequence

        Algorithm: Start with point slope form of a line
        y - y1 = m(x - x1)
        y = y1 + m(x - x1)
        Where x is ch2 magnitude, y is ch1-ch2 color, m is the slope of the
        red sequence, and x1 and y1 are the zeropoint that was returned by
        EzGal.
        """

        return self.color_point + self._slope * (ch2_mag - self.mag_point)


def model_dict(spacing):
    """
    Create a dictionary of model objects, that represent the red sequence at
    different redshifts.

    :param spacing: spacing of the redshifts at which the RS will be
    calculated.
    :return:
    """

    # get the model
    model = _make_model()

    #set the formation redshift and observed redshift
    zf = 3.0
    zs = np.arange(0.7, 1.700001, spacing)

    # normalize to Coma
    model.set_normalization(filter='ks', mag=10.9, apparent=True, vega=True,
                            z=0.023)

    # Calculate the obserables we want in AB mags
    mags = model.get_apparent_mags(zf, filters=["ch1", "ch2"], zs=zs,
                                   ab=True)

    # initialize an empty dictionary that will be filled with RSModel objects
    rs_models = dict()

    # turn redshifts into strings, to avoid floating point errors
    zs = [str(round(z, 5)) for z in zs]

    for z, m in zip(zs, mags):
        ch1, ch2 = m  # split the magnitudes into ch1 and ch2
        this_model = RSModel(z, ch1, ch2)
        rs_models[this_model.z] = this_model

    return rs_models







def _make_model():
    """
    Make an EzGal object.

    This also takes care of making sure that the model exists and is
    pre-evolved in the way we like.

    :return: EzGal object.
    """

    evolved_model = "bc03_exp_0.1_z_0.02_chab_evolved_zf_3.0.model"
    default_model = "bc03_exp_0.1_z_0.02_chab.model"

    try:  # to open the evolved model
        model = ezgal.ezgal(evolved_model)
    except ValueError:  # the default model wasn't found
        try:  # to open the default model
            model = ezgal.model(default_model)
            _evolve_model(model,
                          "bc03_exp_0.1_z_0.02_chab_evolved_zf_3.0.model" )

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