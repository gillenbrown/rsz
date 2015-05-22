import os
import math
import numpy as np

import model
import source
import data
import plotting

class Cluster(object):
    # Keeping the predictions for the red sequence with the cluster object
    # made things a lot easier. And since the predictions are the same for
    # all cluster, this is a class variable rather than an instance variable.
    models = model.model_dict(0.01)

    def __init__(self, filepath):
        self.name = Cluster._name(filepath)

        # I'll initialize emtpy source list, that will be filled as we go
        self.sources_list = []
        self.z = None

        # then read in the objects in the catalog
        with open(filepath) as cat:
            for line in cat:
                if not line.startswith("#"):
                    # Each line has the data of ra, dec, ch1 mag, ch1 error,
                    #  ch2 mag, ch2 error.
                    entries = line.split()
                    ra = float(entries[0])
                    dec = float(entries[1])
                    ch1 = data.Data(float(entries[2]), float(entries[3]))
                    ch2 = data.Data(float(entries[4]), float(entries[5]))

                    this_source = source.Source(ra, dec, ch1, ch2)
                    self.sources_list.append(this_source)

    @staticmethod
    def _name(filepath):
        filename = os.path.split(filepath)[-1]
        #TODO: parse the filename based on the format. I don't know what
        # that will be at the moment.
        return filename

    def __repr__(self):
        return self.name


    def fit_z(self, params):
        """
        Find the redshift of the cluster by matching its red sequence the
        red sequence models produced by EzGal. This is the main
        functionality of the entire code.

        :param params: Dictionary full of the user's parameters for the
                       program from the config file.
        :return: redshift of the cluster. Will be an AsymmetricData
        instance, to account for the possibility of different upper and
        lower limits on the redshift.
        """

        # initialize a list of figures, that will be filled as we go along.
        figures = []

        # do a location cut, to only focus on galaxies near the center of
        # the image.
        self._location_cut(1.5)

        # If the user wants, plot the initial CMD with predictions
        if params["CMD"] == "1":
            fig, ax = plotting.cmd(self)
            plotting.add_all_models(fig, ax)
            figures.append(fig)

        # Do a quick and dirty initial redshift fitting, to get a starting
        # point.
        self.z = self._initial_z()

        # If the user wants to see this initial, fit, plot it.
        if params["fitting_procedure"] == "1":
            fig, ax = plotting.cmd(self)
            plotting.add_one_model(ax, self.models[self.z.value], "k")
            plotting.add_redshift(ax, self.z.value)
            figures.append(fig)

        # set cuts that will be used in the successive iterations to refine
        # the fit of the red sequence
        bluer_color_cut = [0.25, 0.225, 0.20]
        redder_color_cut = [0.4, 0.3, 0.2]
        # the bluer color cut is smaller than the red color cut to try to
        # throw away foregrounds. If too many foregrounds are included,
        # they drag the redshift down, since those foregrounds typically
        # have lower error, which biases the fit. That isn't a problem with
        # objects redder than the RS.
        brighter_mag_cut = 1.4
        dimmer_mag_cut = 0.6

        # do three iterations of fitting, each with a progressively smaller
        # color cut, which is designed to hone in on the red sequence.
        for i in range(0, 3):
            # set red sequence members based on the cuts
            self._set_RS_membership(self.z.value, bluer_color_cut[i],
                                    redder_color_cut[i], brighter_mag_cut,
                                    dimmer_mag_cut)

            self.z = self._chi_square_w_error()

            # if the user wants, plot the procedure
            if params["fitting_procedure"] == "1":
                fig, ax = plotting.cmd(self)
                plotting.add_one_model(ax, self.models[self.z.value], "k")
                plotting.add_redshift(ax, self.z.value)
                figures.append(fig)

        # we now have a final answer for the redshift of the cluster.
        # if the user wants, plot it up
        if params["final_CMD"] == "1":
            fig, ax = plotting.cmd(self)
            plotting.add_one_model(ax, self.models[self.z.value], "k")
            # I want to plot both the low and high models, so get those zs
            high_z = self.z.value + self.z.upper_error
            low_z = self.z.value - self.z.lower_error

            plotting.add_one_model(ax, self.models[high_z], "0.6")
            plotting.add_one_model(ax, self.models[low_z], "0.6")
            plotting.add_redshift(ax, self.z)

            figures.append(fig)






        # now that we are all done, save the figures.
        save_as_one_pdf(figures, params["plot_directory"] +
                                self.name + ".pdf")




    def _location_cut(self, radius):
        """
        Does a location cut on the galaxies in the image.

        Only those within a certain radius are marked as being near_center,
        which is what we will use for the rest of the fitting. This throws
        away foregrounds, leaving the most dense portion of the cluster. A
        good radius to pick isn't necessarily the radius of the cluster.
        Picking a small radius can increase the signal to noise by only
        including the densest part of the cluster, and throwing out most
        foregrounds.

        The radius cut is done from the center of the image, which isn't
        necessarily where the galaxies are most dense.

        :param radius: radius of the cut, in arcminutes.
        :return: None, but galaxies near the center are marked as having
        near_center = True.
        """

        ras = [source.ra for source in self.sources_list]
        decs = [source.dec for source in self.sources_list]

        middle_ra = (max(ras) + min(ras)) / 2.0
        middle_dec = (max(decs) + min(decs)) / 2.0

        for source in self.sources_list:
            # Use pythagorean theorem to find distance in degrees
            dist = math.sqrt((source.ra - middle_ra)**2 +
                             (source.dec - middle_dec)**2)

            if dist < radius/60.0:  # divide by 60 since radius is in arcmin
                source.near_center = True
            else:
                source.near_center = False

    def _initial_z(self):
        """
        Find a decent initial redshift estimate, based on the number of
        galaxies near each model.

        It iterates through all redshifts, and counts the number of galaxies
        that are within XXXXXXX color of the model at that redshift. Then the
        redshift that has the highest number of galaxies is the
        redshift selected. It's a quick a dirty estimate.

        :return: an intial redshift estimate
        """

        # Get models with a large spacing, since we don't need a lot of
        # accuracy here.
        models = model.model_dict(0.05)

        # set placeholder values that will be replaced as we go
        max_nearby = -999
        best_z = -999

        # Iterate through the redshifts
        for z in models:
            nearby = 0  # reset the number of nearby galaxies

            this_model = models[z]
            mag_point = this_model.mag_point
            for source in self.sources_list:
                # get the expected RS color for a galaxy of this magnitude
                color_point = this_model.rs_color(source.ch2.value)
                # see if it passes a color and magnitude cut
                if (mag_point - 1.5 < source.ch2 < mag_point + 0.7) and \
                   (color_point - 0.3 < source.ch1_m_ch2 < color_point + 0.3):
                    # if it did pass, note that it is nearby
                    nearby += 1

            # replace the best values if this is the best
            if nearby > max_nearby:
                max_nearby = nearby
                best_z = z

        # Turn this into a data object, with no errors
        z = data.AsymmetricData(best_z, np.NaN, np.NaN)
        return z

    def _set_RS_membership(self, redshift, bluer, redder, brighter, dimmer):
        """Set some sources to be RS members, based on color and mag cuts.

        :param redshift: redshift of the red sequence that these sources
                         will be identified as belonging to.
        :param bluer: maximum color difference on the blue side from the
                      characteristic color of the red sequence at the
                      magnitude of the galaxy.
        :param redder: maximum color differenece on the red side from the
                       characteristic color of the RS
        :param brighter: maximum magnitude difference (on the bright end)
                         from the characteristic magnitude of the red
                         sequence.
        :param dimmer: maximum magnitude difference (on the faint end)
                         from the characteristic magnitude of the red
                         sequence.

        As an example: At some redshift, the red sequence has a characteristic
        magnitude of 20, and a characteristic color of zero. We pass in
        bluer=0.1, redder=0.2, brighter=2.0, dimmer=1.0. Any galaxies that
        have magnitude 18 to 21 pass the magnitude cut. The color cut is
        trickier, since the RS has some slope. We find the characteristic
        color of the red sequence at the ch2 magnitude of the galaxy (say
        0.05 for a particular galaxy). Then if the color is within the
        bounds set by bluer and redder (in this case -0.05 to 0.25), then it
        passes the color cut. If a galaxy passes both the magnitude and
        color cuts, it is marked as a red sequence galaxy.

        Note: This function marks galaxies both inside and outside the
        location cut as RS members. This is because the loctaion cut may be
        smaller than the cluster (in fact this is often good to do), and we
        don't want to mark galaxies as not RS members just because they are
        outside the (possibly) small location cut.

        :return: None, but galaxies that are RS members are marked as such.
        """

        # get the model, it's characteristic magnitude, and then turn it
        # into magnitude limits bsed on the parameters passed in
        RS_model = self.models[redshift]
        char_mag = RS_model.mag_point
        dim_mag = char_mag + dimmer
        bright_mag = char_mag - brighter

        for source in self.sources_list:
            #get the color correspoinding to the red sequence at the ch2
            # magnitude of this particular source
            char_color = RS_model.rs_color(source.ch2.value)
            # turn it into color limits based on parameters passed in
            red_color = char_color + redder
            blue_color = char_color - bluer
            # a function in the source class does the actual marking
            source.RS_membership(blue_color, red_color, bright_mag, dim_mag)

    def _chi_square_w_error(self):
        """Does chi-squared fitting, and returns the best fit value and the
        1 sigma error.

        The best fit value is simply the value that minimizes the reduced
        chi squared value. The upper and lower errors on this are the value
        at which the chi squared value is 1 greater than the best fit.

        :returns: An AsymmetricData object, where value is the best fit
        values, and upper and lower values are indicated.
        """

        # we only want to do the fitting on those near the center, and those
        #  that are in our tentative RS.
        to_fit = [source for source in self.sources_list if source.RS_member
                  and source.near_center]

        # initialize lists for the chi squared distribution and the redshifts
        chi_sq_values = []
        redshifts = []

        # test each model
        for z in sorted(self.models):
            chi_sq = 0  # reset the chi squared value for this model
            this_model = self.models[z]
            for source in to_fit:
                model_color = this_model.rs_color(source.ch2.value)
                color = source.ch1_m_ch2.value
                error = source.ch1_m_ch2.error
                chi_sq += ((model_color - color)/error)**2

            #reduce the chi square values
            chi_sq /= len(to_fit)
            # put these values into the lists, so we can keep track of them
            chi_sq_values.append(chi_sq)
            redshifts.append(z)

        best_chi = min(chi_sq_values)
        best_chi_index = chi_sq_values.index(best_chi)

        high_idx = best_chi_index
        # find the place where the chi value is one greater than the best
        while high_idx < len(chi_sq_values) - 1 and \
              chi_sq_values[high_idx] - best_chi < 1.0:
            # Note: len(chi_sq_values) - 1 needs that -1 so that we don't
            # end up with an high_idx that is past the end of chi_sq_values
            high_idx += 1
        # high_idx will now point to the first value where the chi squared
        # value is one greater than the best fit value.

        # TODO: should I do some interpolation here to find the place
        # where it is exactly one, or just go with this? I'll go with it
        #  for now, but I need to decide.

        # do the same thing for the low error
        low_idx = best_chi_index
        while low_idx > 0 and chi_sq_values[low_idx] - best_chi < 1.0:
            low_idx -= 1

        # now get the redshifts corresponding to the best fit, the low
        # error, and the high error
        # This works because the indexing will be the same for redshifts as
        # it was for chi_sq_values
        best_z = redshifts[best_chi_index]
        high_z = redshifts[high_idx]
        low_z = redshifts[low_idx]

        # turn these limits into errors, then into an AsymmetricData object
        high_error = high_z - best_z
        low_error = best_z - low_z
        return data.AsymmetricData(best_z, high_error, low_error)
        # TODO: that was a long function. Break it up somehow?










def save_as_one_pdf(figs, filename):
    """
    Save the figures into one long PDF file

    :param figs: list of figures to be saved as PDFs
    :param filename: place where the PDFs will be saved
    :return: none
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    # TODO: move these imports to the top of whatever file this function
    # ends up in

    # Save the pdfs as one file
    pp = PdfPages(filename)
    for fig in figs:
        pp.savefig(fig)
    pp.close()
    for fig in figs:
        plt.close(fig)