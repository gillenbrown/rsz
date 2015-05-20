import os
import math

import model
import source
import data
import plotting

class Cluster(object):
    # Keeping the predictions for the red sequence with the cluster object
    # made things a lot easier. And since the predictions are the same for
    # all cluster, this is a class variable rather than an instance variable.
    predictions_dict = model.model_dict(0.01)

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
            plotting.add_models(fig, ax)
            figures.append(fig)

        # Do a quick and dirty initial redshift fitting, to get a starting
        # point.
        z = self._initial_z()
        print z

        # TODO: Next: plot this intial redshift, and highlight red sequence
        # members






        # now that we are all done, save the figures.
        Cluster.save_as_one_pdf(figures, params["plot_directory"] +
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
                if (mag_point - 1.0 < source.ch2 < mag_point + 1.0) and \
                   (color_point - 0.2 < source.ch1_m_ch2 < color_point + 0.2):
                    # if it did pass, note that it is nearby
                    nearby += 1

            # replace the best values if this is the best
            if nearby > max_nearby:
                max_nearby = nearby
                best_z = z

        return best_z



    @staticmethod
    def save_as_one_pdf(figs, filename):
        """
        Save the figures into one long PDF file

        :param figs: list of figures to be saved as PDFs
        :param filename: place where the PDFs will be saved
        :return: none
        """
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt

        # Save the pdfs as one file
        pp = PdfPages(filename)
        for fig in figs:
            pp.savefig(fig)
        pp.close()
        for fig in figs:
            plt.close(fig)