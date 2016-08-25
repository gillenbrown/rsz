import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from astropy import convolution

import model
from source import Source
import data
import plotting
import config


class Cluster(object):
    """
    Class that holds the data about a cluster, and that does the redshift
    fitting.
    """

    # Keeping the predictions for the red sequence with the cluster object
    # made things a lot easier. And since the predictions are the same for
    # all cluster, this is a class variable rather than an instance variable.
    # The regular one is for fitting, the low res one is for plotting all the
    # models at once
    models = model.model_dict(0.01)
    low_res_models = model.model_dict(0.05)

    def __init__(self, file_path, params):
        """
        Initialize the cluster object. Turns the file path into a nice name,
        and sets other things to emtpy lists or zeros.

        :param file_path: Path to the catalog containing the cluster. The
                         extension will be removed to make the name.
        :param params: Parameters dictionary from the input file.
        """
        self.name = self._name(file_path, params["extension"])

        # I'll initialize empty source list, that will be filled as we go
        self.sources_list = []

        # redshift and flags are dictionaries, since they will have different
        # values for different bancs.
        self.z = dict()
        self.flags = dict()

        self.figures = []

        self.interesting = 0

        self.center_ra = None
        self.center_dec = None

        # then read in the objects in the catalog
        self.read_catalog(file_path, params)

    @staticmethod
    def _name(file_path, extension):
        """
        Turns a file path into a name for the cluster. All this does is remove
        the extension from the filename.

        :param file_path: Path to the catalog containing the data for this
                         cluster.
        :param extension: File extension on the catalog.
        :return: Name for this cluster.
        """
        # get just the filename, ignore the rest of the path.
        filename = os.path.split(file_path)[-1]
        # just remove the extension from the filename
        return filename.rstrip(extension)

    def read_catalog(self, file_path, params):
        """ Read the catalog, parsing things into sources.

        This function calls other functions to do the dirty work.

        :param file_path: path of the catalog to be parsed.
        :param params: parameter dictionary
        :return: none, but the cluster's source list variable is populated.
        """

        with open(file_path) as cat:
            for line_num, line in enumerate(cat, start=1):
                # some catalog formats aren't formatted the way I'd like (the
                # column headers aren't commented out), so ignore that line.
                if not line.startswith("#"):
                    # split the line, to make for easier parsing.
                    split_line = [self.to_float(i) for i in line.split()]

                    # get the various things from the parser functions.
                    ra, dec = self.get_ra_dec(params, split_line, line_num)
                    mags = self.get_mags(params, split_line, line_num)
                    dist = self.get_dist(params, split_line, line_num)

                    # check that the user has specified the appropriate things.
                    if ra is None and dec is None and dist is None:
                        raise TypeError("Specify one of either ra/dec or dist")

                    # turn this info into a source object, then add it.
                    this_source = Source(ra, dec, mags, dist)
                    self.sources_list.append(this_source)

    @staticmethod
    def to_float(item):
        """If an item can be converted to a float, it will be. If not, the item
        will be returned unchanged.

        :param item: Item to potentially convert to a float.
        :return: the item, potentially converted to a float.
        """
        try:
            return float(item)
        except ValueError:  # That's the error from a failed float conversion.
            return item

    @staticmethod
    def _check_valid_int(params, key):
        """
        Checks whether the value of they key in the parameters dictionary is
        an integer. Throwns an error if not.

        :param params: Dictionary with parameters specified by the user.
        :param key: Name of the parameter that we are checking.
        :type key: str
        :return: If the value can be turned into an int, we return the integer
                 value of it. If not, we raise a ValueError.
        """
        try:
            return int(params[key])
        except ValueError:
            raise ValueError("The parameter {} is formatted incorrectly.\n"
                             "\tIf it is supposed to be the index of a\n"
                             "\tband, then it should just be an integer.\n"
                             "\tIf not, then then this parameter is not\n"
                             "\tused by the code. Please remove it."
                             "".format(key))

    def _check_valid_idx(self, split_line, idx, line_number):
        """
        See whether a given index is valid, given the line that is supposed
        to index. If the index is invalid, we raise an error and explain what
        went wront to the user.

        :param split_line: Line that has been split into a list of items.
        :type split_line: list
        :param idx: Index that is supposed to be used on split_line.
        :type idx: int
        :param line_number: Line number in the catalog that split_line
               holds the data for. Only used in error messages if something
               goes wrong.
        :type line_number: int
        :return: The value of split_line[idx]. If that doesn't exist, we raise
                 an error.
        """
        try:
            return split_line[idx]
        except IndexError:
            raise ValueError("The indexing appears to be broken on line {}\n"
                             "\tof the {} catalog.\n "
                             "\tCheck the indexes in the param file, as well\n"
                             "\tas the line in that catalog itself."
                             "".format(line_number, self.name))

    def get_ra_dec(self, params, split_line, line_number):
        """Parses the line to get the ra and dec.

        If the user doesn't specify ra and dec, then return None.

        :param params: Parameter dictionary that is passed all around the code.
        :param split_line: Line of the catalog split into the components.
        :param line_number: Line number in the catalog. Only used for error
                            messages if something doesn't work.
        :return: ra and dec as a tuple
        """

        ra_idx = self._check_valid_int(params, "ra")
        dec_idx = self._check_valid_int(params, "dec")
        ra = self._check_valid_idx(split_line, ra_idx, line_number)
        dec = self._check_valid_idx(split_line, dec_idx, line_number)

        return ra, dec

    def get_mags(self, params, split_line, line_number):
        """ Parses the config file to get the magnitudes.

        :param params: Parameter dictionary that is passed all around the code.
        :param split_line: Line of the catalog split into the components.
        :param line_number: Line number in the catalog. Only used for error
                            messages if something doesn't work.
        :return: dictionary with keys of band names (eg: ch1, r, z), and values
                 of data objects with the mag and mag error in that band.
        """
        # see which bands were specified in the param file
        non_band_params = ["catalog_directory", "extension",
                           "plot_directory", "results_file", "rs_catalog_dir",
                           "type", "mag_zeropoint", "mag_system", "ra", "dec",
                           "dist", "CMD", "fitting_procedure", "final_CMD",
                           "location", "interactive"]

        # we'll iterate through all keys in the user's paramter dictionary,
        # then take the ones that aren't the known non-band keys above.
        bands = []
        for key in params:
            if key in non_band_params:
                continue
            # the key we have is either a flux/mag or an error
            if not key.endswith("_err"):
                bands.append(key)
                # if we have a flux, we need to have an error, too.
                if key + "_err" not in params:
                    raise ValueError("The band {} does not have an error\n"
                                     "\tassociated with it. Please include\n"
                                     "\tthat in the param file. If that \n"
                                     "\tparameter is something other than\n"
                                     "\ta band, please remove it. It is not\n"
                                     "\tneeded by the code.".format(key))

        # then we are ready to get the info for the bands we identified
        mags = dict()
        for band in bands:
            # verify the indicies the user specified, and then get the data.
            band_idx = self._check_valid_int(params, band)
            band_err_idx = self._check_valid_int(params, band + "_err")

            band_data = self._check_valid_idx(split_line, band_idx,
                                              line_number)
            band_data_err = self._check_valid_idx(split_line, band_err_idx,
                                                  line_number)

            # convert fluxes to magnitudes if need be
            if params["type"] == "flux":
                if band_data == 0:
                    # give it a slightly negative flux, which will be
                    # interpreted as a bad mag
                    band_data = -0.1
                # calculate magnitude err first, since we need to use the flux
                # to calculate the magnitude.
                band_data_err = Cluster.percent_flux_errors_to_mag_errors(
                    band_data_err / band_data)
                band_data = Cluster.flux_to_mag(band_data,
                                                params["mag_zeropoint"])

            # Convert to AB mags if needed. If we had flux, they are AB already
            if params["type"] == "mag" and params["mag_system"] == "vega":
                try:
                    band_data -= config.ab_to_vega(band)
                except KeyError:
                    raise KeyError("Please specify the AB/Vega conversion "
                                   "\tfor {} in config.py.".format(band))

            # convert to Data type, and add to dictionary
            mags[band] = data.Data(band_data, band_data_err)

        return mags

    def get_dist(self, params, split_line, line_number):
        """Parse the line to get the distance from center.

        :param params: Parameter dictionary that is passed all around the code.
        :param split_line: Line of the catalog split into the components.
        :param line_number: Line number of this line in the catalog. Only used
                            for error checking.
        :return: The distance of the object from the cluster center.
        """
        if params["dist"] != "-99":  # -99 means its not in the catalog.
            dist_idx = self._check_valid_int(params, "dist")
            return self._check_valid_idx(split_line, dist_idx, line_number)
        else:  # if the user didn't specify
            return None

    @staticmethod
    def flux_to_mag(flux, zeropoint):
        """Convert flux to magnitude with the given zeropoint.

        :param flux: flux in whatever units. Choose your zeropoint correctly
                     to make this work with the units flux is in.
        :param zeropoint: zeropoint of the system, such that
                          m = -2.5 log(flux) + zeropoint
        :return: magnitude that corresponds to the given flux. If flux is
                 negative, -99 is returned.
        """
        if flux <= 0:
            return -99
        return -2.5 * np.log10(flux) + zeropoint

    @staticmethod
    def mag_to_flux(mag, zeropoint):
        """
        Converts a magnitude into a flux, given the zeropoint of the mag system

        :param mag: Magnitude
        :param zeropoint: Zeropoint of the magnitude system, such that
                          m = -2.5 log(flux) + zeropoint
        :return: flux corresponding to the given magnitude. If the mag is less
                 than zero (which while physical, will only happen in this
                 code if there is an error), -99 will be returned.
        """
        if mag < 0:
            return -99
        return 10 ** ((zeropoint - mag) / 2.5)

    @staticmethod
    def mag_errors_to_flux_errors(mag_error, flux):
        """
        Converts magnitude errors into flux errors.

        :param mag_error: Magnitude error
        :param flux: Flux of the object in question. This is needed since
                     magnitude errors correspond to percentage flux errors.
        :return:
        """
        if flux < 0:
            return -99
        return (flux * np.log(10) * mag_error) / 2.5

    @staticmethod
    def percent_flux_errors_to_mag_errors(percent_flux_error):
        """Converts a percentage flux error into a magnitude error.

        m = -2.5 log10(F) + C
        dm = -2.5/(ln(10)) dF/F

        :param percent_flux_error: percentage flux error
        :return: magnitude error corresponding to the percentage flux error.
        """
        if percent_flux_error < 0:
            return 99
        return (2.5 / np.log(10)) * percent_flux_error

    def __repr__(self):
        return self.name

    def fit_z(self, params, cfg):
        """
        Find the redshift of the cluster by matching its red sequence the
        red sequence models produced by EzGal. This is the main
        functionality of the entire code.

        :param params: Dictionary full of the user's parameters for the
                       program from the config file.
        :param cfg: Configuration dictionary for the color combo in question.
        :return: None, but the redshift of the cluster is set.
        """
        # do a location cut, to only focus on galaxies near the center of
        # the image.
        self._location_cut(1.0, params)  # 1.0 is in arcminutes

        # store the color we are working with, to make things cleaner
        this_color = cfg["color"]

        # we need to initialize the flags for the cluster
        self.flags[this_color] = 0

        # If the user wants, plot the initial CMD with predictions
        if params["CMD"] == "1":
            fig, ax = plt.subplots(figsize=(9, 6))
            ax = plotting.cmd(self, ax, cfg)
            vega_color_ax, vega_mag_ax = plotting.add_vega_labels(ax, cfg)
            plotting.add_all_models(fig, ax,
                                    steal_axs=[ax, vega_color_ax, vega_mag_ax],
                                    cfg=cfg, models=self.low_res_models)
            self.figures.append(fig)

        # Do a quick and dirty initial redshift fitting, to get a starting
        # point.
        self.z[this_color] = self._initial_z(cfg)

        # If the user wants to see this initial fit, plot it.
        if params["fitting_procedure"] == "1":
            # set up the plot
            fig, ax = plt.subplots(figsize=(9, 6))
            ax = plotting.cmd(self, ax, cfg)
            plotting.add_vega_labels(ax, cfg)
            # the following line is UGLY, and will unfortunately be
            # repeated. The second parameter in the function call is the model
            # in the correct color at the current redshift of the cluster.
            plotting.add_one_model(ax, self.models[this_color][self.z[this_color].value], "k")
            plotting.add_redshift(ax, self.z[this_color].value)
            self.figures.append(fig)

        # do iterations of fitting, each with a progressively smaller
        # color cut, which is designed to hone in on the red sequence.
        for bluer_cut, redder_cut in zip(cfg["bluer_color_cut"],
                                         cfg["redder_color_cut"]):
            # set red sequence members based on the cuts
            self._set_rs_membership(self.z[this_color].value,
                                    bluer_cut, redder_cut,
                                    cfg["brighter_mag_cut"],
                                    cfg["dimmer_mag_cut"],
                                    cfg)

            # do the chi-squared fitting.
            self.z[this_color] = self._chi_square_w_error(cfg)

            # if the user wants, plot the procedure
            if params["fitting_procedure"] == "1":
                fig, ax = plt.subplots(figsize=(9, 6))
                ax = plotting.cmd(self, ax, cfg)
                plotting.add_vega_labels(ax, cfg)
                plotting.add_one_model(ax, self.models[this_color][self.z[this_color].value], "k")
                plotting.add_redshift(ax, self.z[this_color].value)
                self.figures.append(fig)

        # See if there is a red cloud, rather than a clear red sequence.
        self._clean_rs_check(cfg)

        # Set the final red sequence members
        self._set_rs_membership(self.z[this_color].value,
                                cfg["final_rs_color"][0],
                                cfg["final_rs_color"][1],
                                cfg["final_rs_mag"][0],
                                cfg["final_rs_mag"][1],
                                cfg)

        # and do the location check based on these RS values
        self._location_check(this_color)

        # we now have a final answer for the redshift of the cluster.
        # if the user wants, plot it up
        if params["final_CMD"] == "1":
            # set up the plot
            fig, ax = plt.subplots(figsize=(9, 6))
            ax = plotting.cmd(self, ax, cfg)
            plotting.add_vega_labels(ax, cfg)
            plotting.add_one_model(ax, self.models[this_color][self.z[this_color].value], "k")

            # I want to plot the models that match the 1 sigma errors, so we
            # first need to get those redshifts
            high_z = self.z[this_color].value + self.z[this_color].upper_error
            low_z = self.z[this_color].value - self.z[this_color].lower_error

            # then plot the models in light grey
            light_grey = "#AAAAAA"
            plotting.add_one_model(ax, self.models[this_color][high_z],
                                   light_grey)
            plotting.add_one_model(ax, self.models[this_color][low_z],
                                   light_grey)
            plotting.add_redshift(ax, self.z[this_color])

            self.figures.append(fig)

        # If the user wants, plot the location of the RS members.
        if params["location"] == "1":
            fig, ax = plt.subplots(figsize=(7, 7))
            ax = plotting.location(self, ax, this_color)
            plotting.add_redshift(ax, self.z[this_color])
            self.figures.append(fig)

        # interactive mode requires some more work
        if params["interactive"] == "1":
            # we'll show a two panel plot. The first panel shows the final CMD,
            # the second shows the location plot. This is because they are the
            # most useful for determining if the fitting worked and if the
            # cluster within is interesting
            fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[13, 6],
                                           tight_layout=True)
            # final CMD
            ax1 = plotting.cmd(self, ax1, cfg)
            plotting.add_vega_labels(ax1, cfg)
            plotting.add_redshift(ax1, self.z[this_color])
            plotting.add_one_model(ax1, self.models[this_color][self.z[this_color].value], "k")
            plotting.add_flags(ax1, self.flags[this_color])

            # location plot
            ax2 = plotting.location(self, ax2, this_color)
            plotting.add_redshift(ax2, self.z[this_color])
            plotting.add_flags(ax2, self.flags[this_color])

            # we close all other figures within the cluster object, so they
            # aren't shown to the user. They can still be saved later, though.
            for s_fig in self.figures:
                plt.close(s_fig)

            # We then show the two panel plot to the user, so they can
            # determine what to do with it.
            plt.show(block=False)
            while True:
                flags = raw_input("Enter the flags for this cluster [i/f/enter]: ")
                # validate user input, which can only be i, f, or enter
                if flags not in ["i", "f", ""]:
                    print "That is not a valid choice."
                else:
                    break
            plt.close(fig)

            # we have the user input, so do what we need to with it.
            if flags == "f":  # user flag
                self.flags[this_color] += 8
            elif flags == "i":  # interesting
                self.interesting = 1

    @staticmethod
    def _bin_edges(values, bin_size):
        """ Determines the bin edges given data and a bin size.

        The bin edges will span all the data. The lowest edge will be at the
        minimum of the data, and the highest one will be just above the max
        value of the values.

        :param values: List of data elements to be binned. Note that this
                     method doesn't actually do the binning, it just determines
                     what the appropriate bin edges are.
        :param bin_size: Size of the bins.
        :return: List of bin edges
        """

        # start with the lowest value in the data set as our lowest edge.
        edges = [min(values)]
        # we will add an edge at a time until we have spanned all the values
        while edges[-1] < max(values):
            edges.append(edges[-1] + bin_size)

        return edges

    def _centering(self, ras, decs):
        """
        Finds a guess for the center of the cluster if the user doesn't know.

        This is done in catalog space with a density based-approach. This
        assumes the cluster is the densest thing in the image, which is
        probably will be. Contamination won't affect this too much, hopefully.

        The algorithm starts by making the 2D histogram of the number of
        galaxies in each bin. This is then smoothed on a scale roughly
        equivalent to the size of the core of the cluster. Then the highest
        value is picked as the center of the cluster.

        :param ras: List of RA values we want to find a center for.
        :param decs: List of dec values we want to find a center for.
        :return: Guesses for the RA and dec of the center of the cluster.
        """

        # we want bin sizes of 5 arcseconds. This is small enough that we
        # can still get good resolution, but large enough that future
        # calculations won't take forever
        dec_bin_size = 5.0 / 3600.0  # convert to degrees
        # the ra bin size has to take into account the cos(dec) affect.
        # we can derive a size, assuming we want equally sized bins
        #                delta x = delta y
        #    delta ra * cos(dec) = delta dec
        #               delta ra = delta dec / cos(dec)
        # so since we have our dec bin size, we can get the ra bin size. We
        # will use the declination of the center of the image.
        middle_dec = (max(decs) + min(decs)) / 2.0
        ra_bin_size = dec_bin_size / np.cos(middle_dec * np.pi / 180.0)  # rads

        # we then make the bins themselves
        dec_bin_edges = Cluster._bin_edges(decs, dec_bin_size)
        ra_bin_edges = Cluster._bin_edges(ras, ra_bin_size)

        # we can now make the original density 2D histogram
        hist, _, _ = np.histogram2d(ras, decs,
                                    bins=[ra_bin_edges, dec_bin_edges])

        # we now want to smooth this histogram. I will smooth on 0.5 arcminute
        # scales. This will highlight the center of the cluster, where it is
        # densest. We need to figure out how big this is, but we already know
        # how big our bins are, so this is easy.
        smoothing_scale = 0.5 / (dec_bin_size * 60)  # in arcmin
        kernel = convolution.Gaussian2DKernel(stddev=smoothing_scale)

        # then we can smooth our 2D histogram
        smoothed_hist = convolution.convolve(hist, kernel)

        # then we need to get the max value of that, and turn it into coords
        max_idx = np.argmax(smoothed_hist)  # returns idx of flattened array
        # we turn this idx of flattened array into the appropriate 2D indices
        ra_idx = max_idx // smoothed_hist.shape[1]
        dec_idx = max_idx % smoothed_hist.shape[1]

        # since we know the bin size, we can turn this into a coordinate
        # we add half a bin size to get to the center of each bin
        ra_cen = min(ras) + ra_idx * ra_bin_size + ra_bin_size / 2.0
        dec_cen = min(decs) + dec_idx * dec_bin_size + dec_bin_size / 2.0

        self.center_ra = ra_cen
        self.center_dec = dec_cen

    def _location_cut(self, radius, params):
        """
        Does a location cut on the galaxies in the image.

        Only those within a certain radius are marked as being near_center,
        which is what we will use for the rest of the fitting. This throws
        away foregrounds, leaving the most dense portion of the cluster. A
        good radius to pick isn't necessarily the radius of the cluster.
        Picking a small radius can increase the signal to noise by only
        including the densest part of the cluster, and throwing out most
        foregrounds.

        If the user specified `dist` in their catalog, use that distance. If
        not, calculate our own center.

        :param radius: radius of the cut, in arcminutes.
        :return: None, but galaxies near the center are marked as having
                 near_center = True.
        """

        # first we need to get all the coordinates
        ras = [source.ra for source in self.sources_list]
        decs = [source.dec for source in self.sources_list]

        # we may need to find our own centers
        if params["dist"] == "-99":
            # get the center
            self._centering(ras, decs)

            # then find the distance from the center for each source.
            for source in self.sources_list:
                # Use pythagorean theorem to find distance in degrees, then
                # multiply by 3600 to convert to arcsec
                dec_radians = self.center_dec * np.pi / 180.0
                # we have to account for the cosine(dec) term in the
                # ra separation
                ra_sep = (source.ra - self.center_ra) * np.cos(dec_radians)
                dec_sep = source.dec - self.center_dec
                source.dist = np.sqrt(ra_sep**2 + dec_sep**2) * 3600

        # then set things as near the center if they are indeed near the center
        for source in self.sources_list:
            if source.dist < radius*60.0:  # convert radius to arcsec
                source.near_center = True
            else:
                source.near_center = False

    def _initial_z(self, cfg):
        """
        Find a decent initial redshift estimate, based on the number of
        galaxies near each model.

        It iterates through all redshifts, and counts the number of galaxies
        that are within some color of the model at that redshift. Then the
        redshift that has the highest number of galaxies is the
        redshift selected. It's a quick a dirty estimate.

        :param cfg: The configuration dictionary for the color combination used
        :return: an initial redshift estimate
        """

        # Get models with a large spacing, since we don't need a lot of
        # accuracy here.
        models = self.low_res_models[cfg["color"]]

        # set placeholder values that will be replaced as we go
        max_nearby = -999
        best_z = -999
        z_nearby_pairs = []
        # ^ will be list of tuples, where each tuple is the redshift and
        # number of nearby galaxies

        # Iterate through the redshifts
        for z in sorted(models):
            nearby = 0  # reset the number of nearby galaxies

            # get the m* at this redshift
            this_model = models[z]
            mag_point = this_model.mag_point

            # then iterate through the sources to see which are close to m*
            for source in self.sources_list:
                if source.near_center:
                    # get the mag and color of the source
                    source_mag = source.mags[cfg["red_band"]].value
                    source_color = source.colors[cfg["color"]].value

                    # get the expected RS color for a galaxy of this magnitude
                    rs_color = this_model.rs_color(source.mags[cfg["red_band"]].value)

                    # then determine the limits for a valid RS galaxy
                    bright_mag = mag_point - cfg["initial_mag"][0]
                    faint_mag = mag_point + cfg["initial_mag"][1]
                    blue_color = rs_color - cfg["initial_color"][0]
                    red_color = rs_color + cfg["initial_color"][1]
                    # see if it passes a color and magnitude cut
                    if (bright_mag < source_mag < faint_mag) and \
                       (blue_color < source_color < red_color):
                        # if it did pass, note that it is nearby
                        nearby += 1

            z_nearby_pairs.append((z, nearby))
            # replace the best values if this is the best
            if nearby > max_nearby:
                max_nearby = nearby
                best_z = z

        redshifts, nearbies = zip(*z_nearby_pairs)

        # check for the double red sequence
        self._double_red_sequence(nearbies, cfg["color"])

        # Turn this into a data object, with errors that span the maximum
        # range, since we don't have a good feeling yet
        up_error = sorted(models.keys())[-1] - best_z  # max z - best
        low_error = best_z - sorted(models.keys())[0]  # best - min z
        z = data.AsymmetricData(best_z, up_error, low_error)
        return z

    def _double_red_sequence(self, nearby, color):
        """Checks for the possibility of a double red sequence.

        Is to be used inside the initial_z function, since it uses the list of
        nearby galaxies as a function of redshift. It looks for two or more
        maxima in the number of red sequence galaxies simply by comparing each
        point to its neighbors.

        :param nearby: list of galaxies nearby to the red sequence at a given
                       redshift. Should be sorted in order of increasing
                       redshift. This is taken care of by the initial_z
                       function, though, so as long as it is used there
                       everything will be fine.
        :param color: Color combination being used this time through the
                      fitting. Needed to set the appropriate flag.
        :returns: None, but does set the cluster's flag variable if need be.
        """
        num_local_maxima = 0
        # We will compare each one to its 3 neighbors, so start as the third
        # index.
        idx = 3
        # we can't go past the 4th to last item.
        while 3 <= idx <= len(nearby) - 4:
            item = nearby[idx]
            # compare it to its three neighbors on each side.
            if item > nearby[idx-3] and item > nearby[idx-2] and \
               item >= nearby[idx-1] and item >= nearby[idx+1] and \
               item > nearby[idx+2] and item > nearby[idx+3]:
                # if it is bigger than all 6 neighbors, call it a maxima.
                num_local_maxima += 1
                # move 3 ahead, to avoid the possibility of two maxima being
                # right next to each other, which can happen if the top is
                # flat. That obviously isn't what we are looking for here.
                idx += 3
            else:
                idx += 1
        # If there are 2 or more maxima, increment the flag.
        if num_local_maxima >= 2:
            self.flags[color] += 2

    def _set_rs_membership(self, redshift, bluer, redder,
                           brighter, dimmer, cfg):
        """Set some sources to be RS members, based on color and mag cuts.

        :param redshift: redshift of the red sequence that these sources
                         will be identified as belonging to.
        :param bluer: maximum color difference on the blue side from the
                      characteristic color of the red sequence at the
                      magnitude of the galaxy.
        :param redder: maximum color difference on the red side from the
                       characteristic color of the RS
        :param brighter: maximum magnitude difference (on the bright end)
                         from the characteristic magnitude of the red
                         sequence.
        :param dimmer: maximum magnitude difference (on the faint end)
                         from the characteristic magnitude of the red
                         sequence.
        :param cfg: Configuration dictionary for this color combination.

        As an example: At some redshift, the red sequence has a characteristic
        magnitude of 20, and a characteristic color of zero. We pass in
        bluer=0.1, redder=0.2, brighter=2.0, dimmer=1.0. Any galaxies that
        have magnitude 18 to 21 pass the magnitude cut. The color cut is
        trickier, since the RS has some slope. We find the characteristic
        color of the red sequence at the magnitude of the galaxy (say
        0.05 for a particular galaxy). Then if the color is within the
        bounds set by bluer and redder (in this case -0.05 to 0.25), then it
        passes the color cut. If a galaxy passes both the magnitude and
        color cuts, it is marked as a red sequence galaxy.

        Note: This function marks galaxies both inside and outside the
        location cut as RS members. This is because the location cut may be
        smaller than the cluster (in fact this is often good to do), and we
        don't want to mark galaxies as not RS members just because they are
        outside the (possibly) small location cut.

        :return: None, but galaxies that are RS members are marked as such.
        """

        # get the model, it's characteristic magnitude, and then turn it
        # into magnitude limits based on the parameters passed in
        rs_model = self.models[cfg["color"]][redshift]
        char_mag = rs_model.mag_point
        dim_mag = char_mag + dimmer
        bright_mag = char_mag - brighter

        for source in self.sources_list:
            # get the color correspoinding to the red sequence at the ch2
            # magnitude of this particular source
            char_color = rs_model.rs_color(source.mags[cfg["red_band"]].value)
            # turn it into color limits based on parameters passed in
            red_color = char_color + redder
            blue_color = char_color - bluer
            # a function in the source class does the actual marking
            source.rs_membership(blue_color, red_color, bright_mag, dim_mag,
                                 cfg["color"], cfg["red_band"])

    def _chi_square_w_error(self, cfg):
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
        to_fit = [source for source in self.sources_list
                  if source.RS_member[cfg["color"]] and source.near_center]

        # if there isn't enough to fit to, keep the initial z. We need at
        # least two objects to do the chi squared calculation. This is because
        # we divide by number of objects - 2. More comments explain this when
        # this is done 20 lines below.
        if len(to_fit) <= 2:
            return self.z[cfg["color"]]

        # initialize lists for the chi squared distribution and the redshifts
        chi_sq_values = []
        redshifts = []

        # test each model
        for z in sorted(self.models[cfg["color"]]):
            chi_sq = 0  # reset the chi squared value for this model
            this_model = self.models[cfg["color"]][z]
            for source in to_fit:
                # get the model points, then compare that with the data
                model_color = this_model.rs_color(source.mags[cfg["red_band"]].value)
                color = source.colors[cfg["color"]].value
                error = source.colors[cfg["color"]].error
                chi_sq += ((model_color - color)/error)**2

            # reduce the chi square values. We will divide by the degrees of
            # freedom, which = number of data points - number of parameters = 1
            # here we have only 1 parameter (redshift), so it is just
            # number of data points - 2. We can't divide by zero, and a
            # number wouldn't make any sense, so we need at least 3 objects,
            # which is why we have the check above.
            chi_sq /= (len(to_fit) - 2)
            # put these values into the lists, so we can keep track of them
            chi_sq_values.append(chi_sq)
            redshifts.append(z)

        # get the minimum chi squared value, and it's index in the list.
        best_chi = min(chi_sq_values)
        best_chi_index = chi_sq_values.index(best_chi)

        # Start finding errors, using the chi squared distribution. Where it's
        # one greater than the best chi, that's where the error is.
        high_idx = best_chi_index
        # find the place where the chi value is one greater than the best
        while high_idx < len(chi_sq_values) - 1 and \
                chi_sq_values[high_idx] - best_chi <= 1.0:
            # Note: len(chi_sq_values) - 1 needs that -1 so that we don't
            # end up with an high_idx that is past the end of chi_sq_values
            high_idx += 1
        # high_idx will now point to the first value where the chi squared
        # value is one greater than the best fit value, or to the edge of the
        # redshift range if no limit was found.

        # do the same thing for the low error
        low_idx = best_chi_index
        while low_idx > 0 and chi_sq_values[low_idx] - best_chi <= 1.0:
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

    def _clean_rs_check(self, cfg):
        """ Determine whether or not there is a clean red sequence.

        This works by counting the number of red sequence galaxies of the best
        fit model, then comparing that to the number of red sequence galaxies
        if the red sequence were bluer or redder. If the best fit doesn't have
        a lot more than the non-best fits, then we mark the flag for a
        not clean red sequence.

        :param cfg: Configuration dictionary for the color of interest.
        :return: None, but the flag is incremented if there isn't a clean red
                 sequence.
        """

        # set cuts to be used each time
        bluer = cfg["final_rs_color"][0]
        redder = cfg["final_rs_color"][1]
        brighter = cfg["final_rs_mag"][0]
        dimmer = cfg["final_rs_mag"][1]

        # set the red sequence members for the accepted red sequence.
        self._set_rs_membership(self.z[cfg["color"]].value,
                                bluer=bluer, redder=redder,
                                brighter=brighter, dimmer=dimmer,
                                cfg=cfg)
        # count the RS members in the best fit
        best_rs = self._count_galaxies(cfg["color"])

        # set red sequence members to be in the color range redder then the
        # best fit red sequence. The blue limit is where the red limit used
        # to be, and the red limit is the same distance from the blue limit
        # that is used to be. This essentially creates an adjacent red
        # sequence that we can compare to.
        # Translated into values, the blue limit goes at the negative value of
        # the old red limit (since it's in the opposite direction blue
        # normally is). Red turns into the old red limit plus the distance,
        # which is simply the sum of red and blue limits. Then this reverses
        # for the bluer RS

        red_bluer = -1 * redder
        red_redder = 2 * redder + bluer

        blue_redder = -1 * bluer
        blue_bluer = 2 * bluer + redder

        # red RS cut
        self._set_rs_membership(self.z[cfg["color"]].value,
                                bluer=red_bluer, redder=red_redder,
                                brighter=brighter, dimmer=dimmer,
                                cfg=cfg)
        # again count the galaxies in this red sequence
        red_rs = self._count_galaxies(cfg["color"])

        # blue RS cut
        self._set_rs_membership(self.z[cfg["color"]].value,
                                bluer=blue_bluer, redder=blue_redder,
                                brighter=brighter, dimmer=dimmer,
                                cfg=cfg)
        blue_rs = self._count_galaxies(cfg["color"])

        # Compare the numbers in the 3 red sequences. Set the flag if the
        # number of galaxies in the best red sequence is less than 1.5 times
        # the sum of the two offset red sequences. The 1.5 times is arbitrary.
        # It was chosen to make things I thought looked bad have this flag.
        if ((red_rs + blue_rs) * 1.5) >= best_rs:
            self.flags[cfg["color"]] += 4

    def _count_galaxies(self, color):
        """ Counts the number of red sequence galaxies near the center.

        :return: number of red sequence galaxies that are near the center.
        :rtype: float. I return float so that it can be used in division in
                Python 2.7 smartly.
        """
        count = 0
        for source in self.sources_list:
            if source.near_center and source.RS_member[color]:
                count += 1
        return float(count)

    def _location_check(self, color):
        """Looks to see if the red sequence galaxies are concentrated in the
        middle. If they are not, this will raise a flag.

        It works be getting the percent of galaxies within the location cut
        that are red sequence members, then comparing that to the percentage
        outside the location cut. If the percentage inside isn't 75 percent
        higher, then the flag is raised.

        :param color: Color used for this fitting. Is needed to set the flag.
        :returns: None, but the flag is set if need be.
        """
        # find how many objects are near then center and not near the center
        # that are or are not red sequence members. I converted everything
        # to floats to avoid the rounding that comes from dividing integers
        # in Python 2.x
        total_near_center = float(len([source for source in self.sources_list
                                       if source.near_center]))
        rs_near_center = float(len([source for source in self.sources_list
                                    if source.near_center and
                                    source.RS_member[color]]))

        total_not_near_center = float(len([source for source in
                                           self.sources_list if
                                           not source.near_center]))
        rs_not_near_center = float(len([source for source in self.sources_list
                                        if not source.near_center and
                                        source.RS_member[color]]))

        # Calculate the percent of sources that are red sequence members both
        # near and far from the center.
        try:
            near_rs_percent = rs_near_center / total_near_center
            not_near_rs_percent = rs_not_near_center / total_not_near_center

            # raise the flag if the percent near the center isn't high enough.
            # "high enough" is arbitrary, and can be adjusted.
            if near_rs_percent <= not_near_rs_percent * 1.75:
                self.flags[color] += 1
        except ZeroDivisionError:  # there were zero sources near the center
            self.flags[color] += 1

    def rs_catalog(self, params):
        """ Writes a catalog of all the objects, indicating the RS members.

        :param params: User's parameter dictionary.
        :returns: none, but the catalog is saved.
        """

        # get the path where the file will be saved. It will be of the form
        # name.rs.cat
        filepath = params["rs_catalog_dir"] + os.sep + self.name + ".rs" + \
            params["extension"]

        try:
            cat = open(filepath, "w")
        except IOError:
            raise IOError("The location to save the RS catalogs could not\n"
                          "\tbe located. Make sure you specified the \n"
                          "\t'rs_catalog_dir' parameter appropriately.\n"
                          "\tThe code will create new files, but not new "
                          "directories.")

        # get the type of info, and make a header using that
        d_type = params["type"]

        # first add ra and dec
        header = "# {:<12s} {:<12s}".format("ra", "dec")

        # also make a general formatters that will be used later for ra/dec
        coords_formatter = "  {:<12.7f} {:<12.7f}"

        phot_formatter = " {:<15.3f} {:<15.3f}"
        for band in self.sources_list[0].mags:
            band = band.replace("sloan_", "")  # we don't care about the sloan
            # add these bands to the header
            header += " {:<15s} {:<15s}".format(band + "_" + d_type,
                                                band + "_" + d_type + "_err")

        center_formatter = " {:<7}"
        header += center_formatter.format("center")

        rs_formatter = " {:<10}"
        for color in self.z:
            header += rs_formatter.format("RS_" + color.replace("sloan_", ""))

        # I want the headers to line up over the data
        cat.write(header + "\n")

        # then we can write all the sources to the catalog.
        for source in self.sources_list:
            # start by adding the ra/dec
            line = coords_formatter.format(source.ra, source.dec)

            # then add all the photometry
            for band in source.mags:
                mag = source.mags[band].value
                magerr = source.mags[band].error
                if d_type == "mags":
                    # we don't have to convert to flux
                    if params["mag_system"] == "ab":
                        line += phot_formatter.format(mag, magerr)
                    else:
                        mag += config.ab_to_vega[band]  # convert to Vega
                        line += phot_formatter.format(mag, magerr)
                else:
                    # convert to flux before writing to catalog
                    flux = self.mag_to_flux(mag, params["mag_zeropoint"])
                    fluxerr = self.mag_errors_to_flux_errors(magerr, flux)

                    line += phot_formatter.format(flux, fluxerr)

            # add whether or not it's centered
            line += center_formatter.format(source.near_center)

            # then add RS information
            for color in self.z:
                line += rs_formatter.format(source.RS_member[color])

            line += "\n"
            cat.write(line)
        cat.close()


def save_as_one_pdf(figs, filename):
    """
    Save the figures into one long PDF file

    :param figs: list of figures to be saved as PDFs
    :param filename: place where the PDFs will be saved
    :return: none
    """

    # check if there are actually figures to save.
    if len(figs) == 0:
        return

    # if so, save them.
    try:
        pp = PdfPages(filename)
    except IOError:
        raise IOError("The location to save the plots could not be found.\n"
                      "\tMake sure you specified the 'plot_directory'\n"
                      "\tparameter appropriately. The code will create\n"
                      "\tnew files, but not new directories.")

    for fig in figs:
        pp.savefig(fig)
    pp.close()
    for fig in figs:
        plt.close(fig)
