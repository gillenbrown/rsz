import os

import model
import source
import data

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
                    ra = entries[0]
                    dec = entries[1]
                    ch1 = data.Data(entries[2], entries[3])
                    ch2 = data.Data(entries[4], entries[5])

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


    def fit_z(self):
        """
        Find the redshift of the cluster by matching its red sequence the
        red sequence models produced by EzGal. This is the main
        functionality of the entire code.

        :return: redshift of the cluster. Will be an AsymmetricData
        instance, to account for the possibility of different upper and
        lower limits on the redshift.
        """

        # initialize a list of figures, that will be filled as we go along.
        figures = []

        # do a location cut, to only focus on galaxies near the center of
        # the image.
        self._location_cut(1.5)
