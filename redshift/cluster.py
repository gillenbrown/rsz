import os

import model

class Cluster(object):
    # Keeping the predictions for the red sequence with the cluster object made things a lot easier. And since the
    # predictions are the same for all cluster, this is a class variable rather than an instance variable.
    predictions_dict = model.make_prediction_dictionary(0.01)

    def __init__(self, filepath):
        self.name = Cluster._name(filepath)

        # I'll initialize emtpy source list, that will be filled as we go





    @staticmethod
    def _name(filepath):
        filename = os.path.split(filepath)
        print filename
        return "YELLO"

