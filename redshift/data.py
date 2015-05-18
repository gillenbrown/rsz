class data(object):
    """
    Class that represents a data point. Has a value, as well as errors. Can
    have a single error value, or upper and lower errors.
    """
    def __init__(self, value, error):
        self.value = value
        self.error = error

    def __repr__(self):
        return str(self.value) + u'\u00b1' + str(self.error)

    def __add__(self, other):
        # TODO:

