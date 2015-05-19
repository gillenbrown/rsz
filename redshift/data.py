import math

class Data(object):
    """
    Class that represents a data point. Has a value, as well as errors. Can
    have a single error value, or upper and lower errors.
    """
    def __init__(self, value, error):
        self.value = value
        self.error = error

    def __repr__(self):
        return str(self.value) + "+-" + str(self.error)
        # u'\u00b1' is the code for +- in Unicode, if I can ever get that working

    def __add__(self, other):
        """
        Adds two data objects. The values are added, and the errors are
        added in quadrature. Another data object is returned.
        """
        new_value = self.value + other.value
        new_error = math.sqrt((self.error)**2 + (other.error)**2)
        return Data(new_value, new_error)

    def __sub__(self, other):
        """ Same as add, but subtracting."""
        new_value = self.value - other.value
        new_error = math.sqrt((self.error)**2 + (other.error)**2)
        return Data(new_value, new_error)

    
    # define a bunch of comparison operators, using only > and ==
    # errors are ignored here, only the values are compared.
    # Note: comparisons between different data types can work strangely in
    # Python 2.x, so these may return weird stuff. Be careful comparing
    # to class types that aren't numeric in nature.
    def __lt__(self, other):  # less than:  self < other
        if type(other) is Data:
            return self.value < other.value
        else:
            return self.value < other

    def __eq__(self, other):  # equality self == other
        if type(other) is Data:
            return self.value == other.value
        else:
            return self.value == other

    def __ge__(self, other):  # greater than or equal to: self >= other
        return not self < other

    def __ne__(self, other):  # non-equality: self != other
        return not self == other

    def __gt__(self, other):  # greater than: self > other
        return (not self < other) and (not self == other)

    def __le__(self, other):  # less than or equal to: self <= other
        return (self < other) or (self == other)