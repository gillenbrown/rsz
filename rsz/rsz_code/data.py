import math


class Data(object):
    """
    Class that represents a data point. Has a value, as well as errors.

    Addition and subtraction operators are implemented smartly, too.
    """

    def __init__(self, value, error):
        self.value = value
        self.error = error

    def __repr__(self):
        return str(self.value) + "+-" + str(self.error)
        # u'\u00b1' is the code for +- in Unicode, if I can get that to work.

    def __add__(self, other):
        """
        Adds two data objects. The values are added, and the errors are
        added in quadrature. Another data object is returned.
        """
        if type(other) is Data:
            new_value = self.value + other.value
            new_error = math.sqrt(self.error**2 + other.error**2)
            return Data(new_value, new_error)
        else:  # the other object won't have errors, so we can't do anything
            #  with that.
            new_value = self.value + other
            return Data(new_value, self.error)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """ Same as add, but subtracting."""
        if type(other) is Data:
            new_value = self.value - other.value
            new_error = math.sqrt(self.error**2 + other.error**2)
            return Data(new_value, new_error)
        else:  # the other object won't have errors, so we can't do anything
            # with that.
            new_value = self.value - other
            return Data(new_value, self.error)

    def __rsub__(self, other):
        result = self - other
        result.value *= -1
        return result

    # define a bunch of comparison operators, using only > and ==
    # errors are ignored here, only the values are compared.
    # Note: comparisons between different data types can work strangely in
    # Python 2.x, so these may return weird stuff. Be careful comparing
    # to class types that aren't numeric in nature.
    def __lt__(self, other):  # less than:  self < other
        if isinstance(other, Data):
            return self.value < other.value
        else:
            return self.value < other

    def __eq__(self, other):  # equality self == other
        if isinstance(other, Data):
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


class AsymmetricData(Data):
    """
    Class that represents a data point, with both a value and upper and lower
    errors.

    Addition and subtraction operators are implemented smartly.
    """
    def __init__(self, value, upper_error, lower_error):
        # create a data object, but with None as the error, since it is
        # meaningless with asymmetric errors
        Data.__init__(self, value=value, error=None)
        self.upper_error = upper_error
        self.lower_error = lower_error

    def __repr__(self):
        return str(self.value) + "+" + str(self.upper_error) + \
            "-" + str(self.lower_error)

    def __add__(self, other):
        """
        Addition where upper and lower errors are added in quadrature
        """

        if type(other) is AsymmetricData:
            new_value = self.value + other.value
            new_upper_error = math.sqrt(self.upper_error**2 +
                                        other.upper_error**2)
            new_lower_error = math.sqrt(self.lower_error**2 +
                                        other.lower_error**2)
            return AsymmetricData(new_value, new_upper_error,
                                  new_lower_error)
        elif type(other) is Data:
            new_value = self.value + other.value
            new_upper_error = math.sqrt(self.upper_error**2 +
                                        other.error)**2
            new_lower_error = math.sqrt(self.lower_error**2 +
                                        other.error**2)
            return AsymmetricData(new_value, new_upper_error,
                                  new_lower_error)
        else:   # the other object won't have errors, so we can't do
                # anything with that.
            new_value = self.value + other
            return AsymmetricData(new_value, self.upper_error,
                                  self.lower_error)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """
        Subtraction where upper and lower errors are added in quadrature
        """

        if type(other) is AsymmetricData:
            new_value = self.value - other.value
            new_upper_error = math.sqrt(self.upper_error**2 +
                                        other.upper_error**2)
            new_lower_error = math.sqrt(self.lower_error**2 +
                                        other.lower_error**2)
            return AsymmetricData(new_value, new_upper_error,
                                  new_lower_error)
        elif type(other) is Data:
            new_value = self.value - other.value
            new_upper_error = math.sqrt(self.upper_error**2 +
                                        other.error**2)
            new_lower_error = math.sqrt(self.lower_error**2 +
                                        other.error**2)
            return AsymmetricData(new_value, new_upper_error,
                                  new_lower_error)
        else:   # the other object won't have errors, so we can't do
                # anything with that.
            new_value = self.value - other
            return AsymmetricData(new_value, self.upper_error,
                                  self.lower_error)

    def __rsub__(self, other):
        result = self - other
        result.value *= -1
        return result

    # comparison operators will be unchanged compared to the symmetric data
