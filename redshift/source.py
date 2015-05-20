class Source(object):
    """
    Class representing a source (as detected by something like SExtractor).
    For this implementation, only holds data in ch1 and ch2, since that's
    all we need for this implementation.
    """
    def __init__(self, ra, dec, ch1_mag, ch2_mag):
        """
        Constructor. Pass in Data class objects for the magnitudes if you
        want to include errors.
        :param ra: ra of the object
        :param dec: declination of the object
        :param ch1_mag: ch1 mag of the object. pass in a Data object if you
                        want to include errors.
        :param ch2_mag: ch2 mag of the object. pass in a Data object if you
                        want to include errors.
        """
        self.ra = float(ra)
        self.dec = float(dec)
        # Don't convert the magnitudes to floats, since they will be data
        # objects.
        self.ch1 = ch1_mag
        self.ch2 = ch2_mag
        self.ch1_m_ch2 = self.ch1 - self.ch2

        self.near_center = False
        self.RS_member = False