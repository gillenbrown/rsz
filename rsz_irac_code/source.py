class Source(object):
    """
    Class representing a source (as detected by something like SExtractor).
    For this implementation, only holds data in ch1 and ch2, since that's
    all we need for this implementation.
    """
    def __init__(self, ra, dec, ch1_mag, ch2_mag, dist, ch1mch2=None):
        """
        Constructor. Pass in Data class objects for the magnitudes if you
        want to include errors.
        :param ra: ra of the object
        :param dec: declination of the object
        :param ch1_mag: ch1 mag of the object. pass in a Data object if you
                        want to include errors.
        :param ch2_mag: ch2 mag of the object. pass in a Data object if you
                        want to include errors.
        :param dist: distance of the galaxy from the overdensity center.
        """
        self.ra = ra
        self.dec = dec
        # Don't convert the magnitudes to floats, since they will be data
        # objects.
        self.ch1 = ch1_mag
        self.ch2 = ch2_mag
        if ch1mch2 is not None:
            self.ch1_m_ch2 = ch1mch2
        else:
            self.ch1_m_ch2 = self.ch1 - self.ch2

        self.dist = dist

        self.near_center = False
        self.RS_member = False


    def RS_membership(self, blue, red, bright, faint):
        """Mark sources as red sequence members if they pass the given cuts.

        Sources will be marked as red sequence members if they have a color
        between blue and red, and a ch2 magnitude between bright and faint.

        :param blue: bluest color a source can be to be a RS member
        :param red: reddest color a source can be to be a RS member
        :param bright: brightest ch2 magnitude "" "" "" "" "" "" ""
        :param faint: dimmest ch2 magnitude "" "" "" "" "" "" ""
        :return: None, but some sources will be marked as RS members.
        """

        if blue < self.ch1_m_ch2 < red and bright < self.ch2 < faint \
                and self.ch1_m_ch2.error < 0.2:
            self.RS_member = True
        else:
            self.RS_member = False
