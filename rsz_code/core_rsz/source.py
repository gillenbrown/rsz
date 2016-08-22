class Source(object):
    """
    Class representing a source (as detected by something like SExtractor).
    """
    def __init__(self, ra, dec, mags, dist=None):
        """
        Constructor. Pass in Data class objects for the magnitudes if you
        want to include errors.

        Convert things to numeric data types (ie float, data, int). We will do
        math on them, which won't work if they are strings.

        :param ra: ra of the object
        :type ra: float
        :param dec: declination of the object
        :type dec: float
        :param mags: Dictionary with magnitude information. Keys are bands,
                     values are the data objects with a value and error.
        :type mags: dict
        :param dist: distance of the galaxy from the overdensity center in
                     arcseconds.
        :type dist: float
        """
        self.ra = ra
        self.dec = dec

        self.mags = mags
        self._calculate_colors()

        self.dist = dist

        self.near_center = False
        self.RS_member = dict()

    def _calculate_colors(self):
        """
        Calculates color information from all the mag information. This
        creates all possible colors, even if some of them are meaningless.
        This won't matter in the long run, since those meaningless colors
        won't get called.

        :return: None, but self.colors is initialized appropiately
        """
        self.colors = dict()
        for band_1 in self.mags:
            for band_2 in self.mags:
                color = "{}-{}".format(band_1, band_2)
                self.colors[color] = self.mags[band_1] - self.mags[band_2]

    def rs_membership(self, blue, red, bright, faint, color, red_band):
        """Mark sources as red sequence members if they pass the given cuts.

        Sources will be marked as red sequence members if they have a color
        between blue and red, and a magnitude between bright and faint.

        Also does an error cut on the color.

        :param color: Which color we are using to determin RS membership.
        :param blue: bluest color a source can be to be a RS member
        :param red: reddest color a source can be to be a RS member
        :param bright: brightest magnitude "" "" "" "" "" "" ""
        :param faint: dimmest magnitude "" "" "" "" "" "" ""
        :return: None, but some sources will be marked as RS members.
        """
        if blue < self.colors[color] < red and \
                bright < self.mags[red_band] < faint and \
                self.colors[color].error < 0.2:
            self.RS_member[color] = True
        else:
            self.RS_member[color] = False
