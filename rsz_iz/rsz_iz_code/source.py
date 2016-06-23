class Source(object):
    """
    Class representing a source (as detected by something like SExtractor).
    For this implementation, only holds data in r and z, since that's
    all we need.
    """
    def __init__(self, ra, dec, i_mag, z_mag, dist=None, imz=None):
        """
        Constructor. Pass in Data class objects for the magnitudes if you
        want to include errors.

        Convert things to numeric data types (ie float, data, int). We will do
        math on them, which won't work if they are strings.

        :param ra: ra of the object
        :param dec: declination of the object
        :param i_mag: i band mag of the object. pass in a Data object if you
                      want to include errors.
        :param z_mag: z band mag of the object. pass in a Data object if you
                        want to include errors.
        :param dist: distance of the galaxy from the overdensity center in
                     arcseconds.
        :param rmz: r-z color of the object.
        """
        self.ra = ra
        self.dec = dec
        self.i_mag = i_mag
        self.z_mag = z_mag
        if imz is not None:
            self.imz = imz
        else:
            self.imz = self.i_mag - self.z_mag

        self.dist = dist

        self.near_center = False
        self.RS_member = False

    def rs_membership(self, blue, red, bright, faint):
        """Mark sources as red sequence members if they pass the given cuts.

        Sources will be marked as red sequence members if they have a color
        between blue and red, and a z magnitude between bright and faint.

        Also does an error cut on the color.

        :param blue: bluest color a source can be to be a RS member
        :param red: reddest color a source can be to be a RS member
        :param bright: brightest z magnitude "" "" "" "" "" "" ""
        :param faint: dimmest z magnitude "" "" "" "" "" "" ""
        :return: None, but some sources will be marked as RS members.
        """

        if blue < self.imz < red and bright < self.z_mag < faint:# and \
               # self.rmz.error < 0.2:
            self.RS_member = True
        else:
            self.RS_member = False
