class Model(object):
    """
    Class storing data from the EzGal models.
    """
    @staticmethod
    def make_slopes():
        """
        Find the slope of the red sequence as a function of redshfit.

        Uses the slopes of the red sequence as published in Eisenhardt 2007.

        :return: WHO KNOWS? I DON'T
        """

        # first we need to see at what redshift does ch1-ch2 see the various
        # Eisenhardt colors.
        # I did this in an iPython notebook that is in the
        # calculating_slopes folder, if you want to see the work. I'll just
        # manually add the results here, to avoid unecessary computation
        redshifts = [1, 2, 3]
        slopes = [-0.1, -0.04, -0.08]

        # TODO: do some interpolating here to find the slope at any redshift
        def f(z):
            return 0.05

        return f

    #get slopes for all redshifts
    slope = make_slopes()





    def predict_dict(self, spacing):
         pass