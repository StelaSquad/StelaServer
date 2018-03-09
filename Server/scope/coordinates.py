class Coordinates:
    """Coordinates of the Telescope. Telescope's have the
    following properties:

    Attributes:
        azimuth:  angular distance of an object from the local North, measured along the horizon.
        altitude: angular distance of an object above the local horizon.
    """

    def __init__(self, azimuth = 0.0, altitude=0.0):
        """Return a Customer object whose name is *name* and starting
        balance is *balance*."""
        self.azimuth = azimuth
        self.altitude = altitude

    def updateAzimuth(self, deltaAz):
        """Updates the Azimuth angle after moving a delta amount."""
        #TODO: once we go other 180? reset back to zero... don't want to overflow
        self.azimuth += deltaAz

    def updateAltitude(self, deltaAl):
        """Updates the Altitude angle after moving a delta amount."""
        self.altitude += deltaAl