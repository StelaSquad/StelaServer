from scope import coordinates

class Telescope:
    """Telescope... STELA. Telescope's have the
    following properties:

    Attributes:
        coordinates: coordinates of the telescope
    """

    def __init__(self, azimuth = 0.0, altitude = 0.0):
        """Return a Customer object whose name is *name* and starting
        balance is *balance*."""
        # TODO: update this to better reflect what telescope will look like.
        #   ex: add something for motor configuration.
        self.coordinates = coordinates.Coordinates(azimuth, altitude)

    def rotateAz(self, deltaAz):
        """Return the Azimuth angle after moving a delta amount."""
        self.coordinates.updateAzimuth(deltaAz)

    def rotateAl(self, deltaAl):
        """Return the Altitude angle after moving a delta amount."""
        self.coordinates.updateAzimuth(deltaAl)