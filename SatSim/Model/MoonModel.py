"""

"""

import AbstractObjectModel
import basics.vector_types as vt
import basics.constants as constant
import math


class MoonModel(AbstractObjectModel.AbstractObjectModel):
    """This class models tue location of the sun in J2000 ECI frame.

    :type time: datetime.datetime
    :type direction: vt.Vector
    :type distance: float
    """

    def __init__(self, time):
        """

        :param datetime.datetime time: Initial time to start the Sun model
        :return: None
        """
        self.time = time
        self.direction = vt.Vector([0, 0, 0])
        self.distance = 0.0
        self.update(time)   # Update direction, distance to their correct values

    def update(self, time):
        """Calculates the direction and distance to the moon from the Earth

        :param datetime.datetime time: Update the moon model to the time given.
        :rtype: vt.Vector, float
        :return: Unit moon vector, distance from Earth to Moon
        """

        delta_J2000 = self.time - constant.J2000_DATE
        n_days_J2000 = delta_J2000.days + delta_J2000.seconds/86400

        mean_lon_moon = 218.316 + 13.176396*n_days_J2000
        mean_lon_moon %= 360.0
        mean_lon_moon *= constant.DEG_TO_RAD

        mean_anomaly_moon = 134.963 + 13.064993*n_days_J2000
        mean_anomaly_moon %= 360.0
        mean_anomaly_moon *= constant.DEG_TO_RAD

        mean_dist_moon = 93.272 + 13.229350*n_days_J2000
        mean_dist_moon %= 360.0
        mean_dist_moon *= constant.DEG_TO_RAD

        ecliptic_lon_moon = (mean_lon_moon/constant.DEG_TO_RAD +
                             6.289*math.sin(mean_anomaly_moon) )
        ecliptic_lon_moon *= constant.DEG_TO_RAD

        ecliptic_lat_moon = 5.128*math.sin(mean_dist_moon)
        ecliptic_lat_moon *= constant.DEG_TO_RAD

        self.distance = (385001.0 - 20905.0*math.cos(mean_anomaly_moon))

        obliquity_ecliptic = 23.439 - 0.0000004*n_days_J2000
        obliquity_ecliptic *= constant.DEG_TO_RAD

        x_ec = math.cos(ecliptic_lat_moon)*math.cos(ecliptic_lon_moon)
        y_ec = math.cos(ecliptic_lat_moon)*math.sin(ecliptic_lon_moon)
        z_ec = math.sin(ecliptic_lat_moon)
        rect_ec = vt.Vector([[x_ec],[y_ec],[z_ec]])

        Q_eq_ec = vt.Matrix([[1, 0, 0],
                             [0, math.cos(obliquity_ecliptic), -math.sin(obliquity_ecliptic)],
                             [0, math.sin(obliquity_ecliptic), math.cos(obliquity_ecliptic)]])

        self.direction = Q_eq_ec*rect_ec
        self.time = time
