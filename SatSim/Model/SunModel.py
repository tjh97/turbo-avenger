"""

"""

import AbstractObjectModel
import basics.vector_types as vt
import basics.constants as constant
import math


class SunModel(AbstractObjectModel.AbstractObjectModel):
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
        """Calculates the direction and distance to the sun from the Earth

        :param dt.datetime time: Update the sun model to the time given.
        :rtype: vt.Vector, float
        :return: Unit sun vector, distance from the Earth to the Sun
        """

        delta_J2000 = self.time - constant.J2000_DATE
        n_days_J2000 = delta_J2000.days + delta_J2000.seconds/86400

        mean_lon_sun = 280.460 + 0.9856474*n_days_J2000
        mean_lon_sun %= 360.0
        mean_lon_sun *= constant.DEG_TO_RAD

        mean_anomaly_sun = 357.528 + 0.9856003*n_days_J2000
        mean_anomaly_sun %= 360.0
        mean_anomaly_sun *= constant.DEG_TO_RAD

        ecliptic_lon_sun = ( mean_lon_sun/constant.DEG_TO_RAD +
                             1.915*math.sin(mean_anomaly_sun) +
                             0.020*math.sin(2.0*mean_anomaly_sun) )
        ecliptic_lon_sun *= constant.DEG_TO_RAD

        dist_earth_to_sun = (1.00014 -
                             0.01671*math.cos(mean_anomaly_sun) -
                             0.00014*math.cos(2.0*mean_anomaly_sun) )
        dist_earth_to_sun *= constant.AU_TO_KM

        obliquity_ecliptic = 23.439 - 0.0000004*n_days_J2000
        obliquity_ecliptic *= constant.DEG_TO_RAD

        x_J2000_sun = math.cos(ecliptic_lon_sun)
        y_J2000_sun = math.cos(obliquity_ecliptic)*math.sin(ecliptic_lon_sun)
        z_J2000_sun = math.sin(obliquity_ecliptic)*math.sin(ecliptic_lon_sun)

        self.direction = vt.Vector([x_J2000_sun, y_J2000_sun, z_J2000_sun])
        self.distance = dist_earth_to_sun
        self.time = time
