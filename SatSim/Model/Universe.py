'''
Created on Apr 21, 2015

@author: tjh97
'''
# from Simulator.Simulator import Simulator

import datetime as dt
from math import cos, sin, pi, sqrt, isnan
from docutils.nodes import math
import numpy as np
import scipy.special as sp
import basics.vector_types as vt
import basics.conversions as conv
import SatModel
import AbstractForceModel as afm


class Universe(object):
    """
    classdocs
    """

    # -------------------
    #  Public Constants
    # -------------------
    MU = 398600.4418                # km^3/s^2
    UNIVERSE_TIME_STEP = 0.01       # sec
    RADIUS_OF_EARTH = 6378.1363     # km

    # -------------------
    #  Private Constants
    # -------------------
    __J2000_DATE = dt.datetime(2000, 1, 1, 12, 0, 0)     # J2000 Epoch starts at noon on 1/1/2000
    __DEG_TO_RAD = pi/180.0
    __AU_TO_KM = 149597871.0

    # -------------
    #  Constructor
    # -------------
    def __init__(self, start=None):
        """Initialize the model at start time start.

        Inputs:
            :param dt.datetime start: DateTime object at the desired start

        Parameters:
            time:         datetime object representing the current time of
                            the model
            sun_unit:     unit vector from the center of the Earth to the
                            Sun in J2000
            sun_mag:      distance from the center of the Earth to the Sun
            moon_unit:    unit vector from the center of the Earth to the
                            Moon in J2000
            moon_mag:     distance from the cneter of the Earth to the Moon
        """

        # Set universe time
        self.time = start if start is not None else dt.datetime.today()

        # Set sun vector in J2000
        self.sun_unit, self.sun_mag = self.__update_sun_vector()

        # Set moon vector in J2000
        self.moon_unit, self.moon_mag = self.__update_moon_vector()

        # Initialize force models
        self.force_models = [Gravity()]

        # Initialize satellite
        time_init       = self.time
        mass_prop       = {'mass': 1500.0, 'cog': vt.Vector([0.5, 0.5, 1.5]), 'moi': vt.Matrix(1375, 0, 0, 0, 1375, 0, 0, 0, 275)}
        pos_init        = vt.Vector([7141.9897400, 0.000, 0.000])
        vel_init        = vt.Vector([0.0, 270.96981048, 4306.941804035])
        acc_init        = sum(map(lambda m: m.force(pos_init, time_init), self.force_models)) / mass_prop['mass']
        att_init        = vt.Vector([0.0, 0.0, 1.0])
        ang_vel_init    = vt.Vector([0.0, 0.0, 0.0])
        self.satellite = SatModel.SatModel(time_init, mass_prop, pos_init, vel_init, acc_init, att_init, ang_vel_init)

    # ----------------
    #  Public Methods
    # ----------------
    def update(self, delta_time):
        self.time += delta_time
        self.sun_unit, self.sun_mag = self.__update_sun_vector()
        self.moon_unit, self.moon_mag = self.__update_moon_vector()
        self.satellite.update(self.time)

    # -----------------
    #  Private Methods
    # -----------------
    def __update_sun_vector(self):
        """Calculates the direction and distance to the sun from the Earth

        :rtype: vt.Vector, float
        :return: Unit sun vector, sun magnitude
        """

        delta_J2000 = self.time - Universe.__J2000_DATE
        n_days_J2000 = delta_J2000.days + delta_J2000.seconds/86400

        mean_lon_sun = 280.460 + 0.9856474*n_days_J2000
        mean_lon_sun %= 360.0
        mean_lon_sun *= Universe.__DEG_TO_RAD

        mean_anomaly_sun = 357.528 + 0.9856003*n_days_J2000
        mean_anomaly_sun %= 360.0
        mean_anomaly_sun *= Universe.__DEG_TO_RAD

        ecliptic_lon_sun = ( mean_lon_sun/Universe.__DEG_TO_RAD +
                             1.915*sin(mean_anomaly_sun) +
                             0.020*sin(2.0*mean_anomaly_sun) )
        ecliptic_lon_sun *= Universe.__DEG_TO_RAD

        dist_earth_to_sun = (1.00014 -
                             0.01671*cos(mean_anomaly_sun) -
                             0.00014*cos(2.0*mean_anomaly_sun) )
        dist_earth_to_sun *= Universe.__AU_TO_KM

        obliquity_ecliptic = 23.439 - 0.0000004*n_days_J2000
        obliquity_ecliptic *= Universe.__DEG_TO_RAD

        x_J2000_sun = cos(ecliptic_lon_sun)
        y_J2000_sun = cos(obliquity_ecliptic)*sin(ecliptic_lon_sun)
        z_J2000_sun = sin(obliquity_ecliptic)*sin(ecliptic_lon_sun)

        return vt.Vector([x_J2000_sun, y_J2000_sun, z_J2000_sun]), dist_earth_to_sun

    def __update_moon_vector(self):
        """Calculates the direction and distance to the moon from the Earth

        :rtype: vt.Vector, float
        :return: Unit moon vector, sun magnitude
        """

        delta_J2000 = self.time - Universe.__J2000_DATE
        n_days_J2000 = delta_J2000.days + delta_J2000.seconds/86400

        mean_lon_moon = 218.316 + 13.176396*n_days_J2000
        mean_lon_moon %= 360.0
        mean_lon_moon *= Universe.__DEG_TO_RAD

        mean_anomaly_moon = 134.963 + 13.064993*n_days_J2000
        mean_anomaly_moon %= 360.0
        mean_anomaly_moon *= Universe.__DEG_TO_RAD

        mean_dist_moon = 93.272 + 13.229350*n_days_J2000
        mean_dist_moon %= 360.0
        mean_dist_moon *= Universe.__DEG_TO_RAD

        ecliptic_lon_moon = (mean_lon_moon/Universe.__DEG_TO_RAD +
                             6.289*sin(mean_anomaly_moon) )
        ecliptic_lon_moon *= Universe.__DEG_TO_RAD

        ecliptic_lat_moon = 5.128*sin(mean_dist_moon)
        ecliptic_lat_moon *= Universe.__DEG_TO_RAD

        dist_earth_to_moon = (385001.0 -
                              20905.0*cos(mean_anomaly_moon) )

        obliquity_ecliptic = 23.439 - 0.0000004*n_days_J2000
        obliquity_ecliptic *= Universe.__DEG_TO_RAD

        x_ec = cos(ecliptic_lat_moon)*cos(ecliptic_lon_moon)
        y_ec = cos(ecliptic_lat_moon)*sin(ecliptic_lon_moon)
        z_ec = sin(ecliptic_lat_moon)
        rect_ec = np.matrix([[x_ec],[y_ec],[z_ec]])

        Q_eq_ec = np.matrix([[1, 0, 0],
                             [0, cos(obliquity_ecliptic), -sin(obliquity_ecliptic)],
                             [0, sin(obliquity_ecliptic), cos(obliquity_ecliptic)]])

        rect_eq = Q_eq_ec*rect_ec

        return rect_eq, dist_earth_to_moon

    def __update_orbit(self):
        p = self.satellite.p
        v = self.satellite.v
        a = self.satellite.a




class Gravity(afm.AbstractForceModel):

    __MODEL_ORDER_N = 10
    __COEFFICIENTS = dict()
    __C = np.zeros((__MODEL_ORDER_N+1, __MODEL_ORDER_N+1))
    __S = np.zeros((__MODEL_ORDER_N+1, __MODEL_ORDER_N+1))

    def __init__(self):
        return

    def force(self, r_eci, t):
        """

        :param dt.datetime t: Time at which the force takes place
        :param vt.Vector r_eci: Position in ECI
        :rtype: vt.Vector
        :return: Force vector in ECI
        """

        r_norm = r_eci.norm()
        r_norm3 = r_norm**3

        return -Universe.MU/r_norm3*r_eci

        # ecef = conv.eci_to_ecef(t, r_eci)
        # f_ecef = Gravity.calculate_force(ecef[0])
        # return conv.ecef_to_eci(t, f_ecef)[0]

    @classmethod
    def import_coefficients(cls, filename):
        with open(filename, "r") as fcoeff:
            for line in fcoeff:
                # EGM96 coefficients in the form of:
                #    (n,m,Cnm,Snm,sigmaCnm,sigmaSnm) --> FORMAT(2I4,2E20.12,2E16.8)
                line_entries = line.split()

                n = int(line_entries[0])
                m = int(line_entries[1])
                c_nm = float(line_entries[2].replace('D', 'E'))
                s_nm = float(line_entries[3].replace('D', 'E'))

                cls.__COEFFICIENTS[(n, m)] = (c_nm, s_nm)
                if m <= cls.__MODEL_ORDER_N and n<= cls.__MODEL_ORDER_N:
                    cls.__C[m][n] = c_nm
                    cls.__S[m][n] = s_nm

        return

    @classmethod
    def calculate_geopotential(cls, r, phi, theta):
        if len(cls.__COEFFICIENTS) == 0:
            cls.import_coefficients('D:\Users\\tjh97\Dropbox\Workspace\SatelliteSimulator\Model\egm96_coeff.dat')
        R = Universe.RADIUS_OF_EARTH
        P, P_prime = sp.lpmn(cls.__MODEL_ORDER_N, cls.__MODEL_ORDER_N, sin(theta))

        zonal = 0
        # for n in range(2, cls.__MODEL_ORDER_N):
        #     J_tilde = -1/(Universe.MU*R**n)
        #     P = sp.lpmv(0, n, sin(theta))
        #     den = (r/R)**n
        #
        #     zonal += (J_tilde*P)/den

        tess = 0
        for n in range(2, cls.__MODEL_ORDER_N):
            for m in range(1, n):
                C_nm, S_nm = cls.__C[m][n], cls.__S[m][n]

                C_nm_tilde = -C_nm/(Universe.MU*R**n)
                S_nm_tilde = -S_nm/(Universe.MU*R**n)

                P_mn = P[m][n]
                tg = C_nm_tilde*cos(m*phi) + S_nm_tilde*sin(m*phi)
                den = (r/R)**n

                tess += P_mn*tg/den

        return -(Universe.MU/r)*(1 + zonal + tess)

    @classmethod
    def calculate_force(cls, r_ecef):
        """

        :param dt.datetime t: Time at which the force takes place
        :param vt.Vector r_eci: Position in ECI
        :rtype: vt.Vector
        :return: Force vector in ECI
        """

        r, theta, phi = conv.xyz_to_sph(r_ecef)

        f_r = cls._calculate_radial_force(r, phi, theta)
        f_t = cls._calculate_theta_force(r, phi, theta)
        f_p = cls._calculate_phi_force(r, phi, theta)

        f_ecef = conv.sph_to_xyz_dcm(theta, phi) * vt.Vector([f_r, f_t, f_p])
        return f_ecef

    @classmethod
    def _calculate_radial_force(cls, r, phi, theta):
        if len(cls.__COEFFICIENTS) == 0:
            cls.import_coefficients('D:\Users\\tjh97\Dropbox\Workspace\SatelliteSimulator\Model\egm96_coeff.dat')

        R = Universe.RADIUS_OF_EARTH
        P, P_prime = sp.lpmn(cls.__MODEL_ORDER_N, cls.__MODEL_ORDER_N, sin(theta))

        tess = 0
        for n in range(2, cls.__MODEL_ORDER_N):
            for m in range(1, n):
                C_nm, S_nm = cls.__C[m][n], cls.__S[m][n]

                C_nm_tilde = -C_nm/(Universe.MU*R**n)
                S_nm_tilde = -S_nm/(Universe.MU*R**n)

                P_mn = P[m][n]
                tg = C_nm_tilde*cos(m*phi) + S_nm_tilde*sin(m*phi)
                den = (r/R)**n

                tess += (n+1)*P_mn*tg/den
                # tess += (n+1)*P*(C_nm*cos(m*phi) + S_nm*sin(m*phi))/r**(n+2)
        # tess = 0

        return -(Universe.MU/r**2)  # *(1/Universe.MU - tess)

    @classmethod
    def _calculate_theta_force(cls, r, phi, theta):
        P, P_prime = sp.lpmn(cls.__MODEL_ORDER_N, cls.__MODEL_ORDER_N, sin(theta))

        # top_array = P_prime*cos(theta) * \
        #     (np.dot(np.diag([cos(m*phi) for m in range(cls.__MODEL_ORDER_N+1)]), cls.__C) +
        #      np.dot(np.diag([sin(m*phi) for m in range(cls.__MODEL_ORDER_N+1)]), cls.__S))

        tess = 0
        for n in range(2, cls.__MODEL_ORDER_N):
            for m in range(1, n):
                C_nm, S_nm = cls.__C[m][n], cls.__S[m][n]

                top = P_prime[m][n] * cos(theta) * (C_nm*cos(m*phi) + S_nm*sin(m*phi))
                den = r**(n+2)

                tess += top/den

        return tess if not isnan(tess) else 0.0
        # return np.dot(np.ones((1, cls.__MODEL_ORDER_N+1)), top_array) * np.matrix([[r ** n+2 for n in range(cls.__MODEL_ORDER_N+1)]]).T

    @classmethod
    def _calculate_phi_force(cls, r, phi, theta):
        # if sin(theta) < 1e-16:
        #     return 0.0

        P, P_prime = sp.lpmn(cls.__MODEL_ORDER_N, cls.__MODEL_ORDER_N, sin(theta))

        tess = 0
        for n in range(2, cls.__MODEL_ORDER_N):
            for m in range(1, n):
                C_nm, S_nm, = cls.__C[m][n], cls.__S[m][n]

                top = P[m][n] * (m*S_nm*sin(m*phi) - m*C_nm*cos(m*phi))
                den = r**(n+2) * sin(theta)

                tess += top/den

        return tess if not isnan(tess) else 0.0
