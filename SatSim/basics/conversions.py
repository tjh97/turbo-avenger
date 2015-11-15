'''
Created on Oct 26, 2015

@author: tjh97
'''

import datetime as dt
import math
import basics.vector_types as vt


def sph_to_xyz(r, theta, phi):
    """ Converts from spherical to cartesian coordinates.

    :param float r: the radial component
    :param float theta: the polar angle component, measured from the +Z axis
    :param float phi: the azimuthal angle component, measured form the +X axis
    :rtype: vt.Vector
    :return: xyz vector
    """

    return vt.Vector(r*math.cos(phi)*math.sin(theta),
                     r*math.sin(phi)*math.sin(theta),
                     r*math.cos(theta))


def xyz_to_sph(v):
    """ Converts from cartesian to spherical coordinates.

    :param vt.Vector v: xyz vector input
    :rtype: (float, float, float)
    :return: (r, theta, phi) where r is the radial distance, theta is the polar
             angle, and phi is the azimuthal angle
    """

    r = v.norm()
    theta = math.acos(v.z()/r)
    phi = math.atan2(v.y(), v.x())

    return r, theta, phi


def time_to_mst(time):
    """Convert a time to Greenwhich Mean Sidereal Time.

    :param dt.datetime time: time to be converted
    :rtype: float
    :return: Greenwhich Mean Sidereal Time, in radians
    """

    time_midnight = dt.datetime(time.year, time.month, time.day)
    hrs_after_midnight = (time - time_midnight).total_seconds() / 3600.0

    # Calculate days since J2000 epoch
    dt_j2000 = (time - dt.datetime(2000, 1, 1, 12)).total_seconds() / 86400.0
    dt_m_j2000 = (time_midnight - dt.datetime(2000, 1, 1, 12)).total_seconds() / 86400.0

    # Centuries since J2000
    cent_j2000 = dt_j2000 / 36525.0

    # Calculate GMST in hours and convert to degrees
    gmst_hr = (6.697374558 + 0.06570982441908*dt_m_j2000 +
               1.00273790935*hrs_after_midnight + 0.000026*(cent_j2000**2)) % 24.0
    return math.radians(gmst_hr*15.0)


def time_to_gast(time):
    """Convert a specified time to the Greenwhich Apparent Sidereal Time.

    :param dt.datetime time: time to be converted
    :rtype: float
    :return: Greenwhich Apparent Sidereal Time, in radians
    """

    SECONDS_IN_DAY = 86400

    # Get the mean siderial time in degrees
    th_m = time_to_mst(time)

    # Centuries since J2000
    d_j2000 = time - dt.datetime(2000, 1, 1, 12)
    cent_j2000 = d_j2000.total_seconds()/86400.0/36525.0

    # Mean obliquity of the ecliptic (eps_m)
    # see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 3
    eps_m = 84381448.0 - \
            46.8150*cent_j2000 - \
            .00059*(cent_j2000**2) + \
            .001813*(cent_j2000**3)
    eps_m = math.radians(eps_m)

    # Nutations in obliquity and longitude (degrees)
    # see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 4
    l = math.radians(280.4665 + 36000.7698*cent_j2000)
    dl = math.radians(218.3165 + 481267.8813*cent_j2000)
    omega = math.radians(125.04452 - 1934.136261*cent_j2000)

    # Calculate nutations
    # see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 5
    d_psi = -17.20*math.sin(omega) - 1.32*math.sin(2*l) - \
             0.23*math.sin(2*dl) + 0.21*math.sin(2*omega)
    d_eps = 9.20*math.cos(omega) + 0.57*math.cos(2*l) + \
            0.10*math.cos(2*dl) - 0.09*math.cos(2*omega)

    # Convert from arc-seconds to degrees to radians
    d_psi *= 1/3600.0
    d_psi = math.radians(d_psi)
    d_eps *= 1/3600.0
    d_eps = math.radians(d_eps)

    # GAST in radians
    # see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 1
    return (th_m + d_psi*math.cos(eps_m+d_eps)) % 2*math.pi


def ecef_to_eci(time, r_ecef, v_ecef=None, a_ecef=None):
    """Converts from the WGS 84 ECEF coordinate frame to the J2000 ECI
    coordinate frame.

    :param dt.datetime time:
    :param vt.Vector r_ecef:
    :param vt.Vector v_ecef:
    :param vt.Vector a_ecef:
    :rtype: (vt.Vector, vt.Vector|None, vt.Vector|None)
    :return: (r_eci, v_eci, a_eci) if the corresponding ECEF vector was given as input
    """

    # Average inertial rotation rate of the Earth in rad/sec
    OMEGA_E = 7.29211585275553e-005

    # Calculate the Greenwich Apparent Sidereal Time
    th = time_to_gast(time)

    # Create the transformation matrices
    t = vt.Matrix([[math.cos(th), math.sin(th), 0],
                   [-math.sin(th), math.cos(th), 0],
                   [0, 0, 1]])
    t_dot = vt.Matrix([[-OMEGA_E*math.sin(th), OMEGA_E*math.cos(th), 0],
                       [-OMEGA_E*math.cos(th), -OMEGA_E*math.sin(th), 0],
                       [0, 0, 0]])
    t_ddot = vt.Matrix([[-(OMEGA_E**2)*math.cos(th), -(OMEGA_E**2)*math.sin(th), 0],
                        [(OMEGA_E**2)*math.sin(th), -(OMEGA_E**2)*math.cos(th), 0],
                        [0, 0, 0]])

    # Calculate the ECI vectors
    r_eci = t.T()*r_ecef
    v_eci = t.T()*v_ecef + t_dot.T()*r_ecef if v_ecef else None
    a_eci = t.T()*a_ecef + 2*t_dot.T()*v_ecef + t_ddot.T()*r_ecef if a_ecef else None

    return r_eci, v_eci, a_eci


def eci_to_ecef(time, r_eci, v_eci=None, a_eci=None):
    """Converts from the J2000 ECI coordinate frame to the WGS 84 ECEF
    coordinate frame.

    :param dt.datetime time:
    :param vt.Vector r_eci:
    :param vt.Vector v_eci:
    :param vt.Vector a_eci:
    :rtype: (vt.Vector, vt.Vector|None, vt.Vector|None)
    :return:
    """

    # Average inertial rotation rate of the Earth in rad/sec
    OMEGA_E = 7.29211585275553e-005

    # Calculate the Greenwich Apparent Sidereal Time
    th = time_to_gast(time)

    # Create the transformation matrices
    t = vt.Matrix([[math.cos(th), math.sin(th), 0],
                   [-math.sin(th), math.cos(th), 0],
                   [0, 0, 1]])
    t_dot = vt.Matrix([[-OMEGA_E*math.sin(th), OMEGA_E*math.cos(th), 0],
                       [-OMEGA_E*math.cos(th), -OMEGA_E*math.sin(th), 0],
                       [0, 0, 0]])
    t_ddot = vt.Matrix([[-(OMEGA_E**2)*math.cos(th), -(OMEGA_E**2)*math.sin(th), 0],
                        [(OMEGA_E**2)*math.sin(th), -(OMEGA_E**2)*math.cos(th), 0],
                        [0, 0, 0]])

    # Calculate the ECI vectors
    r_ecef = t*r_eci
    v_ecef = t*v_eci + t_dot*r_eci if v_eci else None
    a_ecef = t*a_eci + 2*t_dot*v_eci + t_ddot*r_eci if a_eci else None

    return r_ecef, v_ecef, a_ecef
