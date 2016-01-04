import datetime as dt
import math

# -------------------
#  Public Constants
# -------------------
MU = 398600.4418                # km^3/s^2
UNIVERSE_TIME_STEP = 0.01       # sec
RADIUS_OF_EARTH = 6378.1363     # km

J2000_DATE = dt.datetime(2000, 1, 1, 12, 0, 0)     # J2000 Epoch starts at noon on 1/1/2000
DEG_TO_RAD = math.pi/180.0
RAD_TO_DEG = 180.0/math.pi
AU_TO_KM = 149597871.0
