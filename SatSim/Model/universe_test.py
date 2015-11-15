from Universe import *
import datetime as dt
import basics.conversions as conv
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

t = dt.datetime.now()
x_plt = np.array([])
y_plt = np.array([])
z_plt = np.array([])

x_plt_eci = np.array([])
y_plt_eci = np.array([])
z_plt_eci = np.array([])

time_step = dt.timedelta(0, 16.6667)

for th in range(12*360):
    th = math.radians(th)
    r_eci = conv.sph_to_xyz(Universe.RADIUS_OF_EARTH+780, th, math.radians(0))
    ecef = conv.eci_to_ecef(t, r_eci)
    r_ecef = ecef[0]
    r_ecef, th_ecef, phi_ecef = conv.xyz_to_sph(r_ecef)
    # print('%f, %f, %f' %(r_ecef, math.degrees(th_ecef), math.degrees(phi_ecef)))
    t += time_step

    x_plt = np.append(x_plt, ecef[0].x())
    y_plt = np.append(y_plt, ecef[0].y())
    z_plt = np.append(z_plt, ecef[0].z())

    x_plt_eci = np.append(x_plt_eci, r_eci.x())
    y_plt_eci = np.append(y_plt_eci, r_eci.y())
    z_plt_eci = np.append(z_plt_eci, r_eci.z())

    # u = Gravity.calculate_phi_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # v = Gravity.calculate_theta_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # w = Gravity.calculate_radial_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # x = Gravity.calculate_geopotential(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # print(str(w) + ' ' + str(v) + '  ' + str(u) + ' ' + str(x))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x_plt,y_plt,z_plt,label='ECEF')
ax.plot(x_plt_eci, y_plt_eci, z_plt_eci, label='ECI')
ax.legend()
plt.show()
