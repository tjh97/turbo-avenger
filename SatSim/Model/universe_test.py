# from Universe import *
import matplotlib
# matplotlib.use('GTKAgg')
from matplotlib import pyplot as plt
import Universe
import SatModel
import datetime as dt
import basics.conversions as conv
import math
import basics.plotting as plotting
from mpl_toolkits.mplot3d import Axes3D
import basics.rkf6 as rkf6
import basics.vector_types as vt
import scipy.integrate as spyi
import numpy as np
import time

# Setup the satellite model
time_init = dt.datetime.now()
mass_prop = {'mass': 1500.0,
             'cog': vt.Vector([0.5, 0.5, 1.5]),
             'moi': vt.Matrix(1375, 0, 0, 0, 1375, 0, 0, 0, 275)}
p_init = vt.Vector([7155.27304548, 0.0, 0.0])
v_init = vt.Vector([0.0, 0.46957258, 7.45043721])
att_init = vt.Vector([0.0, 0.0, 1.0])
ang_vel_init = vt.Vector([0.0, 0.06, 0.0])
satellite = SatModel.SatModel(time_init, mass_prop, p_init, v_init, att_init, ang_vel_init)

# Make universe
universe = Universe.Universe(satellite, time_init)

# Make plotting variables
x_sat = np.array([p_init.x])
y_sat = np.array([p_init.y])
z_sat = np.array([p_init.z])

# Initialize the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plotting.add_earth_sphere(ax)
# ax.plot(x_sat, y_sat, z_sat, label='Orbit')
# ax.legend()

# plt.ion()
# plt.show()
# plt.draw()
# background = fig.canvas.copy_from_bbox(ax.bbox)
# points = ax.plot(x_sat, y_sat, z_sat, label='Orbit')[0]
# plt.draw()
# time.sleep(1)

# Iterate through some time period to test the model
time_step = dt.timedelta(seconds=10.0)
for i in range(600):
    universe.update(time_step)

    x_sat = np.append(x_sat, satellite.p.x)
    y_sat = np.append(y_sat, satellite.p.y)
    z_sat = np.append(z_sat, satellite.p.z)
    # ax.plot(x_sat, y_sat, z_sat, label='Orbit')
    # plt.draw()
    # time.sleep(0.05)

    # points.set_data(satellite.p.x, satellite.p.y)
    # points.set_3d_properties(satellite.p.z)
    # fig.canvas.restore_region(background)
    # ax.draw_artist(points)
    # fig.canvas.blit(ax.bbox)

ax.plot(x_sat, y_sat, z_sat)
plt.show()
