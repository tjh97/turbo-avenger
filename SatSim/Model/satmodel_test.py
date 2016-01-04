from SatModel import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

time_init       = dt.datetime.now()
mass_prop       = {'mass': 1500.0, 'cog': vector_types.Vector([0.5, 0.5, 1.5]), 'moi': vector_types.Matrix(1375, 0, 0, 0, 1375, 0, 0, 0, 275)}
pos_init        = vector_types.Vector([7155.27304548, 0.0, 0.0])
vel_init        = vector_types.Vector([0.0, 0.46957258, 7.45043721])
att_init        = vector_types.Vector([0.0, 0.0, 1.0])
ang_vel_init    = vector_types.Vector([0.0, 0.06, 0.0])

satmodel = SatModel(time_init, mass_prop, pos_init, vel_init, att_init, ang_vel_init)
x_plt, y_plt, z_plt = satmodel.att.get_att_vec()
x_plt = [x_plt]
y_plt = [y_plt]
z_plt = [z_plt]

x_pos, y_pos, z_pos = pos_init.x, pos_init.y, pos_init.z
x_pos = [x_pos]
y_pos = [y_pos]
z_pos = [z_pos]

for i in range(1000):
    satmodel.update(time_init + (i+1)*dt.timedelta(seconds=10.0))
    x, y, z = satmodel.att.get_att_vec()
    print(str(satmodel))

    x_plt = np.append(x_plt, x)
    y_plt = np.append(y_plt, y)
    z_plt = np.append(z_plt, z)

    x_pos = np.append(x_pos, satmodel.p.x)
    y_pos = np.append(y_pos, satmodel.p.y)
    z_pos = np.append(z_pos, satmodel.p.z)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_plt, y_plt, z_plt, label='Attitude')

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_pos, y_pos, z_pos, label='Orbit')
plt.show()

