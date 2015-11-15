from SatModel import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

time_init       = dt.datetime.now()
mass_prop       = {'mass': 1500.0, 'cog': vector_types.Vector([0.5, 0.5, 1.5]), 'moi': vector_types.Matrix(1375, 0, 0, 0, 1375, 0, 0, 0, 275)}
pos_init        = vector_types.Vector([7141.9897400, 0.000, 0.000])
vel_init        = vector_types.Vector([0.0, 270.96981048, 4306.941804035])
att_init        = vector_types.Vector([0.0, 0.0, 1.0])
ang_vel_init    = vector_types.Vector([0.0, 0.06, 0.0])

satmodel = SatModel(time_init, mass_prop, pos_init, vel_init, att_init, ang_vel_init)
x_plt, y_plt, z_plt = satmodel.att.get_att_vec()
x_plt = [x_plt]
y_plt = [y_plt]
z_plt = [z_plt]

for i in range(1000):
    satmodel.update(dt.timedelta(0,1))
    x,y,z = satmodel.att.get_att_vec()
    print(str(satmodel))

    x_plt = np.append(x_plt, x)
    y_plt = np.append(y_plt, y)
    z_plt = np.append(z_plt, z)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x_plt,y_plt,z_plt,label='Attitude')
plt.show()

