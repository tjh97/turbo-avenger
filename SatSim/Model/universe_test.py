from Universe import *
import datetime as dt
import basics.conversions as conv
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import basics.rkf6 as rkf6
import basics.vector_types as vt
import scipy.integrate as spyi

r = vt.Vector([7155.27304548, 0.0, 0.0])
v = vt.Vector([0.0, 0.46957258, 7.45043721])

t = dt.datetime.now()
x_plt = np.array([r.x])
y_plt = np.array([r.y])
z_plt = np.array([r.z])

x_plt_eci = np.array([])
y_plt_eci = np.array([])
z_plt_eci = np.array([])

time_step = dt.timedelta(0, 10.0)
time_now = dt.datetime.now()

force_model = Gravity()
orbit_rk = rkf6.RKF6(force_model.force, time_step, 0.02)

# f = lambda time, x: [x[3], x[4], x[5]].extend(list(force_model.force(vt.Vector(x[0:2]), time_now + dt.timedelta(seconds=time))))
# ode853 = spyi.ode(f)
# ode853.set_integrator("dop853", first_step=16.6667, atol=0.1, dfactor=4, safety=0.84, max_step=1000)
# ode853.set_initial_value([7141.9897400, 0.001, 0.001, 0.000001, 0.27096981048, 7.306941804035])

for th in range(360):
    r, v, t = orbit_rk.step(r, v, t)

    # r = ode853.integrate(th*10.0)

    x_plt = np.append(x_plt, r[0])
    y_plt = np.append(y_plt, r[1])
    z_plt = np.append(z_plt, r[2])

    a, e, i, om, w, nu = conv.rv_to_oe(r, v, Universe.MU)
    print("a: %f, e: %f, i: %f, Om: %f, w: %f, nu: %f" %(a, e, i, om, w, nu))

    # th = math.radians(th)
    # r_eci = conv.sph_to_xyz(Universe.RADIUS_OF_EARTH+780, th, math.radians(0))
    # ecef = conv.eci_to_ecef(t, r_eci)
    # r_ecef = ecef[0]
    # r_ecef, th_ecef, phi_ecef = conv.xyz_to_sph(r_ecef)
    # # print('%f, %f, %f' %(r_ecef, math.degrees(th_ecef), math.degrees(phi_ecef)))
    # t += time_step
    #
    # x_plt = np.append(x_plt, ecef[0].x())
    # y_plt = np.append(y_plt, ecef[0].y())
    # z_plt = np.append(z_plt, ecef[0].z())

    # x_plt_eci = np.append(x_plt_eci, r_eci.x())
    # y_plt_eci = np.append(y_plt_eci, r_eci.y())
    # z_plt_eci = np.append(z_plt_eci, r_eci.z())

    # u = Gravity.calculate_phi_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # v = Gravity.calculate_theta_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # w = Gravity.calculate_radial_force(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # x = Gravity.calculate_geopotential(Universe.RADIUS_OF_EARTH+780, 0, th/57.3)
    # print(str(w) + ' ' + str(v) + '  ' + str(u) + ' ' + str(x))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_plt, y_plt, z_plt, label='ECEF')
# ax.plot(x_plt_eci, y_plt_eci, z_plt_eci, label='ECI')
ax.legend()
plt.show()
