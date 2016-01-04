"""
Created on Apr 25, 2015

@author: tjh97
"""
import AbstractObjectModel
import Universe
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et
from mpl_toolkits.mplot3d import Axes3D

from basics import space_types, vector_types, dop853
from math import sqrt, sin, cos


class SatModel(AbstractObjectModel.AbstractObjectModel):
    """Models the satellite.


    """

    force_model = Universe.Gravity()

    def __init__(self, time_init, mass_prop, pos_init, vel_init, att_init, ang_vel_init):
        """


        :param dt. datetime time_init:
        :param mass_prop:
        :param vt.Vector pos_init:
        :param vt.Vector vel_init:
        :param vt.Vector att_init:
        :param vt.Vector ang_vel_init:
        :return:
        """
        self.time = time_init
        self.mass = mass_prop['mass']
        self.cog = mass_prop['cog']
        self.moi = mass_prop['moi']
        self.att = Attitude(att_init)
        self.w = ang_vel_init
        self.p = pos_init
        self.v = vel_init

        self._ode_orbit = dop853.Dop853(SatModel.__orbit_force_model, [self.p.x, self.p.y, self.p.z,
                                                                       self.v.x, self.v.y, self.v.z],
                                        self.time, params=time_init)
        
    def update(self, time):
        delta_time = time - self.time
        self.att.update(delta_time, self.w)
        self.__update_orbit(time)
        self.time = time

    def __str__(self):
        time_str = 'Time: %s' % self.time.ctime()
        att_str  = 'Att: %s' % self.att
        angv_str = 'AngV: %s' % self.w
        orbit_str = 'Orbit: %s, %s' % (self.p, self.v)
        
        return ' '.join([time_str, att_str, angv_str, orbit_str])

    def __update_orbit(self, time):
        self.p, self.v = self._ode_orbit.step(time)

    @staticmethod
    def __orbit_force_model(t, y, time_init):
        x, y, z, v_x, v_y, v_z = y
        x_dot, y_dot, z_dot = v_x, v_y, v_z
        v_dot = SatModel.force_model.force(vector_types.Vector(x, y, z), time_init + dt.timedelta(seconds=t))

        return [x_dot, y_dot, z_dot, v_dot.x, v_dot.y, v_dot.z]

    @classmethod
    def import_sat(cls, filename):
        single_items = ["name", "shape", "dry_mass", "fuel_mass", "drag_coeff", "rad_coeff"]
        vector_items = ["length", "cog", "cop", "drag_area", "solar_area", "rad_area"]
        matrix_items = ["moi"]
        params = {}

        try:
            tree = et.parse(filename)
        except IOError:
            print("Unable to open satellite file.")
            return
        except et.ParseError:
            print("Invalid XML file.")
            return
        root = tree.getroot()

        # Single items are defined in text
        for item in single_items:
            value = root.find(item).text
            try:
                params[item] = float(value)
            except ValueError:
                params[item] = value

        # Vector items are defined in attributes x, y, and z
        attributes = ['x', 'y', 'z']
        for item in vector_items:
            node = root.find(item)
            value = [float(node.attrib.get(a, 0)) for a in attributes]
            params[item] = vector_types.Vector(value)

        # Matrix items are defined in attributes xx, xy, xz, yx, yy, yz, zx, zy, zz
        attributes = ['xx', 'xy', 'xz',
                      'yx', 'yy', 'yz',
                      'zx', 'zy', 'zz']
        for item in matrix_items:
            node = root.find(item)
            value = [float(node.attrib.get(a, 0)) for a in attributes]
            params[item] = vector_types.Matrix(value)

        return SatModel(params)

        
class Attitude(object):
    """
    classdocs
    """

    def __init__(self, J2000_init):
        self.att_init = J2000_init   # Direction of +z-axis in J2000
        self.q = space_types.Quaternion([1, 0, 0, 0])  # Initialize the quaternion to 1,0,0,0
        
    def update(self, delta_time, w):
        quat = [1, 0, 0, 0]   # Quaternion in form (r, v)
        
        x = w[0] * dt.timedelta.total_seconds(delta_time)
        y = w[1] * dt.timedelta.total_seconds(delta_time)
        z = w[2] * dt.timedelta.total_seconds(delta_time)
        angle = sqrt(x**2 + y**2 + z**2)
        
        quat[0] = cos(angle/2.0)
        quat[1] = x/angle * sin(angle/2.0)
        quat[2] = y/angle * sin(angle/2.0)
        quat[3] = z/angle * sin(angle/2.0)
        
        self.q *= space_types.Quaternion(quat)
        
    def get_att_vec(self):
        return self.q.dcm() * self.att_init
    
    def __str__(self):
        return str(self.get_att_vec())


if __name__ == '__main__':
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

    x_pos, y_pos, z_pos = pos_init.x, pos_init.y, pos_init.z
    x_pos = list(x_pos)
    y_pos = list(y_pos)
    z_pos = list(z_pos)
    
    for i in range(1000):
        satmodel.update(dt.timedelta(0,1))
        x,y,z = satmodel.att.get_att_vec()
        print(str(satmodel))
        
        x_plt = np.append(x_plt, x)
        y_plt = np.append(y_plt, y)
        z_plt = np.append(z_plt, z)

        x_pos = np.append(x_pos, satmodel.p.x)
        y_pos = np.append(y_pos, satmodel.p.y)
        z_pos = np.append(z_pos, satmodel.p.z)
        
    plt.figure(1)
    plt.subplot(111)
    plt.gca(projection='3d')
    plt.plot(x_plt, y_plt, z_plt, label='Attitude')
    
    plt.figure(2)
    plt.subplot(111)
    plt.gca(projections='3d')
    plt.plot(x_pos, y_pos, z_pos, label='Orbit')
    plt.show()