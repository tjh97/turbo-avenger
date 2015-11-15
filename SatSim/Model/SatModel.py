'''
Created on Apr 25, 2015

@author: tjh97
'''
import Universe
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from Simulator.Simulator import Simulator
from basics import space_types, vector_types
from math import sqrt, sin, cos, acos, pi
# from jinja2.nodes import Pos


class SatModel(object):
    """


    """


    def __init__(self, time_init, mass_prop, pos_init, vel_init, acc_init, att_init, ang_vel_init):
        """


        :param time_init:
        :param mass_prop:
        :param pos_init:
        :param vel_init:
        :param acc_init:
        :param att_init:
        :param ang_vel_init:
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
        self.a = acc_init
        
    def update(self, time):
        delta_time = time - self.time
        self.att.update(delta_time, self.w)
        self.time = time

    def __str__(self):
        time_str = 'Time: %s' %self.time.ctime()
        att_str  = 'Att: %s' %self.att
        angv_str = 'AngV: %s' %self.w
        
        return '\n'.join([time_str,att_str,angv_str])

        
class Attitude(SatModel):
    '''
    classdocs
    '''
    
    
    def __init__(self, J2000_init):
        self.att_init = J2000_init   # Direction of +z-axis in J2000
        self.q = space_types.Quaternion([1,0,0,0])  # Initialize the quaternion to 1,0,0,0
        
    def update(self, delta_time, w):
        quat = [0, 0, 0, 0]   # Quaternion in form (r, v)
        
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
        

class Forces(SatModel):
    '''
    classdocs
    '''
    
    
    def __init__(self, params):
        pass
    
    
class Orbit(SatModel):
    '''
    classdocs
    '''
    
    
    def __init__(self, pos, vel):
        '''
        Constructor
        '''
        self.r = pos
        self.v = vel
        
        self.a, self.e, self.i, self.raan, self.w, self.nu = self.__rv2oe()
        
    
    def update(self, pos, vel):
        self.r = pos
        self.v = vel
        
        self.a, self.e, self.i, self.raan, self.w, self.nu = self.__rv2oe()
    
    def __rv2oe(self):
        h = np.cross(self.r.T, self.v.T).T
        n = np.cross(np.array([0,0,1]), h.T)
        n_norm = np.linalg.norm(n)
        
        mu = Universe.Universe.MU
        r = self.r
        v = self.v
        v2 = v.T.dot(v)
        r_norm = np.linalg.norm(r)
        
        e = ( (v2 - mu/r_norm )*r - ( r.T.dot(v) )*v ) / mu
        e_norm = np.linalg.norm(e)
        
        E = v2/2 - mu/r_norm
        
        if e_norm != 1:
            a = mu/(2*E)
            p = a*(1-e_norm**2)
        else:
            p = h.T.dot(h)/mu
            a = np.inf
        
        i    = acos( h[2]/np.linalg.norm(h) )
        raan = acos( n[0]/np.linalg.norm(n) )
        aop  = acos( n.T.dot(r)/(n_norm*e_norm) )
        t_an = acos( e.T.dot(r)/(e_norm*r_norm) )
        
        raan = 2*pi - raan if n[1] < 0 else raan
        aop  = 2*pi - aop  if e[2] < 0 else aop
        t_an = 2*pi - t_an if r.T.dot(v) < 0 else t_an
    
        return a, e_norm, i, raan, aop, t_an
    


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
    
    