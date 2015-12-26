'''
Created on Aug 17, 2015

@author: tjh97
'''


import numpy as np
from basics import vector_types as vt
import math


class Quaternion(vt.Vector):
    """Class representing quaternions for rotations.

    The convention of this quaternion is the scalar component first.
    Q = [w, x, y, z]
    """


    def __init__(self, *args):
        """

        :param args:
        :return:
        """

        if len(args) == 1:
            self.value = np.matrix(args[0])
        elif len(args) == 4:
            self.value = np.matrix(args)
        else:
            pass
        
        if self.value.shape[0] == 1:
            self.value = self.value.T
        
#     def norm(self):
#         return math.sqrt(float(self.q.T * self.q))

    @property
    def w(self):
        return self[0]

    @property
    def x(self):
        return self[1]

    @property
    def y(self):
        return self[2]

    @property
    def z(self):
        return self[3]

    def isunit(self):
        return float(self.value.T * self.value) == 1.0
    
    def unit(self):
        return self.value / self.norm()
    
    def mkunit(self):
        if not self.isunit():
            self.value = self.value / self.norm()
    
    def dcm(self):
        self.mkunit()
        
        w = self[0]
        x = self[1]
        y = self[2]
        z = self[3]
        
        r1 = [1 - 2*(y**2 + z**2), 2*(x*y - w*z), 2*(x*z + w*y)]
        r2 = [2*(x*y + w*z), 1 - 2*(x**2 + z**2), 2*(y*z + w*x)]
        r3 = [2*(x*z - w*y), 2*(y*z - w*x), 1 - 2*(x**2 + y**2)]
        
        return vt.Matrix([r1, r2, r3])
    
    def __add__(self, other):
        return Quaternion(self.q + other.q)
    
    def __mul__(self, other):
        self.mkunit()
        
        w_s = self[0]
        x_s = self[1]
        y_s = self[2]
        z_s = self[3]
        
        w_o = other[0]
        x_o = other[1]
        y_o = other[2]
        z_o = other[3]
        
        return Quaternion(w_s*w_o - x_s*x_o - y_s*y_o - z_s*z_o,
                          w_s*x_o + x_s*w_o + y_s*z_o - z_s*y_o,
                          w_s*y_o - x_s*z_o + y_s*w_o + z_s*x_o,
                          w_s*z_o + x_s*y_o - y_s*x_o + z_s*w_o)
        
    def __imul__(self, other):
        self.mkunit()
         
        w_s = other[0]
        x_s = other[1]
        y_s = other[2]
        z_s = other[3]
         
        w_o = self[0]
        x_o = self[1]
        y_o = self[2]
        z_o = self[3]
        
        self.value = np.matrix([w_s*w_o - x_s*x_o - y_s*y_o - z_s*z_o,
                            w_s*x_o + x_s*w_o + y_s*z_o - z_s*y_o,
                            w_s*y_o - x_s*z_o + y_s*w_o + z_s*x_o,
                            w_s*z_o + x_s*y_o - y_s*x_o + z_s*w_o]).T
         
        return self
        
        
    @staticmethod
    def rot2quat(rot):
        """

        :param vt.Matrix rot:
        :rtype: Quaternion
        :return:
        """

        r00 = rot[0,0]
        r01 = rot[0,1]
        r02 = rot[0,2]
        r10 = rot[1,0]
        r11 = rot[1,1]
        r12 = rot[1,2]
        r20 = rot[2,0]
        r21 = rot[2,1]
        r22 = rot[2,2]

        w = math.sqrt(max(0, (r00 + r11 + r22 + 1.0)/4.0))
        x = math.sqrt(max(0, (r00 - r11 - r22 + 1.0)/4.0))
        y = math.sqrt(max(0, (-r00 + r11 - r22 + 1.0)/4.0))
        z = math.sqrt(max(0, (-r00 - r11 + r22 + 1.0)/4.0))
        
        max_el = max(w, x, y, z)

        if (max_el == w):
            w *= 1.0
            x = math.copysign(x, r21 - r12)
            y = math.copysign(y, r02 - r20)
            z = math.copysign(z, r10 - r01)
        elif (max_el == x):
            w = math.copysign(w, r21 - r12)
            x *= 1.0
            y = math.copysign(y, r10 + r01)
            z = math.copysign(z, r02 + r20)
        elif (max_el == y):
            w = math.copysign(w, r02 - r20)
            x = math.copysign(x, r10 + r01)
            y *= 1.0
            z = math.copysign(z, r21 + r12)
        elif (max_el == z):
            w = math.copysign(w, r10 - r01)
            x = math.copysign(x, r20 + r02)
            y = math.copysign(y, r21 + r12)
            z *= 1.0
        else:
            print("Error!")

        q = Quaternion([w, x, y, z])
        q.mkunit()

        return q

    @staticmethod
    def makeQ(ax, ang):
        """

        :param vt.Vector axis: Axis of the rotation.
        :param float angle: Angle of rotation in radians.
        :rtype: Quaternion
        :return: quaternion representing the rotation
        """

        # Get copy of axis
        axis = vt.Vector(ax.x, ax.y, ax.z)

        # Make unit vector
        axis /= axis.norm()

        sangle = math.sin(ang/2.0)
        cangle = math.cos(ang/2.0)

        axis *= sangle

        return Quaternion([cangle, axis.x, axis.y, axis.z])