'''
Created on Aug 17, 2015

@author: tjh97
'''


import numpy as np
from basics import vector_types


class Quaternion(vector_types.Vector):
    '''
    classdocs
    '''


    def __init__(self, *args):
        '''
        Constructor
        '''
        
        
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
        
    def isunit(self):
        return float(self.value.T * self.value) == 1
    
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
        
        return vector_types.Matrix([r1, r2, r3])
    
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
        
        
        
        
        
        
        