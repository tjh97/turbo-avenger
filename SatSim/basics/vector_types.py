'''
Created on Aug 16, 2015

@author: tjh97
'''

import numpy as np
import math
from __builtin__ import list


class Matrix(object):
    """
    classdocs
    """

    def __init__(self, *args):
        """
        Constructor
        """
        
        if len(args) == 1:
            self.value = np.matrix(args[0])
        elif len(args) == 3:
            if all(isinstance(x, (list, tuple)) for x in args):
                self.value = np.matrix(args)
            elif all(isinstance(x, (float, int)) for x in args):
                self.value = np.matrix(args)
            else:
                raise TypeError('Arguement types are not consistent')
        elif len(args) == 9:
            if all(isinstance(x, (float, int)) for x in args):
                self.value = np.matrix([args[0:3], args[3:6], args[6:]])
            else:
                raise TypeError('Arguements are not all number types')
        else:
            raise IndexError('Number of arguements must be equal to 1, 3, or 9.')

        # Save any vector as a column
        if self.value.shape[0] == 1:
            self.value = self.value.T
        
    def __str__(self):
        return str(self.value)
    
    def __repr__(self):
        return 'Matrix(' + str(self) + ')'
    
    def __add__(self, other):
        if isinstance(self, Vector) and isinstance(other, Vector):
            return Vector(self.value + other.value)
        return Matrix(self.value + other.value)
    
    def __sub__(self, other):
        if isinstance(self, Vector) and isinstance(other, Vector):
            return Vector(self.value - other.value)
        return Matrix(self.value - other.value)
    
    def __mul__(self, other):
        result = self.value * other.value
        return Vector(result) if any(filter(lambda n: n == 1, result.shape)) else Matrix(result)
    
    def __getitem__(self, k):
        return self.value.item(k)

    def T(self):
        """Returns the transpose of the matrix
        """
        return Matrix(self.value.T)

    
class Vector(Matrix):
    """
    """
    
    def __init__(self, *args):
        # super(Vector, self).__init__(args)
        if len(args) == 1:
            self.value = np.matrix(args[0])
        elif len(args) == 3:
            self.value = np.matrix(args)
        else:
            raise IndexError('Number of arguements must be equal to 1, 3, or 9.')

        # Save any vector as a column
        if self.value.shape[0] == 1:
            self.value = self.value.T
    
    def norm(self):
        return math.sqrt(float(self.value.T * self.value))
    
    def __repr__(self):
        return 'Vector(' + str(self) + ')'

    def x(self):
        return self[0]

    def y(self):
        return self[1]

    def z(self):
        return self[2]