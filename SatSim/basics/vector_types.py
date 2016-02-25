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
        if isinstance(other, (float, int)):
            result = self.value * other
        # elif isinstance(self, (float, int)):
        #     result = self * other.value
        else:
            result = self.value * other.value
        result = Vector(result) if any(filter(lambda n: n == 1, result.shape)) else Matrix(result)
        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        if isinstance(other, (int, float)):
            result = self.value / other
        else:
            raise NotImplementedError
        result = Vector(result) if any(filter(lambda n: n == 1, result.shape)) else Matrix(result)
        return result

    def __rdiv__(self, other):
        return self.__div__(other)

    def __getitem__(self, k):
        return self.value.item(k)

    def __eq__(self, other):
        return (self.value == other.value).all()

    def T(self):
        """Returns the transpose of the matrix
        """
        return Matrix(self.value.T)

    def trace(self):
        return self[0, 0] + self[1, 1] + self[2, 2]

    
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

    def norm2(self):
        return float(self.value.T * self.value)

    def cross(self, other):
        """Return self x other

        :param Vector other:
        :rtype: Vector
        :return: self x other
        """

        return Vector([self.y*other.z - self.z*other.y,
                       self.z*other.x - self.x*other.z,
                       self.x*other.y - self.y-other.x])

    def __mul__(self, other):
        if isinstance(other, Vector):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return super(Vector, self).__mul__(other)
    
    def __repr__(self):
        return 'Vector(' + str(self) + ')'

    @property
    def x(self):
        return self[0]

    @property
    def y(self):
        return self[1]

    @property
    def z(self):
        return self[2]
