
import abc
import datetime as dt
import basics.vector_types as vt


class AbstractForceModel(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def force(self, r_eci, t):
        """

        :param dt.datetime t: Time at which the force takes place
        :param vt.Vector r_eci: Position in ECI
        :rtype: vt.Vector
        :return: Force vector in ECI
        """

        return