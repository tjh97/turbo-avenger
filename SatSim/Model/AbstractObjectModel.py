"""

"""

import abc


class AbstractObjectModel(object):
    """

    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, time, *args):
        self.time = time

    @abc.abstractmethod
    def update(self, time):
        pass
