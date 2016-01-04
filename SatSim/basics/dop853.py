import scipy.integrate as sp
import datetime as dt
import vector_types as vt


class Dop853(object):
    """
    :type f: (float, list[float]) -> list[float]
    :type t0: dt.datetime
    :type t: dt.datetime
    :type y: list[float]
    """

    def __init__(self, f, y0, t0, params=None):
        """

        :param (float, list[float]) -> list[float] f:
        :param  list[float] y0:
        :param float t0:
        :return:
        """

        self.f = f
        self.t0 = t0
        self.t = t0
        self.y = y0

        self._ode = sp.ode(f)
        self._ode.set_integrator("dop853")
        self._ode.set_initial_value(y0)

        if params is not None:
            self._ode.set_f_params(params)

    def step(self, t):
        """

        :param datetime.datetime t:
        :rtype: (vt.Vector, vt.Vector)
        :return:
        """

        delta_t = t - self.t
        self.t = t

        y = self._ode.integrate(self._ode.t + delta_t.total_seconds())
        self.y = y

        r = vt.Vector(y[0:3])
        v = vt.Vector(y[3:])

        return r, v
