import math
import datetime as dt


class RKF6(object):
    """A Runge-Kutta-Fehlberg ODE numerical integration
    method. Used for a second order, 3-D ODE, such as for
    orbit propagation.

    """

    import vector_types as vt

    def __init__(self, f, h, eps):
        """

        :param function f:
        :param dt.timedelta h:
        :return:
        """

        self.h = h
        self.f = f
        self.eps = eps

    def step(self, r_init, v_init, t):
        """

        :param vector_types.Vector r_init:
        :param vector_types.Vector v_init:
        :param dt.datetime t:
        :rtype: vector_types.Vector, vector_types.Vector
        :return: r_new, v_new
        """

        h = self.h
        h_sec = h.total_seconds()
        f = self.f

        k11 = h_sec*v_init
        k12 = h_sec*f(r_init, t)

        k21 = h_sec*(v_init + k12/4.0)
        k22 = h_sec*f(r_init + k11/4.0, t + h/4)

        k31 = h_sec*(v_init + 3.0/32*k12 + 9.0/32*k22)
        k32 = h_sec*f(r_init + 3.0/32*k11 + 9.0/32*k21, t + h*3/8)

        k41 = h_sec*(v_init + 1932.0/2197*k12 - 7200.0/2197*k22 + 7296.0/2197*k32)
        k42 = h_sec*f(r_init + 1932.0/2197*k11 - 7200.0/2197*k21 + 7296.0/2197*k31, t + h*12/13)

        k51 = h_sec*(v_init + 439.0/216*k12 - 8.0*k22 + 3680.0/513*k32 - 845.0/4104*k42)
        k52 = h_sec*f(r_init + 439.0/216*k11 - 8.0*k21 + 3680.0/513*k31 - 845.0/4104*k41, t + h*23/24)

        k61 = h_sec*(v_init - 8.0/27*k12 + 2.0*k22 - 3544.0/2565*k32 + 1859.0/4104*k42 - 11.0/40*k52)
        k62 = h_sec*f(r_init - 8.0/27*k11 + 2.0*k21 - 3544.0/2565*k31 + 1859.0/4104*k41 - 11.0/40*k51, t + h/2)

        r_new = r_init + 25.0/216*k11 + 1408.0/2565*k31 + 2197.0/4104*k41 - k51/5.0
        v_new = v_init + 25.0/216*k12 + 1408.0/2565*k32 + 2197.0/4104/k42 - k52/5.0

        r_tilde = r_init + 16.0/135*k11 + 6656.0/12825*k31 + 28561.0/56430*k41 - 9.0/50*k51 + 2.0/55*k61
        v_tilde = v_init + 16.0/135*k12 + 6656.0/12825*k32 + 28561.0/56430*k42 - 9.0/50*k52 + 2.0/55*k62

        r_state = (r_tilde-r_new) / h_sec
        v_state = (v_tilde-v_new) / h_sec

        norm = math.sqrt(r_state.norm2() + v_state.norm2())
        delta = 0.84 * (self.eps/norm)**0.25

        # Adapt step size
        # if delta <= 0.1:
        #     self.h /= 10
        # elif delta >= 4.0:
        #     self.h *= 4
        # elif 1.0 < delta < 4.0:
        #     self.h = dt.timedelta(seconds=(h_sec*delta))

        # Accept or reject the new step value
        # if norm <= self.eps:
        return r_new, v_new
        # else:
        #     return self.step(r_init, v_init, t)
