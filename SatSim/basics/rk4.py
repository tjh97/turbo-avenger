def rk4(y, dydt, t, h, f):
    """

    :param y:
    :param dydt:
    :param t:
    :param h:
    :param f:
    :return:
    """

    hh = h*0.5
    h6 = h/6.0
    xh = t+hh

    yt = y + hh*dydt
    dyt = f(xh, yt)

    yt = y + hh*dyt
    dym = f(xh, yt)

    yt = y + h*dym
    dym += dyt
    dyt = f(t+h, yt)

    return y + h6*(dydt + dyt + 2.0*dym)


def rk4_sample():
    import numpy as np
    import matplotlib.pyplot as plt
    import math

    x = 0
    y = 1
    deriv = lambda x, y: 5*y-3
    dxy = deriv(x, y)
    h = 0.0001

    x_plt = np.array([x])
    y_plt = np.array([y])

    for ii in range(10001):
        y = rk4(y, dxy, x, h, deriv)
        x += h
        x_plt = np.append(x_plt, x)
        y_plt = np.append(y_plt, y)
        dxy = deriv(x, y)

    x_tru = np.array(range(10001), dtype=np.float64)*h
    y_tru = map(lambda _: (2/5.0)*math.exp(5*_) + 3/5.0, x_tru)

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(x_plt, y_plt, label="RK4")
    ax.plot(x_tru, y_tru, label="Truth")
    plt.show()

if __name__ == "__main__":
    rk4_sample()