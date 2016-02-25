from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cbook as mpl_cbook
# import matplotlib._png as mpl_png
import matplotlib.pyplot as plt
import mayavi.mlab as mlab
import basics.constants as constant
import basics.vector_types as vt
import numpy as np


def plot_earth_map():
    # set up orthographic map projection with
    # perspective of satellite looking down at 50N, 100W.
    # use low resolution coastlines.
    map = Basemap(projection='ortho', lat_0=45, lon_0=-100, resolution='l')

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral', lake_color='aqua')

    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')

    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0, 360, 30))
    map.drawparallels(np.arange(-90, 90, 30))

    plt.show()


def add_earth_sphere(center=None):
    """Add a sphere the radius of the Earth to the current figure.

    :param [float] center: Coordinates to locate the center of the Earth in the plot
    :return:
    """
    phi = np.linspace(0, 2*np.pi, 100)
    theta = np.linspace(0, np.pi, 100)

    if center is None:
        center = vt.Vector([0, 0, 0])

    x = constant.RADIUS_OF_EARTH * np.outer(np.cos(phi), np.sin(theta)) + center[0]
    y = constant.RADIUS_OF_EARTH * np.outer(np.sin(phi), np.sin(theta)) + center[1]
    z = constant.RADIUS_OF_EARTH * np.outer(np.ones(np.size(phi)), np.cos(theta)) + center[2]

    # earth_image_file_path = os.path.join(os.environ.get('SATSIMDIR'), 'basics', 'earth_map.png')
    # fn = mpl_cbook.get_sample_data(earth_image_file_path, asfileobj=False)
    # img = mpl_png.read_png(fn).transpose(1, 0, 2)

    # ax.plot_surface(x, y, z, color='r', shade=True)# rstride=1, cstride=1, facecolors=img)
    mlab.mesh(x, y, z, color=(0, 1, 0))
