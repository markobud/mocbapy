import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull


def draw3d(polygon):
    """ Returns a 3D figure of the Pareto Front. Input: A mo-fba Polygon (eg. sol_mofba.Primal) """
    points = polygon.vertex_value[[x == 1 for x in polygon.vertex_type]]
    hull = ConvexHull(points)
    pd = Poly3DCollection([hull.points[s] for s in hull.simplices])
    fig = plt.figure(figsize=(9, 10))
    ax = fig.add_subplot(111, projection='3d')
    pd.set_facecolor('yellow')
    pd.set_alpha(0.4)
    pd.set_edgecolor('black')

    ax.add_collection3d(pd)
    return fig, ax


def draw2d(polygon):
    """"""


def drawNd(polygon):
    """"""