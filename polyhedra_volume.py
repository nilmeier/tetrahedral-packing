# from http://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy
# from pyhull.convex_hull import ConvexHull 


from scipy.spatial import ConvexHull
import numpy as np

def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

def convex_hull_volume_bis(pts):
    ch = ConvexHull(pts)

    simplices = np.column_stack((np.repeat(ch.vertices[0], ch.nsimplex),
                                 ch.simplices))
    tets = ch.points[simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))

#pts = np.random.rand(1000, 3)
