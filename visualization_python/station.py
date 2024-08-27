from meshloader import load_meshes
from utils import obs_trisurf, set_aspect_equal_3d
import numpy as np
import matplotlib.pyplot as plt

def obs_normals(ax, meshes, show=False):
    print('normal magnitudes:')
    for struct in meshes:
        for i in range(struct.stlmesh.normals.shape[0]):
            normal = struct.stlmesh.normals[i]
            print(np.linalg.norm(normal))
            v0 = struct.stlmesh.v0[i]
            v1 = struct.stlmesh.v1[i]
            v2 = struct.stlmesh.v2[i]
            centroid = (v0 + v1 + v2) / 3
            ax.quiver(centroid[0], centroid[1], centroid[2], normal[0], normal[1], normal[2], color='r')
    return ax

def visualize_station(show=False, show_start=False):
    x,y,z = 1.8, 4.7, 2.7
    start = np.array([2, 3, 4.2])
    meshes = load_meshes()
    ax = obs_trisurf(meshes, show=False)
    # ax = plot_path_direct('station_viewpoints_coverage.csv', ax=ax)
    # ax = obs_normals(ax, meshes, show=False)

    # show axes also
    # ax.plot(x,y,z, 'rx')
    if show_start: ax.plot(start[0], start[1], start[2], 'rx') # start point
    set_aspect_equal_3d(ax)
    if show: plt.show()
    else: return ax

if __name__ == "__main__":
    visualize_station(show=True, show_start=False)