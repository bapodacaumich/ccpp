from meshloader import load_meshes
from utils import obs_trisurf, set_aspect_equal_3d
import os
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

def visualize_station(show=False, coverage_file='coverage_2m_coverage.csv', show_start=False):
    x,y,z = 1.8, 4.7, 2.7
    start = np.array([2, 3, 4.2])
    meshes = load_meshes()
    coverage = None
    if coverage_file is not None:
        # coverage = np.loadtxt(os.path.join(os.getcwd(), '..', 'data', 'coverage_viewpoint_sets', coverage_file)).flatten()
        coverage = np.loadtxt(os.path.join(os.getcwd(), 'coverage', coverage_file), delimiter='').flatten()
        print('Coverage: ', np.sum(coverage)/len(coverage) * 100, '%')
    ax = obs_trisurf(meshes, coverage, show=False, alph=0.5)
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