from meshloader import load_meshes, load_vx_station, load_original_stl
from utils import obs_trisurf, set_aspect_equal_3d
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def obs_normals(ax, meshes, normals, show=False):
    for i in range(len(meshes)):
        mesh = meshes[i]
        for j in range(len(mesh)):
            face = np.array(mesh[j])
            normal = np.array(normals[i])
            centroid = np.mean(face, axis=0)
            ax.quiver(centroid[0], centroid[1], centroid[2], normal[0], normal[1], normal[2], color='r')
    return ax

def visualize_station(show=False, coverage_file='coverage_2m_coverage.csv', show_start=False, vx=False, convex=False, ax=None, final=False, original=False):
    x,y,z = 1.8, 4.7, 2.7
    start = np.array([2, 3, 4.2])
    if vx: meshes, normals = load_vx_station()
    elif original: meshes = load_original_stl()
    else : meshes = load_meshes(convex=convex)
    coverage = None
    if coverage_file is not None:
        if vx: coverage = np.loadtxt(os.path.join(os.getcwd(), '..', 'data_vx', 'coverage_viewpoint_sets', coverage_file)).flatten()
        elif final: 
            coverage = np.genfromtxt(os.path.join(os.getcwd(), 'coverage', coverage_file.split('/')[0], coverage_file.split('/')[1]), delimiter=1,dtype=int)
        else: coverage = np.loadtxt(os.path.join(os.getcwd(), '..', 'data', 'coverage_viewpoint_sets', coverage_file)).flatten()
        # uncoverable = [ 148, 152, 156, 158, 186, 190, 194, 196, 230, 231, 234, 235, 260, 272, 305, 316, 318, 333, 334, 380, 381, 392, 393, 396, 397, 422, 467, 478, 480, 495, 496, 728, 730, 733, 735, 744, 745, 746, 747, 749, 753, 758, 760, 763, 765, 774, 775, 779, 780, 781, 783]
        # for i in range(coverage.shape[0]):
        #     if i in uncoverable: coverage[i] = 1
        # print('Coverage: ', np.sum(coverage)/len(coverage) * 100, '%')
        # idxs = np.arange(len(coverage))
        # print(''.join([str(x) + '\n' for x in idxs[coverage==0]]))

    if original:
        n_coverage = np.sum(np.array([len(mesh) for mesh in meshes]))
        coverage = np.zeros(n_coverage)
    if convex:
        n_coverage = np.sum(np.array([len(mesh) for mesh in meshes]))
        coverage = np.ones(n_coverage)

    if ax is None:
        fig = plt.figure(figsize=(10,8))
        ax = Axes3D(fig, auto_add_to_figure=False, computed_zorder=False)
        fig.add_axes(ax)
    # ax = obs_normals(ax, meshes, normals, show=False)
    # if convex: coverage = 0 * coverage
    ax = obs_trisurf(meshes, ax, coverage, show=False, alph=0.3)
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