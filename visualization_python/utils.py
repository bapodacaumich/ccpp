import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def obs_trisurf(meshes, coverage=None, show=True, surface=True, wireframe=True, lw=0.1, alph=0.2):
    """
    plots a list of OBS objects by face
    """
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    coverage_idx = 0
    num_faces = np.sum(np.array([len(mesh) for mesh in meshes]))
    n_coverage = len(coverage) if coverage is not None else 0
    highlight_faces = [756,1032] # for debugging purposes
    highlight_faces = []
    assert num_faces == n_coverage, f'coverage length {n_coverage} does not match number of faces {num_faces}'
    for mesh in meshes:
        for face in mesh:
            x = [x[0] for x in face]
            y = [x[1] for x in face]
            z = [x[2] for x in face]
            verts = [list(zip(x,y,z))]

            if coverage is not None:
                if coverage[coverage_idx] == 1:
                    if coverage_idx in highlight_faces:
                        pc = Poly3DCollection(verts, fc='tab:purple', alpha=1.0)
                    else:
                        pc = Poly3DCollection(verts, fc='tab:blue', alpha=alph)
                else:
                    if coverage_idx in highlight_faces:
                        pc = Poly3DCollection(verts, fc='tab:green', alpha=1.0)
                    else:
                        pc = Poly3DCollection(verts, fc='red', alpha=alph)
            else:
                pc = Poly3DCollection(verts, fc='tab:blue', alpha=alph)
            if surface: ax.add_collection3d(pc)
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            if wireframe: ax.plot(x,y,z,lw=0.01, c='k')
            coverage_idx += 1

    # ax.set_aspect('equal')
    if show == True: plt.show()
    return ax

def coverage_area(self, pose, d=0.5, hfov=60, vfov=60):
    """
    get bounds of the coverage area for a given pose
    pose - np.array(x, y, z, pan(ccw around z), tilt(positive up)) assume no roll/swing angle
        > Note: pan is bounded by -pi to pi, tilt goes 0 to pi
    """
    # compute rotation matrices
    R_z = np.array([[np.cos(pose[3]), -np.sin(pose[3]), 0.],
                    [np.sin(pose[3]),  np.cos(pose[3]), 0.],
                    [0.,               0.,              1.]])
    R_x = np.array([[1.,              0.,               0.],
                    [0., np.cos(pose[4]), -np.sin(pose[4])],
                    [0., np.sin(pose[4]),  np.cos(pose[4])]])

    dx = d*np.sin(hfov/2)
    dy = d*np.sin(vfov/2)
    dz = -np.sqrt(d**2 - dx**2 - dy**2)
    tl = R_z @ R_x @ np.array([-dx, dy, dz]) + pose[:3] # top left
    tr = R_z @ R_x @ np.array([dx, dy , dz]) + pose[:3] # top right
    br = R_z @ R_x @ np.array([dx, -dy , dz]) + pose[:3] # bottom right
    bl = R_z @ R_x @ np.array([-dx, -dy , dz]) + pose[:3] # bottom left
    ct = R_z @ R_x @ np.array([0, 0, -d]) + pose[:3] # center
    return tl, tr, br, bl, ct


# def plot_path_direct(cam_dist, local, folder='always_start_best', station=True, ax=None):
def plot_path_direct(folder, file, ax=None):
    filepath = os.path.join(os.getcwd(), '..', 'data', folder, file)
    # filepath = os.getcwd() + '\\data\\coverage_viewpoint_sets\\' + file
    # for file in os.listdir(dir):
    #     if file.endswith(".npy"): continue
    #     if (cam_dist == float(file[:3])):
    #         if not ((file[5] == 'l') ^ local):
    #             filepath=os.path.join(dir, file)
    # print('found file: ', filepath)
    viewpoints = np.loadtxt(filepath, delimiter=',')
    print(viewpoints.shape)
    if ax is None: ax = plt.figure(figsize=(8, 8)).add_subplot(projection='3d')
    ax.plot(viewpoints[:,0], viewpoints[:,1], viewpoints[:,2], 'k-')
    ax.scatter(viewpoints[:,0], viewpoints[:,1], viewpoints[:,2], c='k', alpha=0.2)
    # ax.quiver(viewpoints[:,0], viewpoints[:,1], viewpoints[:,2], viewpoints[:,3], viewpoints[:,4], viewpoints[:,5], length=1)
    weights = np.linspace(0,1,viewpoints.shape[0])
    # weights = np.ones(viewpoints.shape[0])
    sc = ax.scatter(*[viewpoints[:,i] for i in range(3)], c=weights, cmap='YlOrBr', alpha=0.9)
    ax.quiver(*[viewpoints[:,i] for i in range(viewpoints.shape[1]-1)], length=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # plt.colorbar(sc, shrink=0.5, pad=-0.04, label='Path Progression')
    return ax

def set_aspect_equal_3d(axes):
    # get aspect ratios
    xlim = axes.get_xlim()
    ylim = axes.get_ylim()
    zlim = axes.get_zlim()
    axes.set_box_aspect([xlim[1]-xlim[0], ylim[1]-ylim[0], zlim[1]-zlim[0]])