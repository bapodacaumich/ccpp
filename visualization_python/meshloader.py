from os import getcwd
from os.path import join
import numpy as np

def load_single_mesh(mesh_file, new=False):
    """
    """
    mesh_dir = join(getcwd(), mesh_file)
    faces_normals = np.loadtxt(mesh_dir, delimiter=' ')
    f = []
    n = []
    for face in faces_normals:
        f.append([(face[0], face[1], face[2]), (face[3], face[4], face[5]), (face[6], face[7], face[8])])
        n.append([face[9], face[10], face[11]])
    return f,n

def load_meshes():
    """
    load meshes from /models/remeshed_parts folder using OBS
    """
    meshes = []
    normals = []
    num_meshes = 10
    mesh_dir = '../data/model_remeshed/'
    for i in range(num_meshes):
        meshfilename = mesh_dir + str(i) + '_faces_normals.csv'
        faces, normals = load_single_mesh(meshfilename)
        meshes.append(faces)
        normals.append(normals)

    return meshes