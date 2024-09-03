import os
from os import getcwd
from os.path import join
import numpy as np

def load_single_mesh(mesh_file, new=False, offset=[0,0,0]):
    """
    """
    faces_normals = np.loadtxt(mesh_file, delimiter=' ')
    f = []
    n = []
    offset = np.array(offset)
    for face in faces_normals:
        f.append([tuple(np.array([face[0], face[1], face[2]])+offset), 
                  tuple(np.array([face[3], face[4], face[5]])+offset), 
                  tuple(np.array([face[6], face[7], face[8]])+offset)])
        n.append([face[9], face[10], face[11]])
    return f,n

def load_meshes(convex=False):
    """
    load meshes from /models/remeshed_parts folder using OBS
    """
    meshes = []
    normals = []
    if convex:
        num_meshes = 15
        offset = [2.529, 4.821, 2.591]
    else:
        num_meshes = 10
        offset = [0,0,0]
    mesh_dir = os.path.abspath(join(getcwd(),'..', 'data', 'model_remeshed'))
    for i in range(num_meshes):
        meshfilename = join(mesh_dir, str(i) + '_faces_normals.csv')
        faces, normals = load_single_mesh(meshfilename, offset)
        meshes.append(faces)
        normals.append(normals)

    return meshes