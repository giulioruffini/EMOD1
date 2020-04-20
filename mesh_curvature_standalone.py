# 1. Install miniconda (python 3.7 version, windows 64 bits)
# 2. In a miniconda terminal run conda install igl
# 3. Swap to this python compiler in pycharm


import igl
import numpy as np
import os

#Change this to point to your subject folder
root_folder = r'D:\Neuroelectrics\Other Research topics\Ephaptic matrix project\Code to share\EMOD calculations'
subject_id_usr = os.path.join(root_folder, 'subjects folder','subject_01')

#Make sure you've at least simplified the mesh of the subject
vertices = np.load(subject_id_usr + '/vertices_simple.npy')
faces = np.load(subject_id_usr + '/faces_simple.npy')

ret = igl.write_triangle_mesh(os.path.join(subject_id_usr, "bunny_out.obj"), vertices, faces)

v, f = igl.read_triangle_mesh(os.path.join(subject_id_usr, "bunny_out.obj"))
[pd1, pd2, pv1, pv2] = igl.principal_curvature(v, f, radius = 5, use_k_ring = True)
#Mean curvature is (pv1+pv2)/2 and Gaussian curvature is pv1*pv2
# k = igl.gaussian_curvature(v, f)
np.save(subject_id_usr + '/mean_curv_simple.npy',(pv1+pv2)/2)

#Deletes temporary object
os.remove(os.path.join(subject_id_usr, "bunny_out.obj"))

