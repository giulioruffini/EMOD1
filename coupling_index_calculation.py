import os
import MeshProcessingPipeline as mpp
import numpy as np
import tkinter as tk
from tkinter import filedialog
from mayavi import mlab
from progress.bar import Bar

os.environ['ETS_TOOLKIT'] = 'qt4'
# from pyface.qt import QtGui, QtCore
# from traits.api import HasTraits, Instance, on_trait_change, Str, Float, Range, Trait, Enum
# from traitsui.api import View, Item, HSplit, Group, VSplit
# from mayavi.core.api import PipelineBase, Engine
# from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor

#Simpler matrix to linear...Check notebook Sketch 1 on 31-10-2018
def m2l(i,j):
    positions = (i-1)*(i-1+3)/2+1+j

    positions = np.matrix(positions)

    return positions.getA1()


def coupling_index_visualizer(l3_vis = 20):
    #First load DM and mesh
    initial_subj_folder = r'\\antares.starlab.es\Head_Models'
    root = tk.Tk()
    root.withdraw()
    subject_folder_path = filedialog.askdirectory(initialdir = initial_subj_folder, title = "Select subject folder")
    root.destroy()
    nodes = np.load(subject_folder_path + r'\Communication\NIC_interface\Simplified Geometry\verticesGM.npy')
    faces = np.load(subject_folder_path + r'\Communication\NIC_interface\Simplified Geometry\facesGM.npy')

    #Calculates EMOD
    print('Calculating EMOD for this subject!')
    EMOD_hm = mpp.EMOD1Calc(subj_folder = subject_folder_path, l3 = l3_vis,mode = 'mesh',sigma = 0.33)

    EMOD_sum = np.sum(EMOD_hm)*1e3

    print('EMOD calculated!')
    print('Subject EMOD1: {0:3e} microV'.format(EMOD_sum*2))

    #Number of mesh points
    n_points = len(nodes)

    #Creates position matrices...
    x = np.arange(n_points)

    yv, xv = np.meshgrid(x, x)

    places = np.where(yv > xv)

    yv_tmp = yv.copy()

    yv[places] = xv[places]
    xv[places] = yv_tmp[places]

    #Variable to plot
    clim_2_plot0 = np.zeros(n_points)

    bar = Bar('Processing', max=n_points)
    for i in np.arange(0,n_points,1):
        positions = m2l(xv[i,:], yv[i,:])
        clim_2_plot0[i] = np.sum(EMOD_hm[positions.astype(int)])*1e3
        bar.next()

    bar.finish()

    # max_val = np.max(clim_2_plot0)
    #
    #
    # faces_new = faces.copy()
    # other_faces = faces.copy()
    #
    # for i in  np.arange(0,n_points,1):
    #     if clim_2_plot0[i]/max_val<=0.001:
    #         faces_new = np.where(faces_new==i,-1,faces_new)
    #     elif clim_2_plot0[i]/max_val>=0.001:
    #         other_faces = np.where(other_faces==i,-1,other_faces)
    #
    # faces_new = faces_new[~np.any(faces_new == -1, axis=1)]
    # other_faces = other_faces[~np.any(other_faces == -1, axis=1)]

    fig1 = mlab.figure(1)
    #src = mlab.pipeline.triangular_mesh_source(nodes[:,0],nodes[:,1],nodes[:,2],faces,scalars = DM[0,:])
    #src = mlab.pipeline.triangular_mesh_source(nodes[:, 0], nodes[:, 1], nodes[:, 2], faces_new)
    src = mlab.pipeline.triangular_mesh_source(nodes[:, 0], nodes[:, 1], nodes[:, 2], faces)
    #src.representation = 'wireframe'
    src.mlab_source.scalars = clim_2_plot0
    surf = mlab.pipeline.surface(src)
    #mlab.pipeline.surface(src, vmin=0, vmax=np.max(clim_2_plot0))
    mlab.pipeline.surface(src, vmin=0, vmax=0.01)
    mlab.colorbar(title='Ephaptic coupling index on surface (microV)', orientation='vertical')

    #, representation='wireframe'
    #src_wireframe = mlab.triangular_mesh(nodes[:, 0], nodes[:, 1], nodes[:, 2], faces,scalars=None,line_width=0.1,color=(0.5, 0.5, 0.5),opacity = 0.55)
    mlab.triangular_mesh(nodes[:, 0], nodes[:, 1], nodes[:, 2], faces,representation='wireframe',scalars=None,line_width=0.1,color=(0.5, 0.5, 0.5),opacity = 0.55)

    mlab.view(azimuth=360, elevation=0, distance=None, focalpoint=None, roll=None, reset_roll=True, figure=None)
    # cursor3d = mlab.points3d(nodes[1,0], nodes[1,1], nodes[1,2], mode='axes',color=(0, 0, 0),scale_factor=1.0)
    #
    # def picker_callback(picker_obj):
    #     picked = picker_obj.actors
    #     if surf.actor.actor._vtk_obj in [o._vtk_obj for o in picked]:
    #         print('Node coordinates (index: {0:3d}): ({1:3f},{2:3f},{3:3f}) mm'.format(picker_obj.point_id,nodes[picker_obj.point_id,0],nodes[picker_obj.point_id,1],nodes[picker_obj.point_id,2]))
    #         #mlab.points3d(nodes[picker_obj.point_id,0],nodes[picker_obj.point_id,1],nodes[picker_obj.point_id,2], mode='2dthick_cross',color=(0, 0, 0),scale_factor=10)
    #         cursor3d.mlab_source.reset(x=nodes[picker_obj.point_id,0], y=nodes[picker_obj.point_id,1], z=nodes[picker_obj.point_id,2])
    #         #affected_indexes = matrix_2_linear(picker_obj.point_id,':',n_points)
    #
    #
    #         print('EMOD1 node value: {0:3e} microV'.format(clim_2_plot0[picker_obj.point_id]))
    #
    #         src.mlab_source.scalars = clim_2_plot0
    #         src.mlab_source.update()
    #         #mlab.pipeline.surface(src, vmin=0, vmax=np.max(clim_2_plot0))
    #         mlab.pipeline.surface(src, vmin=0, vmax=0.01)
    #         mlab.colorbar(title='Ephaptic coupling index on surface (microV)', orientation='vertical')
    #         mlab.view(azimuth=360, elevation=0, distance=None, focalpoint=None, roll=None, reset_roll=True,figure=None)
    #
    # fig1.on_mouse_pick(picker_callback)
    mlab.show()

coupling_index_visualizer(l3_vis = 20)