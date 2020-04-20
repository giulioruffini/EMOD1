"""
NE_calculate_normals.py
This script is the .py version of the omonimous script in Matlab.  That script uses surfacenorm to reorganize the faces so that
they all point outwards.
This script uses a similar function by trimesh package. Trimesh and networks must be installed.
It starts using the simplified meshes

Modified by Maria Chiara Biagi (September 2019)

"""

import pygmsh

import os
import sys
import argparse
# import nibabel
import meshio
from sys import platform

if platform == "linux" or platform == "linux2":
    from xvfbwrapper import Xvfb

    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

from mayavi import mlab
import scipy.io as sio

from src.NE_simplify_mesh import simplify_mesh_subj
from src.NE_giu_processing import reduce_giu, fix_giu
from src.NE_convert_mat_npy import mesh_mat2npy, dm_npy2mat
from src.NE_dist_matrix import calculateDM_gm, calculateDM_wm, calculate_parcelDM
import easygui as gui
import trimesh
import networkx
import numpy as np

def ismember_rows(a, b):
    '''Equivalent of 'ismember' from Matlab
    a.shape = (nRows_a, nCol)
    b.shape = (nRows_b, nCol)
    return the idx where b[idx] == a
    '''
#    idx = np.nonzero(np.all(b == a[:, np.newaxis], axis=2))[1]
    # voida, voidb = map(asvoid,(a[:, 0],b[:, 0]))
    # idx = np.where(np.in1d(voida, voidb))[0]
    idx = []
    for j in range(a[:,1].size):
        x = np.argwhere((b[j,:]==a[:]))
        idx.append(x[0, 0])
    return np.asarray(idx)

def calculate_normals(subject_fldr, **kwargs):
    " described above"
    # Check which arguments are present:
    for key in kwargs:
        if key == 'no_cereb':
            no_cb = kwargs[key]
        if key == 'surf':
            surf = kwargs[key]

    #######################################################################################
    #  load the geometry
    if surf == 'gm':
        if no_cb == 1:
            nodes_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                             'verticesGM_nocb.mat'))['verticesGM']
            faces_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                             'facesGM_nocb.mat'))['facesGM']
        elif no_cb == 0:
            nodes_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'verticesGM.mat'))[
                'verticesGM']
            faces_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesGM.mat'))[
                'facesGM']
        NodeStruct = sio.loadmat(os.path.join(subject_fldr, 'Optimization', 'Giu_files', 'Nodes_struct.mat'))['Nodes']
        nodes_orig = NodeStruct[0][0]['nodes_g']
        faces_orig = NodeStruct[0][0]['faces_g']
    elif surf == 'wm':
        if no_cb == 1:
            nodes_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                             'verticesWM_nocb.mat'))['verticesWM']
            faces_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                             'facesWM_nocb.mat'))['facesWM']
        elif no_cb == 0:
            nodes_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'verticesWM.mat'))[
                'verticesWM']
            faces_simpl = sio.loadmat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesWM.mat'))[
                'facesWM']
        NodeStruct = sio.loadmat(os.path.join(subject_fldr, 'Optimization', 'Giu_files', 'Nodes_struct.mat'))['Nodes']
        nodes_orig = NodeStruct[0][0]['nodes_w']
        faces_orig = NodeStruct[0][0]['faces_w']
        print('Surface mesh loaded')
    faces_simpl = np.array(faces_simpl) - 1  # -1 is because it's faces coming from matlab
    faces_orig = np.array(faces_orig) - 1
    print('calculate normals')
    #######################################################################################
    #  re-mesh with Trimesh,  calculate normals on faces and ensure outward orientation,
    # here for simplified geometry
    meshTRIM_simply = trimesh.Trimesh(nodes_simpl, faces_simpl)
    trimesh.repair.fix_inversion(meshTRIM_simply, multibody=True)
    trimesh.repair.fix_normals(meshTRIM_simply, multibody=True)
    nodesTrimesh_simply = meshTRIM_simply.vertices
    facesTrimesh_simply = meshTRIM_simply.faces
    s_snorm = meshTRIM_simply.face_normals
    # here for original geometry
    meshTRIM_orig = trimesh.Trimesh(nodes_orig, faces_orig)
    trimesh.repair.fix_inversion(meshTRIM_orig, multibody=True)
    trimesh.repair.fix_normals(meshTRIM_orig, multibody=True)
    nodesTrimesh_orig = meshTRIM_orig.vertices
    facesTrimesh_orig = meshTRIM_orig.faces
    o_snorm = meshTRIM_orig.face_normals
    #######################################################################################

    # normals in nodes. trimesh alternative to triangles2nodes
    s_snorm_nodes = triangles2nodes(nodesTrimesh_simply, facesTrimesh_simply, s_snorm)
    o_snorm_nodes_tmp = triangles2nodes(nodesTrimesh_orig, facesTrimesh_orig, o_snorm)
    #######################################################################################

    # Checks if normals point outwards or inwards
    if np.mean(s_snorm_nodes[np.where(nodesTrimesh_simply[:, 0] > 0), 0]) > 0:
        s_snorm_nodes = -s_snorm_nodes

    if np.mean(o_snorm_nodes_tmp[np.where(nodesTrimesh_orig[:, 0] > 0), 0]) > 0:
        o_snorm_nodes_tmp = -o_snorm_nodes_tmp

    indices = ismember_rows(np.asarray(nodesTrimesh_orig),nodes_orig)
    o_snorm_nodes = o_snorm_nodes_tmp[indices, :]
    #######################################################################################
    # SAVES the normals and the new meshes into the subject folder
    os.mkdir(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals'))
    if surf == 'gm':

        sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                 'o_norm_gm.mat'), {'o_snorm_nodes': o_snorm_nodes})
        if no_cb == 1:
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                     'simpl_norm_gm_nocb.mat'), {'s_snorm_nodes': s_snorm_nodes})
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                                     'verticesGM_nocb.mat'), {'verticesGM': nodesTrimesh_simply})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesGM_nocb.mat'),
                {'facesGM': np.int32(facesTrimesh_simply)+1})
        elif no_cb == 0:
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                     'simpl_norm_gm.mat'), {'s_snorm_nodes': s_snorm_nodes})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'verticesGM.mat'),
                {'verticesGM': nodesTrimesh_simply})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesGM.mat'),
                {'facesGM': np.int32(facesTrimesh_simply)+1})

    elif surf == 'wm':
        sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                 'o_norm_wm.mat'), {'o_snorm_nodes': o_snorm_nodes})
        if no_cb == 1:
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                     'simpl_norm_wm_nocb.mat'), {'s_snorm_nodes': s_snorm_nodes})
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry',
                                     'verticesWM_nocb.mat'), {'verticesWM': nodesTrimesh_simply})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesWM_nocb.mat'),
                {'facesWM': np.int32(facesTrimesh_simply)+1})
        elif no_cb == 0:
            sio.savemat(os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'Normals',
                                     'simpl_norm_wm.mat'), {'s_snorm_nodes': s_snorm_nodes})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'verticesWM.mat'),
                {'verticesWM': nodesTrimesh_simply})
            sio.savemat(
                os.path.join(subject_fldr, 'Communication', 'NIC_interface', 'Simplified Geometry', 'facesWM.mat'),
                {'facesWM': np.int32(facesTrimesh_simply)+1})

    print('Normals and new geometry saved!')

    #######################################################################################


def triangles2nodes(nodes, faces, norms):
    " This function maps from triangle-defined-data to node-defined-data using weights given by the areas of the triangles."
    norms_in_nodes = np.zeros((nodes.shape[0], 3))
    area = triangle_area(nodes, faces)

    Lia_a = np.isin(np.arange(0, nodes.shape[0]), np.array(faces)[:, 0])
    Lia_b = np.isin(np.arange(0, nodes.shape[0]), np.array(faces)[:, 0])
    Lia_c = np.isin(np.arange(0, nodes.shape[0]), np.array(faces)[:, 0])

    for i in range(nodes.shape[0]):
        Lia_a = np.where(faces[:, 0] == i)
        Lia_b = np.where(faces[:, 1] == i)
        Lia_c = np.where(faces[:, 2] == i)
        Lia_abc = np.concatenate((np.asarray(Lia_a[0]), np.asarray(Lia_b[0]), np.asarray(Lia_c[0])), axis=0)
        Lia_all = np.unique([Lia_abc])

        area_repmat = np.tile(np.asarray(area), (3, 1))
        area_reshaped = area_repmat.transpose()
        norms_in_nodes[i, :] = np.divide(np.sum(norms[Lia_all, :] * area_reshaped[Lia_all, :], axis=0),
                                         np.sum(area_reshaped[Lia_all, :]))

        # normalize data
    norm_of_normals = np.tile(np.sqrt(
        np.power(norms_in_nodes[:, 0], 2) + np.power(norms_in_nodes[:, 1], 2) + np.power(norms_in_nodes[:, 2], 2)),
                              (3, 1))
    norms_in_nodes = np.divide(norms_in_nodes, norm_of_normals.transpose())
    return norms_in_nodes


def triangle_area(nodes, faces):
    "This function calculates triangle area using Heron's formula"
    a = np.power((np.power((nodes[faces[:, 0], 0] - nodes[faces[:, 1], 0]), 2) + np.power(
        (nodes[faces[:, 0], 1] - nodes[faces[:, 1], 1]), 2) + np.power((nodes[faces[:, 0], 2] - nodes[faces[:, 1], 2]),
                                                                       2)), 1 / 2)
    b = np.power((np.power((nodes[faces[:, 0], 0] - nodes[faces[:, 2], 0]), 2) + np.power(
        (nodes[faces[:, 0], 1] - nodes[faces[:, 2], 1]), 2) + np.power((nodes[faces[:, 0], 2] - nodes[faces[:, 2], 2]),
                                                                       2)), 1 / 2)
    c = np.power((np.power((nodes[faces[:, 2], 0] - nodes[faces[:, 1], 0]), 2) + np.power(
        (nodes[faces[:, 2], 1] - nodes[faces[:, 1], 1]), 2) + np.power((nodes[faces[:, 2], 2] - nodes[faces[:, 1], 2]),
                                                                       2)), 1 / 2)
    s = (a + b + c) / 2

    area = np.power(np.multiply(np.multiply(np.multiply(s, (s - a)), (s - b)), (s - c)), 1 / 2)
    return area



