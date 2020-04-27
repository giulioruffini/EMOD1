'''EMOD calculation pipeline
2020 Neuroelectrics Corporation
@autor: Ricardo Salvador (NE)
        Giulio Ruffini (NE)


This module currently includes the following functionalities:
-Simplification of gm and wm surface meshes
-Calculation of distance matrices (geodesic distance) from meshes
-Reconstruction of distance matrices
-Reduction and fixing of giu files
-Converting simplified meshes in a subject folder from .npy to .mat (and vice versa)
 to allow working in both pipelines easily
-Converting distance matrix files in .mat to .npy to use the data with the python pipeline

Required libraries:
-pymeshfix
-mayavi
-scipy
-easygui
-numpy
-tvb-gdist
-h5py
'''

from multiprocessing import Pool
import pymeshfix as fix
import numpy as np
import itertools
import gdist
import time
import gc
import os
import trimesh
import progressbar
import math
import warnings
from mayavi import mlab

# Class Subject
##########Multiprocess test function for DM calculation
class subject:
    def __init__(self, subject_id, reductions_ratio, quiet):
        # import numpy as np
        # import easygui as gui
        self.Subj_folder = subject_id

        # Loads original mesh
        faces_full = np.loadtxt(os.path.join(self.Subj_folder, 'faces_GM_full.txt'))
        self.faces_full = faces_full[:, 0:3] - 1  # Removes 1 because faces index start from 1 in Matlab

        vertices_full = np.loadtxt(os.path.join(self.Subj_folder,'nodes_GM_full.txt'))
        self.vertices_full = vertices_full[:, 0:3]

        self.faces_simple = []  # faces of the simplified mesh
        self.vertices_simple = []  # vertices of the simplified mesh
        self.normals_simple = []  # normals of each node
        self.area_faces_simple = []  # area of each triangle in the simplified mesh
        self.area_nodes_simple = []  # area of each node in the simplified mesh

        self.reductions_ratio = reductions_ratio  # reduction ratio.
        self.quiet = quiet

        # Euclidean distance half-matrix
        # self.euclid_dist = []
        # self.geodesic_dist = []

    def SimplifyMeshSubj(self):
        '''Function to simplify a mesh from a subject folder, the simplification also cleans problems and
        suppresses small disconnected objects. The meshes created are displayed for the user to check that
        there aren´t any problems, close the window of the displayed mesh for the script to continue.
        The function doesn´t return outputs, it saves the new faces
        and vertices in .npy files in the subject folder.

        Function call:
        SimplifyMeshSubj()'''
        ###### Mesh simplification

        warnings.filterwarnings("ignore")

        proceed = True

        if os.path.isfile(os.path.join(self.Subj_folder, 'vertices_simple.npy')) and\
                os.path.isfile(os.path.join(self.Subj_folder,'faces_simple.npy')):
            print('Simplified mesh already exists, do you want to re-calculate?')
            proceed_ = input("Yes(y)/No(n): ")

            if proceed_ == 'y' or proceed_ == 'Y':
                proceed = True
            else:
                proceed = False

        if proceed:
            # Plots the original surf
            # Plot the final reduced mesh

            if self.quiet == 0:
                print('Loaded mesh. Close window to proceed!')
                mlab.triangular_mesh(self.vertices_full[:, 0], self.vertices_full[:, 1], self.vertices_full[:, 2],
                                     self.faces_full)
                mlab.show()

            ######## GM ########
            print('Simplifying cortical surface')
            ####### Initial Mesh cleaning (fixes problems and deletes disconnected components)
            repaired = fix.MeshFix(self.vertices_full, self.faces_full)
            repaired.repair()
            vertices_simple = repaired.v
            faces_simpl = repaired.f

            ######## Mesh decimation ########
            x = vertices_simple[:, 0]
            y = vertices_simple[:, 1]
            z = vertices_simple[:, 2]
            mesh = mlab.pipeline.triangular_mesh_source(x, y, z, faces_simpl)  # creates mesh structure
            dec = mlab.pipeline.quadric_decimation(mesh)  # Creates the decimation filter
            dec.filter.target_reduction = self.reductions_ratio  # Set the decimation value (reduces to 25%,changed due to post-fixing)
            new_nodes = np.array(dec.outputs[0]._get_output()._get_points())  # Get the nodes
            new_faces = np.array(dec.outputs[0]._get_output()._get_polys()._get_data()).astype(int)  # get the faces
            bad_index = range(0, len(new_faces),4)  # values to discard (probably reference to the number of nodes in the face)
            new_faces = np.delete(new_faces, bad_index)  # obtain the real face values
            new_faces = np.reshape(new_faces, (-1, 3))  # generate the faces matrix

            ####### Final Mesh cleaning (fixes problems and deletes disconnected components)
            repaired = fix.MeshFix(new_nodes, new_faces)
            repaired.repair()
            vertices_simple = repaired.v
            faces_simple = repaired.f
            if self.quiet == 0:
                # Plot the final reduced mesh
                print('Check the mesh! Close the window to proceed!')
                mlab.triangular_mesh(vertices_simple[:, 0], vertices_simple[:, 1], vertices_simple[:, 2], faces_simple)
                mlab.show()

            print('Done!')

            # Calculates normals!
            print('Calculating faces normals!')
            # Calculates faces normals with trimesh
            meshTRIM_simply = trimesh.Trimesh(vertices_simple, faces_simple)
            trimesh.repair.fix_inversion(meshTRIM_simply, multibody=True)
            trimesh.repair.fix_normals(meshTRIM_simply, multibody=True)

            #Updates attributes of the class!
            self.vertices_simple = meshTRIM_simply.vertices
            self.faces_simple = np.int32(meshTRIM_simply.faces)
            s_snorm = meshTRIM_simply.face_normals

            # Saves nodes and faces!
            np.save(os.path.join(self.Subj_folder, 'vertices_simple.npy'), self.vertices_simple)

            np.save(os.path.join(self.Subj_folder, 'faces_simple.npy'), self.faces_simple)

            print('Done!')

            print('Calculating faces/normals areas!')
            self.Area_calc()
            print('Done!')

            print('Calculates normals in each node!')
            self.normals_simple = self.triangles2nodes(data_2_remap=s_snorm, normalize=True)
            np.save(os.path.join(self.Subj_folder, 'normals_simple.npy'), self.normals_simple)

            print('Done!')

            print('Calculates dot products of normals (convenient for EMOD calculation)!')

            n = self.normals_simple.shape[0]
            normals_vector = np.zeros(np.int(0.5 * n * (n - 1) + n))

            finish_vector = 1
            init_vector = 0
            final_node = 1

            for i in self.normals_simple:
                nodes2use = self.normals_simple[0:final_node, ]
                dot_norm = np.dot(nodes2use, i)

                normals_vector[init_vector:finish_vector] = dot_norm
                init_vector = finish_vector
                finish_vector = init_vector + dot_norm.shape[0] + 1
                final_node = final_node + 1

            np.save(os.path.join(self.Subj_folder, 'dot_normals_half.npy'), np.float32(normals_vector))


            print('Done!')

            print('Calculates products of areas (convenient for EMOD calculation)!')

            n = self.area_nodes_simple.shape[0]
            areas_vector = np.zeros(np.int(0.5 * n * (n - 1) + n))
            single_areas_vector = np.zeros(np.int(0.5 * n * (n - 1) + n))

            finish_vector = 1
            init_vector = 0
            final_node = 1

            for i in progressbar.progressbar(self.area_nodes_simple):
                nodes2use = self.area_nodes_simple[0:final_node, ]
                dot_norm = nodes2use * i
                dot_norm_single = i * np.ones(np.shape(nodes2use))

                areas_vector[init_vector:finish_vector] = dot_norm
                single_areas_vector[init_vector:finish_vector] = dot_norm_single
                init_vector = finish_vector
                finish_vector = init_vector + dot_norm.shape[0] + 1
                final_node = final_node + 1

            np.save(os.path.join(self.Subj_folder, 'dot_areas_half.npy'), np.float32(areas_vector))
            np.save(os.path.join(self.Subj_folder, 'dot_areas_single_half.npy'), np.float32(single_areas_vector))


            print('Done!')

            print('Saved simplified mesh and its attributes in subject folder!')

        else:
            self.vertices_simple = np.load(os.path.join(self.Subj_folder, 'vertices_simple.npy'))
            self.faces_simple = np.load(os.path.join(self.Subj_folder, 'faces_simple.npy'))
            self.normals_simple = np.load(os.path.join(self.Subj_folder, 'normals_simple.npy'))
            self.area_nodes_simple = np.load(os.path.join(self.Subj_folder, 'area_per_node.npy'))
            self.area_faces_simple = np.load(os.path.join(self.Subj_folder, 'area_per_face.npy'))


    def _DM_operation(self, iteration):

        OneSourceDistance = gdist.compute_gdist(self.vertices_simple, self.faces_simple,
                                                source_indices=np.array([iteration]),
                                                target_indices=np.array(range(
                                                    iteration + 1)))  # Distance from node 0 to the rest (around 0.8 s for 43350 nodes
        vector = np.float32(OneSourceDistance)
        return vector

    def EuclidDistCalcSubj(self):
        '''Function to calculate euclidean distance in simplified mesh. Assumes that a simplified mesh has been calculated
            and saved in the subject's folder.

            Function call:
            EuclidDistCalcSubj(subject = 'jane_doe')'''

        proceed = True

        if os.path.isfile((self.Subj_folder + r'\dm_euclid_half.npy')):
            print('Euclidean distance matrix already exists, do you want to re-calculate?')
            proceed_ = input("Yes(y)/No(n): ")

            if proceed_ == 'y' or proceed_ == 'Y':
                proceed = True
            else:
                proceed = False

        if proceed:

            print('Calculating Euclidean distance matrix!')

            # Calculates Euclidean matrix
            n = self.vertices_simple.shape[0]
            dm_vector = np.zeros(np.int(0.5 * n * (n - 1) + n))
            finish_vector = 1
            init_vector = 0
            final_node = 1
            for i in self.vertices_simple:
                nodes2use = self.vertices_simple[0:final_node, ]
                distances = i - nodes2use
                distances = distances ** 2
                distances = np.sum(distances, 1)
                distances = np.sqrt(distances)
                dm_vector[init_vector:finish_vector] = distances
                init_vector = finish_vector
                finish_vector = init_vector + distances.shape[0] + 1
                final_node = final_node + 1

            # Saves matrix
            np.save(os.path.join(self.Subj_folder, 'dm_euclid_half.npy'), np.float32(dm_vector))


            print('Euclidean distance matrix saved!')
        else:
            print('Euclidean distance matrix already exists. Skipping calculation!')

    def GeodesicDistCalcSubj(self, divisions=120):
        '''Function to calculate geodesic distance in simplified mesh. Assumes that a simplified mesh has been calculated
            and saved in the subject's folder.

            Function call:
            EuclidDistCalcSubj(subject = 'jane_doe')'''

        proceed = True

        if os.path.isfile((self.Subj_folder + r'\dm_geodesic_half.npy')):
            print('Geodesic distance matrix already exists, do you want to re-calculate?')
            proceed_ = input("Yes(y)/No(n): ")

            if proceed_ == 'y' or proceed_ == 'Y':
                proceed = True
            else:
                proceed = False

        if proceed:

            n = self.vertices_simple.shape[0]
            val = n / divisions
            val_pieces = np.int(np.floor(val))
            dm_vector = np.zeros(np.int(0.5 * n * (n - 1) + n))

            print(
                'Calculating the Geodesic distance matrix. This may take a few hours! You will need a computer with enough RAM!')
            gc.collect()
            pool = Pool()
            vector = pool.map(self._DM_operation, range(val_pieces))
            pool.close()
            pool.join()

            merged = list(itertools.chain.from_iterable(vector))
            l = len(merged)
            dm_vector[0:l] = merged

            t = time.time()
            for i in progressbar.progressbar(range(1, divisions)):
                gc.collect()
                pool = Pool()
                vector = pool.map(self._DM_operation, range(val_pieces * i, val_pieces * (i + 1)))
                merged = list(itertools.chain.from_iterable(vector))
                dm_vector[l:l + len(merged)] = merged
                l = l + len(merged)
                pool.close()
                pool.join()

            gc.collect()
            pool = Pool()
            finalvector = pool.map(self._DM_operation, range(divisions * val_pieces, n))
            pool.close()

            merged = list(itertools.chain.from_iterable(finalvector))
            dm_vector[l:] = merged
            t = time.time() - t
            print('Geodesic distance Matrix calculated after ' + str(t) + ' s')
            print('Saving Half Distance Matrix in subject folder, it may take some time...')

            np.save(os.path.join(self.Subj_folder, 'dm_geodesic_half.npy'), np.float32(dm_vector))


        else:
            print('Geodesic distance matrix already exists. Skipping calculation!')

    def Area_calc(self):

        proceed = True

        if os.path.isfile((self.Subj_folder + r'\area_per_node.npy')):
            print('Areas per node are already calculated, do you want to re-calculate?')
            proceed_ = input("Yes(y)/No(n): ")

            if proceed_ == 'y' or proceed_ == 'Y':
                proceed = True
            else:
                proceed = False

        if proceed:
            "This function calculates triangle area using Heron's formula"
            a = np.power((np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 0] - self.vertices_simple[self.faces_simple[:, 1], 0]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 1] - self.vertices_simple[self.faces_simple[:, 1], 1]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 2] - self.vertices_simple[self.faces_simple[:, 1], 2]),
                2)), 1 / 2)
            b = np.power((np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 0] - self.vertices_simple[self.faces_simple[:, 2], 0]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 1] - self.vertices_simple[self.faces_simple[:, 2], 1]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 0], 2] - self.vertices_simple[self.faces_simple[:, 2], 2]),
                2)), 1 / 2)
            c = np.power((np.power(
                (self.vertices_simple[self.faces_simple[:, 2], 0] - self.vertices_simple[self.faces_simple[:, 1], 0]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 2], 1] - self.vertices_simple[self.faces_simple[:, 1], 1]),
                2) + np.power(
                (self.vertices_simple[self.faces_simple[:, 2], 2] - self.vertices_simple[self.faces_simple[:, 1], 2]),
                2)), 1 / 2)
            s = (a + b + c) / 2

            # Area per triangle/face
            self.area_faces_simple = np.power(np.multiply(np.multiply(np.multiply(s, (s - a)), (s - b)), (s - c)),1 / 2)

            # Saves face area
            np.save(os.path.join(self.Subj_folder, 'area_per_face.npy'), self.area_faces_simple)



            # Area per node
            self.area_nodes_simple = np.zeros(np.shape(self.vertices_simple[:, 0]))
            area_tmp = []
            for i in progressbar.progressbar(range(self.vertices_simple.shape[0])):
                test = np.where(np.in1d(self.faces_simple[:, 0], i) | np.in1d(self.faces_simple[:, 1], i) | np.in1d(self.faces_simple[:, 2], i))
                area_tmp.append(self.area_faces_simple[test])
                sum_areas = np.sum(np.array(area_tmp))
                self.area_nodes_simple[i] = sum_areas / 3
                area_tmp = []

            np.save(os.path.join(self.Subj_folder, 'area_per_node.npy'), self.area_nodes_simple)


        else:  # Loads the area info and updates the appropriate attribute
            self.area_nodes_simple = np.load(os.path.join(self.Subj_folder, 'area_per_node.npy'))


            self.area_faces_simple = np.load(os.path.join(self.Subj_folder, 'area_per_face.npy'))

    def triangles2nodes(self, data_2_remap=[], normalize=False):
        " This function maps from triangle-defined-data to node-defined-data using weights given by the areas of the triangles."

        data_in_nodes = np.zeros((self.vertices_simple.shape[0], 3))

        if len(data_2_remap.shape) == 2:
            print('Data has ' + str(data_2_remap.shape[0]) + ' rows and ' + str(data_2_remap.shape[1]) + ' columns!')
        elif len(data_2_remap.shape) == 1:
            print('Data has ' + str(data_2_remap.shape[0]) + ' rows!')

        for i in progressbar.progressbar(range(self.vertices_simple.shape[0])):
            Lia_a = np.where(self.faces_simple[:, 0] == i)
            Lia_b = np.where(self.faces_simple[:, 1] == i)
            Lia_c = np.where(self.faces_simple[:, 2] == i)
            Lia_abc = np.concatenate((np.asarray(Lia_a[0]), np.asarray(Lia_b[0]), np.asarray(Lia_c[0])), axis=0)
            Lia_all = np.unique([Lia_abc])

            if len(data_2_remap.shape) == 2:
                for j in range(data_2_remap.shape[1]):
                    data_in_nodes[i, j] = np.divide(
                        np.sum(data_2_remap[Lia_all, j] * self.area_faces_simple[Lia_all], axis=0),
                        np.sum(self.area_faces_simple[Lia_all]))
            else:
                data_in_nodes[i, j] = np.divide(np.sum(data_2_remap[Lia_all] * self.area_faces_simple[Lia_all], axis=0),np.sum(self.area_faces_simple[Lia_all]))

        # normalize data
        if normalize: #only makes sense if it's a vector!
            if len(data_2_remap.shape) == 2 and data_2_remap.shape[1] == 3:
                norm_of_data = np.tile(np.sqrt(
                    np.power(data_in_nodes[:, 0], 2) + np.power(data_in_nodes[:, 1], 2) + np.power(data_in_nodes[:, 2], 2)),
                    (3, 1))
                data_in_nodes = np.divide(data_in_nodes, norm_of_data.transpose())
            else:
                print('Can only normalize vector data (nx3 array)!')

        return data_in_nodes

    def plot_data(self, data_2_plot=[], node_id=[]):

        warnings.filterwarnings("ignore")

        if data_2_plot == 'normals_x':
            data2plot = self.normals_simple[:,0]
            title2plot = 'x component of the vector normal to the surface'
            lim2plot = np.array([-1, 1])
            print('Plotting x component of normal vector!')
        elif data_2_plot == 'normals_y':
            data2plot = self.normals_simple[:,0]
            title2plot = 'y component of the vector normal to the surface'
            lim2plot = np.array([-1, 1])
            print('Plotting y component of normal vector!')
        elif data_2_plot == 'normals_z':
            data2plot = self.normals_simple[:,0]
            title2plot = 'z component of the vector normal to the surface'
            lim2plot = np.array([-1, 1])
            print('Plotting z component of normal vector!')
        elif data_2_plot == 'area':
            data2plot = self.area_nodes_simple
            title2plot = 'Area of each node (mm^2)'
            lim2plot = np.array([np.min(self.area_nodes_simple), np.max(self.area_nodes_simple)])
            print('Plotting area of mesh elements (per node)!')
        elif data_2_plot == 'EMOD':
            EMOD_half = np.load(os.path.join(self.Subj_folder, 'EMOD_half.npy'))


            # Number of mesh points
            n_points = len(self.vertices_simple)

            # Creates position matrices...
            x = np.arange(n_points)

            yv, xv = np.meshgrid(x, x)

            places = np.where(yv > xv)

            yv_tmp = yv.copy()

            yv[places] = xv[places]
            xv[places] = yv_tmp[places]

            # Variable to plot
            data2plot = np.zeros(n_points)

            for i in progressbar.progressbar(np.arange(0, n_points, 1)):
                positions = m2l(xv[i, :], yv[i, :])
                data2plot[i] = np.sum(EMOD_half[positions.astype(int)]) * 1e3

            title2plot = 'EMOD index ($\micro V$)'
            lim2plot = np.array([0, self.EMOD_sum*1000])

            print('Plotting EMOD!')
        elif data_2_plot == 'mean_curv':
            data2plot = self.mean_curv
            title2plot = 'Mean curvature (mm^-1)'
            lim2plot = np.array([-0.5/math.sqrt(2), 0.5/math.sqrt(2)])

        mlab.figure(1)

        src = mlab.pipeline.triangular_mesh_source(self.vertices_simple[:, 0], \
                                                   self.vertices_simple[:, 1], \
                                                   self.vertices_simple[:, 2], \
                                                   self.faces_simple)

        src.mlab_source.scalars = data2plot
        mlab.pipeline.surface(src)
        mlab.pipeline.surface(src, vmin=lim2plot[0], vmax=lim2plot[1])
        mlab.colorbar(title=title2plot, orientation='vertical')

        mlab.triangular_mesh(self.vertices_simple[:, 0], self.vertices_simple[:, 1], self.vertices_simple[:, 2], \
                             self.faces_simple, representation='wireframe', scalars=None,
                             line_width=0.1, color=(0.5, 0.5, 0.5), opacity=0.55)

        mlab.view(azimuth=360, elevation=0, distance=None, focalpoint=None, roll=None, reset_roll=True, figure=None)

        print('Close window to proceed!')

        mlab.show()

    def calcEMOD(self,lambda_sc = 1, p0 = 0.5, sigma = 0.4, l0 = 1, version = '1'):
        '''Calculates EMOD for the subject and prints the value!'''

        print('Calculating EMOD!')

        if version == '1' or version == '1a':
            print('Version of EMOD: ' + version)
            print('Lambda (mm): ' + str(lambda_sc))
            print('Dipole moment surface density (nA*m/mm^2): ' + str(p0))
            print('Electrical conductivity of GM (S/m): ' + str(sigma))
            print('Nearest neighbors distance threshold (mm): ' + str(l0))

            k_const = lambda_sc * p0 / (2 * np.pi * sigma)  # It will calculate everything in mV

            if version == '1a':
                # Areas in half matrix representation
                areas_single_hm = np.load(os.path.join(self.Subj_folder, 'dot_areas_single_half.npy'))

        elif version == '2':
            print('Version of EMOD: ' + version)
            print('Lambda (mm/mm^2): ' + str(lambda_sc))
            print('Dipole moment surface density (nA*m/mm^2): ' + str(p0))
            print('Electrical conductivity of GM (S/m): ' + str(sigma))
            print('Nearest neighbors distance threshold (mm): ' + str(l0))

            k_const = lambda_sc * p0 / (2 * np.pi * sigma)  # It will calculate everything in mV

            #Loads half matrices only required for this version of emod

            # Product of areas
            areas_hm = np.load(os.path.join(self.Subj_folder, 'dot_areas_half.npy'))


            # Sum of mean curvatures
            mean_curv_hm = np.load(os.path.join(self.Subj_folder, 'mean_curv_sum.npy'))



        #Loads Distance matrices
        #DM_geodesic = np.load(self.Subj_folder + r'\dm_geodesic_half.npy')
        DM_euclid = np.load(os.path.join(self.Subj_folder, 'dm_euclid_half.npy'))


        #Loads normals dot product matrix
        normals_hm = np.load(os.path.join(self.Subj_folder, 'dot_normals_half.npy'))


        # Calculates EMOD
        # Auxiliary terms
        diag_zeros_euclid = np.where(DM_euclid == 0)

        # Number of mesh points
        n_points = len(self.vertices_simple)

        if version == '1':
            cim = k_const * (1 / DM_euclid ** 3) * np.abs(normals_hm)
            diag_zeros = np.where(np.logical_or(DM_euclid > l0, normals_hm > 0))
            cim[diag_zeros] = 0
            cim[diag_zeros_euclid] = 0
            self.EMOD_sum = 2*np.sum(cim)/n_points
            print('EMOD 1 value: ' + "{0:.4g}".format(self.EMOD_sum*1000) + ' microV')

            #Saves EMOD mesh file (EMOD for each node of the simplified mesh)
            np.save(os.path.join(self.Subj_folder, 'EMOD_half.npy'),np.float32(cim))



        elif version == '1a':
            cim = k_const * (1 / DM_euclid ** 3) * np.abs(normals_hm) * areas_single_hm
            diag_zeros = np.where(np.logical_or(DM_euclid > l0, normals_hm > 0))
            cim[diag_zeros] = 0
            cim[diag_zeros_euclid] = 0
            self.EMOD_sum = 2*np.sum(cim)/n_points
            print('EMOD 1a value: ' + "{0:.4g}".format(self.EMOD_sum*1000) + ' microV')
            #Saves EMOD mesh file (EMOD for each node of the simplified mesh)
            np.save(os.path.join(self.Subj_folder, 'EMOD_half.npy'),np.float32(cim))



        elif version == '2':
            cim = k_const * (1 / DM_euclid ** 3) * np.abs(normals_hm) * areas_hm * np.exp(2*mean_curv_hm*DM_euclid)
            diag_zeros = np.where(np.logical_or(DM_euclid > l0, normals_hm > 0))
            cim[diag_zeros] = 0
            cim[diag_zeros_euclid] = 0
            self.EMOD_sum = 2 * np.sum(cim) / n_points
            print('EMOD 2 value: ' + "{0:.4g}".format(self.EMOD_sum * 1000) + ' microV')

            # Saves EMOD mesh file (EMOD for each node of the simplified mesh)
            np.save(os.path.join(self.Subj_folder, 'EMOD_half.npy'), np.float32(cim))


    def curvatures(self):

        #Mean curvature loaded if available!
        self.mean_curv = np.load(os.path.join(self.Subj_folder, 'mean_curv_simple.npy'))


        proceed = True

        if os.path.isfile((self.Subj_folder + r'\mean_curv_sum.npy')):
            print('Sums of mean curvature terms are already calculated, do you want to re-calculate?')
            proceed_ = input("Yes(y)/No(n): ")

            if proceed_ == 'y' or proceed_ == 'Y':
                proceed = True
            else:
                proceed = False

        if proceed:
            n = self.mean_curv.shape[0]
            sum_curvs_hm = np.zeros(np.int(0.5 * n * (n - 1) + n))

            finish_vector = 1
            init_vector = 0
            final_node = 1

            for i in self.mean_curv:
                nodes2use = self.mean_curv[0:final_node, ]
                sum_curvs = i + nodes2use

                sum_curvs_hm[init_vector:finish_vector] = sum_curvs
                init_vector = finish_vector
                finish_vector = init_vector + sum_curvs.shape[0] + 1
                final_node = final_node + 1

            np.save(os.path.join(self.Subj_folder, 'mean_curv_sum.npy'), np.float32(sum_curvs_hm))




#Auxiliary function
def m2l(i,j):
    '''Converts from the position in the matrix to the linear index in the half matrix representation'''
    positions = (i-1)*(i-1+3)/2+1+j

    positions = np.matrix(positions)

    return positions.getA1()