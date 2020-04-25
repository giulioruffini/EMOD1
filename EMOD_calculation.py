#Imports original mesh (Use Freesurfer and reference code to begin with)

from mayavi import mlab
import EMODCalculationPipeline as EMOD
import os

if __name__ == '__main__':
    #Change your root folder and path as needed
    root_folder = r'C:\Users\rnsal\Documents\Documents\Neuroelectrics\EMOD calculation standalone'
    subject_id_usr = os.path.join(root_folder, 'subjects folder','template_subject')

    #Reduces mesh to (1-reductions_ratio_usr)% of initial size
    reductions_ratio_usr = 0.73

    #If 1, won't show the plots of the imported and simplified mesh
    quiet_usr = 0

    subject_obj = EMOD.subject(subject_id = subject_id_usr, reductions_ratio = reductions_ratio_usr, quiet = quiet_usr)

    #Simplifies mesh and calculates areas and normals.
    subject_obj.SimplifyMeshSubj()

    # For EMOD2 calculations you'll need to run the included mesh_curvature_standalone.py script.
    # You'll need to set-up a miniconda environment for that.

    #Calculates Euclidean distance matrix (the saved npy file is quite large, >1GB)
    subject_obj.EuclidDistCalcSubj()

    # No need for DM calculations, but may be useful for other metrics
    #subject_obj.GeodesicDistCalcSubj( divisions=120)

    # Reads calculated curvature (calculation is done separately)
    #subject_obj.curvatures()

    # Calculates EMOD. Check documentation (https://www.biorxiv.org/content/10.1101/688101v1.full.pdf) for the meaning of the different inputs.
    subject_obj.calcEMOD(lambda_sc=1, p0=0.5, sigma=0.4, l0=5, version='1a')

    #Plots desired quantity on the mesh
    subject_obj.plot_data(data_2_plot='area')

