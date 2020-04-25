# EMOD calculation (Python 3)

Calculation of EMOD index for a head model. The only required inputs are two lists: a list of vertices 
(n_vertices x 3 txt file with the x, y and z coordinates of each vertex of the mesh) and a list of faces (m_faces x 3, indexing the positions in the vertices list of each point comprising the triangle).

  - Simplifies the mesh
  - Calculates areas of triangles, normals to the surface (pointing outwards) and euclidean and geodesic distance matrices
  - Calculates EMOD on the simplified mesh
  - Visualizes distribution of certain quantities on the mesh

# Getting started
This section will explain how to configure the software to use it locally in your computer.

## Pre-requisites (Python 3):
You will need to install the following python toolboxes:

 
 - pymeshfix (v0.13.14)
 - mayavi (v4.7.1)
 - numpy (v1.18.1)
 - tvb-gdist (v2.0.0)
 - trimesh (v3.5.25)
 - PyQt5 (v5.11.2) (for Mayavi)
 - progressbar2 (v3.5.0)
 
 For curvature calculation (stand-alone script) you will need:
 - igl (installed via conda, v0.4.1)
 - numpy (v1.18.1)
 
The generation of the original cortical surface is not done within these scripts. Surfaces generated with typical segmentation sofware (Freesurfer, Simnibs) are adequate provided they are converted to text files.
Place the list of vertices/triangles inside a subject folder and rename them to: nodes_GM_full/faces_GM_full.txt.

# Running the code
 - the main script is EMOD_calculation.py
 - run this script in Python 3 from the EMOD1 directory (working directory) 
 - Close the Mayavi windows after each step to proceed

# Todos
 - Implement other versions of EMOD
 - Implement options for batch calculation of EMODs
 - Implement original cortical surface generation within the script

# Known-bugs
 - Calculations require a substantial amount of RAM memory (around 32 GB for the mesh sizes used by default) and may fail in some computers. 
 - Let us know

# License
----

This software belongs to Neuroelectrics.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
