# HalbachMRIDesigner
This tool creates an OpenSCAD geometry of a halbach cylinder which generates a homogeneous magnetic field inside. The design is fully parameterized so that it can easily be adapted. The parameters of the magnet are encapsulated in a json file. In addition to parameters, custom SCAD code can be specified in the json file, so that structures which are not parameterized, ie coil mounting planes, can be realized. An example can be found in examples/mri1.json.

To frabricate the individual slices by milling or 3d printing, you can import the OpenSCAD model into FreeCAD and then export a separate file for every slice.

# Running the script
- pip install -U pip
- pip install -r requirements.txt
- run 'python HalbachMRIDesigner.py examples/mri1.conf --contour --quiver' with python 3.8 or higher
- the generated .scad file should look like examples/cylinder.scad
- in OpenSCAD it should look like this

![OpenSCAD screenshot](https://github.com/menkueclab/HalbachMRIDesigner/blob/master/examples/cylinder.png?raw=true)

# Converting OpenSCAD files to FreeCAD
- install OpenSCAD (https://openscad.org/)
- change Edit->Preferences->Advanced
  - Turn off rendering at 1000000
  - CGAL Cache size 1000
  - PolySet Cache size 1000
- open cylinder.scad
- File->export->CSG
- Import .csg file into FreeCAD

# Contributions
- Benjamin Menkuec (code and ideas)
- Lukas Winter (ideas)
- Thomas O'Reilly (ideas, (https://github.com/LUMC-LowFieldMRI/HalbachOptimisation))

# License
GNU Version 3
