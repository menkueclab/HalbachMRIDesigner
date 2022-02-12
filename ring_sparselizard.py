import pickle
import numpy as np
import matplotlib as plt
from spylizard import *

mymesh = mesh()
mymesh.load("ring.msh")

settingsDict = []
with open("ring.py", "rb") as python_file:
    settingsDict = pickle.load(python_file)
numMagnets = len(settingsDict['dataDict'])

pi = getpi()
mu0 = 4*pi*1e-7
mu = parameter()
br = parameter(3,1)
magnets = []
magnetSurfaces = []
for dataItem in settingsDict['dataDict']:
    mu.setvalue(dataItem['id'], mu0*dataItem['mur'])
    temp = array3x1(dataItem['BR'], 0, 0)
    temp.rotate(0, 0, dataItem['angle']/pi*180)
    br.setvalue(dataItem['id'], temp)
    magnets.append(dataItem['id'])
    magnetSurfaces.append(dataItem['id'] + settingsDict['SurfaceRegionOffset'])
air = numMagnets + 1
inf = numMagnets + 2
dsv = numMagnets + 3
#dsvBoundary = numMagnets + 4
wholedomain = selectunion([*magnets, air, inf, *magnetSurfaces])
mu.setvalue(air, mu0)
br.setvalue(air, array3x1(0, 0, 0))

phi = field("h1")
phi.setorder(wholedomain, 1)

# The magnetic scalar potential (just like the electric potential)
# needs to be fixed at least at one node to have a well-posed problem.
# Here it is forced to 0 at boundary of air box.
phi.setconstraint(inf)

magnetostatics = formulation()
magnetostatics += integral(wholedomain, -grad(dof(phi)) * mu * grad(tf(phi)) )
for magnet in magnets:
    magnetostatics += integral(magnet, -br * grad(tf(phi)))

magnetostatics.generate()

sol = solve(magnetostatics.A(), magnetostatics.b())
phi.setdata(wholedomain, sol)
h = (-grad(phi))
b = mu*h - br
h.write(wholedomain, "ring_sparselizard_h.pos", 1)
b.write(wholedomain, "ring_sparselizard_b.pos", 1)

# generate dsv sphere
resolution = 0.005
dsv = 0.2
simDimensions = (dsv, dsv, dsv)

# sphere
x = np.linspace(-simDimensions[0]/2, simDimensions[0]/2, int(simDimensions[0]/resolution)+1, dtype=np.float32)
y = np.linspace(-simDimensions[1]/2, simDimensions[1]/2, int(simDimensions[1]/resolution)+1, dtype=np.float32)
z = np.linspace(-simDimensions[2]/2, simDimensions[2]/2, int(simDimensions[2]/resolution)+1, dtype=np.float32)
grid = np.meshgrid(x,y,z)
mask = np.zeros(np.shape(grid[0]))
mask[np.square(grid[0]) + np.square(grid[1]) + np.square(grid[2]) <= (dsv/2)**2] = 1   
evalPoints = [g[mask==1] for g in grid]    
interpolated = []
isFound = []
b.interpolate(wholedomain, list(np.array(evalPoints).reshape(len(evalPoints[0])*3, order='F')), interpolated, isFound)
interp3 = np.array(interpolated).reshape(3,len(evalPoints[0]), order='F')
scatterwrite("interpolated_b_sphere.pos", evalPoints[0], evalPoints[1], evalPoints[2], interp3[0], interp3[1], interp3[2])

# slice at z=0
x = np.linspace(-simDimensions[0]/2, simDimensions[0]/2, int(simDimensions[0]/resolution)+1, dtype=np.float32)
y = np.linspace(-simDimensions[1]/2, simDimensions[1]/2, int(simDimensions[1]/resolution)+1, dtype=np.float32)
z = [0]
grid = np.meshgrid(x,y,z)
mask = np.zeros(np.shape(grid[0]))
mask[np.square(grid[0]) + np.square(grid[1]) <= (dsv/2)**2] = 1   
evalPoints = [g[mask==1] for g in grid]    
interpolated = []
isFound = []
b.interpolate(wholedomain, list(np.array(evalPoints).reshape(len(evalPoints[0])*3, order='F')), interpolated, isFound)
interp3 = np.array(interpolated).reshape(3,len(evalPoints[0]), order='F')
scatterwrite("interpolated_b_z0.pos", evalPoints[0], evalPoints[1], evalPoints[2], interp3[0], interp3[1], interp3[2])
interpolated = []
isFound = []
h.interpolate(wholedomain, list(np.array(evalPoints).reshape(len(evalPoints[0])*3, order='F')), interpolated, isFound)
interp3 = np.array(interpolated).reshape(3,len(evalPoints[0]), order='F')
scatterwrite("interpolated_h_z0.pos", evalPoints[0], evalPoints[1], evalPoints[2], interp3[0], interp3[1], interp3[2])

