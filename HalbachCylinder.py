import numpy as np
import HalbachRing as HalbachRing
import CubeMagnet as CubeMagnet
import HalbachSlice as HalbachSlice
import matplotlib.pyplot as plt
from solid import *
from solid.utils import *

class HalbachCylinder:
    def __init__(self):
        self.slices = []

    def addSlice(self, slice):
        self.slices.append(slice)

    def calculateB(self, grid):
        B = np.zeros((grid[0].size, 3))
        for slice in self.slices:
            for ring in slice.rings:
                B += ring.calculateB(grid)
        return B
    
    def sliceAtPosition(self, position):
        for slice in self.slices:
            if slice.position == position:
                return slice
        return []

    def ringsAtPosition(self, position):
        returnRings = []
        for slice in self.slices:
            for ring in slice.rings:
                if ring.position == position:
                    returnRings.append(ring)
        return returnRings

    def generateSCADObject(self):
        for index, slice in enumerate(self.slices):
            if index == 0:
                cylinderSCADObject = slice.generateSCADObject()
            else:
                cylinderSCADObject += slice.generateSCADObject()
        return cylinderSCADObject

    def generateSCADFile(self, filename):        
        scad_render_to_file(self.generateSCADObject(), filename)

    def generateRingRadii(self, magnetSpacing, startRadius, radiusIncrement, numRadii):
        numMagnets = []
        for n in np.arange(numRadii):
            numMagnets.append()
        
if __name__ == '__main__':
    # same settings as in https://github.com/LUMC-LowFieldMRI/HalbachOptimisation/blob/master/homogeneityOptimisation.py
    innerRingRadii = np.array([148, 151, 154, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 201])*1e-3
    innerNumMagnets = np.array([50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67, 68])

    index = 0
    print('radius(numMagnets): ', end='')
    for r in np.nditer(innerRingRadii.T):
        print(str(r) + "(" + str(innerNumMagnets[index]) + ")", end=' ')
        index = index + 1
    print('')
    
    outerRingRadii = innerRingRadii + 21*1e-3
    outerNumMagnets = innerNumMagnets + 7
    
    numRings = 23
    ringSep = .022 
    magnetLength = (numRings - 1) * ringSep
    ringPositions = np.linspace(-magnetLength/2, magnetLength/2, numRings)   
    ringPositionsPositiveOnly = ringPositions[ringPositions >= 0]    

    resolution = 0.005
    dsv = 0.2
    simDimensions = (dsv, dsv, dsv)
    x = np.linspace(-simDimensions[0]/2, simDimensions[0]/2, int(simDimensions[0]/resolution)+1, dtype=np.float32)
    y = np.linspace(-simDimensions[1]/2, simDimensions[1]/2, int(simDimensions[1]/resolution)+1, dtype=np.float32)
    z = np.linspace(-simDimensions[2]/2, simDimensions[2]/2, int(simDimensions[2]/resolution)+1, dtype=np.float32)
    grid = np.meshgrid(x,y,z)
    mask = np.zeros(np.shape(grid[0]))
    mask[np.square(grid[0]) + np.square(grid[1]) + np.square(grid[2]) <= (dsv/2)**2] = 1   
    evalPoints = [g[mask==1] for g in grid]

    maxOuterRadius = np.max(outerRingRadii)
    standHeight = maxOuterRadius + 0.05
    standWidth = maxOuterRadius*1.6
    halbachCylinder = HalbachCylinder()
    for positionIndex, position in enumerate(ringPositionsPositiveOnly):
        sizeIndex = innerRingRadii.size - positionIndex
        if sizeIndex >= innerRingRadii.size:
            sizeIndex = innerRingRadii.size-1
        halbachSlice = HalbachSlice.HalbachSlice(position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, standHeight, standWidth)
        halbachSlice.addRing(HalbachRing.HalbachRing(position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachSlice.addRing(HalbachRing.HalbachRing(position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachCylinder.addSlice(halbachSlice)
        if positionIndex != 0:
            halbachSlice = HalbachSlice.HalbachSlice(-position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, standHeight, standWidth)
            halbachSlice.addRing(HalbachRing.HalbachRing(-position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachSlice.addRing(HalbachRing.HalbachRing(-position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachCylinder.addSlice(halbachSlice)
    B0 = halbachCylinder.calculateB(evalPoints)

    halbachCylinder.sliceAtPosition(0).generateSCADFile('slice.scad')
    halbachCylinder.generateSCADFile('cylinder.scad')

    print("Max B0 amplitude is " + str(np.amax(B0)) + " T")
    fig = plt.figure(figsize=(16,12))
    B0z0 = B0[evalPoints[2]==0,:]
    B_abs = np.linalg.norm(B0z0[:,0:1], axis=1)
    print("Homogeneity: " + str(((np.max(B_abs)-np.min(B_abs))/np.mean(B_abs))*1e6) + " ppm")
    evalPointsz0=(evalPoints[0][evalPoints[2]==0], evalPoints[1][evalPoints[2]==0], evalPoints[2][evalPoints[2]==0])
    qq = plt.quiver(evalPointsz0[0], evalPointsz0[1], B0z0[:,0], B0z0[:,1], B_abs, cmap=plt.cm.jet)
    plt.colorbar(qq)
    rings = halbachCylinder.ringsAtPosition(0)
    for ring in rings:
        ring.plotMagnets()
    plt.show(block=False)
    fig2 = plt.figure(figsize=(16,12))
    qq = plt.tricontour(evalPointsz0[0], evalPointsz0[1], B_abs, cmap=plt.cm.jet)
    plt.colorbar(qq)
    plt.show()    
    