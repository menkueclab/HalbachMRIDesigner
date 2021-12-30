import numpy as np
import HalbachRing as HalbachRing
import CubeMagnet as CubeMagnet
import HalbachSlice as HalbachSlice
import matplotlib.pyplot as plt
import json as json
import copy
import argparse
import gmsh
from solid import *
from solid.utils import *
from shutil import copyfile

class HalbachCylinder:
    def __init__(self):
        self.slices = []
        self.params = []
        gmsh.initialize()


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
    
    def loadJSON(self, filename):
        f = open(filename)
        params = json.load(f)
        self.params = params
        mirrorSlices = params['mirrorSlices']
        if 'standHeight' in params:
            standHeightGlobal = float(params['standHeight'])/1e3
        else:
            standHeightGlobal = 0
        if 'standWidth' in params:
            standWidthGlobal = float(params['standWidth'])/1e3
        else:
            standWidthGlobal = 0
        numConnectionRods = int(params['numConnectionRods'])
        connectionRodsDiameter = float(params['connectionRodsDiameter'])/1e3
        connectionRodsArcRadius = float(params['connectionRodsArcRadius'])/1e3
        rings = []
        for ring in params['rings']:
            rings.append(HalbachRing.HalbachRing(0, float(ring['radius'])/1e3, ring['numMagnets'], CubeMagnet))
        for index, slice in enumerate(params['slices']):
            if 'standHeight' in slice:
                standHeight = float(slice['standHeight'])/1e3
            else:
                standHeight = standHeightGlobal
            if 'standWidth' in slice:
                standWidth = float(slice['standWidth'])/1e3
            else:
                standWidth = standWidthGlobal
            halbachSlice = HalbachSlice.HalbachSlice(float(slice['position'])/1e3, float(slice['innerRadius'])/1e3, float(slice['outerRadius'])/1e3, 
                numConnectionRods, connectionRodsArcRadius, connectionRodsDiameter, standHeight, standWidth)
            for ring in slice['rings']:
                halbachSlice.addRing(copy.deepcopy(rings[ring['id']]), float(slice['position'])/1e3)
            halbachSlice.setParams(params)
            halbachSlice.setId(index)
            self.addSlice(halbachSlice)
            if mirrorSlices and float(slice['position']) != 0:
                halbachSlice = HalbachSlice.HalbachSlice(-float(slice['position'])/1e3, float(slice['innerRadius'])/1e3, float(slice['outerRadius'])/1e3, 
                    numConnectionRods, connectionRodsArcRadius, connectionRodsDiameter, standHeight, standWidth)
                for ring in slice['rings']:
                    halbachSlice.addRing(copy.deepcopy(rings[ring['id']]), -float(slice['position'])/1e3)
                halbachSlice.setParams(params)
                halbachSlice.setId(index)
                self.addSlice(halbachSlice)
        f.close()
def generateExampleGeometry():
    innerRingRadii = np.array([148, 151, 154, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 201])*1e-3
    innerNumMagnets = np.array([50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67, 68])    
    outerRingRadii = innerRingRadii + 21*1e-3
    outerNumMagnets = innerNumMagnets + 7
    
    numRings = 23
    ringSep = .022 
    magnetLength = (numRings - 1) * ringSep
    ringPositions = np.linspace(-magnetLength/2, magnetLength/2, numRings)   
    ringPositionsPositiveOnly = ringPositions[ringPositions >= 0]    
    maxOuterRadius = np.max(outerRingRadii)
    standHeight = maxOuterRadius + 0.05
    standWidth = maxOuterRadius*1.6
    halbachCylinder = HalbachCylinder()
    for positionIndex, position in enumerate(ringPositionsPositiveOnly):
        sizeIndex = innerRingRadii.size - positionIndex
        if sizeIndex >= innerRingRadii.size:
            sizeIndex = innerRingRadii.size-1
        halbachSlice = HalbachSlice.HalbachSlice(position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, 12, maxOuterRadius + 0.025, 0.005, standHeight, standWidth)
        halbachSlice.addRing(HalbachRing.HalbachRing(position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachSlice.addRing(HalbachRing.HalbachRing(position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachCylinder.addSlice(halbachSlice)
        if positionIndex != 0:
            halbachSlice = HalbachSlice.HalbachSlice(-position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, 12, maxOuterRadius + 0.025, 0.005, standHeight, standWidth)
            halbachSlice.addRing(HalbachRing.HalbachRing(-position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachSlice.addRing(HalbachRing.HalbachRing(-position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachCylinder.addSlice(halbachSlice)
    return halbachCylinder

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create Halbach-MRI geometry')   
    parser.add_argument('filename', nargs='+', help='a .json file that describes the geometry')
    parser.add_argument('--contour', action='store_true', help='creates a contour plot at z=0')
    parser.add_argument('--quiver', action='store_true', help='creates a quiver plot at z=0')
    parser.add_argument('--fem', action='store_true', help='generate a .geo and .pro file for simulation with GMSH=GetDP')
    parser.add_argument('-o', nargs='?', default='out.scad', help='name of output scad file')
    args = parser.parse_args()     
    halbachCylinder = HalbachCylinder()

    print("loading file...")
    halbachCylinder.loadJSON(args.filename[0])
    print(args.filename[0] + " loaded")

    print("writing " + args.o)
    halbachCylinder.generateSCADFile(args.o)

    # alternative to json file
    # halbachCylinder = generateExampleGeometry()

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
    B0 = halbachCylinder.calculateB(evalPoints)
    print("Max B0 amplitude is " + str(np.amax(B0)) + " T")
    B0z0 = B0[evalPoints[2]==0,:]
    B_abs = np.linalg.norm(B0z0[:,0:1], axis=1)
    print("Homogeneity: " + str(((np.max(B_abs)-np.min(B_abs))/np.mean(B_abs))*1e6) + " ppm")
    evalPointsz0=(evalPoints[0][evalPoints[2]==0], evalPoints[1][evalPoints[2]==0], evalPoints[2][evalPoints[2]==0])

    if args.quiver:
        fig = plt.figure(figsize=(16,12))
        qq = plt.quiver(evalPointsz0[0], evalPointsz0[1], B0z0[:,0], B0z0[:,1], B_abs, cmap=plt.cm.jet)
        plt.colorbar(qq)
        rings = halbachCylinder.ringsAtPosition(0)
        for ring in rings:
            ring.plotMagnets()
        plt.show(block=False)
    if args.contour:
        fig2 = plt.figure(figsize=(16,12))
        qq = plt.tricontour(evalPointsz0[0], evalPointsz0[1], B_abs, cmap=plt.cm.jet)
        plt.colorbar(qq)
    if args.fem:
        gmsh.initialize()        
        gmsh.model.add("cylinder")
        meshResolution = 0.024
        boxDimensions = (.300, .300, .300)
        gmsh.model.occ.synchronize()    
        gmsh.option.setNumber("Mesh.Optimize",1)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        gmsh.option.setNumber("Solver.AutoMesh",2)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        numMagnets = 0
        magnetData = "DefineConstant[\n"
        for numSlice, slice in enumerate(halbachCylinder.slices):
            print("Slice " + str(numSlice))
            for numRing, ring in enumerate(slice.rings):
                numMagnetsProcessed, magnetAngles = ring.generateGeometry(numMagnets)
                numMagnets += numMagnetsProcessed
                magnetData += magnetAngles
                print("   Ring " + str(numRing))
        gmsh.model.occ.synchronize()
        magnetData += "NumMagnets = " + str(numMagnets) + "\n"
        magnetData += "SurfaceOffset = 10000\n"
        magnetData += "];"
        with open("cylinder_magnets_data.pro", "w") as text_file:
            text_file.write(magnetData)            

        # add bounding box 
        #airVol, airSL = addBox(*tuple(x*(-1) for x in boxDimensions), *tuple(x*2 for x in boxDimensions))             
        airVol = gmsh.model.occ.addBox(*tuple(x*(-1) for x in boxDimensions), *tuple(x*2 for x in boxDimensions))             
        gmsh.model.occ.synchronize()
        airVolBoundary = [x[1] for x in gmsh.model.getBoundary([[3,airVol]], oriented=False)]
        gmsh.model.occ.fragment(gmsh.model.occ.getEntities(3), [])
        gmsh.model.occ.synchronize()
        physicalTag = gmsh.model.addPhysicalGroup(3, [airVol], numMagnets+1) 
        physicalTag = gmsh.model.addPhysicalGroup(2, airVolBoundary, numMagnets+2)       
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)    
        gmsh.write("cylinder.geo_unrolled")
        copyfile("cylinder.geo_unrolled", "cylinder.geo") # opening the .pro file in gmsh GUI searches for a .geo file
        gmsh.write("cylinder.msh")
        #gmsh.fltk.run()
        gmsh.finalize()
    plt.show()    
    