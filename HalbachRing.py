import numpy as np
from CubeMagnet import CubeMagnet
from Magnet import Magnet
import matplotlib.pyplot as plt
import gmsh
from shutil import copyfile
import os
import pickle

def addBox(x, y, z, dx, dy, dz, meshSize = 0):
    # needed becaues occ.addBox does not return the surface loop tag
    p1 = gmsh.model.occ.addPoint(x, y, z, meshSize)
    p2 = gmsh.model.occ.addPoint(x+dx, y, z, meshSize)
    p3 = gmsh.model.occ.addPoint(x+dx, y+dy, z, meshSize)
    p4 = gmsh.model.occ.addPoint(x, y+dy, z, meshSize)
    p5 = gmsh.model.occ.addPoint(x, y, z+dz, meshSize)
    p6 = gmsh.model.occ.addPoint(x+dx, y, z+dz, meshSize)
    p7 = gmsh.model.occ.addPoint(x+dx, y+dy, z+dz, meshSize)
    p8 = gmsh.model.occ.addPoint(x, y+dy, z+dz, meshSize)
    l12 = gmsh.model.occ.addLine(p1, p2)
    l23 = gmsh.model.occ.addLine(p2, p3)
    l34 = gmsh.model.occ.addLine(p3, p4)
    l14 = gmsh.model.occ.addLine(p1, p4)
    l56 = gmsh.model.occ.addLine(p5, p6)
    l67 = gmsh.model.occ.addLine(p6, p7)
    l78 = gmsh.model.occ.addLine(p7, p8)
    l58 = gmsh.model.occ.addLine(p5, p8)
    l48 = gmsh.model.occ.addLine(p4, p8)
    l37 = gmsh.model.occ.addLine(p3, p7)
    l26 = gmsh.model.occ.addLine(p2, p6)
    l56 = gmsh.model.occ.addLine(p5, p6)
    l15 = gmsh.model.occ.addLine(p1, p5)
    l26 = gmsh.model.occ.addLine(p2, p6)
    l58 = gmsh.model.occ.addLine(p5, p8)
    l37 = gmsh.model.occ.addLine(p3, p7)
    front = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l12, l23, l34, -l14])])
    back = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l58, -l78, -l67, -l56])])
    top = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l34, l48, -l78, -l37])])
    bottom = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l12, l26, -l56, -l15])])
    left = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l15, l58, -l48, -l14])])
    right = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l23, l37, -l67, -l26])])
    slTag = gmsh.model.occ.addSurfaceLoop([front, back, top, bottom, left, right])
    volTag = gmsh.model.occ.addVolume([slTag])
    return volTag, slTag

class HalbachRing:
    # all units are SI units
    def __init__(self, ringId, position, radius, numMagnets, magnetType):
        self.ringId = ringId
        self.radius = radius
        self.numMagnets = numMagnets
        self.position = position

        magnetAngles = np.linspace(0, 2*np.pi, self.numMagnets, endpoint=False)
        self.magnets = []
        k = 2
        for angle in magnetAngles:
            self.magnets.append(Magnet((radius*np.cos(angle), radius*np.sin(angle), 0), angle*k, magnetType))

    def setPosition(self, position):
        self.position = position

    def calculateB(self, evalPoints):
        x, y, z = evalPoints
        B0 = np.zeros( (x.size, 3) )
        for magnet in self.magnets:
            magnetAbsoluteZPosition = magnet.position[2] + self.position
            magneticMoment_2 = np.multiply(1/(4*np.pi), magnet.magneticMoment_withoutMu0())
            r = [(x+magnet.position[0]), (y+magnet.position[1]), (z+magnetAbsoluteZPosition)]
            r_dot_magneticMoment_2 = 3*((r[0])*magneticMoment_2[0] + (r[1])*magneticMoment_2[1])
            # magnetic moment has no z-component
            r_abs_p2 = np.square(r[0]) + np.square(r[1]) + np.square(r[2])
            r_abs_p3 = np.power(r_abs_p2,1.5)
            r_abs_p5 = np.power(r_abs_p2,2.5)
            B0[:,0] += np.divide(np.multiply(r[0], r_dot_magneticMoment_2),r_abs_p5) - np.divide(magneticMoment_2[0],r_abs_p3)
            B0[:,1] += np.divide(np.multiply(r[1], r_dot_magneticMoment_2),r_abs_p5) - np.divide(magneticMoment_2[1],r_abs_p3)
            B0[:,2] += np.divide(np.multiply(r[2], r_dot_magneticMoment_2),r_abs_p5)
        return B0

    def plotMagnets(self):
        for magnet in self.magnets:
            magnet.draw2DMatplotlib()
        plt.gca().add_patch(plt.Circle((0,0), self.radius, fill=False, color='green'))
        #circle=plt.Circle((magnet.position[0],magnet.position[1]),0.001,color='green')
        #plt.gca().add_patch(circle)

    def getDict(self, startIndex):
        dataDict = []
        for index, magnet in enumerate(self.magnets, start=1):
            dictItem = {
                "angle" : magnet.angle,
                "BR"  : magnet.magnetType.BR, 
                "mur" : magnet.magnetType.mur,
                "id" : index + startIndex
                }
            dataDict.append(dictItem)
        return dataDict

    def generateGeometry(self, startIndex=0):        
        meshResolution = 0.006
        magnetData = ""
        for index, magnet in enumerate(self.magnets, start=1):
            boxPos = np.array(magnet.position) - magnet.magnetType.dimension/2 + (0,0,self.position)
            magnetTag = gmsh.model.occ.addBox(boxPos[0], boxPos[1], boxPos[2], magnet.magnetType.dimension, magnet.magnetType.dimension, magnet.magnetType.dimension)
            gmsh.model.occ.rotate([[3, magnetTag]], magnet.position[0], magnet.position[1], magnet.position[2] + self.position, 0, 0, 1, magnet.angle)
            gmsh.model.occ.synchronize()            
            entities = gmsh.model.getBoundary([[3, magnetTag]], oriented=False)
            physicalTag = gmsh.model.addPhysicalGroup(2, [x[1] for x in entities], index + startIndex + 10000)
            physicalTag = gmsh.model.addPhysicalGroup(3, [magnetTag], index + startIndex)
            gmsh.model.occ.mesh.setSize(gmsh.model.getBoundary([[3, magnetTag]], oriented = False, recursive = True) ,meshResolution)
            magnetData += "angle_" + str(index + startIndex) + " = " + str(magnet.angle) + "\n"
            magnetData += "BR_" + str(index + startIndex) + " = " + str(magnet.magnetType.BR) + "\n"
            magnetData += "mur_" + str(index + startIndex) + " = " + str(magnet.magnetType.mur) + "\n"
        return len(self.magnets), magnetData


def threadedGmshGUI():
    gmsh.fltk.run()


if __name__ == '__main__':
    '''
    testcode
    '''
    gmsh.initialize()
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
    
    cubeMagnet = CubeMagnet(0.012, 1.3, 1)
    testRing1 = HalbachRing(0, 0.195, 66, cubeMagnet)
    testRing2 = HalbachRing(0, 0.195+0.021, 66+7, cubeMagnet)
    B0 = testRing1.calculateB(evalPoints)
    B0 += testRing2.calculateB(evalPoints)
    
    print("Max B0 amplitude is " + str(np.amax(B0)) + " T")
    B0z0 = B0[evalPoints[2]==0,:]
    B_abs = np.linalg.norm(B0z0[:,0:1], axis=1)
    print("Homogeneity: " + str(((np.max(B_abs)-np.min(B_abs))/np.mean(B_abs))*1e6) + " ppm")
    if False:
        fig = plt.figure(figsize=(16,12))
        evalPointsz0=(evalPoints[0][evalPoints[2]==0], evalPoints[1][evalPoints[2]==0], evalPoints[2][evalPoints[2]==0])
        qq = plt.quiver(evalPointsz0[0], evalPointsz0[1], B0z0[:,0], B0z0[:,1], B_abs, cmap=plt.cm.jet)
        plt.colorbar(qq)
        rings = [testRing1, testRing2]
        for ring in rings:
            ring.plotMagnets()
        plt.show(block=False)
        fig2 = plt.figure(figsize=(16,12))
        qq = plt.tricontour(evalPointsz0[0], evalPointsz0[1], B_abs, cmap=plt.cm.jet)
        plt.colorbar(qq)
        plt.show()    
    if True:
        gmsh.model.add("ring")
        meshResolution = 0.05
        meshResolutionDSV = 0.006
        meshResolutionDSVMagnet = 0.05
        meshResolutionMagnet = 0.003
        DSV = 0.2
        BoundingBoxDiameter = 0.25
        boxDimensions = (BoundingBoxDiameter, BoundingBoxDiameter, BoundingBoxDiameter)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2)  # for sparselizard
        gmsh.option.setNumber("Geometry.ExactExtrusion", 0)
        gmsh.option.setNumber("Solver.AutoMesh", 2)
        gmsh.option.setNumber("Geometry.ExactExtrusion", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFactor", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.003)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1)
        index, magnetData1 = testRing1.generateGeometry()
        index, magnetData2 = testRing2.generateGeometry(index)
        gmsh.model.occ.synchronize()
        numMagnets = len(testRing1.magnets)+len(testRing2.magnets)
        magnetData = "DefineConstant[\n"
        magnetData += magnetData1
        magnetData += magnetData2
        magnetData += "NumMagnets = " + str(numMagnets) + "\n"
        magnetData += "SurfaceRegionOffset = 10000\n"
        magnetData += "DSV = " + str(DSV) + "\n"
        magnetData += "outputFilename = " + "\"ring\"" + "\n" 
        magnetData += "];\n"
           

        # add bounding box
        #airVol, airSL = addBox(*tuple(x*(-1) for x in boxDimensions), *tuple(x*2 for x in boxDimensions))             
        airVol = gmsh.model.occ.addBox(*tuple(x*(-1) for x in boxDimensions), *tuple(x*2 for x in boxDimensions))   
        gmsh.model.occ.synchronize()
        airVolBoundary = [x[1] for x in gmsh.model.getBoundary([[3,airVol]], oriented=False)]
        gmsh.model.occ.fragment(gmsh.model.occ.getEntities(3), [])
        # add after fragment() so that dsv is still included in the airVol
        dsvVol = gmsh.model.occ.addSphere(0, 0, 0, dsv/2)   
        gmsh.model.occ.synchronize()
        dsvVolBoundary = [x[1] for x in gmsh.model.getBoundary([[3,dsvVol]], oriented=False)]
        gmsh.model.occ.synchronize()
        physicalTag = gmsh.model.addPhysicalGroup(3, [airVol, dsvVol], numMagnets+1) 
        physicalTag = gmsh.model.addPhysicalGroup(2, airVolBoundary, numMagnets+2)        
        physicalTag = gmsh.model.addPhysicalGroup(3, [dsvVol], numMagnets+3) 
        physicalTag = gmsh.model.addPhysicalGroup(2, dsvVolBoundary, numMagnets+4) 
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.field.add("Ball", 1)
        gmsh.model.mesh.field.setNumber(1, "Radius", dsv/2)
        gmsh.model.mesh.field.setNumber(1, "Thickness", 0.02)
        gmsh.model.mesh.field.setNumber(1, "VIn", meshResolutionDSV)
        gmsh.model.mesh.field.setNumber(1, "VOut", meshResolutionDSVMagnet)
        gmsh.model.mesh.field.setAsBackgroundMesh(1)
        gmsh.model.mesh.generate(3)
        # with open("ring.pro", "w") as text_file:
        #     text_file.write(magnetData)    
        #     text_file.write("Include \"templates/ring_template.pro\"\n")
        dataDict = testRing1.getDict(0)
        dataDict += testRing2.getDict(len(testRing1.magnets))
        settingsDict = {
                "dataDict" : dataDict,
                "SurfaceRegionOffset" : 10000
            }
        with open("ring.pickle", "wb") as python_file:
            pickle.dump(settingsDict, python_file)
        gmsh.write("ring.geo_unrolled")
        copyfile("ring.geo_unrolled", "ring.geo") # opening the .pro file in gmsh GUI searches for a .geo file
        os.remove("ring.geo_unrolled")
        gmsh.write("ring.msh")
        gmsh.write("ring.geo.opt")
        # gmsh.merge("ring.pro")
        #gmsh.onelab.run()
        #gmsh.fltk.run()
        gmsh.finalize()

