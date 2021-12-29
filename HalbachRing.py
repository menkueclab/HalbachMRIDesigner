import numpy as np
import CubeMagnet as CubeMagnet
import matplotlib.pyplot as plt
import gmsh
import math
import pygetdp

# all units are SI units
# TODO: make magnet type a template paramter

class HalbachRing:
    def __init__(self, position, radius, numMagnets, TMagnet):
        self.radius = radius
        self.numMagnets = numMagnets
        self.position = position

        magnetAngles = np.linspace(0, 2*np.pi, self.numMagnets, endpoint=False)
        self.magnets = []
        k = 2
        for angle in magnetAngles:
            self.magnets.append(TMagnet.Magnet((radius*np.cos(angle), radius*np.sin(angle), 0), angle*k))

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

    def generateGeometry(self):        
        magnetTags = []   
        meshResolution = 0.003
        magnetData = "DefineConstant[\n"
        magnetData += "NumMagnets = " + str(len(self.magnets)) + ",\n"
        for index, magnet in enumerate(self.magnets):
            index += 1
            boxPos = np.array(magnet.position) - magnet.dimension/2
            magnetTag = gmsh.model.occ.addBox(boxPos[0], boxPos[1], boxPos[2], magnet.dimension, magnet.dimension, magnet.dimension)
            magnetTags.append(magnetTag+1)
            gmsh.model.occ.rotate([[3, magnetTag]], magnet.position[0], magnet.position[1], magnet.position[2], 0, 0, 1, magnet.angle)
            gmsh.model.occ.synchronize()            
            entities = gmsh.model.occ.getEntitiesInBoundingBox(*gmsh.model.occ.getBoundingBox(3,magnetTag), dim=2)
            physicalTag = gmsh.model.addPhysicalGroup(2, [x[1] for x in entities], index + 1000)
            #gmsh.model.setPhysicalName(2, physicalTag, "MagnetSurf" + str(index))
            gmsh.model.occ.synchronize()              
            physicalTag = gmsh.model.addPhysicalGroup(3, [magnetTag], index)
            #gmsh.model.setPhysicalName(3, physicalTag, "MagnetVol" + str(index))
            magnetData += "angle_" + str(index) + " = " + str(magnet.angle) + "\n"
            entities = gmsh.model.occ.getEntitiesInBoundingBox(*gmsh.model.occ.getBoundingBox(3,magnetTag), dim=0)            
            gmsh.model.occ.mesh.setSize(entities,meshResolution)
        magnetData += "];"
        with open("magnets_data.pro", "w") as text_file:
            text_file.write(magnetData)        
        return magnetTags


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
    
    testRing1 = HalbachRing(0, 0.195, 66, CubeMagnet)
    testRing2 = HalbachRing(0, 0.195+0.021, 66+7, CubeMagnet)
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
        meshResolution = 0.024
        numMagnets = len(testRing1.magnets)
        boxDimensions = (.250, .250, .250)
        boxTag = gmsh.model.occ.addBox(*tuple(x*(-1) for x in boxDimensions), *tuple(x*2 for x in boxDimensions))     
        entities = gmsh.model.occ.getEntitiesInBoundingBox(*gmsh.model.occ.getBoundingBox(3,boxTag), dim=2)   
        gmsh.model.occ.mesh.setSize(entities,meshResolution)        
        gmsh.model.occ.synchronize()    
        # remove bounding box volume to prevent overlapping volumes
        gmsh.model.occ.remove(gmsh.model.getEntities(3))
        # but keep bounding box surface Loop
        gmsh.model.occ.synchronize()    
        gmsh.option.setNumber("Mesh.Optimize",1)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        gmsh.option.setNumber("Solver.AutoMesh",2)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        surfaceLoops = testRing1.generateGeometry()
        gmsh.model.occ.synchronize()
        # add bounding box 
        airVol = gmsh.model.occ.addVolume([boxTag] + surfaceLoops)
        gmsh.model.occ.synchronize()
        physicalTag = gmsh.model.addPhysicalGroup(3, [airVol], numMagnets+1) 
        physicalTag = gmsh.model.addPhysicalGroup(2, [boxTag], numMagnets+2)        
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)    
        gmsh.write("ring.geo_unrolled")
        gmsh.write("ring.msh")
        gmsh.finalize()

