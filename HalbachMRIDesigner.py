import numpy as np
import HalbachRing
import CubeMagnet
import HalbachSlice
import HalbachCylinder
import matplotlib.pyplot as plt
import json as json
import copy
import argparse
import gmsh
from solid import *
from solid.utils import *
from shutil import copyfile
import os

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
        halbachSlice = HalbachSlice(position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, 12, maxOuterRadius + 0.025, 0.005, standHeight, standWidth)
        halbachSlice.addRing(HalbachRing(position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachSlice.addRing(HalbachRing(position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), position)
        halbachCylinder.addSlice(halbachSlice)
        if positionIndex != 0:
            halbachSlice = HalbachSlice(-position, innerRingRadii[sizeIndex]-0.020, maxOuterRadius+0.040, 12, maxOuterRadius + 0.025, 0.005, standHeight, standWidth)
            halbachSlice.addRing(HalbachRing(-position, innerRingRadii[sizeIndex], innerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachSlice.addRing(HalbachRing(-position, outerRingRadii[sizeIndex], outerNumMagnets[sizeIndex], CubeMagnet), -position)
            halbachCylinder.addSlice(halbachSlice)
    return halbachCylinder

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create Halbach-MRI geometry')   
    parser.add_argument('filename', nargs='+', help='a .json file that describes the geometry')
    parser.add_argument('--contour', action='store_true', help='creates a contour plot at z=0')
    parser.add_argument('--quiver', action='store_true', help='creates a quiver plot at z=0')
    parser.add_argument('--fem', action='store_true', help='generate a .geo and .pro file for simulation with GMSH=GetDP')
    parser.add_argument('-o', nargs='?', default='input filename without extension', help='output filename without extension')
    args = parser.parse_args()     
    halbachCylinder = HalbachCylinder.HalbachCylinder()

    print("loading file...")
    halbachCylinder.loadJSON(args.filename[0])
    print(args.filename[0] + " loaded")

    outputFilename = args.o
    if args.o == "input filename without extension":
        outputFilename = os.path.splitext(os.path.basename(args.filename[0]))[0]
    print("writing " + outputFilename + ".scad")
    halbachCylinder.generateSCADFile(outputFilename + ".scad")

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
        BoundingBoxDiameter = 0.3        
        DSV = 0.2
        boxDimensions = (BoundingBoxDiameter, BoundingBoxDiameter, BoundingBoxDiameter)
        gmsh.model.occ.synchronize()    
        gmsh.option.setNumber("Mesh.Optimize",1)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        gmsh.option.setNumber("Solver.AutoMesh",2)
        gmsh.option.setNumber("Geometry.ExactExtrusion",0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.003)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFactor", 1)        
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
        gmsh.model.occ.synchronize()
        physicalTag = gmsh.model.addPhysicalGroup(3, [airVol], numMagnets+1) 
        physicalTag = gmsh.model.addPhysicalGroup(2, airVolBoundary, numMagnets+2)       
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)    

        with open(outputFilename + ".pro", "w") as text_file:
            text_file.write(magnetData)       
            text_file.write("Include \"templates/cylinder_template.pro\"\n")       
        gmsh.write(outputFilename + ".geo_unrolled")
        copyfile(outputFilename + ".geo_unrolled", outputFilename + ".geo") # opening the .pro file in gmsh GUI searches for a .geo file
        os.remove(outputFilename + ".geo_unrolled")
        gmsh.write(outputFilename + ".msh")
        gmsh.write(outputFilename + ".geo.opt")
        #gmsh.fltk.run()
        gmsh.finalize()
    plt.show()    
    