import numpy as np
from HalbachRing import HalbachRing
from CubeMagnet import CubeMagnet
from HalbachSlice import HalbachSlice
import json
import copy
import gmsh
import solid
from shutil import copyfile
import os

class HalbachCylinder:
    def __init__(self):
        self.slices = []
        self.params = []
        gmsh.initialize()


    def addSlice(self, sliceObj):
        self.slices.append(sliceObj)

    def calculateB(self, grid):
        B = np.zeros((grid[0].size, 3))
        for sliceObj in self.slices:
            for ring in sliceObj.rings:
                B += ring.calculateB(grid)
        return B

    def sliceAtPosition(self, position):
        for sliceObj in self.slices:
            if sliceObj.position == position:
                return sliceObj
        return []

    def ringsAtPosition(self, position):
        returnRings = []
        for sliceObj in self.slices:
            for ring in sliceObj.rings:
                if ring.position == position:
                    returnRings.append(ring)
        return returnRings

    def generateSCADObject(self):
        for index, sliceObj in enumerate(self.slices):
            if index == 0:
                cylinderSCADObject = sliceObj.generateSCADObject()
            else:
                cylinderSCADObject += sliceObj.generateSCADObject()
        return cylinderSCADObject

    def generateSCADFile(self, filename):        
        solid.scad_render_to_file(self.generateSCADObject(), filename)

    def generateRingRadii(self, magnetSpacing, startRadius, radiusIncrement, numRadii):
        numMagnets = []
        for _ in np.arange(numRadii):
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
        magnetTypes = []
        for magnet in params['magnets']:
            if magnet['shape'] == "cube":
                dimension = float(magnet['dimension'])*1e-3
                BR = float(magnet['BR'])
                mur = float(magnet['mur'])
                cubeMagnet = CubeMagnet(dimension, BR, mur)
                magnetTypes.append(cubeMagnet)
            else:
                print("Magnet shape \"" + magnet['shape'] + "\" is not supported!")
        defaultMagnetType = params['defaultMagnetType']
        rings = []
        for ring in params['rings']:
            rings.append(HalbachRing(0, float(ring['radius'])/1e3, ring['numMagnets'], magnetTypes[defaultMagnetType]))
        for index, slice in enumerate(params['slices']):
            if 'standHeight' in slice:
                standHeight = float(slice['standHeight'])/1e3
            else:
                standHeight = standHeightGlobal
            if 'standWidth' in slice:
                standWidth = float(slice['standWidth'])/1e3
            else:
                standWidth = standWidthGlobal
            halbachSlice = HalbachSlice(float(slice['position'])/1e3, float(slice['innerRadius'])/1e3,
                float(slice['outerRadius'])/1e3, numConnectionRods, connectionRodsArcRadius, connectionRodsDiameter,
                standHeight, standWidth)
            for ring in slice['rings']:
                halbachSlice.addRing(copy.deepcopy(rings[ring['id']]), float(slice['position'])/1e3)
            halbachSlice.setParams(params)
            halbachSlice.setId(index)
            self.addSlice(halbachSlice)
            if mirrorSlices and float(slice['position']) != 0:
                halbachSlice = HalbachSlice(-float(slice['position'])/1e3, float(slice['innerRadius'])/1e3,
                    float(slice['outerRadius'])/1e3, numConnectionRods, connectionRodsArcRadius, connectionRodsDiameter,
                    standHeight, standWidth)
                for ring in slice['rings']:
                    halbachSlice.addRing(copy.deepcopy(rings[ring['id']]), -float(slice['position'])/1e3)
                halbachSlice.setParams(params)
                halbachSlice.setId(index)
                self.addSlice(halbachSlice)
        f.close()
