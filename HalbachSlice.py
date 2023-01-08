import numpy as np
import solid

def arc( height, depth, radius, degrees ):
    retVal = solid.render()(
            solid.difference()(
            # Outer ring
            solid.rotate_extrude(segments = 100)(
                solid.translate([radius - height, 0, 0])(
                    solid.square([height,depth])
                )
            ),
            # Cut half off
            solid.translate([0,-(radius+1),-.5]) (
                solid.cube([radius+1,(radius+1)*2,depth+1])
            ),
            # Cover the other half as necessary
            solid.rotate([0,0,180-degrees])(
                solid.translate([0,-(radius+1),-.5]) (
                    solid.cube ([radius+1,(radius+1)*2,depth+1])
                )
            )
        )
    )
    return retVal

class HalbachSlice:
    def __init__(self, position, innerRadius, outerRadius, numConnectionRods, connectionRodsArcRadius,
        connectionRodsDiameter, standHeight = 0, standWidth = 0):
        self.position = position
        self.innerRadius = innerRadius
        self.outerRadius = outerRadius
        self.standHeight = standHeight
        self.standWidth = standWidth
        self.numConnectionRods = numConnectionRods
        self.connectionRodDiameter = connectionRodsDiameter
        self.shimTrayHeight = 0.014
        self.shimTrayAngle = 25
        self.shimTrayRadius = self.outerRadius*1e3 - 10
        self.rodRadius = connectionRodsArcRadius*1e3
        self.rings = []
        self.id = 0

    def setParams(self, params):
        self.params = params

    def setId(self, id):
        self.id = id

    def addRing(self, ring, position):
        ring.setPosition(position)
        self.rings.append(ring)

    def setPosition(self, position):
        self.position = position
        for ring in self.rings:
            ring.setPosition(position)

    def generateSCADObject(self):
        shimTray = arc(16,2000,self.shimTrayRadius,self.shimTrayAngle)
        stand = solid.linear_extrude(height = 12, center = True)(
            solid.polygon(points = [
                [-self.standWidth/2*1e3,-self.standHeight*1e3],
                [self.standWidth/2*1e3,-self.standHeight*1e3],
                [self.outerRadius*1e3,0],
                [-self.outerRadius*1e3,0]])
        )
        disc = solid.translate([0, 0, self.position*1e3])(
            stand + solid.cylinder(h=12, d=self.outerRadius*2*1e3, center=True) 
            - solid.cylinder(h=13, d=self.innerRadius*2*1e3, center=True)
            )
        cubes = []
        for ring in self.rings:
            for magnet in ring.magnets:
                cubes.append(solid.translate([magnet.position[0]*1e3, magnet.position[1]*1e3,
                    (magnet.position[2] + ring.position)*1e3])(
                    solid.rotate([0, 0, magnet.angle*180/np.pi])(
                        solid.cube([magnet.magnetType.dimension*1e3, magnet.magnetType.dimension*1e3,
                            magnet.magnetType.dimension*1e3+2], center=True)
                    )
                ))
        rods = []
        for n in np.arange(self.numConnectionRods):
            rods.append(
                solid.rotate([0, 0, 360/self.numConnectionRods*n + 360/self.numConnectionRods*0.5])(
                    solid.translate([self.rodRadius, 0, 0])(
                        solid.cylinder(h=2000, d=self.connectionRodDiameter*1e3, center=True)
                    )
                ))
        shimTrays = []
        for n in np.arange(self.numConnectionRods):
            shimTrays.append(
                solid.translate([0,0,-1000])(
                    solid.rotate([0, 0, 360/self.numConnectionRods*n + self.shimTrayAngle/2])(
                        shimTray
                    )
                ))
        addSCAD = []
        if 'addSCAD' in self.params['slices'][self.id]:
            addSCAD.append(
                solid.translate([0,0,self.position*1e3])(solid.code(self.params['slices'][self.id]['addSCAD'])))
        subSCAD = []
        if 'subSCAD' in self.params['slices'][self.id]:
            addSCAD.append(
                solid.translate([0,0,self.position*1e3])(solid.code(self.params['slices'][self.id]['subSCAD'])))
        return disc - cubes - rods - shimTrays + addSCAD - subSCAD

    def generateSCADFile(self, filename):
        solid.scad_render_to_file(self.generateSCADObject(), filename, include_orig_code=False)
