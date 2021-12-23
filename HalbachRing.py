import numpy as np
import CubeMagnet as CubeMagnet
import matplotlib.pyplot as plt

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



if __name__ == '__main__':
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
    fig = plt.figure(figsize=(16,12))
    B0z0 = B0[evalPoints[2]==0,:]
    B_abs = np.linalg.norm(B0z0[:,0:1], axis=1)
    print("Homogeneity: " + str(((np.max(B_abs)-np.min(B_abs))/np.mean(B_abs))*1e6) + " ppm")
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
