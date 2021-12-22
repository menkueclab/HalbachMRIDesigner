import numpy as np
import matplotlib.pyplot as plt

class Magnet:
    def __init__(self, position, angle):
        self.position = position
        self.angle = angle
        self.dimension = 0.012
        self.B_rem = 1.3
    
    def setZPosition(self, position):
        self.position = (self.position[0], self.position[1], position)

    def magneticMoment_withoutMu0(self):
        magneticMoment_abs = self.B_rem * self.dimension**3
        magneticMoment = [magneticMoment_abs*np.cos(self.angle), magneticMoment_abs*np.sin(self.angle)]
        return magneticMoment        

    def draw2DMatplotlib(self):
        # rectangle rotates around corner, but we want rotation around center, so some math is needed here
        x = self.position[0] + self.dimension/np.sqrt(2)*np.sin(self.angle - np.pi/4)
        y = self.position[1] - self.dimension/np.sqrt(2)*np.cos(self.angle - np.pi/4)
        rectangle = plt.Rectangle((x, y), self.dimension, self.dimension, self.angle*360/(2*np.pi), fc='white',ec="red")
        plt.gca().add_patch(rectangle)
        x = self.position[0] + self.dimension/np.sqrt(2)*np.sin(self.angle - np.pi/4) - self.dimension/2*np.sin(self.angle)
        y = self.position[1] - self.dimension/np.sqrt(2)*np.cos(self.angle - np.pi/4) + self.dimension/2*np.cos(self.angle)            
        plt.arrow(x, y, self.dimension/3*np.cos(self.angle), self.dimension/3*np.sin(self.angle), color='red')

    def drawGMSH(self):
        pass
