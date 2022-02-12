import numpy as np
import matplotlib.pyplot as plt

class Magnet:
    def __init__(self, position, angle, magnetType):
        self.position = position
        self.angle = angle
        self.magnetType = magnetType
    
    def magneticMoment_withoutMu0(self):
        magneticMoment_abs = self.magnetType.BR * self.magnetType.dimension**3
        magneticMoment = [magneticMoment_abs*np.cos(self.angle), magneticMoment_abs*np.sin(self.angle)]
        return magneticMoment        

    def draw2DMatplotlib(self):
        dimension = self.magnetType.dimension
        # rectangle rotates around corner, but we want rotation around center, so some math is needed here
        x = self.position[0] + dimension/np.sqrt(2)*np.sin(self.angle - np.pi/4)
        y = self.position[1] - dimension/np.sqrt(2)*np.cos(self.angle - np.pi/4)
        rectangle = plt.Rectangle((x, y), dimension, dimension, self.angle*360/(2*np.pi), fc='white',ec="red")
        plt.gca().add_patch(rectangle)
        x = self.position[0] + dimension/np.sqrt(2)*np.sin(self.angle - np.pi/4) - dimension/2*np.sin(self.angle)
        y = self.position[1] - dimension/np.sqrt(2)*np.cos(self.angle - np.pi/4) + dimension/2*np.cos(self.angle)            
        plt.arrow(x, y, dimension/3*np.cos(self.angle), dimension/3*np.sin(self.angle), color='red')

    def drawGMSH(self):
        pass
