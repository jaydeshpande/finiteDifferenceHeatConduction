import matplotlib.pyplot as plt 
from matplotlib import pylab 
from matplotlib import cm
import numpy as np
class createPlots:
    def __init__(self, nx, ny):
        self.x = np.linspace(0, 1, nx)
        y = np.linspace(0, 1, ny)
        self.X, self.Y = np.meshgrid(self.x, y)
        plt.ion()
        fig = plt.figure(1)
        fig2 = plt.figure(2)
        self.ax2 = fig2.gca()
        self.ax = fig.gca(projection='3d')

    def on(self,u):
        self.ax.clear()
        self.ax2.clear()
        self.ax.set_xlabel('X Dimension [m]')
        self.ax.set_ylabel('Y Dimension [m]')
        self.ax.set_zlabel('Temperature [K]')
        self.ax2.set_xlabel('Y Dimension [m]')
        self.ax2.set_ylabel('Temperature [K]')
        self.ax.plot_surface(self.X, self.Y, u, cmap=cm.rainbow,linewidth=0, antialiased=False)
        self.ax2.plot(self.x,u[:,5])
        plt.show()
        plt.pause(0.04) 
        
    def end(self):    
        pylab.show(block=True)