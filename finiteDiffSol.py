import scipy as sp
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from math import pow
import createPlots as tplot
import sys
class Heat_Equation(object):
    """
    Class which implements a numerical solution of the 2d heat equation
    """
    def __init__(self, dx, dy, a):
                self.dx = dx # Interval size in x-direction.
                self.dy = dy # Interval size in y-direction.
                self.a = a # Diffusion constant.
                self.dx2 = dx**2
                self.dy2 = dy**2
                self.nx = int(1/dx)
                self.ny = int(1/dy)
                self.boundaries = {}
                self.bccoords = {}
                self.isBoundaryApplied = {}
                # For stability, this is the largest interval possible
                # for the size of the time-step:
                self.dt = self.dx2*self.dy2/( 2*a*(self.dx2+self.dy2) )
                self.set_initial_conditions()
                self.showPlots = 0

    def set_heat_source(self, Q):
        self.Q = Q

    def set_initial_conditions(self):
        # Start u and ui off as zero matrices:
        self.ui = np.zeros((self.nx,self.ny))
        self.u = np.zeros((self.nx,self.ny))
        self.u.fill(300)
        self.ui.fill(300)
        self.boundaries['up'] = self.ui[1:-1,-1]
        self.boundaries['down'] = self.ui[1:-1,0]
        self.boundaries['left'] = self.ui[0,1:-1]
        self.boundaries['right'] = self.ui[-1,1:-1]
        self.bccoords['up'] = (slice(1, -1), -1)
        self.bccoords['down'] = (slice(1, -1), 0)
        self.bccoords['left'] = (0, slice(1,-1))
        self.bccoords['right'] = (-1, slice(1, -1))
        self.isBoundaryApplied['up'] = 0
        self.isBoundaryApplied['down'] = 0
        self.isBoundaryApplied['left'] = 0
        self.isBoundaryApplied['right'] = 0  
        # Now, set the initial conditions (ui).
    
    def set_boundary_condition(self, name, val, boundary = 'all'):
        if name == 'constant' and boundary != 'all':
            try:
                self.ui[self.bccoords[boundary]]= val
                self.isBoundaryApplied[boundary] = 1
            except KeyError:
                print ("Invalid boundary name: {}".format(boundary))
                sys.exit(1)
        elif name == 'constant' and boundary == 'all':
            self.ui[self.bccoords['up']]= val
            self.ui[self.bccoords['down']]= val
            self.ui[self.bccoords['left']]= val
            self.ui[self.bccoords['right']]= val
            self.isBoundaryApplied['up'] = 1
            self.isBoundaryApplied['down'] = 1
            self.isBoundaryApplied['left'] = 1
            self.isBoundaryApplied['right'] = 1            
            #self.u = self.ui
        elif name == 'convection':
            Q = self.Q
            try:
                self.ui[self.bccoords[boundary]]= val
                self.isBoundaryApplied[boundary] = 1
            except KeyError:
                print ("Invalid boundary name: {}".format(boundary))
                sys.exit(1)
            if boundary == 'left':
                self.u[0,1:-1] = self.ui[0,1:-1] + self.a*self.dt*( (2*self.ui[1, 1:-1] - 2*(self.ui[0, 1:-1] + self.dx*((self.ui[1, 1:-1]-300)*3/1.9)))/self.dx2 + (self.ui[1, 2:] - 2*self.ui[1, 1:-1] + self.ui[1, :-2])/self.dy2 ) + self.dt*(Q/(753*2400))
            elif boundary == 'down':
                self.u[1:-1,0] = self.ui[1:-1,0] + self.a*self.dt*( (2*self.ui[1:-1, 1] - 2*(self.ui[1:-1, 0] + self.dx*((self.ui[1:-1, 1]-300)*3/1.9)))/self.dx2 + (self.ui[2:, 1] - 2*self.ui[1:-1, 1] + self.ui[:-2, 1])/self.dy2 ) + self.dt*(Q/(753*2400))
            else:
                print ("Condition {} not implemented for boundary: {}".format(name, boundary))
                sys.exit(1)
        elif name == 'mixed':
            try:
                self.ui[self.bccoords[boundary]]= val
                self.isBoundaryApplied[boundary] = 1
            except KeyError:
                print ("Invalid boundary name: {}".format(boundary))
                sys.exit(1)
            Q = self.Q
            sigma =  5.6703e-8
            radPowi = []
            for i in range (0, self.nx):
                radPowi.append(0.9*sigma*(pow(self.ui[0, i],4) - pow(300, 4)/1.9))

            #((0.9*sigma/1.9)*(0))
            self.u[0,1:-1] = self.ui[0,1:-1] + self.a*self.dt*( (2*self.ui[1, 1:-1] - 2*(self.ui[0, 1:-1] + self.dx*( (radPowi[1:-1]) + (self.ui[1, 1:-1]-300)*3/1.9)))/self.dx2 + (self.ui[1, 2:] - 2*self.ui[1, 1:-1] + self.ui[1, :-2])/self.dy2 ) + self.dt*(Q/(753*2400))

        elif name == 'mixed-test':
            try:
                self.ui[self.bccoords[boundary]]= val
                self.isBoundaryApplied[boundary] = 1
            except KeyError:
                print ("Invalid boundary name: {}".format(boundary))
                sys.exit(1)

            Q = self.Q
            sigma =  5.6703e-8
            radPowi = []
            for i in range (0, self.nx):
                radPowi.append(0.9*sigma*(pow(self.ui[0, i],4) - pow(300, 4))/1.9)
            
            #((0.9*sigma/1.9)*(0))
            self.u[0,1:-1] = self.ui[0,1:-1] + self.a*self.dt*( (2*self.ui[1, 1:-1] - 2*(self.ui[0, 1:-1] + self.dx*( radPowi[1:-1] + (self.ui[1, 1:-1]-300)*3/1.9)))/self.dx2 + (self.ui[1, 2:] - 2*self.ui[1, 1:-1] + self.ui[1, :-2])/self.dy2 ) + self.dt*(Q/(753*2400))

    def evolve_ts(self):
        Q = self.Q
        self.u[1:-1, 1:-1] = self.ui[1:-1, 1:-1] + self.a*self.dt*( (self.ui[2:, 1:-1] - 2*self.ui[1:-1, 1:-1] + self.ui[:-2, 1:-1])/self.dx2 + (self.ui[1:-1, 2:] - 2*self.ui[1:-1, 1:-1] + self.ui[1:-1, :-2])/self.dy2 ) + self.dt*(Q/(753*2400))
    
    def nextTime(self):
        self.ui = self.u.copy()

    def handleCorners(self):
        self.u[0,0] = np.mean([self.u[0,1], self.u[1,0]])
        self.u[-1,0] = np.mean([self.u[-1,1], self.u[-2,0]])
        self.u[0,-1] = np.mean([self.u[1,-1], self.u[0,-2]])
        self.u[-1,-1] = np.mean([self.u[-1,-2], self.u[-2,-1]])

    def run(self, timesteps):
        if 0 in self.isBoundaryApplied.values():
            print ("Not enough boundary conditions specified...")
            sys.exit(1)

        if self.showPlots == 1:
            plotter = tplot.createPlots(self.nx, self.ny)
        
        for i in range(1, timesteps, 1):
            if self.showPlots == 1:
                plotter.on(self.u)
            self.evolve_ts()
            self.set_boundary_condition('mixed-test',300, boundary = 'left')
            self.set_boundary_condition('mixed-test',300, boundary = 'down')
            self.handleCorners()
            self.nextTime()
        plotter.end()

if __name__=="__main__":
    plate = Heat_Equation(0.05, 0.05, 3.5e-7)
    plate.showPlots = 1
    plate.set_heat_source(300)
    plate.set_boundary_condition('mixed-test', 300,  boundary = 'left')
    plate.set_boundary_condition('mixed-test', 300,  boundary = 'down')
    plate.set_boundary_condition('constant', 300, boundary = 'right')
    plate.set_boundary_condition('constant', 300, boundary = 'up')
    plate.run(100)
