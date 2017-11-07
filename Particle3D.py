"""
CMod Ex. 3: Particle3D.py

Author: D Harper & M Galli
Version: 11/2017

"""
import numpy as np
import math

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    position(array) - position as an array
    velocity(array) - velocity as an array
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, pos, vel, mass):
        """
        Initialise a Particle1D instance
        
        :param pos: position as float
        :param vel: velocity as float
        :param mass: mass as float
        """
 
        self.position = pos
        self.velocity = vel
        self.mass = mass
    
    def __str__(self):
        """
        Define output format.
        For particle p=(2.0, 0.5, 1.0) this will print as
        "x = 2.0, v = 0.5, m = 1.0"
        """

        return "x = " + str(self.position[0]) + "y= " + str(self.position[1]) + "z= " + str(self.position[2]) +  ", m = " + str(self.mass)
    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        
        return 0.5*self.mass*math.sqrt(self.velocity[0]**2+self.velocity[1]**2=self.velocity[2]**2)

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        for i in len(velocity):
            self.velocity[i] = self.velocity[i] + dt*force/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        for i in len(velocity):
            self.position[i] = self.position[i] + dt*self.velocity[i]

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        for i in len(velocity):
            self.position[i] = self.position[i] + dt*self.velocity[i] + 0.5*dt**2*force/self.mass

    def init_from_file(filename,n):
        file_handle = open(filename, "r")
        for i in n:
            particle[i] = Particle3D(
        
