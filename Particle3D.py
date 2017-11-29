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
    position(array) - position as an array (floats)
    velocity(array) - velocity as an array (floats)
    mass(float) - particle mass
    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    * initialize particles from input file
    """

    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle3D instance
        
        :param pos: position as an array of floats
        :param vel: velocity as an array of floats
        :param mass: mass as float
        """
 
        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label
    
    def __str__(self):
        """
        Define output format.
        """

        return str(self.label) + "x = " + str(self.position[0]) + "y= " + str(self.position[1]) + "z= " + str(self.position[2]) +  ", m = " + str(self.mass)
    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*||vel||^2
        """
        
        return 0.5*self.mass*np.inner(self.velocity,self.velocity)

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t).

        :param dt: timestep as float
        :param force: force on particle as array of floats
        """
        self.velocity = self.velocity + dt*force/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        r(t+dt) = r(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position = self.position + dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as array of floats
        """
        self.position = self.position + dt*self.velocity + 0.5*(dt**2)*force/self.mass
    
    @staticmethod
    def init_from_file(filename,l):
       """
       Initialize particles from file.

       :param filename: name of the file (should be .txt )
       :param l: number of line that stores particle parameters 
       : return : a Particle3D instance
       """
       # lines in file should be written as: Re,De,alpha,label,mass,pos[0 -> 2], vel[0 -> 2]
       file_handle = open(filename, "r") 
       lines = file_handle.readlines()  # initialise list of lists, all the lines
       choose = lines[l].split(",")  #split line "l" 
       pos = np.array([float(0) for x in xrange(0,3)]) #initialize empty arrays (floats)
       vel = np.array([float(0) for x in xrange(0,3)])
       label = choose[3]    # read mass, label
       mass = float(choose[4])
       for i in xrange(5,8):   #read pos and vel
           pos[i-5] = float(choose[i])
       for j in xrange(8,11):
           vel[j-8] = float(choose[j])
    
       return Particle3D(pos,vel,mass,label)

    @staticmethod
    def relative_pos(particleA,particleB):
        """
        Method to return relative position 
        of 2 particles.

        :param particleA, particleB : Particle3D instances
        :return : relative position as an array (floats)
        """
        return np.subtract(particleA.position,particleB.position)
