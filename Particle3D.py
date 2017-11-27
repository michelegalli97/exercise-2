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

    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle1D instance
        
        :param pos: position as float
        :param vel: velocity as float
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
        1/2*mass*vel^2
        """
        
        return 0.5*self.mass*np.inner(self.velocity,self.velocity)

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        self.velocity = self.velocity + dt*force/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position = self.position + dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        self.position = self.position + dt*self.velocity + 0.5*(dt**2)*force/self.mass
    
    @staticmethod
    def init_from_file(filename):
        file_handle = open(filename,"r")
        pos = np.array([0.1,0.1,0.1])
        vel = np.array([0.1,0.1,0.1])
        alpha = float(file_handle.readline())
        De = float(file_handle.readline())
        Re = float(file_handle.readline())
        pos[0] = float(file_handle.readline())
        pos[1] = float(file_handle.readline())
        pos[2] = float(file_handle.readline())
        vel[0] = float(file_handle.readline())
        vel[1] = float(file_handle.readline())
        vel[2] = float(file_handle.readline())
        mass = float(file_handle.readline())
        label = file_handle.readline()
        file_handle.close()
        fin = open(filename,"r")
        data_list = fin.readlines()
        fin.close()
        del data_list[3:10+1]
        print (data_list)
        fout = open("newfile.txt", "w")
        fout.writelines(data_list)
        return Particle3D(pos,vel,mass,label)

    @staticmethod
    def relative_pos(particleA,particleB):
        return np.subtract(particleA.position,particleB.position)
