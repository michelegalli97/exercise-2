"""
CMod Ex3: velocity Verlet time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

def force_dw(particle1, particle2, alpha, De, Re):
    """
    Method to return the force on a particle
    in a double well potential.
    :return: force acting on particle as Numpy array
    """
    r12 = Particle3D.relative_pos(particle1,particle2)
    modr12 = np.linalg.norm(r12)
    force = -2*alpha*De*(1-math.exp(-alpha*(modr12-Re)))*math.exp(-alpha*(modr12-Re))*r12
    return force

def pot_energy_dw(particle1, particle2, alpha, De, Re):
    """
    Method to return potential energy 
    of particle in double-well potential
    :param a: parameter a from potential
    :param b: parameter b from potential
    :return: potential energy of particle as float
    """
    modr12 = np.linalg.norm(Particle3D.relative_pos(particle1, particle2))
    potential = De*((1-math.exp(-alpha*(modr12-Re)))**2-1)
    return potential

# Begin main code
def main():

    # Read name of output file from command line
    if len(sys.argv)!=3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[1]
        infile_name = sys.argv[2]

    # Open output file
    outfile = open(outfile_name, "w")
    infile = open(infile_name, "r")

    # Set up simulation parameters
    dt = 0.001
    numstep = 10000
    time = 0.0
    alpha = float(infile.readline())
    De =  float(infile.readline())
    Re = float(infile.readline())
    infile.close()

    # Set up particle initial conditions:
    #  position x0 = 0.0
    #  velocity v0 = 1.0
    #  mass      m = 1.0
    p1 = Particle3D.init_from_file(infile_name)
    p2 = Particle3D.init_from_file("newfile.txt")

    # Write out initial conditions
    energy = p1.kinetic_energy() +p2.kinetic_energy() + pot_energy_dw(p1, p2, alpha, De, Re)
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))

    # Get initial force
    force = force_dw(p1, p2, alpha, De, Re)

    # Initialise data lists for plotting later
    time_list = [time]
    relpos_list = [Re]
    energy_list = [energy]

    # Start the time integration loop

    for i in range(numstep):
        # Update particle position
        p1.leap_pos2nd(dt, force)
        p2.leap_pos2nd(dt, -force)
        
        # Update force
        force_new = force_dw(p1, p2, alpha, De, Re)
        # Update particle velocity by averaging
        # current and new forces
        p1.leap_velocity(dt, 0.5*(force+force_new))
        p2.leap_velocity(dt, -0.5*(force+force_new))
        # Re-define force value
        force = force_new

        # Increase time
        time = time + dt
        
        # Output particle information
        energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_dw(p1, p2, alpha, De, Re)
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))

        # Append information to data lists
        time_list.append(time)
        relpos_list.append(np.linalg.norm(Particle3D.relative_pos(p1,p2)))
        energy_list.append(energy)

    # Post-simulation:
    
    # Close output file
    outfile.close()

   # Plot particle trajectory to screen
    pyplot.title('Velocity Verlet: relative position vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Relative Position')
    pyplot.plot(time_list, relpos_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method:
main()

