"""
CMod Ex3: symplectic Euler time integration of
a particle moving in a Morse potential.

D Harper & M Galli

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

def force_dw(particle1, particle2, alpha, De, Re):
    """

    :param particle: Particle3D instance
    :return: force acting on particle as Numpy array
    """
    r12 = Particle3D.relative_pos(particle1,particle2)
    modr12 = np.linalg.norm(r12)
    r12hat = (1/modr12)*r12
    force = -2*alpha*De*(1-math.exp(-alpha*(modr12-Re)))*math.exp(-alpha*(modr12-Re))*r12hat
    return force

def pot_energy_dw(particle1, particle2, alpha, De, Re):
    """

    :param particle: Particle3D instance
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
    dt = 0.0092
    numstep = 1000
    time = 0.0
    alpha = float(infile.readline())
    De = float(infile.readline())
    Re = float(infile.readline())
    infile.close()
    # Set up particle initial conditions:
    
    p1 = Particle3D.init_from_file(infile_name)
    p2 = Particle3D.init_from_file("newfile.txt")

    # Write out initial conditions
    energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_dw(p1 ,p2 , alpha, De,Re)
    initial_energy = energy
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))


    # Initialise data lists for plotting later
    time_list = [time]
    relpos_list = [Re]
    energy_list = [energy]
    
    # Start the time integration loop

    for i in range(numstep):
        # Update particle position
        p1.leap_pos1st(dt)
        p2.leap_pos1st(dt)
        
        # Calculate force
        force = force_dw(p1,p2, alpha, De, Re)
        # Update particle velocity 
        p1.leap_velocity(dt, force)
        p2.leap_velocity(dt, -force)
        
        # Increase time
        time = time + dt
        
        # Output particle information
        energy = p1.kinetic_energy() +p2.kinetic_energy() + pot_energy_dw(p1, p2, alpha, De, Re)
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))

        # Append information to data lists
        time_list.append(time)
        relpos_list.append(np.linalg.norm(Particle3D.relative_pos(p1,p2)))
        energy_list.append(energy)

    # Post-simulation:

    # Close output file
    delta_E = max(energy_list)-min(energy_list)
    avg_E = sum(energy_list)/len(energy_list)
    if -delta_E/avg_E <= 1e-3:
        print ("Average OK " + str(-delta_E/avg_E))
    else:
        print ("Average Fail " + str(-delta_E/avg_E))
    if -delta_E/initial_energy <= 1e-3:
        print ("Initial OK " + str(-delta_E/initial_energy))
    else:
        print ("Initial Fail " + str(-delta_E/initial_energy))
    outfile.close()

    # Plot particle trajectory to screen
    pyplot.title('Symplectic Euler: relative position vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Relative Position')
    pyplot.plot(time_list, relpos_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Symplectic Euler: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method:
main()

