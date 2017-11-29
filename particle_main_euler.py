"""
CMod Ex3: symplectic Euler time integration of
2 particles moving in a Morse pair potential.
D Harper & M Galli
Produces plots of the relative distance of the particles
and the total energy, both as function of time. Also
saves both to file.
The Morse potential needs parameters Re, De, alpha. these are taken
from the input file that sets up initial conditions. (should be included by user)
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

def force_dw(particle1, particle2, alpha, De, Re):
    """
    Method to define force on a particle.

    :param particle1, particle2 : Particle3D instances
    :param alpha, De, Re: floats (potential parameters read from a file)
    :return: force acting on  (by convention) particle1 as array (floats)
    """
    r12 = Particle3D.relative_pos(particle1,particle2)  # define r12 (vector), ||r12||, and unit vector r12hat to make expression for force easier 
    modr12 = np.linalg.norm(r12)
    r12hat = (1/modr12)*r12
    force = -2*alpha*De*(1-math.exp(-alpha*(modr12-Re)))*math.exp(-alpha*(modr12-Re))*r12hat
    return force

def pot_energy_dw(particle1, particle2, alpha, De, Re):
    """
    Method to define total potential energy 

    of a pair of particles
    :param particle1, particle2: Particle3D instances
    :param alpha, De, Re: floats
    :return: total potential energy as float
    """

    modr12 = np.linalg.norm(Particle3D.relative_pos(particle1, particle2)) #define ||r12||
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
    lines = infile.readlines()  # initialise list of lists, similar to Particle3D.init_from_file()
    choose = lines[1].split(",")  #split line "1" to get potential parameters

    # Set up simulation parameters
    dt = 0.014               ######### dt is here ######## 
    numstep = 1000           ######### numsteps is here #######
    time = 0.0
    alpha = float(choose[2])  # take potential parameters from line[1]
    De = float(choose[1])     
    Re = float(choose[0])     
    infile.close()
    # Set up particles initial conditions:
    
    p1 = Particle3D.init_from_file(infile_name,1)  
    p2 = Particle3D.init_from_file(infile_name,2)  

    # Write out initial conditions
    energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_dw(p1 ,p2 , alpha, De, Re)
    initial_energy = energy
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))


    # Initialise data lists for plotting later
    time_list = [time]
    relpos_list = [Re]
    energy_list = [energy]
    
    # Start the time integration loop

    for i in range(numstep):

        # Update particles relative  position
        p1.leap_pos1st(dt)
        p2.leap_pos1st(dt)
        
        # Calculate force vector
        force = force_dw(p1,p2, alpha, De, Re)

        # Update particle velocity 
        p1.leap_velocity(dt, force) # use convention defined in force_dw 
        p2.leap_velocity(dt, -force)
        
        # Increase time
        time = time + dt
        
        # Output particles information
        energy = p1.kinetic_energy() +p2.kinetic_energy() + pot_energy_dw(p1, p2, alpha, De, Re)
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,np.linalg.norm(Particle3D.relative_pos(p1,p2)),energy))

        # Append information to data lists
        time_list.append(time)
        relpos_list.append(np.linalg.norm(Particle3D.relative_pos(p1,p2))) # list of distances (scalars) not positions (for plotting)
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

    # Plot relative distance to screen
    pyplot.title('Symplectic Euler: relative distance vs time')
    pyplot.xlabel('Time (natural unit ~ 10.18e-15 s)')
    pyplot.ylabel('Relative Distance (Angstrom)')
    pyplot.plot(time_list, relpos_list)
    pyplot.show()

    # Plot total energy to screen
    pyplot.title('Symplectic Euler: total energy vs time')
    pyplot.xlabel('Time (natural unit ~ 10.18e-15 s)')
    pyplot.ylabel('Energy (eV)')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method:
main()
