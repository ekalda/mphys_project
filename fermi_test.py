from sim import *
from system import *
import matplotlib.pyplot as plt
import numpy as np
import constants

eV = constants.eV

energies = np.arange(0, 0.01, 0.0001)
print(energies)

sys = System()

def find_dist(temp, energies, E_f):
    dist = []
    for e in energies:
        sim = Simulation(sys)
        sim.T = temp
        dist.append(sim.boltzmann_dist(e*eV, 0, E_f*eV))
    return dist

#dist0 = find_dist(0.1, energies)
#dist0_1 = find_dist(3, energies)
#dist20 = find_dist(20, energies)
dist1 = find_dist(10, energies, 0.005)
dist2 = find_dist(10, energies, -0.01)

#plt.plot(energies, dist0, label='0K')
#plt.plot(energies, dist0_1, label='3K')
#plt.plot(energies, dist20, label='20K')
plt.plot(energies, dist1, label='77K, E_f = 0.005')
plt.plot(energies, dist2, label='77K, E_f = -0.01')
plt.ylabel('probability density')
plt.xlabel('E')
plt.legend()
plt.show()