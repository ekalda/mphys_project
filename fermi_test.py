from sim import *
from system import *
import matplotlib.pyplot as plt
import numpy as np
import constants

eV = constants.eV

energies = np.arange(0, 0.01, 0.0001)
print(energies)

sys = System()

def find_dist(temp, energies):
    dist = []
    for e in energies:
        sim = Simulation(sys)
        sim.T = temp
        dist.append(sim.boltzmann_dist(e*eV, 0, 0.005*eV))
    return dist

dist0 = find_dist(0.1, energies)
dist0_1 = find_dist(3, energies)
dist20 = find_dist(20, energies)
dist77 = find_dist(77, energies)

plt.plot(energies, dist0, label='0K')
plt.plot(energies, dist0_1, label='3K')
plt.plot(energies, dist20, label='20K')
plt.plot(energies, dist77, label='77K')
plt.ylabel('probability density')
plt.xlabel('E')
plt.legend()
plt.show()