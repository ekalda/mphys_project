import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('sim_data\\test_runs\IV_77K.txt')
data2 = np.loadtxt('sim_data\\test_runs\IV_77K_ver3.txt')
#data3 = np.loadtxt('IV_300K.txt')
#print(data1)

plt.plot(data1[:, 0], data1[:, 1], label='77K, energy range [0, 0.005]')
plt.plot(data2[:, 0], data2[:, 1], label='77K, energy range [0, 0.6]')
#plt.plot(data3[:, 0], data3[:, 1], label='300K')
plt.xlabel('V(eV)')
plt.ylabel('I(arbitrary units)')
plt.legend()
plt.yscale('log')
plt.show()