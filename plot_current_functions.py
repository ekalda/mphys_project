import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('current_functions\\sim1.txt')
data2 = np.loadtxt('current_functions\\sim4.txt')
data3 = np.loadtxt('current_functions\\sim6.txt')
data4 = np.loadtxt('current_functions\\sim8.txt')

plt.plot(data1[:, 0], data1[:, 1], 'black', linestyle='-', label='E_f = 0.005, T = 1K')
plt.plot(data1[:, 0], data1[:, 2], 'gray', linestyle='-', label='E_f = 0.005, T = 1K')
plt.plot(data2[:, 0], data2[:, 1], 'black', linestyle='--', label='E_f = 0.5, T = 1K')
plt.plot(data2[:, 0], data2[:, 2], 'gray', linestyle='--', label='E_f = 0.5, T = 1K')
plt.plot(data3[:, 0], data3[:, 1], 'black', linestyle='-.', label='E_f = 0.1, T = 77K')
plt.plot(data3[:, 0], data3[:, 2], 'gray', linestyle='-.', label='E_f = 0.1, T = 77K')
plt.plot(data4[:, 0], data4[:, 1], 'black', linestyle=':', label='E_f = 0.5, T = 300K')
plt.plot(data4[:, 0], data4[:, 2], 'gray', linestyle=':', label='E_f = 0.5, T = 300K')
plt.legend()
plt.xlabel('E(eV)')
plt.yscale('log')
plt.show()