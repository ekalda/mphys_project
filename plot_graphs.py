import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('GaAs\\20_50_20_current.csv')
data2 = np.loadtxt('GaAs\\20_50_20_long_barriers.txt')
#data3 = np.loadtxt('sim_data\\current\\sim10.txt')
#data4 = np.loadtxt('sim_data\\current\\sim11.txt')
#data5 = np.loadtxt('sim_data\\test_runs\\test5.txt')
#data3 = np.loadtxt('IV_300K.txt')
#print(data1)

plt.plot(data1[:, 0], data1[:, 1], label='emitter/collector widths 1A')
#plt.plot(data2[:, 0], data2[:, 1], linestyle='--', label='emitter/collector widths 10A')
#plt.plot(data3[:, 0], data3[:, 1], linestyle=':', label='300K')
#plt.plot(data4[:, 0], data4[:, 1], linestyle='-', label='400K')
#plt.plot(data5[:, 0], data5[:, 1], label='1K, E_f = 0.199 eV')
#plt.plot(data3[:, 0], data3[:, 1], label='300K')
plt.xlabel('V(eV)')
plt.ylabel('I')
plt.legend()
plt.yscale('log')
plt.show()