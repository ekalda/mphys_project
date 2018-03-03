import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('sim_data\\f2_data\\sim4.txt')
data2 = np.loadtxt('sim_data\\f2_data\\sim5.txt')
data3 = np.loadtxt('sim_data\\after_bug_fix\\fermi3.txt')
data4 = np.loadtxt('sim_data\\after_bug_fix\\fermi4.txt')
#data5 = np.loadtxt('sim_data\\test_runs\\test5.txt')
#data3 = np.loadtxt('IV_300K.txt')
#print(data1)

plt.plot(data1[:, 0], data1[:, 1], linestyle='-', label='sim4')
plt.plot(data2[:, 0], data2[:, 1], linestyle='--', label='sim5')
#plt.plot(data3[:, 0], data3[:, 1], linestyle=':', label='77K, E_f = 0.05')
#plt.plot(data4[:, 0], data4[:, 1], linestyle='-.', label='77K, E_f = 0.1')
#plt.plot(data5[:, 0], data5[:, 1], label='1K, E_f = 0.199 eV')
#plt.plot(data3[:, 0], data3[:, 1], label='300K')
plt.xlabel('V(eV)')
plt.ylabel('I(arbitrary units)')
plt.legend()
plt.yscale('log')
plt.show()