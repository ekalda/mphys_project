from system import *
from sys_objects import *
from mat import *
from sim import *
import constants
import matplotlib.pyplot as plt
import numpy as np

me = constants.me
A = constants.angstrom
eV = constants.eV
print(me, A, eV)


e_num = 100
v_num = 100
runs = 1

my_sys = System()
print(my_sys.app_volt)
my_sys.add_system_object(SysObject(type='emitter', width_array=[2*A], height_array=[0.0], m_effective=0.066*me))
my_sys.add_system_object(SysObject(type='barrier', width_array=[20*A], height_array=[0.341*eV], m_effective=0.101*me))
my_sys.add_system_object(SysObject(type='well', width_array=[50*A], height_array=[0.0], m_effective=0.066*me))
my_sys.add_system_object(SysObject(type='barrier', width_array=[20*A], height_array=[0.341*eV], m_effective=0.101*me))
my_sys.add_system_object(SysObject(type='collector', width_array=[2*A], height_array=[0.0], m_effective=0.066*me))

def run(out_file):

    simulation = Simulation(my_sys)

    e_range = [0.0001, 1.0]
    e_inc = (e_range[1] - e_range[0])/e_num
    e_array = np.arange(e_range[0], e_range[1], e_inc)

    v_range = [0.00001, 0.75]
    v_inc = (v_range[1] - v_range[0]) / v_num
    v_array = np.arange(v_range[0], v_range[1], v_inc)

    results_arr = np.ndarray(shape=(runs+1, e_num))
    results_arr[0, :] = e_array
    print(results_arr)

    for i in range(runs):
        print(i)
        #results = simulation.integrate_energy(e_range, e_inc, surface_roughness=True, dev=2*A)
        results = simulation.find_current(v_range, e_range, 0.01, Ef_right=0.01, e_inc=e_inc, v_inc=v_inc, sr=True, dev=2*A)
        results_arr[i+1] = results[1]

    print(results_arr)
    with open(out_file, 'wb') as f:
        np.savetxt(f, results_arr, fmt='%.4e', delimiter=' ')
        #results_arr.tofile(f)


def plot_data(data_file, data_file2=None):
    results_arr = np.loadtxt(data_file)
    if data_file2: results_arr2 = np.loadtxt(data_file2)
    print(results_arr.shape)
    plt.plot(results_arr[0, :], np.mean(results_arr[1:, :], axis=0), label='mean')
    #plt.plot(results_arr[0, :], results_arr[1, :], label='dev=2A')
    #plt.plot(results_arr[0, :], np.std(results_arr[1:, :], axis=0), label='std')
    #plt.plot(results_arr2[0, :], results_arr2[1, :], label='dev=0A', linestyle='--')
    plt.xlabel('E ')
    plt.ylabel('T(E)')
    plt.legend()
    #plt.yscale('log')
    plt.show()

out_file = 'surface_roughness\\sf_20_50_20_current_rough.txt'
#out_file2 = 'surface_roughness\\one_run_5_20_5_no_sf_test.txt'
run(out_file)
plot_data(out_file, data_file2=None)
