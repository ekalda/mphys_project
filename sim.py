import numpy as np
import copy as cp
import constants
from types import *
import matplotlib.pyplot as plt

eV = constants.eV
A = constants.angstrom
k = constants.k
#m = 0.066 * constants.me
hbar = constants.hbar
i_const = constants.i_const

class Simulation(object):
    def __init__(self, system, write_data=False, filename=None):
        self.system = system
        self.T = 1 #K
        self.write_data = write_data
        self.f = filename

    def integrate_energy(self, e_interval, inc, surface_roughness=False, dev=2*A):
        assert isinstance(e_interval, list), 'e_interval is not a list: %r' % e_interval
        e_array = np.arange(e_interval[0], e_interval[1], inc) #not in eV
        t_array = []
        for e in e_array:
            #print(e)
            #print(e*eV)
            sys_copy = cp.deepcopy(self.system)
            if surface_roughness:
                assert dev is not None, 'no max dev specified'
                sys_copy.adjust_widths(dev)
            trans_coeff = sys_copy.find_transmission_coefficient(e*eV)
            print(trans_coeff)
            t_array.append(trans_coeff)
            if self.write_data: self.write_to_file(e, trans_coeff)
        return e_array, t_array

    def integrate_voltage(self, v_interval, v_inc=0.01, obj_inc=1*A, const_E=0.0):
        assert isinstance(v_interval, list), 'v_interval is not a list: %r' % v_interval
        v_array = np.arange(v_interval[0], v_interval[1], v_inc)
        t_array = []
        for v in v_array:
            sys_copy = cp.deepcopy(self.system)
            sys_copy.apply_linear_voltage(v*eV, obj_inc)
            trans_coeff = sys_copy.find_transmission_coefficient(const_E*eV)
            t_array.append(trans_coeff)
            if self.write_data: self.write_to_file(v, trans_coeff)
        return v_array, t_array

    def boltzmann_dist(self, energy, voltage, E_f):
        #return 1/(np.exp((energy-voltage-E_f)/(k*self.T)) + 1)
        return 1/(np.exp((E_f - energy - voltage)/(k*self.T)) + 1)

    #finding the current density through the system. For every v_inc integrate over e_interval
    def find_current(self, v_interval, e_interval, Ef_left, Ef_right=0.005*eV, e_inc=None, v_inc=0.01, obj_inc=1*A):
        assert e_inc is not None, "need to specify the energy increment!"
        v_array = np.arange(v_interval[0], v_interval[1], v_inc)
        i_array = []
        e_array = np.arange(e_interval[0], e_interval[1], e_inc)
        for v in v_array:
            sys_copy = cp.deepcopy(self.system)
            sys_copy.apply_linear_voltage(v*eV, obj_inc)
            #   print(sys_copy.sys[-1].height_array[0])
            # array to hold the values of discrete approximation for the integral
            integral = []
            for e in e_array:
                t = sys_copy.find_transmission_coefficient(e*eV)
                #assuming that the applied potential at the collector is zero!
                #integral.append(t * (self.boltzmann_dist(e*eV, v*eV, Ef_left) - self.boltzmann_dist(e*eV, 0.0, Ef_right)) * np.sqrt((v+e)*eV))
                #integral.append(eV * m * k * self.T * t / (2*np.pi*hbar**3) * (np.log(self.boltzmann_dist(e*eV, v*eV, Ef_left)) - np.log(self.boltzmann_dist(e*eV, 0.0, Ef_right))))
                #that's going to be messy
                boltz_left = np.log(self.boltzmann_dist(e*eV, sys_copy.sys[0].height_array[-1], Ef_left))
                boltz_right = np.log(self.boltzmann_dist(e*eV, sys_copy.sys[-1].height_array[0], Ef_right))
                i = -i_const * self.system.sys[0].mass * self.T * t * (boltz_right - boltz_left) #/ (2*np.pi**2*hbar**3) * eV
                integral.append(i)
                #print(boltz_left, boltz_right)
                #print(self.boltzmann_dist(e * eV, sys_copy.sys[0].height_array[-1], Ef_left))
                #print(i_array)
            i_array.append(sum(integral))
            if self.write_data: self.write_to_file(v, sum(integral))
        return v_array, i_array

    #function for plotting the graph
    def plot_graph(self, param1, param2, x_label=None, y_label=None, label=None):
        plt.plot(param1, param2, label=label)
        #plt.yscale('log')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        plt.show()

    def write_to_file(self, param1, param2):
        with open(self.f, 'a+') as f:
            f.write(str(param1)+' '+str(param2)+'\n')
