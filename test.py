from system import *
#import constants
import matplotlib.pyplot as plt
import copy
import sim


eV = constants.eV
A = constants.angstrom
me = constants.me

my_system = System()

#width_arr = [2., 2., 2.]
#height_arr = [0.5, 0.5, 0.5]
#my_system.add_system_object(SysObjSpecified('barrier', width_arr, height_arr, 0.067))

widths_e = np.array([1.])*A
heights_e = np.array([0.])*eV
widths1 = np.array([2, 2, 2, 2])*A
heights1 = np.array([0.5, 0.5, 0.5, 0.5])*eV
widths2 = np.array([12, 12, 12, 12])*A
heights2 = np.array([0., 0., 0., 0.])*eV
widths3 = np.array([2, 2, 2, 2])*A
heights3 = np.array([0.5, 0.5, 0.5, 0.5])*eV
widths_c = np.array([1.])*A
heights_c = np.array([0.])*eV

my_system.add_system_object(SysObjSpecified('source', widths_e, heights_e, 0.067*me))
my_system.add_system_object(SquareBarrier('barrier', 10.1*A, 0.5*eV, 0.067*me))
my_system.add_system_object(SquareBarrier('well', 20.*A, 0.0, 0.067*me))
my_system.add_system_object(SquareBarrier('barrier', 10.*A, 0.5*eV, 0.067*me))
my_system.add_system_object(SysObjSpecified('drain', widths_c, heights_c, 0.067*me))

#my_system.find_transmission_coefficient(0.05)
#print(my_system.sys[1].width_array)
#print(my_system.sys[1].height_array)
print('sys length', my_system.find_sys_length())
print(my_system.sys[1].width_array)
#my_system.apply_linear_voltage(0.3*eV, 1.*A)
my_system.adjust_widths(1.*A)
print(my_system.sys[1].width_array)

#print(my_system.sys[1].height_array)
#print(my_system.app_volt)
#print('final length', my_system.find_sys_length())
#print(my_system.sys[0].height_array)
#print(my_system.sys[1].height_array)
#print(my_system.sys[2].height_array)
#print(my_system.sys[3].height_array)

print(my_system.find_transmission_coefficient(0.0))
#print(my_system.sys[1].width_array)
#my_sim = sim.create_voltage_array(my_system, [0.03, 0.3], 0.01, 1.*A)

#plt.plot(my_sim[0], my_sim[1])
#plt.yscale('log')
#plt.show()

'''v_array = np.arange(0.01, 0.35, 0.001)
print(v_array)
#e1_array = np.arange(0.05*eV, 1.5*eV, 0.145*eV)
# print(e_array)
T_array = []
for v in v_array:
    print('my_system', len(my_system.sys[1].width_array))
    my_sys = copy.deepcopy(my_system)
    print('my_sys', len(my_sys.sys[1].width_array))
    my_sys.apply_linear_voltage(v*eV, 0.1*A)
    print('my_sys after', len(my_sys.sys[1].width_array))
    T = my_sys.find_transmission_coefficient(0.0)
    T_array.append(T)

plt.plot(v_array, T_array, label='square barriers')
plt.xlabel('voltage (eV)')
plt.ylabel('transmission coefficient')
#plt.yscale('log')
#plt.legend()
plt.show()'''

#NEW CHANGES