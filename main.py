import numpy as np
import simplejson as json
from collections import OrderedDict
import ast
from system import *
from sim import *
import os.path

eV = constants.eV
A = constants.angstrom
me = constants.me


def main():
    print('reading the input data...')
    with open('input.dat', 'r') as f:
        config = json.load(f, object_pairs_hook=OrderedDict)
    # effective masses
    m_b = config["barrier_effective_mass"]*me
    m_w = config["well_effective_mass"]*me
    objects = config["barriers_and_wells"]
    shield_barrier = None
    if config['shield barrier'] != 'None':
        shield_barrier = config['shield barrier']

    print('creating the system object...')
    system = System()

    print('adding the emitter...')
    system.add_system_object(SysObject(type='emitter', width_array=np.array([config["emitter_width"]])*A, height_array=np.array([0]), m_effective=m_w))
    #print(system.sys[0].width_array)

    if shield_barrier:
        print('adding the first shielding barrier...')
        system.add_system_object(SysObject(type='shield barrier', width_array=np.array([shield_barrier[1]])*A, height_array=np.array([shield_barrier[0]])*eV, m_effective=m_b))

    print('adding the wells and barriers...')
    for obj in objects:
        if obj[0] == 0:
            system.add_system_object(SysObject(type='well', width_array=np.array([obj[1]])*A, height_array=np.array([obj[0]])*eV, m_effective=m_w)) #it's a well!
        else:
            system.add_system_object(SysObject(type='barrier', width_array=np.array([obj[1]])*A, height_array=np.array([obj[0]])*eV, m_effective=m_b)) #it's a barrier!

    if shield_barrier:
        print('adding the second shielding barrier...')
        system.add_system_object(SysObject(type='shield barrier', width_array=np.array([shield_barrier[1]])*A, height_array=np.array([shield_barrier[0]])*eV, m_effective=m_b))

    print('adding the collector...')
    system.add_system_object(SysObject(type='collector', width_array=np.array([config["collector_width"]])*A, height_array=np.array([0]), m_effective=m_w))

    smooth_corners = ast.literal_eval(config["smooth_corners"])
    # smooth the corners if input boolean is True
    if smooth_corners:
        system.smooth_corners()

    apply_constant_voltage = ast.literal_eval(config["apply_voltage"])
    if apply_constant_voltage:
        voltage = config["voltage_drop"]
        print('applying constant voltage of %s to the system...' %voltage)
        system.apply_linear_voltage(voltage * eV, 1*A)
        print(system.app_volt)

    print('setting up the simulation...')
    write_data = ast.literal_eval(config["write_data"])
    filename = os.path.normpath(config["filename"])

    sim = Simulation(system, write_data=write_data, filename=filename)
    print(sim.system.app_volt)

    integrate_energy = ast.literal_eval(config["integrate_energy"])
    if integrate_energy:
        print('integrating the energy...')
        e_range = config["energy_range"]
        # into how many increments we want to divide the range
        inc = (e_range[1] - e_range[0])/100
        surface_rough = ast.literal_eval(config["surface_roughness"])
        dev = config["deviation"] * A
        results = sim.integrate_energy(e_range, inc, surface_roughness=surface_rough, dev=dev)
        sim.plot_graph(results[0], results[1], 'E', 'T')

    integrate_voltage = ast.literal_eval(config["integrate_voltage"])
    if integrate_voltage:
        print('integrating the voltage...')
        v_range = config["voltage_range"] #not in eV
        inc = (v_range[1] - v_range[0])/100 #not in eV
        const_e = config["constant_energy"] * eV
        results = sim.integrate_voltage(v_range, v_inc=inc, const_E=const_e)
        sim.plot_graph(results[0], results[1], 'V', 'T')

    find_current = ast.literal_eval(config["find_current"])
    if find_current:
        print('finding the current...')
        sim.T = config["temperature"]

        v_range = config["voltage_range"]  #not in eV
        v_inc = (v_range[1] - v_range[0])/200  #not in eV
        e_range = config["energy_range"]  #not in eV
        e_inc = (e_range[1] - e_range[0])/150  #not in eV
        E_f = config["E_f"] * eV
        results = sim.find_current(v_range, e_range, E_f, e_inc=e_inc, v_inc=v_inc)
        sim.plot_graph(results[0], results[1], 'V', 'I')


if __name__ == '__main__':
    main()