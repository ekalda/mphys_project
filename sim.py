import numpy as np
import copy as cp
import constants
from types import *

eV = constants.eV

def create_energy_array(system, e_interval, inc, surface_roughness=False, dev=None, smooth_corners=False):
    assert isinstance(e_interval, list), 'e_interval is not a list: %r' % e_interval
    e_array = np.arange(e_interval[0], e_interval[1], inc) #in eV
    t_array = []
    for e in e_array:
        sys_copy = cp.deepcopy(system)
        if surface_roughness:
            assert dev is not None, 'no max dev specified'
            sys_copy.adjust_widths(dev)
        if smooth_corners:
            sys_copy.smooth_corners()
        trans_coeff = sys_copy.find_transmission_coefficient(e*eV)
        t_array.append(trans_coeff)
    return e_array, t_array

def create_voltage_array(system, v_interval, v_inc, obj_inc):
    assert isinstance(v_interval, list), 'v_interval is not a list: %r' % v_interval
    v_array = np.arange(v_interval[0], v_interval[1], v_inc)
    t_array = []
    for v in v_array:
        sys_copy = cp.deepcopy(system)
        sys_copy.apply_linear_voltage(v*eV, obj_inc)
        trans_coeff = sys_copy.find_transmission_coefficient(0.0)
        t_array.append(trans_coeff)
    return v_array, t_array
