import math as mt
import constants
from mat import *

hbar = constants.hbar
eV = constants.eV
A = constants.angstrom

# class representing system object (either well or barrier)
class SysObject(object):
    def __init__(self, type=None, width_array=None, height_array=None, m_effective=None):
        self.total_width = sum(width_array)
        self.width_array = width_array
        self.height_array = height_array
        self.mass = m_effective
        self.type = type

    def find_wave_vector(self, E, V):
        assert E != V, 'E and V must not equal to each other'
        if E > V:
            return complex(mt.sqrt(2 * self.mass * abs(E - V)) / hbar, 0)
        elif E < V:
            return complex(0, mt.sqrt(2 * self.mass * abs(E - V)) / hbar)

    # function for finding the total propagation matrix through the object
    def find_total_obj_matrix(self, width_array, height_array, E):
        # initialising system matrix
        sys_mat_local = np.matlib.identity(2)

        # looping over the system object
        for i in range(len(width_array) - 1):
            # discontinuity matrix
            k_a = self.find_wave_vector(E, height_array[i])
            k_b = self.find_wave_vector(E, height_array[i + 1])
            disc_mat_local = DiscontinuityMatrix(k_a, k_b, self.mass,
                                                 self.mass).disc_mat
            # update system matrix
            sys_mat_local = np.dot(sys_mat_local, disc_mat_local)
            # propagation matrix
            prop_mat_local = PropagationMatrix(k_b, width_array[i + 1]).prop_mat
            # update system matrix
            sys_mat_local = np.dot(sys_mat_local, prop_mat_local)
        return sys_mat_local

    def smooth_obj_corners(self, mod_width=2*A, mod_height=0.05*eV, incs=5):
        assert 2 * mod_width < self.total_width, 'object is too narrow for smoothing the corners'
        obj_height = self.height_array[0]
        new_w_arr = []
        new_h_arr = []
        new_w_arr.append(self.total_width - 2 * mod_width)
        new_h_arr.append(obj_height)
        pot_inc = mod_height/incs
        width_inc = mod_width/incs
        n = 1 #variable to count the pot incs
        for i in range(incs):
            new_w_arr.insert(0, width_inc)
            new_w_arr.append(width_inc)
            if self.type == 'well':
                new_h_arr.insert(0, n * pot_inc)
                new_h_arr.append(n * pot_inc)
            elif self.type == 'barrier':
                new_h_arr.insert(0, obj_height - n* pot_inc)
                new_h_arr.append(obj_height - n * pot_inc)
            else:
                print('BAD OBJECT TYPE')
            n += 1
        #updating the height and width arrays
        assert len(new_h_arr) == len(new_w_arr), 'width and height arrays are not of a same length'
        self.width_array = new_w_arr
        self.height_array = new_h_arr

class SysObjSpecified(SysObject):
    def __init__(self, obj_type, width_array, height_array, m_effective):
        super(SysObjSpecified, self).__init__(obj_type, sum(width_array),
                                              width_array, height_array,
                                              m_effective)
        assert len(width_array) == len(
            height_array), 'width_array and height_array need to have same lengths'
        # self.width_array = width_array
        # self.height_array = height_array


class RectObject(SysObject):
    def __init__(self, width=None, height=None, m_effective=None):
        # dividing barrier into bits for functions to work
        self.width_array = width
        self.height_array = height
        super(RectObject, self).__init__(width, self.width_array,
                                         self.height_array, m_effective)

'''
class RoundBarrier(SysObject):
    def __init__(self, obj_type, width, min_height, max_height, divisions,
                 m_effective):
        self.min_height = min_height
        self.max_height = max_height
        self.barrier_width = width
        self.divs = divisions
        self.width_array = self.construct_width_array()
        self.height_array = self.construct_height_array()
        super(RoundBarrier, self).__init__(obj_type, width, self.width_array,
                                           self.height_array, m_effective)

    def construct_width_array(self):
        return [self.barrier_width / self.divs for i in range(self.divs)]

    def construct_height_array(self):
        h_max = self.max_height
        h_min = self.min_height
        diff = h_max - h_min
        divs = self.divs
        h_arr = []
        if divs % 2 == 1:
            print('odd')
            inc = diff / ((divs - 1) / 2)
            h_arr.append(h_max)
            for i in reversed(range(int((divs - 1) / 2))):
                h_arr.insert(0, h_min + i * inc)
                h_arr.append(h_min + i * inc)
        else:
            print('even')
            inc = diff / ((divs - 2) / 2)
            h_arr = [h_max, h_max]
            for i in reversed(range(int(divs / 2) - 1)):
                h_arr.insert(0, h_min + i * inc)
                h_arr.append(h_min + i * inc)
        return h_arr
        '''