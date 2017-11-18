import math as mt
import constants
from mat import *

hbar = constants.hbar

# class representing system object (either well or barrier)
class SysObject(object):
    def __init__(self, obj_type, width, width_array, height_array, m_effective):
        self.obj_type = obj_type
        self.total_width = width
        self.width_array = width_array
        self.height_array = height_array
        self.mass = m_effective

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


class SysObjSpecified(SysObject):
    def __init__(self, obj_type, width_array, height_array, m_effective):
        super(SysObjSpecified, self).__init__(obj_type, sum(width_array),
                                              width_array, height_array,
                                              m_effective)
        assert len(width_array) == len(
            height_array), 'width_array and height_array need to have same lengths'
        # self.width_array = width_array
        # self.height_array = height_array


class SquareBarrier(SysObject):
    def __init__(self, obj_type, width, height, m_effective):
        # dividing barrier into bits for functions to work
        self.width_array = [width]
        self.height_array = [height]
        super(SquareBarrier, self).__init__(obj_type, width, self.width_array,
                                            self.height_array, m_effective)


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