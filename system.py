from mat import *
import numpy as np
from numpy import matlib
import math as mt
from constants import Constants
import constants

c = Constants()
#hbar = c.hbar
hbar = constants.hbar
A = constants.angstrom
eV = constants.eV


# class for the system
class System(object):
    def __init__(self):
        self.sys = []
        self.app_volt = 0.0

    def add_system_object(self, new_object):
        self.sys.append(new_object)

    def find_sys_length(self):
        # sum of the length arrays
        l = 0
        for i in range(len(self.sys)):
            #print(i)
            l += np.sum(self.sys[i].width_array)
            #print(l)
        return l

    def find_transmission_coefficient(self, E):
        E += self.app_volt + 0.00000001*eV #to avoid division by 0
        #print('energy', E)
        # initialising system matrix
        sys_mat = np.matlib.identity(2)
        # loop over all the system elements
        for i in range(1,len(self.sys)+1):
            obj = self.sys[i-1]
            #print(obj.width_array, obj.height_array, E)
            obj_mat = obj.find_total_obj_matrix(obj.width_array, obj.height_array, E)
            #print(obj.width_array, obj.height_array, E)
            sys_mat = np.dot(sys_mat, obj_mat)

            if i == len(self.sys):
                break
            # enters the next object
            obj_next = self.sys[i]
            # finding wave vectors on the boundary of two objects
            k_obj = obj.find_wave_vector(E, obj.height_array[len(obj.height_array)-1])
            k_next = obj_next.find_wave_vector(E, obj_next.height_array[0])
            # connecting to sys objects with discontinuity matrix
            boundary_disc_mat = DiscontinuityMatrix(k_obj, k_next, obj.mass, obj_next.mass).disc_mat

            sys_mat = np.dot(sys_mat, boundary_disc_mat)
            boundary_prop_mat = PropagationMatrix(k_next, obj_next.width_array[0]).prop_mat
            sys_mat = np.dot(sys_mat, boundary_prop_mat)
        return 1 / ((abs(sys_mat[0, 0])) ** 2)

    def apply_linear_voltage(self, applied_voltage, inc):
        print('applying voltage...')
        self.app_volt = applied_voltage
        l = self.find_sys_length()
        #print('sys length in function', l)
        divs = l/inc
        #print('divs', divs)
        height_inc = applied_voltage/divs
        #print('height inc', height_inc)
        #variable that keeps track of the height_incs
        inc_num = divs
        #print('h_inc', h_inc)
        # splitting the arrays into increments of width inc
        for obj in self.sys: #looping over sys objects
            # fraction this object's lenght is of the total system length
            obj_frac = sum(obj.width_array)/l
            # potential drop in the object
            obj_pot_drop = applied_voltage * obj_frac
            inc = self.determine_inc(inc, obj.width_array)
            #print(obj)
            w_array_new = []
            h_array_new = []
            for i in range(len(obj.width_array)): #looping over width array elements
                #print(i)
                #number of incs in the array element
                inc_num = obj.width_array[i]/inc
                #print('inc_num', inc_num)
                #print(inc, inc_num, int(inc_num))
                for j in range(int(inc_num)): #looping over new incs
                    #print(j)
                    w_array_new.append(inc)
                    h_array_new.append(obj.height_array[i]+inc_num*height_inc)
                    inc_num -= 1
                    #print(h_inc)
                    #print(w_array_new, h_array_new)
            obj.width_array = np.array(w_array_new)
            obj.height_array = np.array(h_array_new)
            #print(obj.height_array)

    def determine_inc(self, inc_given, obj):
        print('modifying the increment...')
        inc = inc_given
        old_w_array = obj.width_array
        old_h_array = obj.height_array
        # new width and height arrays
        new_w_arr = []
        new_h_arr = []
        # finding the largest inc that is smaller than the smallest increment
        for i in old_w_array:
            while i < inc:
                inc /= 2.
        # counting the old array divisions
        n = 0
        # looping over the widths in width array
        for w in old_w_array:
            # how many incs fit into the division?
            total_inc = int(w//inc)
            # how much to add to every increment (remainder per inc)
            r_per_inc = (w - inc * total_inc)/total_inc
            for i in range(total_inc):
                new_w_arr.append(inc + r_per_inc)
                new_h_arr.append(old_h_array[n])
            n += 1
        assert sum(new_w_arr)==sum(old_w_array), 'lengths are not the same'
        assert len(new_h_arr)== sum(new_w_arr), 'width and height arrays do not have the same lengths'
        obj.width_array = new_w_arr
        obj.height_array = new_h_arr

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