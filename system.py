from mat import *
import numpy as np
from numpy import matlib
import math as mt
from constants import Constants
from sys_objects import *
import constants
import random
import prof

c = Constants()
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
        #print('finding transmission coefficient...')
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
        #print('applying voltage...')
        self.app_volt = applied_voltage
        l = self.find_sys_length()
        # find the array of potential drops
        v_arr = self.find_pot_drops(l)
        for n in range(len(self.sys)): #looping over sys objects
            obj = self.sys[n]
            bounds = v_arr[n]
            #adjust the increments. Calling this function should change the h and w arrays of the object
            self.determine_inc(inc, obj)
            # number of divisions in the object
            obj_divs = len(obj.width_array)
            #potential drop per division
            pot_inc = (bounds[0]-bounds[1])/obj_divs
            # new height array
            h_array_new = []
            #potential added to the first division
            add_pot = bounds[0] - pot_inc
            # looping over width array elements
            for i in range(obj_divs):
                h_array_new.append(obj.height_array[i] + add_pot)
                add_pot -= pot_inc
            obj.height_array = h_array_new

    def determine_inc(self, inc_given, obj):
        #print('modifying the increment...')
        inc = inc_given
        old_w_array = obj.width_array
        old_h_array = obj.height_array
        # new width and height arrays
        new_w_arr = []
        new_h_arr = []
        # finding the largest inc that is smaller than the smallest increment
        for i in old_w_array:
            while i < inc:
                #print('WHILE')
                inc /= 2.
        # counting the old array divisions
        n = 0
        # looping over the widths in width array
        for w in old_w_array:
            # how many incs fit into the division?
            total_incs = int(w//inc)
            # how much to add to every increment (remainder per inc)
            r_per_inc = (w - inc * total_incs)/total_incs
            for i in range(total_incs):
                new_w_arr.append(inc + r_per_inc)
                new_h_arr.append(old_h_array[n])
            n += 1
        #print(sum(new_w_arr), sum(old_w_array))
        assert abs(sum(new_w_arr)-sum(old_w_array)) < 0.000001, 'lengths are not the same'
        #assert len(new_h_arr) == sum(new_w_arr), 'width and height arrays do not have the same lengths' DOESNT WORK
        obj.width_array = new_w_arr
        obj.height_array = new_h_arr

    #function that returns the potential drops
    def find_pot_drops(self, sys_length):
        v = self.app_volt
        l = sys_length
        #array for voltages
        v_arr = []
        for obj in self.sys:
            frac = obj.total_width/l
            #potential drop per object
            pot_drop = self.app_volt*frac
            #append bounadaries for the object to the v_arr
            v_arr.append([v, v-pot_drop])
            v -= pot_drop
        assert len(v_arr)==len(self.sys), 'voltage drop array does not have a length which is same as the elements in a system'
        return v_arr

    #@prof.do_cprofile
    def adjust_widths(self, max_deviation):
        assert max_deviation >= 0., 'max_deviation has to be positive real number'
        init_sys_l = self.find_sys_length()
        for n in range(len(self.sys)-1):
            obj = self.sys[n]
            obj_next = self.sys[n+1]
            dev = random.uniform(-max_deviation, max_deviation)
            obj.width_array[-1] += dev
            obj_next.width_array[0] -= dev
        adj_sys_l = self.find_sys_length()
        assert abs(init_sys_l-adj_sys_l) < 0.000001, 'adjusting object widths changed the system length'

    def smooth_corners(self, objects=None):
        # check if there is a list of object for which to do the smoothing and if not, do it for all the objects
        if objects is None:
            i = 0
            objects = []
            for i in range(len(self.sys)):
                objects.append(i)
                i += 1
        assert isinstance(objects, list), 'argument needs to be a list of object indices'
        for n in objects:
            if self.sys[n].type == 'barrier' or self.sys[n].type == 'well':
                self.sys[n].smooth_obj_corners()


