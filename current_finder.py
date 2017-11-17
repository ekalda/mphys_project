import numpy as np
import math as mt
import cmath
import matplotlib.pyplot as plt
import pylab

# constants
angstrom = 10**(-10) # m
eV = 1.602 * 10**(-19) # J
hbar = 1.0546 * 10**(-34) # J
#hbar = 6.582 * 10**(-16) # eV*s
q = 1.602 * 10**(-19) # electron charge
me = 9.109 * 10**(-31)


# class representing barrier
class Barrier():
    def __init__(self, width, height, divisions, m_effective):
        self.width = width
        self.height = height
        self.div = divisions
        self.mass = m_effective


# class representing the well
class Well():
    def __init__(self, width, height, divisions, m_effective):
        self.width = width
        self.height = height
        self.div = divisions
        self.mass = m_effective



def wave_vector(m_effective, E, V):
    if E > V:
        return complex(mt.sqrt(2*m_effective*abs(E-V))/hbar, 0)
    elif E < V:
        return complex(0, mt.sqrt(2*m_effective*abs(E-V))/hbar)


class Matrix2():
    def __init__(self):
        # creating a 2x2 matrix with zeros as all the entries
        self.M = np.complex64(np.zeros((2, 2)))

    # method for creating propagation matrix in a system with constant effective mass
    def mat_prop(self, k_i, d):
        self.k_i = k_i
        self.d = d
        # assigning values to the matrix elements
        self.M[0, 0] = cmath.exp(complex(0, -1) * complex(self.d, 0) * self.k_i)
        self.M[1, 1] = cmath.exp(complex(0, 1) * complex(self.d, 0) * self.k_i)
        return self.M


    # method for creating discontinuity matrix which takes into account the effective masses in the two regions
    def mat_disc_effective(self, k_i, k_f, m_i_effective, m_f_effective):
        self.k_i = k_i
        self.k_f = k_f
        self.m_i_effective = m_i_effective
        self.m_f_effective = m_f_effective
        self.rho = self.m_i_effective / self.m_f_effective * self.k_f / self.k_i
        # assigning values to the matrix elements
        self.M[0, 0] = 1. / 2. * (1 + self.rho)
        self.M[0, 1] = 1. / 2. * (1 - self.rho)
        self.M[1, 0] = 1. / 2. * (1 - self.rho)
        self.M[1, 1] = 1. / 2. * (1 + self.rho)
        return self.M

    # method that returns M as a complex matrix
    def get_matrix(self):
        return self.M


# function that finds the whole transmission coefficient
def find_T(energy, pot_diff, system, m_i, m_f):
    # TODO: figure out how to represent non-zero applied voltage
    # wave vector before entering the first barrier, using the effective mass specified in the arguments
    k_i = wave_vector(m_i, energy, pot_diff)
    # wave vector in the first barrier
    k_1 = wave_vector(system[0].mass, energy,
                      pot_diff + system[0].height)
    # creating the discontinuity matrix between emitter and first barrier
    M_disc1 = Matrix2()
    mat_disc1 = M_disc1.mat_disc_effective(k_i, k_1, m_i, system[0].mass)
    # creating the system matrix and assigning the first discontinuity matrix to it
    sys_mat = mat_disc1

    # looping over rest of the system
    for n in range(len(system)):
        # making the propagation matrix
        k_prop = wave_vector(system[n].mass, energy,
                             pot_diff + system[n].height)
        M_prop = Matrix2()
        mat_prop = M_prop.mat_prop(k_prop, system[n].width)
        # adding this matrix to the system matrix
        sys_mat = np.dot(sys_mat, mat_prop)

        # stop after the propagation matrix for the last system element is done. Need to add final discontinuity matrix "manually"
        if n == len(system) - 1:
            break

        # wave number for the next item, needed for discontinuity matrix
        k_a = k_prop
        k_b = wave_vector(system[n + 1].mass, energy,
                          pot_diff + system[n + 1].height)
        M_disc = Matrix2()
        mat_disc = M_disc.mat_disc_effective(k_a, k_b, system[n].mass, system[n+1].mass)
        #mat_disc = M_disc.mat_disc_effective(k_a, k_b)

        # multiplying into the system matrix
        sys_mat = np.dot(sys_mat, mat_disc)

    # the wave vector in the last barrier
    k_f_bar = wave_vector(system[len(system) - 1].mass, energy,
                          pot_diff + system[
                              len(system) - 1].height)
    # wave vector on the right of the barrier-tunnel system
    k_r = wave_vector(m_f, energy, pot_diff)
    # discontinuity matrix connecting the last barrier and the region on the right of the barrier
    M_disc_f = Matrix2()
    mat_disc_f = M_disc_f.mat_disc_effective(k_f_bar, k_r, system[
        len(system) - 1].mass, m_f)

    # multiplying M_disc_f into the system matrix
    sys_mat = np.dot(sys_mat, mat_disc_f)

    # transmission coefficient from the system matrix
    T = 1 / ((abs(sys_mat[0, 0])) ** 2)
    return T

# wells and barriers
b1_1 = Barrier(10 * angstrom, 0.5 * eV, 5, 0.067 * me)
w1_1 = Well(50 * angstrom, 0.0, 5, 0.067 * me)
b1_2 = Barrier(10 * angstrom, 0.5 * eV, 5, 0.067 * me)

# w1_2=Well(50*angstrom, 0.0, 5, 0.067*me)
# b1_3=Barrier(20*angstrom, 0.5*eV, 5, 0.067*me)
# array for the system
system1 = [b1_1, w1_1, b1_2]

e1_array = np.arange(0.05, 1.5, 5.5 * 10 ** (-4))
# print(e_array)
T1_array = []
for e in e1_array:
    T = find_T(e * eV, 0.0, system1, 0.067 * me, 0.067 * me)
    T1_array.append(T)
print (find_T(0.08*eV, 0, system1, 0.067*me, 0.067*me))

# print (T_array)
plt.plot(e1_array, T1_array)
# plt.yscale('log')
plt.show()

