import numpy as np
import cmath


class PropagationMatrix(object):
    def __init__(self, k_i, d):
        self.k_i = k_i
        self.width = d
        self.prop_mat = self.find_propagation_matrix()

    def find_propagation_matrix(self):
        # creating numpy matrix
        self.prop_mat = np.complex64(np.zeros((2, 2)))
        # assigning values to the matrix elements
        self.prop_mat[0, 0] = cmath.exp(
            complex(0, -1) * complex(self.width, 0) * self.k_i)
        self.prop_mat[1, 1] = cmath.exp(
            complex(0, 1) * complex(self.width, 0) * self.k_i)
        return self.prop_mat


class DiscontinuityMatrix(object):
    def __init__(self, k_i, k_f, m_i_effective, m_f_effective):
        self.rho = m_i_effective / m_f_effective * k_f / k_i
        self.disc_mat = self.find_discontinuity_matrix()

    def find_discontinuity_matrix(self):
        # creating numpy matrix
        self.disc_mat = np.complex64(np.zeros((2, 2)))
        # assigning values to the matrix elements
        self.disc_mat[0, 0] = 1. / 2. * (1 + self.rho)
        self.disc_mat[0, 1] = 1. / 2. * (1 - self.rho)
        self.disc_mat[1, 0] = 1. / 2. * (1 - self.rho)
        self.disc_mat[1, 1] = 1. / 2. * (1 + self.rho)
        return self.disc_mat
