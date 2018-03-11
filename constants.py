
angstrom = 10 ** (-10)  # meters
eV = 1.602 * 10 ** (-19)  # Joules
#eV = 1 # in eV
hbar = 1.0546 * 10 ** (-34)  # in J
#hbar = 6.582 * 10**(-16) # eV*s
q = 1.602 * 10 ** (-19)  # in J
#q = -1 # in eV
me = 9.109 * 10 ** (-31)
#me = 0.511 * 10**6 # in eV
k = 1.380 * 10 **(-23) # J/K
#k = 8.617 * 10**(-5) # in eV/K

class Constants(object):
    def __init__(self):
        self.angstrom = 10 ** (-10)  # meters
        self.eV = 1.602 * 10 ** (-19)  # Joules
        self.hbar = 1.0546 * 10 ** (-34)  # in J
        # hbar = 6.582 * 10**(-16) # eV*s
        self.q = 1.602 * 10 ** (-19)  # electron charge
        self.me = 9.109 * 10 ** (-31)