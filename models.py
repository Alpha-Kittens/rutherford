import numpy as np

def const_cross_section(c, theta):
    return c / (np.sin(theta))**4
    
def cross_section(Z, E, theta, Z_alpha = 2, e = 1):
    return (Z * Z_alpha * e**2 / 4*E)**2 / (np.sin(theta))**4

def convolve(x, f, g):
    return np.convolve(f(x), g(x))
