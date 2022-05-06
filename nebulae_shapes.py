"""
Contains some sample nebula shapes.
"""


from math import *
from nebula_generator import Nebula, sph_coord, r_spherical

class Shell(Nebula):
    """
    Shell with exponential falloff around r=1
    """
    @sph_coord
    def rho(self, r, th, phi):
        return 1 if 1 < r < 1.5 else 1e-8


class GaussianShell(Nebula):
    """
    Shell with normal distribution around r=1
    """
    @sph_coord
    def rho(self, r, th, phi):
        return exp(-10 * (r - 1)**2)

class Quadrants(Nebula):
    """
    Constant density in each of the 4 quadrants in x and y
    """
    def rho(self, x, y, z):
        return 1 if x >= 0 and y >= 0 else .7 if x < 0 and y >= 0 else .5 if x < 0 and y < 0 else .2


class Butterfly(Nebula):
    """
    Exponential falloff around the equation z = r**2 in cylindrical coordinates
    """

    white_point = 0.3, 0.3
    black_point = 0, 1.4

    def rho(self, x, y, z):
        return exp(-10 * abs((x**2 + y**2)**2 - abs(z)))


class Disk(Nebula):
    """
    Normal distribution around r=1 and z=0 in cylindrical coordinates
    """
    def rho(self, x, y, z):
        return exp(-10 * (sqrt(x**2 + y**2) - 1)**2 - 10 * z**2)

class Grapher(Nebula):
    """
    Takes a string representing an equation and creates a normally distributed density around that equation.
    The equation can reference r_sph, r_cyl, x, y, and z
    """
    def __init__(self, *args, equation, falloff=1):
        super().__init__(*args)

        left, right = equation.split('=')
        self.eqn = f'({left}) - ({right})'
        self.falloff = falloff

    def rho(self, x, y, z):
        r_sph = r_spherical(x, y, z)
        r_cyl = r_spherical(x, y)
        return exp(-10 * self.falloff * eval(self.eqn)**2)