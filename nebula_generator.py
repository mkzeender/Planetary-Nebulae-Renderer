from PIL import Image
import numpy as np
from numpy import array, ndarray
from functools import wraps
from math import sqrt, sin, cos, acos, atan, pi, exp
from typing import Union, Iterable
from tkinter.filedialog import asksaveasfilename, askdirectory
from os import  path

BOUND = 1.5

def r_spherical(x, y, z=0):
    """Gives the radius of a given x, y, and optionally z coordinate"""
    return sqrt(x**2 + y**2 + z**2)

def sph_coord(method):
    """
    method decorator to convert its parameters to spherical coordinates. used as following:

    @sph_coord
    def method(self, r, th, phi):
        ...

    """
    @wraps(method)
    def wrapper(self, x, y, z):
        phi = (
                atan(y / x) if x > 0 else
                atan(y / x) + pi if x < 0 <= y else
                atan(y / x) - pi if x < 0 and y < 0 else
                pi / 2 if x == 0 < y else
                -pi / 2 if y < 0 == x else
                0 if x == y == 0 else
                None
        )
        return method(self,
                      r = r_spherical(x, y, z),
                      th = acos(z / r_spherical(x, y, z)),
                      phi = phi
                      )
    return wrapper


def constrain(point, data_shape):
    """
    constrain a value into a valid position in the data
    """

    return min(point[0], data_shape.shape[0]), min(point[1], data_shape.shape[1])

def mean(args):
    """Mean of a sequence"""
    return sum(args) / len(args)


class RenderedNebula:
    """
    A 2D rendering of a planetary nebula
    """
    def __init__(self, data: ndarray):
        """
        :param data: the raw intensity data for each pixel
        """
        self.data = data

        self.brightness = 0
        self.contrast = 1

        self.black = 0
        self.white = 255


    def to_image(self) -> Image.Image:
        """
        Converts the data into an image
        :return:
        """
        return Image.fromarray(self.magnitude()).rotate(90, expand=True)

    def magnitude(self) -> ndarray:
        """
        Retrieve the image data in a "usable" form. It will be adjusted according to the brightness, contrast, and
        normalization.

        :return: array of relative magnitudes
        """
        a = np.log(self.data)

        # normalize it
        a -= (self.white + self.black) /2

        a *= 255 / (self.white - self.black)

        # apply brightness and contrast
        a += self.brightness
        a = a * self.contrast + 127
        return a

    def show(self):
        """
        Show the image.
        :return:
        """
        self.to_image().show()

    def save(self, filename = None):
        """
        Save the image as a PNG. Prompts the user if the file is not specified
        :param filename: (str) the name of the file.
        :return:
        """
        if filename is None:
            filename = asksaveasfilename(filetypes=[('PNG files', '*.png')], confirmoverwrite=True)
            if filename == '':
                return

        if filename[-4:].lower() != '.png':
            filename += '.png'

        img = self.to_image().convert('L')
        with open(filename, 'wb') as f:
            img.save(f, format='PNG')


class Nebula:
    white_point = None
    black_point = None

    def __init__(self, width_pix, height_pix, iters=200):

        self.width_pix = width_pix
        self.height_pix = height_pix
        self.iters = iters

        self.white = 255
        self.black = 0

        self._last_progress = None

        self._renders = {}

    def __getitem__(self, item):
        return self._renders[item]

    def __iter__(self):
        for item in self._renders.values():
            item: RenderedNebula
            yield item

    def __str__(self):
        return f'<{self.__class__.__name__} Nebula>'

    def __repr__(self):
        return str(self)

    def show_progress(self, prog):
        prog = int(100 * prog)
        if prog != self._last_progress:
            print(f'{prog}%')
            self._last_progress = prog

    def rho(self, x, y, z):
        return 1

    def j(self, x, y, z):
        return 1 / r_spherical(x, y, z) ** 2

    def _render_ray(self, x0, y0, z0, ds, dx, dy, dz):
        intensity = 0
        x = x0
        y = y0
        z = z0
        s = -1

        while s < 1:
            dI = self.j(x, y, z) * self.rho(x, y, z) * ds
            intensity += dI

            x += dx
            y += dy
            z += dz
            s += ds

        return intensity

    def render(self, th_obs: Union[float, Iterable]=0.0):

        width_pix, height_pix, iters = self.width_pix, self.height_pix, self.iters
        # if th_obs is a list, render for all those values
        try:
            iter(th_obs)
        except TypeError: pass
        else:
            for th_obs_ in th_obs:
                self.render(th_obs_)
            return self

        # iterate from -BOUND to BOUND
        i0 = width_pix  / min(width_pix, height_pix) * BOUND
        j0 = height_pix / min(width_pix, height_pix) * BOUND

        di = dj = 3 / min(width_pix, height_pix)

        # position of the pixels in plane
        dx = di
        dy = dj * sin(th_obs)
        dz = dj * cos(th_obs)

        # direction of the ray
        ds = 2 * BOUND / iters
        dx_ray = 0
        dy_ray = ds * cos(th_obs)
        dz_ray = ds * -sin(th_obs)


        # i, j refer to coordinates of 2d image plane, not to be confused with xyz of 3d plane
        i = -i0
        x = -i0

        m = []
        while i < i0:
            self.show_progress((i+i0) / (2*i0))

            # return to bottom of image
            j = -j0
            y = -BOUND * cos(th_obs) - j0 * sin(th_obs)
            z = BOUND * sin(th_obs) - j0 * cos(th_obs)

            row = []
            while j < j0:
                row.append(self._render_ray(x, y, z, ds, dx_ray, dy_ray, dz_ray))

                j += dj
                y += dy
                z += dz

            #put in reverse order to prevent image flips
            m.append(row)

            i += di
            x += dx

        self.show_progress(1)

        rend = RenderedNebula(array(m))

        rend.black = self.black
        rend.white = self.white

        self._renders[th_obs] = rend
        return self

    def normalize(self, black_point, white_point, th_obs=None):
        if th_obs is None:
            th_obs = next(iter(self._renders))

        rend = self[th_obs]

        a = np.log(rend.data)

        self.black_point = constrain(black_point, a)
        self.white_point = constrain(white_point, a)

        bpx = int(a.shape[0] / 3 * self.black_point[0] + a.shape[0] / 2)
        bpy = int(a.shape[1] / 3 * self.black_point[1] + a.shape[1] / 2)
        wpx = int(a.shape[0] / 3 * self.white_point[0] + a.shape[0] / 2)
        wpy = int(a.shape[1] / 3 * self.white_point[1] + a.shape[1] / 2)

        for render in self._renders.values():
            render.black = a[bpx, bpy]
            render.white = a[wpx, wpy]

        self.white = a[bpx, bpy]
        self.black = a[wpx, wpy]

    @property
    def brightness(self):
        return mean([th.brightness for th in self])

    @brightness.setter
    def brightness(self, val):
        for item in self:
            item.brightness = val

    @property
    def contrast(self):
        return mean([th.brightness for th in self])

    @contrast.setter
    def contrast(self, val):
        for item in self:
            item.contrast = val

    def show(self):
        for item in self:
            item.show()

    def save(self, folder=None, name_prefix=None):
        if folder is None:
            folder = askdirectory(mustexist=True, title='Choose Output Folder')
            if folder == '':
                return

        if name_prefix is None:
            name_prefix = self.__class__.__name__


        for k, v in self._renders.items():
            v: RenderedNebula
            if k == 0:
                current_name = f'{name_prefix} 0'
            else:
                denom = str(round(pi/k, ndigits=2)).strip('0').strip('.')
                current_name = f'{name_prefix} pi_over_{denom}'

            v.save(filename=path.join(folder, current_name))

class Shell(Nebula):
    white_point = 1, 0
    black_point = 1.49, 0

    @sph_coord
    def rho(self, r, th, phi):
        return 1 if 1 < r < 1.5 else 1e-8


class GaussianShell(Nebula):
    @sph_coord
    def rho(self, r, th, phi):
        return exp(-10 * (r - 1)**2)

class Quadrants(Nebula):

    def rho(self, x, y, z):
        return 1 if x >= 0 and y >= 0 else .7 if x < 0 and y >= 0 else .5 if x < 0 and y < 0 else .2


class Butterfly(Nebula):

    white_point = 0.3, 0.3
    black_point = 0, 1.4

    def rho(self, x, y, z):
        return exp(-10 * abs((x**2 + y**2)**2 - abs(z)))


class Disk(Nebula):

    def rho(self, x, y, z):
        return exp(-10 * (sqrt(x**2 + y**2) - 1)**2 - 10 * z**2)

class Grapher(Nebula):
    def __init__(self, *args, equation, falloff=1):
        super().__init__(*args)

        left, right = equation.split('=')
        self.eqn = f'{left} - {right}'
        self.falloff = falloff

    def rho(self, x, y, z):
        r_sph = r_spherical(x, y, z)
        r_cyl = r_spherical(x, y)
        return exp(-10 * self.falloff * eval(self.eqn)**2)


if __name__ == '__main__':
    ...