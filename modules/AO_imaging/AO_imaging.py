import numpy as np
import pyfftw
from scipy import ndimage


class Pupil(object):
    def __init__(self, wl, nfrac, NA, f, ds = 170):
        '''
        wl: wavelength
        nfrac: refractive index
        NA: numerical aperture
        ds: the thickness of the coverslip
        '''
        self.wave_len = wl
        self.nfrac = nfrac
        self.NA = NA
        self.foc_len = f
        self.k_max = NA/l



class Simulation(Pupil):
    def __init__(self, nr = 256, dx = 0.102, wl = 0.515, nfrac = 1.33,NA = 1.00, f = 9000, nwl = 10, w_step = 0.005):
       self.nx = nr
       self.ny = nr
       Pupil.__init__(self, wl, nfrac, NA, f)

