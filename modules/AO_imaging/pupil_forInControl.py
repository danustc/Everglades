#!/usr/bin/python

import numpy as _np
from scipy import fftpack as _fftpack
from scipy import ndimage
from numpy.lib.scimath import sqrt as _msqrt
import tempfile as _tempfile
#from Utilities import zernike as _zernike
from matplotlib.pylab import *
import pyfftw


class Pupil(object):

    def __init__(self, l, n, NA, f, d=170):

        self.l = float(l)
        self.n = float(n)
        self.f = float(f)
        self.NA = NA
        self.s_max = f*NA
        self.k_max = NA/l
        self.d = float(d)


    def unit_disk_to_spatial_radial_coordinate(self, unit_disk):

        return self.s_max * unit_disk


    def spatial_radial_coordinate_to_optical_angle(self, radial_coordinate):

        alphac = _np.arcsin(radial_coordinate/(self.n*self.f) + 1j*0)
        return _np.real(alphac) - _np.imag(alphac)


    def spatial_radial_coordinate_to_xy_slope(self, s):

        _mc = _msqrt((self.n*self.f/self.s)**2-1)
        return _np.real(_mc) - _np.imag(_mc)


    def get_pupil_function(self, z, n_photons=1000, coverslip_tilt=0, coverslip_tilt_direction=0):

        '''
        Computes the complex pupil function of single point source.

        Parameters
        ----------
        z:
            Axial coordinate of the point source, where the origin of the
            coordinate system is at focus. Positive values are towards the
            objective.
	n_photons: int
	    The number of collected photons.

        Returns
        -------
        PF: 2D array or list of 2D arrays
            If parameter *z* was a single number, a 2D array of the complex
            pupil function is returned. If *z* was an iterable of floats, a list
            of 2D pupil functions is returned.
        '''

        f = self.f
        m = self.m
        l = self.l
        n = self.n
        a = self.alpha
        t = self.theta
        d = self.d
        e = coverslip_tilt
        g = 2*_np.pi*coverslip_tilt_direction/360.
        ng = 1.5255 # Coverslip refractive index


    def apply_NA_restriction(self, PF):
        '''
        Sets all values of a given pupil function to zero that are out of NA range.
        '''
        PF[self.r>1] = 0
        return PF
    
    def compute_Fisher_Information(self, PSF, poisson_noise, voxel_size,
            mask=None):

        '''
        Computes the Fisher information matrix as described in Ober et al.,
        Biophysical Journal, 2004, taking into account Poisson background noise.

        Parameters
        ----------
        PSF: array
            The 3D Point Spread Function, where the intensity values are the
	    number of photons.
        poisson_noise: float
            The poisson background noise per pixel in photons
        voxel_size: float
            The side length of the PSF voxel in micrometer.
        mask: 2D boolean numpy.array
            Optional: A two dimensional array of the same shape as one PSF
            slice. If mask is set, the Fisher Information is only calculated for
            regions where mask is True.
        '''

        dPSF = _np.gradient(PSF, voxel_size[0], voxel_size[1], voxel_size[2])
        if mask is not None:
            dPSF = [_*mask for _ in dPSF]
            noisy = PSF + poisson_noise

        FI = [[0,0,0] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                FI[i][j] = _np.sum(dPSF[i]*dPSF[j]/noisy,-1).sum(-1)

        return FI


    def compute_CRLB(self, PSF, poisson_noise, voxel_size, mask=None):

        FI = self.compute_Fisher_Information(PSF, poisson_noise, voxel_size,
                mask)
        return [_np.sqrt(1.0/FI[_][_]) for _ in range(3)]



class Simulation(Pupil):

    '''
    Simulates the behaviour of a microscope based on Fourier optics.

    Parameters
    ----------
    nx: int
        The side length of the pupil function or microscope image in pixels.
    dx: float
        The pixel size in the image plane.
    l: float
        Light wavelength in micrometer.
    n: float
        The refractive index of immersion and sample medium. Mismatching media
        are currently not supported.
    NA: float
        The numerical aperture of the microscope objective.
    f: float
        The objective focal length in micrometer.
    '''

    def __init__(self, nx=256, dx=0.1, l=0.68, n=1.33, NA=1.27, f=9000, wavelengths=10, wavelength_halfmax=0.005):

        dx = float(dx)
        self.dx = dx
        l = float(l)
        n = float(n)
        NA = float(NA)
        f = float(f)
        self.nx = nx
        self.ny = nx

        Pupil.__init__(self, l, n, NA, f)

        self.numWavelengths = wavelengths

        # Frequency sampling:
        dk = 1/(nx*dx)
        # Pupil function pixel grid:
        x,y = _np.mgrid[-nx/2.:nx/2.,-nx/2.:nx/2.]+0.5
        self.x_pxl = x
        self.y_pxl = y
        self.r_pxl = _msqrt(x**2+y**2)

        # Pupil function frequency space:
        kx = dk*x
        ky = dk*y
        self.k = _msqrt(kx**2+ky**2)

        # Axial Fourier space coordinate:
        self.kz = _msqrt((n/l)**2-self.k**2)
        self.kzs = _np.zeros((self.numWavelengths,self.kz.shape[0],self.kz.shape[1]),dtype=self.kz.dtype)
        ls = _np.linspace(l-wavelength_halfmax,l+wavelength_halfmax,self.kzs.shape[0])
        for i in range(0,self.kzs.shape[0]):
            self.kzs[i] = _msqrt((n/ls[i])**2-self.k**2)

        # Scaled pupil function radial coordinate:
        self.r = self.k/self.k_max

        self.s = self.unit_disk_to_spatial_radial_coordinate(self.r)
        self.alpha = self.spatial_radial_coordinate_to_optical_angle(self.s)
        self.m = self.spatial_radial_coordinate_to_xy_slope(self.alpha)

        # Plane wave:
        self.plane = _np.ones((nx,nx))+1j*_np.zeros((nx,nx))
        self.plane[self.k>self.k_max] = 0
        self.pupil_npxl = abs(self.plane.sum())

        self.kx = kx
        self.theta = _np.arctan2(y,x)


    def pf2psf(self, PF, zs, intensity=True, verbose=False, use_pyfftw=True):

        """
        Computes the point spread function for a given pupil function.

        Parameters
        ----------
        PF: array
            The complex pupil function.
        zs: number or iteratable
            The axial position or a list of axial positions which should be
            computed. Focus is at z=0.
        intensity: bool
            Specifies if the intensity or the complex field should be returned.

        Returns
        -------
        PSF: array or memmap
            The complex PSF. If the memory is to small, a memmap will be
            returned instead of an array.
        """

        nx = self.nx

        if _np.isscalar(zs):
            zs = [zs]
        nz = len(zs)
        kz = self.kz

        # The normalization for ifft2:
        N = _np.sqrt(nx*self.ny)

        # Preallocating memory for PSF:
        try:
            if intensity:
                PSF = _np.zeros((nz,nx,nx))
            else:
                PSF = _np.zeros((nz,nx,nx))+1j*_np.zeros((nz,nx,nx))
        except MemoryError:
            print('Not enough memory for PSF, \
                    using memory map in a temporary file.')
            temp_file = _tempfile.TemporaryFile()
            if intensity:
                temp_type = float
            else:
                temp_type = complex
            PSF = _np.memmap(temp_file, dtype=temp_type, mode='w+',
                shape=(nz,nx,nx))

        for i in range(nz):
            if verbose: print('Calculating PSF slice for z={0}um.'.format(zs[i]))
            if use_pyfftw:
                aligned = pyfftw.n_byte_align(_np.exp(2*_np.pi*1j*kz*zs[i])*PF,16)
                U = N * pyfftw.interfaces.numpy_fft.ifft2(aligned)
            else:
                U = N*_fftpack.ifft2(_np.exp(2*_np.pi*1j*kz*zs[i])*PF)
            for j in range(0,self.kzs.shape[0]):
                if use_pyfftw:
                    aligned = pyfftw.n_byte_align(_np.exp(2*_np.pi*1j*self.kzs[j]*zs[i])*PF,16)
                    U = U + N*pyfftw.interfaces.numpy_fft.ifft2(aligned)
                else:
                    U = U + N*_fftpack.ifft2(_np.exp(2*_np.pi*1j*self.kzs[j]*zs[i])*PF)
            U = U/(1+self.kzs.shape[0])
            _slice_ = _fftpack.ifftshift(U)
            if intensity:
                _slice_ = _np.abs(_slice_)**2 # intensity
            PSF[i] = _slice_

        if nz == 1:
            PSF = PSF[0]

        return PSF




    def psf2pf(self, PSF, dz, mu, A, nIterations=10, z_offset=0, use_pyfftw=True, resetAmp=True,
               symmeterize=False):

        '''
        Retrieves the complex pupil function from an intensity-only
        PSF stack by relative entropy minimization. The algorithm is
        based on Kner et al., 2010, doi:10.1117/12.840943, which in turn
        is based on Deming, 2007, J Opt Soc Am A, Vol 24, No 11, p.3666.

        Parameters
        ---------
        PSF: 3D numpy.array
            An intensity PSF stack. PSF.shape has to be
            (nz, Simulation.nx, Simulation.nx), where nz is the arbitrary
            number of z slices.
        dz: float
            The distance between two PSF slices.
        mu: float
            The noise level of the PSF.
        A: 2D numpy.array
            The initial guess for the complex pupil function with shape
            (Simulation.nx, Simulation.nx).
        '''

        # z spacing:
        dz = float(dz)
        # Number of z slices:
        nz = PSF.shape[0]
        # Noise level:
        mu = float(mu)

        kz = self.kz
        k = self.k
        k_max = self.k_max

        # Z position of slices:
        upper = 0.5*(nz-1)*dz
        zs = _np.linspace(-upper,upper,nz) - z_offset

        # Normalization for fft2:
        N = _np.sqrt(self.nx*self.ny)

        if use_pyfftw:
            pyfftw.interfaces.cache.enable()

        Ue = _np.ones_like(PSF).astype(_np.complex128)
        U = _np.ones_like(PSF).astype(_np.complex128)
        Uconj = _np.ones_like(PSF).astype(_np.complex128)
        Ic = _np.ones_like(PSF).astype(_np.complex128)

        expr1 = "Ue = (PSF/Ic)*U"
        expr2 = "Ic = mu + (U * Uconj)"


        for ii in range(nIterations):

            print('Iteration',ii+1)
            # Calculate PSF field from given PF:
            U = self.pf2psf(A, zs, intensity=False)
            # Calculated PSF intensity with noise:
            Uconj = _np.conj(U)
            #weave.blitz(expr2)
            Ic = mu + (U * Uconj)

            minFunc = _np.mean(PSF*_np.log(PSF/Ic))
            print('Relative entropy per pixel:', minFunc)
            redChiSq = mean((PSF-Ic)**2)
            print('Reduced Chi square:', redChiSq)

            # Comparing measured with calculated PSF by entropy minimization:
            Ue = (PSF/Ic)*U
            #weave.blitz(expr1)
            # New PF guess:
            A = _np.zeros_like(Ue) + 1j*_np.zeros_like(Ue)
            for i in range(len(zs)):
                #Ue[i] = _fftpack.fftshift(Ue[i])
                if use_pyfftw:
                    Ue_aligned = pyfftw.n_byte_align(_fftpack.fftshift(Ue[i]),16)
                    fted_ue = pyfftw.interfaces.numpy_fft.fft2(Ue_aligned)
                    A[i] = fted_ue/_np.exp(2*_np.pi*1j*kz*zs[i])/N
                else:
                    fted_ue = _fftpack.fft2(_fftpack.fftshift(Ue[i]))
                    A[i] = fted_ue/_np.exp(2*_np.pi*1j*kz*zs[i])/N
                for j in range(0,self.kzs.shape[0]):
                    A[i] = A[i] + fted_ue/_np.exp(2*_np.pi*1j*self.kzs[j]*zs[i])/N
                A[i] = A[i]/(1+self.kzs.shape[0])
            A = _np.mean(A,axis=0)

            #mean(abs(A))*_np.exp(1j*_np.angle(A))

            # NA restriction:
            A[k>k_max] = 0
            if resetAmp:
                amp = ndimage.gaussian_filter(np.abs(A),15)
                A = amp*np.nan_to_num(A/np.abs(A))

            if symmeterize:
                if ii>(nIterations/2):
                    A = 0.5*(A+_np.flipud(A))

                #counts = sum(abs(A))/self.pupil_npxl
                #A = counts*_np.exp(1j*angle(A))
                #A[k>k_max] = 0

        return A
