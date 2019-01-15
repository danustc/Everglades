#!/usr/bin/python


import inLib
import numpy as np
from myWidget import pupil
from scipy.ndimage import gaussian_filter as gf
from scipy import optimize
import os
import time
import libtim.zern
from . import signalForAO
from skimage.restoration import unwrap_phase
from .deskew import deskew_stack
# Constant parameter

PXL = 0.097
WL_fluor = 0.550
NA_new = 1.0
N_Refrac = 1.33
Focal = 9000
GS = 'plane'




class Control(inLib.Module):

    def __init__(self, control, settings):
        print('Initializing AO imaging module.') 
        inLib.Module.__init__(self, control, settings)
        dim = self._control.camera.getDimensions()
        # A list to store the indices of the 'Other' modulations, given by the SLM API:
        self._modulations = []
        print('AO imaging module initialized.')

        self.hasSLM = settings['hasSLM']
        self.hasMirror = settings['hasMirror']
        self.scan_device = settings['scan_device']

        self._PSF = None
        self._sharpness = None

        self.sharpnessList = []

        self.zernModesToFit = 20 #Number of zernike modes to fit unwrapped pupil functions

        self._zernFitUnwrapped = None
        self._zernFitUnwrappedModes = None
        self.varyAOactive = True


        self._center = []


    def updateImSize(self):
        '''
        Gets the image size from the camera.
        _ui uses this for sharpness calculations...
        '''
        return self._control.camera.getDimensions()
    
    
    def preview(self):
        self._control.camera.beginPreview()

    def getImageForPreview(self):
        return self._control.camera.getImageForPreview()


    def acquirePSF(self, range_, nSlices, nFrames, center_xy=True, filename=None,
                   mask_size = 40, mask_center = (-1,-1), deskew = True):
        '''
        Acquires a PSF stack. The PSF is returned but also stored internally.

        :Parameters:
            *range_*: float
                The range of the scan around the current axial position in micrometers.
            *nSlices*: int
                The number of PSF slices to acquire.
            *nFrames*: int
                The number of frames to be averaged for each PSF slice.
            *filename*: str
                The file name into which the PSF will be saved.

        :Returns:
            *PSF*: numpy.array
                An array of shape (k,l,m), where k are the number of PSF slices and
                (l,m) the lateral slice dimensions.
        '''

        # Logging
        self._settings['range'] = range_
        self._settings['nSlices'] = nSlices
        self._settings['nFrames'] = nFrames
        self._settings['filename'] = filename

        self.filename = filename

        # Some parameters
        start = range_/2.0
        end = -range_/2.0
        if self._settings['scan_device'] == 'marzhauser':
            start,end = end, start # the direction of the thorlab motor and the marzhauser stage are opposite
        self._dz = abs(range_/(nSlices-1.0))

        # Scan the PSF:
        origin_x, origin_y, origin_z = self._control.stage.position() 
        print("Positions:", origin_x, origin_y)
        scan = self._control.piezoscan.scan(start, end, nSlices, nFrames, filename)

        if deskew:
            scan = deskew_stack(scan, range_)

        nz, ny, nx = scan.shape
        # An empty PSF
        PSF = np.zeros_like(scan)
        g = pupil.Geometry((ny,nx), ny/2.-0.5, nx/2.-0.5, 16)
        # cyl = np.array(nz*[g.r_pxl<16])
        new_cyl = np.array(nz*[g.r_pxl<mask_size])
        # Hollow cylinder
        hcyl = np.array(nz*[np.logical_and(g.r_pxl>=50, g.r_pxl<61)])
        if mask_center[0] > -1:
            g2 = pupil.Geometry((ny,nx), mask_center[0], mask_center[1], 16)
            mask = g2.r_pxl<mask_size
        else:
            mask = g.r_pxl<mask_size
        if filename:
            np.save(filename+"pre-cut", scan)
        if center_xy:
            # The coordinates of the brightest pixel
            cz, cy, cx = np.unravel_index(np.argmax(gf(scan*mask,2)), scan.shape)
            print("Center found at: ", (cz,cy,cx))
            print("Mask size:", mask_size)
            self._center = [cz, cy, cx] # save this for one-run procedure
            # We laterally center the scan at the brightest pixel

            mid_y = int(ny/2)
            mid_x = int(nx/2)
            cut = scan[:,cy-mask_size:cy+mask_size,cx-mask_size:cx+mask_size]
            PSF = np.zeros((nz,ny,nx))
            PSF[:,mid_y-mask_size:mid_y+mask_size,mid_x-mask_size:mid_x+mask_size] = cut
            # Background estimation
            self._background = np.mean(scan[hcyl])
            print("Background guess: ", self._background)
        else:
            self._background = np.mean(scan[hcyl])
        PSF[np.logical_not(new_cyl)] = self._background

        self._PSF = PSF
        if filename:
            np.save(filename, PSF)
        return PSF

    def getSharpness(self, pixelSize=163, diffLimit=800):
        '''
        Finds the sharpnes of self._PSF
        '''
        if self._PSF is None:
            time.sleep(2)
        sharpness = signalForAO.secondMomentOnStack(self._PSF,pixelSize,diffLimit)
        self._sharpness = sharpness
        print("Maximum sharpness = ", sharpness.max())
        return sharpness


    def retrievePF(self, dx, l, n, NA, f, guess, nIt, neglect_defocus=True,
                   invert=False, wavelengths=1, resetAmp=True,
                   symmeterize=False):
        '''
        Retrieves the phase of the most recent PSF. The retrieved pupil function is both
        returned and stored internally.

        :Parameters:
            *dx*: float
                The lateral pixel size in micrometers.
            *l*: float
                Emission wavelength in micrometers.
            *n*: float
                Immersion medium refractive index.
            *NA*: float
                Objective numerical aperture.
            *f*: float
                Objective focal length in micrometers.
            *guess*: [('plane',) | ('mirror', z0) | ('file', path)]
                If *('plane',)*, the initial guess is a plane wave.
                If *('mirror',z0)*, the initial guess is the pupil function of an emitter
                with distance *z0* (in micrometers) to a mirror.
                If *('file', path)*, the initial guess is loaded from a file given by *path*.
                The file at *path* has to be a numpy array stored in the .npy format and contain
                the guess in radiants.
            *nIt*: int
                Number of iterations of the phase retrieval algorithms to be performed.
            *gaussian*: float
                If not *None*, a Gaussian filter is applied to the phase after each iteration.
                The standard deviation of the Gaussian kernel is given with this parameter.

        :Returns:
            *PF*: numpy.array
                The complex pupil function. The shape is the same as the shape of the acquired
                PSF slices.
        '''

        # Logging
        self._settings['nIterations'] = nIt

        z_offset = 0
        if neglect_defocus:
            # We try to estimate the axial position of the emitter
            # The coordinates of the brightest pixel
            cz, cx, cy = np.unravel_index(self._PSF.argmax(), self._PSF.shape)
            # Intensity trace along z
            i = self._PSF[:,cx,cy]
            # Get z positions
            nz = self._PSF.shape[0]
            upper = 0.5*(nz-1)*self._dz
            z = np.linspace(-upper, upper, nz)
            # Initial fit parameters
            b = np.mean((i[0],i[-1]))
            a = i.max() - b
            w = l/3.2
            p0 = (a,0,w,b)
            def gaussian(z, a, z0, w, b):
                return a * np.exp(-(z-z0)**2/w) + b
            # Fit gaussian to axial intensity trace
            popt, pcov = optimize.curve_fit(gaussian, z, i, p0)
            # Where we think the emitter is axially located:
#             z_offset = -1.0*popt[1] #Added on April 3, 2015
            z_offset = popt[1] # edited on Aug01, 2016

        nx,ny = self._PSF.shape[1:3]
        self._pupil = pupil.Simulation(nx,dx,l,n,NA,f,wavelengths=wavelengths)
        if guess[0] == 'plane':
            A = self._pupil.plane
        elif guess[0] == 'mirror':
            z0 = guess[1]
            A = np.angle(self._pupil.get_sli_pupil_function(z0,0.0))
            A = np.nan_to_num(A)
        elif guess[0] == 'file':
            A = np.load(guess[1])
        else:
            A = self._pupil.plane

        print("Finding PF...")
        print("   Using parameters:")
        print("   dz = ", self._dz)
        print("   background = ", self._background)
        print("   z_offset = ", z_offset)
        complex_PF = self._pupil.psf2pf(self._PSF, self._dz, self._background, A, nIt, z_offset,
                                        resetAmp=resetAmp,symmeterize=symmeterize)

        if invert:
            complex_PF = abs(complex_PF) * np.exp(-1*1j*np.angle(complex_PF))

        self._PF = _PupilFunction(complex_PF, self._pupil)
        self._PF.phase = self.unwrap()
        nz = self._PSF.shape[0]
        upper = 0.5*(nz-1)*self._dz
        z = np.linspace(-upper, upper, nz)
        psft = self._pupil.pf2psf(complex_PF, z)
        np.save('retrieved_psf', psft)

        # added on 08/25, automatically save the PF
        self.savePF(self.filename)
        Pupil_final = self._PF
        self.pf_complex = Pupil_final.complex
        self.pf_phase = unwrap_phase(Pupil_final.phase)
        self.pf_ampli = Pupil_final.amplitude
        sth = self.Strehl_ratio()
        print(("Strehl ratio:", sth))

        return self._PF.phase


    def fit(self):
        '''
        Fits the most recent pupil function to Zernike modes up to the fourth order.

        :Returns:
            *p*: numpy.array
                A 1-dimensional array of Zernike coefficients. Array indices correspond to
                the Noll coefficient *j* of Zernike modes.
            *PF*: numpy.array
                The phase pupil function as produced by the phase retrieval algorithm with Zernike
                coefficients.
        '''
        ''' Old fitting procedure
        p, success, error = zernike.fit_to_basic_set(self._PF.phase, np.zeros(15),
                self._pupil.r, self._pupil.theta)
        self._PF.zernike_coefficients = p
        '''
        nx,ny = self._PF.phase.shape
        radius = np.sum(self._pupil.r[nx/2,:]<1)/2
        print("Radius of fitted zernike: ", radius)
        #fitResults = libtim.zern.fit_zernike(self._PF.phase, rad=radius, nmodes=15)
        unwrapped=self.unwrap()
        print("Unwrapped!")

#        fitResults_wrong = libtim.zern.fit_zernike(self._PF.phase, rad=radius, nmodes=25)
        fitResults = libtim.zern.fit_zernike(unwrapped, rad=radius, nmodes=25)
        self._PF.zernike_coefficients = fitResults[0]

        return self._PF

    def unwrap(self):
        unwrapped = unwrap_phase(self._PF.phase)
        print("The unwrapped phase has a size of :")
        print(unwrapped.shape)

        return unwrapped

    def setZernModesToFit(self, nmodes):
        self.zernModesToFit = nmodes

    

    def removePTTD(self):
        p = self._PF.zernike_coefficients
        p[0:4] = 0
        self._PF.zernike_coefficients = p
        #self._PF.zernike_coefficients[0:4] = 0
        return self._PF


    def savePF(self, filename):
        '''
        Saves the most recent pupil function in InControl's working directory.

        :Parameters:
            *filename*: str
        '''
        working_dir = self._control.getWorkingDir()
        np.save(os.path.join(working_dir, filename + '_complex.npy'), self._PF.complex)
        np.save(os.path.join(working_dir, filename + '_amplitude.npy'), self._PF.amplitude)
        np.save(os.path.join(working_dir, filename + '_phase.npy'), self._PF.phase)


    
    def findSharpnessEachFrame(self, pixelSize, diffLimit, mask=None):
        frame_length = 1.0/self._control.camera.getFrameRate()
        #np_image = self._control._control.camera.getImageForPreview()
        im = self._control.camera.getMostRecentImageNumpy()
        if mask is not None:
            background_est = np.median(im * np.logical_not(mask))
            im = (im * mask) + (background_est * np.logical_not(mask))
        if im is not None:
            sharpness = signalForAO.secondMoment(im,pixelSize,diffLimit)
            if len(self.sharpnessList)<200:
                self.sharpnessList.append(sharpness)
            else:
                self.sharpnessList = self.sharpnessList[1:] + [sharpness]
        else:
            sharpness = None
        time.sleep(frame_length)
        return sharpness, self.sharpnessList





class _PupilFunction(object):
    '''
    A pupil function that keeps track when when either complex or amplitude/phase
    representation is changed.
    '''
    def __init__(self, complex, geometry):
        self.complex = complex
        self._geometry = geometry

    @property
    def complex(self):
        return self._complex

    @complex.setter
    def complex(self, new):
        self._complex = new
        self._amplitude = abs(new)
        self._phase = np.angle(new)

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, new):
        self._amplitude = new
        self._complex = new * np.exp(1j*self._phase)

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, new):
        self._phase = new
        self._complex = self._amplitude * np.exp(1j*new)

    @property
    def zernike_coefficients(self):
        return self._zernike_coefficients

    @zernike_coefficients.setter
    def zernike_coefficients(self, new):
        self._zernike_coefficients = new
        #self._zernike = zernike.basic_set(new, self._geometry.r, self._geometry.theta)
        self._zernike = libtim.zern.calc_zernike(new, self._geometry.nx/2.0)

    @property
    def zernike(self):
        return self._zernike
