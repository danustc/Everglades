#!/usr/bin/python


import inLib
import numpy as np
from . import signalForAO
from .deskew import deskew_stack
from .PR_core import Core
import os

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
        self._modulations = []
        print('AO imaging module initialized.')

        self.hasSLM = settings['hasSLM']
        self.hasMirror = settings['hasMirror']
        self.scan_device = settings['scan_device']

        self._PSF = None
        self._sharpness = None


        self._center = []
        self._PRcore = Core()

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
        self._PRcore.PSF = PSF # reset the psf
        if filename:
            np.save(filename, PSF)
        return PSF
    
    

    def getSharpness(self, pixelSize=103, diffLimit=800):
        '''
        Finds the sharpnes of self._PSF
        '''
        if self._PSF is None:
            print("No psf acquired.")
            return
        sharpness = signalForAO.secondMomentOnStack(self._PSF,pixelSize,diffLimit)
        self._sharpness = sharpness
        print("Maximum sharpness = ", sharpness.max())
        return sharpness


    def retrievePF(self, dx, l, n, NA, f, guess, nIt, neglect_defocus=True,
                   invert=False, wavelengths=1, wavestep = 0.005, resetAmp=True,
                   symmeterize=False):
        '''
        Updated by Dan on 01/15/2019, replacing the computational part with PR_core
        '''
        self._PRcore.lcenter = l
        self._PRcore.pxl = dx
        self._PRcore.nfrac = n
        self._PRcore.NA = NA
        self._PRcore.objf = f # set all the properties
        self._PRcore.n_wave = wavelengths
        self._PRcore.d_wave = wavestep
        self._PRcore.pupil_Simulation() # simulate a pupil
        self._PRcore.retrievePF()
        
       

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


    def modulateMirror(self):
        '''
        Apply to mirror
        1. Synthesize the pattern 
        
        '''
        pass
