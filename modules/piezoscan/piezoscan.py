#!/usr/bin/python


import numpy as np
import time
import inLib


class Control(inLib.Module):

    def __init__(self, control, settings):
        print('Initializing Piezoscan.')
        inLib.Module.__init__(self, control, settings)
        print('Piezoscan initialized.')
        self.active = False

        if settings['ThorlabsMotor'] == True:
            self.useThorlabs = True
        elif settings['Marzhauser'] == True:
            self.useThorlabs = False
            self.useMarzhauser = True


    def calcScanParams(self, start, end, nSteps):
        '''
        Calculates the step size of a scan.
        This function seems useless. ---Dan
        :Parameters:
            *start*: float
                Start of scan in micrometers, relative to current position.
            *end*: float
                End of scan in micrometers, relative to current position.
            *nSteps*: int
                The number of steps.

        :Returns:
            *stepSize*: float
                The step size in nanometers.
        '''
        stepSize = abs(1000.0*(end-start)/(nSteps-1.0))
        return stepSize

    def scan(self, start, end, nSteps, nFrames, filename=None, later_direction = 0.0):
        '''
        Performs a scan
        later_direction: the lateral direction in the translation stage (marzhauser. Default: 0.0, scan in the x direction only.)
        '''
        if self.useThorlabs:
            if start<end:
                up=True
            else:
                up=False
            data = self.scan_thorlabs(nSteps, nFrames, up=up, filename=filename)
        elif self.useMarzhauser:
            print("Use marzhauser.")
            data = self.scan_stage(start, end, nSteps, nFrames, filename=filename, later_direction=later_direction)
        else:
            data = self.scan_piezo(start, end, nSteps, nFrames, filename=filename)
        return data


    def scan_thorlabs(self, nSteps, nFrames, up=True, filename=None):
        self.active = True
        dim = self._control.camera.getDimensions()
        data = np.zeros((nSteps,) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        for i in range(nSteps):
            if self.active:
                if up:
                    self._control.servo.jogUp()
                else:
                    self._control.servo.jogDown()
                time.sleep(4*frame_length)
                for j in range(nFrames):
                    im = self._control.camera.getMostRecentImageNumpy()
                    if im is None:
                        time.sleep(frame_length)
                        im = self._control.camera.getMostRecentImageNumpy()
                    slicesFrames[j] = im
                    time.sleep(frame_length)
                data[i] = np.mean(slicesFrames, axis=0)
            else:
                break
        if self.active and filename:
            print('piezoscan: Saving scan to', filename)
            np.save(filename, data)
            self.active = False
        return data


    def scan_stage(self, start, end, nSteps, nFrames, filename = None, later_direction = 0.0):
        '''
        Performs a motor stage scan. Added by Dan on 11/13/2017.
        '''
        print("Step 1, stage scan: Scanning with parameters:", start, end, nSteps)
        self.active = True
        origin_x, origin_y, origin_z = self._control.stage.position() 
        print("Step 2, Original position:", origin_x, origin_y, origin_z) 
        ds = (end-start)/(nSteps-1.0) # the step size, pos or neg
        dx = ds*np.cos(later_direction)
        dy = ds*np.sin(later_direction)
        rs_x = start*np.cos(later_direction)
        rs_y = start*np.sin(later_direction)
        self._control.stage.goRelative(rs_x, rs_y)
        dim = self._control.camera.getDimensions()
        print("Step 3, Dimension:", dim)
        data = np.zeros((nSteps,) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        print("frame length:", frame_length)
        for ii in range(nSteps):
            if self.active:
                print("step:", ii)
                self._control.stage.goRelative(dx, dy)
                time.sleep(2*frame_length)
                for jj in range(nFrames):
                    slicesFrames[jj] = self._control.camera.getMostRecentImageNumpy()
                    time.sleep(frame_length)
                data[ii] = np.mean(slicesFrames, axis=0)
            else:
                break
        self._control.stage.goAbsolute(origin_x, origin_y)
        if self.active and filename:
            print('piezoscan: Saving scan to', filename)
            np.save(filename, data)
            self.active = False
        return data


    def scan_piezo(self, start, end, nSteps, nFrames, filename=None):
        '''
        Performs a piezo scan.

        :Parameters:
            *start*: float
                The starting point relative to the current position, in micrometers.
            *end*: float
                The end point, relative to the current position, in micrometers.
            *nSteps*: int
                The number of steps.
            *nFrames*: int
                The number of frames to be averaged in each step.
            *filename*: str
                If not *None*, the image data will be saved in a .npy file with this name.
        '''
        print('piezoscan: Scanning with params:', start, end, nSteps)
        self.active = True
        dz = abs((end-start)/(nSteps-1.0))
        orig_z = self._control.piezo.getPosition(3)
        start += orig_z
        end += orig_z
        zs = np.linspace(start, end, nSteps)
        dim = self._control.camera.getDimensions()
        data = np.zeros((nSteps,) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        for i in range(nSteps):
            if self.active:
                self._control.piezo.moveTo(3, zs[i])
                time.sleep(2*frame_length)
                for j in range(nFrames):
                    slicesFrames[j] = self._control.camera.getMostRecentImageNumpy()
                    time.sleep(frame_length)
                data[i] = np.mean(slicesFrames, axis=0)
            else:
                break
        self._control.piezo.moveTo(3, orig_z)
        if self.active and filename:
            print('piezoscan: Saving scan to', filename)
            np.save(filename, data)
            self.active = False
        return data

    def bl_correct(self, nsteps = 31, stepsize = 0.3, z_correct = 3.0, z_start = None):
        '''
        Backlash correct function. Updated on 12/09/2016.
        nsteps: the planned scan steps in the forward direction
        stepsize: the input stepsize
        z_start: supposed start point (in the forward scan)
        '''
        N_correct = int(z_correct/stepsize)
        if z_start is None:
            '''
            If the start point coordinate is not specified
            '''
            for ii in range(nsteps+N_correct):
                self._control.servo.jogUp()
            time.sleep(3)

        else:
            self._control.moveTo(z_start+z_correct)
            time.sleep(2)
        # now, move the piezo back
        for ii in range(N_correct):
            self._control.servo.jogDown()
            time.sleep(0.25)
        # done with bl_correct


    def stop(self):
        self.active = False
