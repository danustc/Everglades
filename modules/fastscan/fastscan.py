#!/usr/bin/python

import inLib
import numpy as np

class Control(inLib.Module):

    def __init__(self, control, settings):
        print('Initializing fastscan module.')
        inLib.Module.__init__(self, control, settings)

        self.cam_thread = None
        self.stage_thread = None


    def updateImSize(self):
        '''
        Gets the image size from the camera.
        _ui uses this for sharpness calculations...
        '''
        return self._control.camera.getDimensions()


    def moveTo(self, dest, speed = None):
        '''
        move the stage to a certain position
        '''
        dx, dy = dest
        print("Stage moved!")
        if speed is not None:
            self._control.stage.setSpeed('x', speed)

        self._control.stage.goAbsolute(dx, dy)


