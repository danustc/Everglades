#!/usr/bin/python

from PyQt5 import QtWidgets,QtCore
import inLib
from functools import partial
from Utilities import QExtensions as qext

class UI(inLib.DeviceUI):
    def __init__(self, control):
        design_path = 'arduino.parallax.parallax_design'
        inLib.DeviceUI.__init__(self, control, design_path)
        self._ui.radioButton_488OD0.toggled.connect(partial(self.set_488, 0))
        self._ui.radioButton_488OD1.toggled.connect(partial(self.set_488, 1))
        self._ui.radioButton_488OD2.toggled.connect(partial(self.set_488, 2))
        self._ui.radioButton_488OD3.toggled.connect(partial(self.set_488, 3))

        self._ui.radioButton_561OD0.toggled.connect(partial(self.set_561, 0))
        self._ui.radioButton_561OD1.toggled.connect(partial(self.set_561, 1))
        self._ui.radioButton_561OD2.toggled.connect(partial(self.set_561, 2))
        self._ui.radioButton_561OD3.toggled.connect(partial(self.set_561, 3))

    def set_488(self, n_OD):
        print("Optical density:", 0)
        self._control.set_OD(488, n_OD)


    def set_561(self, n_OD):
        print("Optical density:", 0)
        self._control.set_OD(561, n_OD)

