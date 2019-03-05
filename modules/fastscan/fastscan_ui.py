#!/usr/bin/python

import inLib
from PyQt5 import QtCore, QtGui
import sys
from functools import partial
import time
import numpy as np
import os
from .Dev_threads import stage_mover, TL_acquisition


def clickable(widget):
    class Filter(QtCore.QObject):
        clicked = QtCore.pyqtSignal(int,int)
        def eventFilter(self,obj,event):
            if obj == widget:
                if event.type() == QtCore.QEvent.MouseButtonRelease:
                    buttonState = event.button()
                    self.clicked.emit(event.x(),event.y())
                    return True
            return False
    filter=Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked

class UI(inLib.ModuleUI):

    def __init__(self, control, ui_control):
        design_path = 'modules.fastscan.fastscan_design'
        inLib.ModuleUI__init__(self, control, ui_control, design_path)

        self.stage_thread = None
        self.cam_thread = None
        self.exp_time = 50.0 # miliseconds


    def quick_scan(self):
        self.stage_thread = stage_mover(self._control, self.destination)
        self.cam_thread = TL_acquisition(self._control, self.exp_time)

        self._ui.pushButton_startmv.setEnabled(False)
        self.stage_thread.start()
        self.cam_thread.start()


    def shutDown(self):
        if self.stage_thread:
            self.stage_thread.wait()
        if self.cam_thread:
            self.cam_thread.wait()
