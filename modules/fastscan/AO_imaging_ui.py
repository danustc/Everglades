#!/usr/bin/python


from PyQt5 import QtWidgets,QtCore, QtGui
import inLib
import numpy as np
import time
from myWidget import QExtensions as qext

def clickable(widget):
    class Filter(QtCore.QObject):
        clicked = QtCore.pyqtSignal(int,int)
        def eventFilter(self,obj,event):
            if obj == widget:
                if event.type() == QtCore.QEvent.MouseButtonRelease:
                    #buttonState = event.button()
                    self.clicked.emit(event.x(),event.y())
                    return True
            return False
    filter=Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked



class UI(inLib.ModuleUI):

    def __init__(self, control, ui_control):
        design_path = 'modules.fastscan.fastscan_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)


        self._ui.doubleSpinBoxRange.setValue(self._control._settings['range'])
        self._ui.spinBoxSlices.setValue(self._control._settings['nSlices'])
        self._ui.spinBoxFrames.setValue(self._control._settings['nFrames'])
        self._ui.lineEdit_imagedest.setText(self._control._settings['filename'])

        self._ui.spinBoxIterations.setValue(self._control._settings['nIterations'])
        self._ui.pushButton_Acquire.clicked.connect(self.acquirePSF)
        self._ui.pushButtonPhase.clicked.connect(self.retrievePF)
        self._ui.pushButton_save.clicked.connect(self.savePF)
        self._ui.pushButton_addmod.clicked.connect(self.addModulate)
        self._ui.pushButton_delmod.clicked.connect(self.removeModulate)
        self._ui.pushButton_ADM.clicked.connect(self.applyToMirror)
        self._ui.pushButton_saveconfig.clicked.connect(self.saveConfig)


        self._ui.labelDisplay.paintEvent = self._labelDisplay_paintEvent
        clickable(self._ui.labelDisplay).connect(self._mousePressEvent)



        self.zernRadius = 0

        self._modulations = []
        self.use_zernike = False
        self.remove_PTTD = True
        self._pixmap = None
        self._scanner = None

        self._autoscale = True
        self.imsize = (256,256)
        self.pixelSize = 103
        self.diffLimit = 800
        self._cmap = None
        self._xres, self._yres = self.imsize


        self._ui.tabWidget_viewer.setEnabled(True)

        self.hasSLM = self._control.hasSLM
        self.hasMirror = self._control.hasMirror
        self._control.preview()
        self._updater = QtCore.QTimer()
        self._updater.timeout.connect(self._update)
        self._updater.start(50)

    def mirror_segs(self):
        '''
        pass the pattern to the DM and convert it into segments.
        save the mirror segs as a temporary file
        '''
        self._control.mirror_segs()


    def displayPhase(self, phase):
        self._ui.mplwidget_phase.figure.axes[0].get_xaxis().set_visible(False)
        self._ui.mplwidget_phase.figure.axes[0].get_yaxis().set_visible(False)
        self._ui.mplwidget_phase.figure.axes[0].matshow(phase, cmap='RdBu')
        self._ui.mplwidget_phase.draw()

    def displayAmpli(self, ampli):
        self._ui.mplwidget_ampli.figure.axes[0].get_xaxis().set_visible(False)
        self._ui.mplwidget_ampli.figure.axes[0].get_yaxis().set_visible(False)
        self._ui.mplwidget_ampli.figure.axes[0].matshow(ampli, cmap='RdBu')
        self._ui.mplwidget_ampli.draw()


    def _updateImSize(self):
        self.imsize = self._control.updateImSize()


    def acquirePSF(self):
        self._updater.stop()
        range_ = self._ui.doubleSpinBoxRange.value()
        nSlices = self._ui.spinBoxSlices.value()
        nFrames = self._ui.spinBoxFrames.value()
        center_xy = self._ui.checkBoxCenterLateral.isChecked()
        maskRadius = self._ui.spinBox_maskRadius.value()
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        imagedest = str(self._ui.lineEdit_imagedest.text())
        self._scanner = Scanner(self._control, range_, nSlices, nFrames, center_xy,
                                imagedest, maskRadius, (cX,cY))
        self._scanner.finished.connect(self._on_scan_done)
        self._ui.pushButton_Acquire.setEnabled(False)
        self._scanner.start()
        # if self._scanner.finish():



    def _on_scan_done(self):
        self._ui.pushButton_Acquire.setEnabled(True)
        self._ui.groupBoxPhase.setEnabled(True)
        time.sleep(1)
        self._updater.start(50)
 
    def retrievePF(self):
        '''
        perform phase retrieval
        '''
        pxlSize = self._ui.doubleSpinBoxPixel.value()
        l = self._ui.doubleSpinBoxWavelength.value()
        n = self._ui.doubleSpinBoxIndex.value()
        NA = self._ui.doubleSpinBoxNA.value()
        f = self._ui.doubleSpinBoxFocal.value()
        numWaves = self._ui.spinBox_numWavelengths.value()
        dWave = self._ui.doubleSpinBox_dwave.value()
        neglect_defocus = self._ui.checkBoxNeglectDefocus.isChecked()
        nIt = self._ui.spinBoxIterations.value()
        invertPF = self._ui.checkBox_invertPF.isChecked()
        resetAmp = self._ui.checkBox_resetAmp.isChecked()
        symmeterize = self._ui.checkBox_symmeterize.isChecked()

        checked = self._ui.buttonGroupGuess.checkedButton()
        if checked == self._ui.radioButtonPlane:
            guess = ('plane',)
        elif checked == self._ui.radioButtonMirror:
            z0 = self._ui.doubleSpinBoxMirrorDistance.value()
            guess = ('mirror',z0)
        elif checked == self._ui.radioButtonFromFile:
            filename = QtWidgets.QFileDialog.getOpenFileName(None,'Open initial guess',
                                                  '','*.npy')
            if filename:
                guess = ('file', str(filename))
        if (guess[0] != 'file') or (guess[0] == 'file' and filename):
            PF, Amp = self._control.retrievePF(pxlSize, l, n, NA, f, guess, nIt, neglect_defocus, invert=invertPF,wavelengths=numWaves, wavestep = dWave, resetAmp=resetAmp,symmeterize=symmeterize)
            self._ui.tabWidget_viewer.setEnabled(True)
        self.displayPhase(PF)
        self.displayAmpli(Amp)
        self.use_zernike = False

    def addModulate(self):
        new_modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, new_modulation.checkbox)
        self._modulations.append(new_modulation)
        self._control.addMOD()


    def removeModulate(self):
        '''
        Check the number of checked modulations
        '''
        mod_list = self.selectModulation()
        if len(mod_list) > 0:
            for idx in mod_list:
                item = self._modulations[idx].checkbox
                item.setParent(None)
                #self._ui.verticalLayoutModulations.removeWidget(widget)
            self._control.removeMOD(mod_list)
        else:
            print("Please select modulation first.")



    def applyToMirror(self):
        '''
        1. Find the number of checked modulations
        2. Synthesize the pattern
        '''
        #print("The function has not been defined. Please give it a definition.")
        mod_list = self.selectModulation()
        if len(mod_list) > 0:
            self._control.modulateMirror(mod_list)
        else:
            print("Please select modulation first.")

    def selectModulation(self):
        '''
        select modulation from current list
        '''
        ic = 0
        mod_list = []
        for mod_item in self._modulations:
            if mod_item.checkbox.isChecked():
                mod_list.append(ic)
            ic +=1
        return mod_list

    def savePF(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(None,'Save to file',
                                                  '','*.npy')
        if filename:
            self._control.savePF(str(filename))

    def saveConfig(self):
        '''
        save configuration into a file
        '''
        conf_dest = self._ui.lineEdit_conf.text()
        if conf_dest == '':
            exp_destination = QtWidgets.QFileDialog.getOpenFileName(None, 'Open Configuration file:', '', '*.*')[0]


    def unwrap(self):
        '''
        This is redundant because we do not need to unwrap.
        '''
        unwrappedPhase = self._control.unwrap()
        self.displayPhase(unwrappedPhase)




    # ----------------------------------Below is a set of events, copied from orcaflase ----------------
    def _update(self):
        np_image = self._control.getImageForPreview()
        if np_image is not None:
            #self._plotCrossSection(np_image)
            np_min = np_image[2:-2,2:-2].min()
            np_max = np_image[2:-2,2:-2].max()
            if self._autoscale:
                self.vmin = np_min
                self.vmax = np_max

            if self.vmin==self.vmax or self.vmin>self.vmax:
                self.vmax = 1+self.vmin
            qt_image = qext.numpy_to_qimage8(np_image, self.vmin, self.vmax, self._cmap)
            self._pixmap = QtGui.QPixmap.fromImage(qt_image)
            #xdim = np.minimum(512, self._control._props["dimensions"][0])
            #ydim = np.minimum(512, self._control._props["dimensions"][1])
            #self._pixmap = self._pixmap.scaled(xdim, ydim, QtCore.Qt.KeepAspectRatio)
            self._pixmap = self._pixmap.scaled(512, 512, QtCore.Qt.KeepAspectRatio)
            #if self._ui.checkBox_showTarget.isChecked():
            #    self._drawTarget()
            self._ui.labelDisplay.update()


    def _labelDisplay_paintEvent(self, event):
        qp = QtGui.QPainter(self._ui.labelDisplay)
        if self._pixmap != None:
            qp.drawPixmap(0,0,self._pixmap)
        qp.setPen(QtCore.Qt.red)

    def _mousePressEvent(self,x,y):
        self.x = int(x*(self._yres/512.))
        self.y = int(y*(self._xres/512.))



    def shutDown(self):
        if self._scanner:
            self._scanner.wait()


class Modulation:
    def __init__(self, index, ui):
        self.index = index
        self.checkbox = QtWidgets.QCheckBox(str(self.index))
        self.checkbox.stateChanged.connect(ui.selectModulation)
        self.checkbox.toggle()

    def index_update(self, new_idx):
        '''
        if some checkboxes are removed, update the indices of the rest accordingly.
        '''
        pass



class Scanner(QtCore.QThread):

    def __init__(self, control, range_, nSlices, nFrames, center_xy, imagedest, maskRadius, maskCenter):
        QtCore.QThread.__init__(self)

        self.control = control
        self.range_ = range_
        self.nSlices = nSlices
        self.nFrames = nFrames
        self.center_xy = center_xy
        self.imagedest = imagedest
        self.maskRadius = maskRadius
        self.maskCenter = maskCenter


    def run(self):
        self.control.acquirePSF(self.range_, self.nSlices, self.nFrames,
                                self.center_xy, self.imagedest,
                                self.maskRadius, self.maskCenter)
