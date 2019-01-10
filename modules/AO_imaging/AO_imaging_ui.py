#!/usr/bin/python


from PyQt5 import QtWidgets,QtCore, QtGui
import inLib
import numpy as np
from . import fit_results_design
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
        design_path = 'modules.AO_imaging.AO_imaging_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        self._ui.buttonGroupGuess = QtWidgets.QButtonGroup()
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonPlane)
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonMirror)
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonFromFile)
        self._ui.pushButton_zfit.clicked.connect(self.fitPF)
        self._ui.pushButtonPhase.clicked.connect(self.retrievePF)

        self._ui.doubleSpinBoxRange.setValue(self._control._settings['range'])
        self._ui.spinBoxSlices.setValue(self._control._settings['nSlices'])
        self._ui.spinBoxFrames.setValue(self._control._settings['nFrames'])
        self._ui.lineEdit_imagedest.setText(self._control._settings['filename'])

        self._ui.spinBoxIterations.setValue(self._control._settings['nIterations'])
        self._ui.pushButton_Acquire.clicked.connect(self.acquirePSF)
        self._ui.pushButton_modulate.clicked.connect(self.modulate)
        self._ui.pushButton_save.clicked.connect(self.savePF)

        self._ui.labelDisplay.paintEvent = self._labelDisplay_paintEvent
        clickable(self._ui.labelDisplay).connect(self._mousePressEvent)


        self._ui.spinBox_zernModesToFit.setValue(self._control.zernModesToFit)
        self._ui.spinBox_zernModesToFit.valueChanged.connect(self.setZernModesToFit)

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
        

    def _displayPhase(self, phase):
        self._ui.mplwidget_PF.figure.axes[0].get_xaxis().set_visible(False)
        self._ui.mplwidget_PF.figure.axes[0].get_yaxis().set_visible(False)
        self._ui.mplwidget_PF.figure.axes[0].matshow(phase, cmap='RdBu')
        self._ui.mplwidget_PF.draw()


    def _updateImSize(self):
        self.imsize = self._control.updateImSize()


    def _modulation_toggled(self, state):
        pass
        '''
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setModulationActive(m.index, state)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()
        '''


    def _modulations_toggled(self, state):
        pass
        '''
        if state == False:
            for m in self._modulations:
                self._control.setModulationActive(m.index, state)
        else:
            self._modulation_toggled(state)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()
        '''

    def set_modulations(self):
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setModulationActive(m.index, state)



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
        self._updater.start()


    def _on_scan_done(self):
        self._ui.pushButton_Acquire.setEnabled(True)
        self._ui.groupBoxPhase.setEnabled(True)
        time.sleep(1)
        
      

    def foundMaxArgSharp(self, argmax):
        '''
        Connects to signal 'maxArgSharpness(int)'
        '''
        self._ui.label_sharpnessArgMax.setText("Arg. max: %i" % argmax)

    def foundMaxFitSharp(self, fit):
        self._ui.label_sharpnessFitMax.setText("Fit max: %.2f" % fit)

    def advanceModulation(self, mod_index):
        '''
        Connected to signal 'nextModulation(int)' which is called by RunningSharpness

        This is triggered each time a new pattern is to be displayed on the adaptive optics device.
        '''
        coeff = self._control.advanceModulation()
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay() #updates display
        self._ui.label_mod_index.setText("Index: %i" % mod_index)
        self._ui.label_mod_value.setText("Value: %.2f" % coeff)

    def fitPF(self):
        PF = self._control.fit()
        fit_result_dialog = FitResultsDialog(PF)
        if fit_result_dialog.exec_():
            print('AO_imaging: Fit accepted.')
            self.use_zernike = True
            remove = fit_result_dialog.getRemove()
            if remove:
                print('AO_imaging: Remove pistion, tip, tilt and defocus.')
                PF = self._control.removePTTD()
            self._displayPhase(PF.zernike)
        else:
            self.use_zernike = False


    def retrievePF(self, crop = True):
        '''
        perform phase retrieval, crop if necessary.
        '''
        pxlSize = self._ui.doubleSpinBoxPixel.value()
        l = self._ui.doubleSpinBoxWavelength.value()
        n = self._ui.doubleSpinBoxIndex.value()
        NA = self._ui.doubleSpinBoxNA.value()
        f = self._ui.doubleSpinBoxFocal.value()
        numWaves = self._ui.spinBox_numWavelengths.value()
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
            PF = self._control.retrievePF(pxlSize, l, n, NA, f, guess, nIt, neglect_defocus, invert=invertPF,wavelengths=numWaves, resetAmp=resetAmp,symmeterize=symmeterize)
            self._ui.tabWidget_viewer.setEnabled(True)
            if crop:
                PF_crop = pupil_crop(PF)
            else:
                self._displayPhase(PF)
        self.use_zernike = False


    def modulate(self): # pushButtonModulate triggering
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modulatePF(self.use_zernike)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()

    def modulateUnwrapped(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modulatePF_unwrapped()
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()

    def savePF(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(None,'Save to file',
                                                  '','*.npy')
        if filename:
            self._control.savePF(str(filename))

    def unwrap(self):
        unwrappedPhase = self._control.unwrap()
        self._displayPhase(unwrappedPhase)

    def setZernModesToFit(self):
        nmodes = self._ui.spinBox_zernModesToFit.value()
        self._control.setZernModesToFit(nmodes)



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
        self.checkbox.stateChanged.connect(ui._modulation_toggled)
        self.checkbox.toggle()



class FitResultsDialog(QtWidgets.QDialog):

    def __init__(self, PF, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        self.PF = PF
        self.ui = fit_results_design.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.lineEditCoefficients.setText(str(PF.zernike_coefficients))
        self.ui.mplwidget.figure.delaxes(self.ui.mplwidget.figure.axes[0])
        axes_raw = self.ui.mplwidget.figure.add_subplot(131)
        axes_raw.matshow(PF.phase, cmap='RdBu')
        axes_raw.set_title('Raw data')
        axes_fit = self.ui.mplwidget.figure.add_subplot(132)
        axes_fit.matshow(PF.zernike, cmap='RdBu', vmin=PF.phase.min(), vmax=PF.phase.max())
        axes_fit.set_title('Fit')
        axes_coefficients = self.ui.mplwidget.figure.add_subplot(133)
        axes_coefficients.bar(np.arange(25), PF.zernike_coefficients)
        axes_coefficients.set_title('Zernike coefficients')


    def getRemove(self):
        return self.ui.checkBoxRemovePTTD.isChecked()


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
