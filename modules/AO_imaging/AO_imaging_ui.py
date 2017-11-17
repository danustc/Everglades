#!/usr/bin/python


from PyQt5 import QtWidgets,QtCore
import inLib
from Utilities import QExtensions as qext
import numpy as np
from numpy.lib.scimath import sqrt as _msqrt
from . import fit_results_design
import copy
import time
from libs import scipy_gaussfitter

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
        self._ui.lineEdit_fname.setText(self._control._settings['filename'])

        self._ui.spinBoxIterations.setValue(self._control._settings['nIterations'])
        self._ui.pushButton_Acquire.clicked.connect(self.acquirePSF)
        self._ui.pushButton_modulate.clicked.connect(self.modulate)
        self._ui.pushButton_save.clicked.connect(self.savePF)

        self._ui.pushButton_modUnwrapped.clicked.connect(self.modulateUnwrapped)

        self._ui.groupBoxModulations.toggled.connect(self._modulations_toggled)

        self._ui.pushButton_unwrap.clicked.connect(self.unwrap)

        self._ui.pushButton_zernFitUnwrapped.clicked.connect(self.fitUnwrapped)

        self._ui.spinBox_zernModesToFit.setValue(self._control.zernModesToFit)
        self._ui.spinBox_zernModesToFit.valueChanged.connect(self.setZernModesToFit)

        self._ui.pushButton_sync.clicked.connect(self.set_modulations)

        self.zernRadius = 0

        self._modulations = []
        self.use_zernike = False
        self.remove_PTTD = True


        self._scanner = None

        self._sharpnessPlot = None

        self.imsize = (256,256)
        self.pixelSize = 163
        self.diffLimit = 800


        self._ui.tabWidget_viewer.setEnabled(True)

        self.hasSLM = self._control.hasSLM
        self.hasMirror = self._control.hasMirror

    def _displayPhase(self, phase):
        self._ui.mplwidget_PF.figure.axes[0].get_xaxis().set_visible(False)
        self._ui.mplwidget_PF.figure.axes[0].get_yaxis().set_visible(False)
        self._ui.mplwidget_PF.figure.axes[0].matshow(phase, cmap='RdBu')
        self._ui.mplwidget_PF.draw()

    def _plotSharpness(self, sharpness):
        self._ui.mplwidget_sharpness.figure.axes[0].plot(sharpness)
        self._ui.mplwidget_sharpness.draw()

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
        range_ = self._ui.doubleSpinBoxRange.value()
        nSlices = self._ui.spinBoxSlices.value()
        nFrames = self._ui.spinBoxFrames.value()
        center_xy = self._ui.checkBoxCenterLateral.isChecked()
        maskRadius = self._ui.spinBox_maskRadius.value()
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        fname = None
        fname = str(self._ui.lineEdit_fname.text())
        self._scanner = Scanner(self._control, range_, nSlices, nFrames, center_xy,
                                fname, maskRadius, (cX,cY))
        self._scanner.finished.connect(self._on_scan_done)
        self._ui.pushButton_Acquire.setEnabled(False)
        self._scanner.start()


    def _on_scan_done(self):
        self._ui.pushButton_Acquire.setEnabled(True)
        self._ui.groupBoxPhase.setEnabled(True)
        time.sleep(2)
        sharpness = self._control.getSharpness()
        self._plotSharpness(sharpness)

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

    def fitUnwrapped(self):
        ignore4 = self._ui.checkBox_ignore4.isChecked()
        resultFit = self._control.zernFitUnwrapped(skip4orders=ignore4)
        self._ui.mplwidget_PF_2.figure.axes[0].matshow(resultFit, cmap='RdBu')
        self._ui.mplwidget_PF_2.draw()

    def modulateUnwrappedZernike(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        mask = self._ui.checkBox_useMask.isChecked()
        self._control.modZernFitUnwrapped(useMask=mask, radius=self.zernRadius)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()

    def oneRun(self):
        # added by Dan to perform one-Run experiment
        self._control.one_Run(4)



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

    def __init__(self, control, range_, nSlices, nFrames, center_xy, fname, maskRadius, maskCenter):
        QtCore.QThread.__init__(self)

        self.control = control
        self.range_ = range_
        self.nSlices = nSlices
        self.nFrames = nFrames
        self.center_xy = center_xy
        self.fname = fname
        self.maskRadius = maskRadius
        self.maskCenter = maskCenter

    def run(self):
        self.control.acquirePSF(self.range_, self.nSlices, self.nFrames,
                                self.center_xy, self.fname,
                                self.maskRadius, self.maskCenter)
