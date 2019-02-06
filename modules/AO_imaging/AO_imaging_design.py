# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AO_imaging_design.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from myWidget.matplotlibwidget import MatplotlibWidget
class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1053, 656)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(912, 0))
        Form.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.horizontalLayout = QtWidgets.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox = QtWidgets.QGroupBox(Form)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setMaximumSize(QtCore.QSize(150, 16777215))
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout_4.addWidget(self.lineEdit, 12, 0, 1, 2)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem, 16, 0, 1, 1)
        self.checkBoxCenterLateral = QtWidgets.QCheckBox(self.groupBox)
        self.checkBoxCenterLateral.setChecked(True)
        self.checkBoxCenterLateral.setObjectName("checkBoxCenterLateral")
        self.gridLayout_4.addWidget(self.checkBoxCenterLateral, 7, 0, 1, 2)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.gridLayout_4.addWidget(self.label_3, 2, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout_4.addWidget(self.label_4, 3, 0, 1, 1)
        self.pushButton_loadconfig = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_loadconfig.setObjectName("pushButton_loadconfig")
        self.gridLayout_4.addWidget(self.pushButton_loadconfig, 14, 0, 1, 2)
        self.label_CY = QtWidgets.QLabel(self.groupBox)
        self.label_CY.setObjectName("label_CY")
        self.gridLayout_4.addWidget(self.label_CY, 6, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 1, 0, 1, 1)
        self.label_fname = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_fname.sizePolicy().hasHeightForWidth())
        self.label_fname.setSizePolicy(sizePolicy)
        self.label_fname.setObjectName("label_fname")
        self.gridLayout_4.addWidget(self.label_fname, 9, 0, 1, 2)
        self.pushButton_saveconfig = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_saveconfig.setObjectName("pushButton_saveconfig")
        self.gridLayout_4.addWidget(self.pushButton_saveconfig, 13, 0, 1, 2)
        self.label_mask = QtWidgets.QLabel(self.groupBox)
        self.label_mask.setObjectName("label_mask")
        self.gridLayout_4.addWidget(self.label_mask, 4, 0, 1, 1)
        self.doubleSpinBoxRange = QtWidgets.QDoubleSpinBox(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.doubleSpinBoxRange.sizePolicy().hasHeightForWidth())
        self.doubleSpinBoxRange.setSizePolicy(sizePolicy)
        self.doubleSpinBoxRange.setDecimals(1)
        self.doubleSpinBoxRange.setMinimum(-100.0)
        self.doubleSpinBoxRange.setSingleStep(0.1)
        self.doubleSpinBoxRange.setProperty("value", 10.0)
        self.doubleSpinBoxRange.setObjectName("doubleSpinBoxRange")
        self.gridLayout_4.addWidget(self.doubleSpinBoxRange, 1, 1, 1, 1)
        self.spinBox_maskRadius = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_maskRadius.setMinimum(16)
        self.spinBox_maskRadius.setMaximum(512)
        self.spinBox_maskRadius.setSingleStep(10)
        self.spinBox_maskRadius.setProperty("value", 80)
        self.spinBox_maskRadius.setObjectName("spinBox_maskRadius")
        self.gridLayout_4.addWidget(self.spinBox_maskRadius, 4, 1, 1, 1)
        self.pushButton_Acquire = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_Acquire.sizePolicy().hasHeightForWidth())
        self.pushButton_Acquire.setSizePolicy(sizePolicy)
        self.pushButton_Acquire.setObjectName("pushButton_Acquire")
        self.gridLayout_4.addWidget(self.pushButton_Acquire, 15, 0, 1, 2)
        self.lineEdit_cX = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_cX.setMaximumSize(QtCore.QSize(50, 16777215))
        self.lineEdit_cX.setObjectName("lineEdit_cX")
        self.gridLayout_4.addWidget(self.lineEdit_cX, 5, 1, 1, 1)
        self.label_CX = QtWidgets.QLabel(self.groupBox)
        self.label_CX.setObjectName("label_CX")
        self.gridLayout_4.addWidget(self.label_CX, 5, 0, 1, 1)
        self.lineEdit_imagedest = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_imagedest.setMinimumSize(QtCore.QSize(130, 0))
        self.lineEdit_imagedest.setObjectName("lineEdit_imagedest")
        self.gridLayout_4.addWidget(self.lineEdit_imagedest, 10, 0, 1, 1)
        self.lineEdit_cY = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_cY.setMaximumSize(QtCore.QSize(50, 16777215))
        self.lineEdit_cY.setObjectName("lineEdit_cY")
        self.gridLayout_4.addWidget(self.lineEdit_cY, 6, 1, 1, 1)
        self.spinBoxSlices = QtWidgets.QSpinBox(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spinBoxSlices.sizePolicy().hasHeightForWidth())
        self.spinBoxSlices.setSizePolicy(sizePolicy)
        self.spinBoxSlices.setReadOnly(False)
        self.spinBoxSlices.setMinimum(3)
        self.spinBoxSlices.setSingleStep(2)
        self.spinBoxSlices.setProperty("value", 21)
        self.spinBoxSlices.setObjectName("spinBoxSlices")
        self.gridLayout_4.addWidget(self.spinBoxSlices, 2, 1, 1, 1)
        self.spinBoxFrames = QtWidgets.QSpinBox(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spinBoxFrames.sizePolicy().hasHeightForWidth())
        self.spinBoxFrames.setSizePolicy(sizePolicy)
        self.spinBoxFrames.setMinimum(1)
        self.spinBoxFrames.setProperty("value", 5)
        self.spinBoxFrames.setObjectName("spinBoxFrames")
        self.gridLayout_4.addWidget(self.spinBoxFrames, 3, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout_4.addWidget(self.label, 11, 0, 1, 2)
        self.checkBox_deskew = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_deskew.setObjectName("checkBox_deskew")
        self.gridLayout_4.addWidget(self.checkBox_deskew, 8, 0, 1, 1)
        self.horizontalLayout.addWidget(self.groupBox)
        self.groupBoxPhase = QtWidgets.QGroupBox(Form)
        self.groupBoxPhase.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBoxPhase.sizePolicy().hasHeightForWidth())
        self.groupBoxPhase.setSizePolicy(sizePolicy)
        self.groupBoxPhase.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.groupBoxPhase.setObjectName("groupBoxPhase")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBoxPhase)
        self.gridLayout_3.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.checkBox_ignore4 = QtWidgets.QCheckBox(self.groupBoxPhase)
        self.checkBox_ignore4.setObjectName("checkBox_ignore4")
        self.gridLayout_3.addWidget(self.checkBox_ignore4, 17, 1, 1, 1)
        self.spinBox_numWavelengths = QtWidgets.QSpinBox(self.groupBoxPhase)
        self.spinBox_numWavelengths.setMinimum(1)
        self.spinBox_numWavelengths.setProperty("value", 1)
        self.spinBox_numWavelengths.setObjectName("spinBox_numWavelengths")
        self.gridLayout_3.addWidget(self.spinBox_numWavelengths, 5, 1, 1, 1)
        self.pushButton_zernFitUnwrapped = QtWidgets.QPushButton(self.groupBoxPhase)
        self.pushButton_zernFitUnwrapped.setObjectName("pushButton_zernFitUnwrapped")
        self.gridLayout_3.addWidget(self.pushButton_zernFitUnwrapped, 13, 1, 1, 1)
        self.pushButton_unwrap = QtWidgets.QPushButton(self.groupBoxPhase)
        self.pushButton_unwrap.setObjectName("pushButton_unwrap")
        self.gridLayout_3.addWidget(self.pushButton_unwrap, 11, 1, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_13.setObjectName("label_13")
        self.gridLayout_3.addWidget(self.label_13, 5, 0, 1, 1)
        self.spinBox_zernModesToFit = QtWidgets.QSpinBox(self.groupBoxPhase)
        self.spinBox_zernModesToFit.setObjectName("spinBox_zernModesToFit")
        self.gridLayout_3.addWidget(self.spinBox_zernModesToFit, 14, 1, 1, 1)
        self.checkBox_resetAmp = QtWidgets.QCheckBox(self.groupBoxPhase)
        self.checkBox_resetAmp.setObjectName("checkBox_resetAmp")
        self.gridLayout_3.addWidget(self.checkBox_resetAmp, 13, 0, 1, 1)
        self.checkBox_invertPF = QtWidgets.QCheckBox(self.groupBoxPhase)
        self.checkBox_invertPF.setObjectName("checkBox_invertPF")
        self.gridLayout_3.addWidget(self.checkBox_invertPF, 12, 0, 1, 1)
        self.checkBox_symmeterize = QtWidgets.QCheckBox(self.groupBoxPhase)
        self.checkBox_symmeterize.setObjectName("checkBox_symmeterize")
        self.gridLayout_3.addWidget(self.checkBox_symmeterize, 10, 0, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem1, 17, 0, 1, 1)
        self.doubleSpinBoxPixel = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBoxPixel.setDecimals(3)
        self.doubleSpinBoxPixel.setMaximum(9.0)
        self.doubleSpinBoxPixel.setSingleStep(0.01)
        self.doubleSpinBoxPixel.setProperty("value", 0.102)
        self.doubleSpinBoxPixel.setObjectName("doubleSpinBoxPixel")
        self.gridLayout_3.addWidget(self.doubleSpinBoxPixel, 0, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 0, 0, 1, 1)
        self.doubleSpinBoxFocal = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBoxFocal.setMaximum(9999.99)
        self.doubleSpinBoxFocal.setProperty("value", 5000.0)
        self.doubleSpinBoxFocal.setObjectName("doubleSpinBoxFocal")
        self.gridLayout_3.addWidget(self.doubleSpinBoxFocal, 4, 1, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_14.setObjectName("label_14")
        self.gridLayout_3.addWidget(self.label_14, 14, 0, 1, 1)
        self.doubleSpinBoxWavelength = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBoxWavelength.setDecimals(3)
        self.doubleSpinBoxWavelength.setProperty("value", 0.515)
        self.doubleSpinBoxWavelength.setObjectName("doubleSpinBoxWavelength")
        self.gridLayout_3.addWidget(self.doubleSpinBoxWavelength, 1, 1, 1, 1)
        self.doubleSpinBoxIndex = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBoxIndex.setSingleStep(0.01)
        self.doubleSpinBoxIndex.setProperty("value", 1.33)
        self.doubleSpinBoxIndex.setObjectName("doubleSpinBoxIndex")
        self.gridLayout_3.addWidget(self.doubleSpinBoxIndex, 2, 1, 1, 1)
        self.spinBoxIterations = QtWidgets.QSpinBox(self.groupBoxPhase)
        self.spinBoxIterations.setMinimum(1)
        self.spinBoxIterations.setMaximum(999)
        self.spinBoxIterations.setProperty("value", 20)
        self.spinBoxIterations.setObjectName("spinBoxIterations")
        self.gridLayout_3.addWidget(self.spinBoxIterations, 9, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_10.setObjectName("label_10")
        self.gridLayout_3.addWidget(self.label_10, 9, 0, 1, 1)
        self.pushButtonPhase = QtWidgets.QPushButton(self.groupBoxPhase)
        self.pushButtonPhase.setObjectName("pushButtonPhase")
        self.gridLayout_3.addWidget(self.pushButtonPhase, 11, 0, 1, 1)
        self.pushButton_modUnwrapped = QtWidgets.QPushButton(self.groupBoxPhase)
        self.pushButton_modUnwrapped.setObjectName("pushButton_modUnwrapped")
        self.gridLayout_3.addWidget(self.pushButton_modUnwrapped, 12, 1, 1, 1)
        self.checkBoxNeglectDefocus = QtWidgets.QCheckBox(self.groupBoxPhase)
        self.checkBoxNeglectDefocus.setChecked(True)
        self.checkBoxNeglectDefocus.setObjectName("checkBoxNeglectDefocus")
        self.gridLayout_3.addWidget(self.checkBoxNeglectDefocus, 7, 0, 1, 2)
        self.label_9 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_9.setObjectName("label_9")
        self.gridLayout_3.addWidget(self.label_9, 4, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_7.setObjectName("label_7")
        self.gridLayout_3.addWidget(self.label_7, 2, 0, 1, 1)
        self.doubleSpinBoxNA = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBoxNA.setSingleStep(0.01)
        self.doubleSpinBoxNA.setProperty("value", 1.0)
        self.doubleSpinBoxNA.setObjectName("doubleSpinBoxNA")
        self.gridLayout_3.addWidget(self.doubleSpinBoxNA, 3, 1, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_8.setObjectName("label_8")
        self.gridLayout_3.addWidget(self.label_8, 3, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 1, 0, 1, 1)
        self.groupBox_4 = QtWidgets.QGroupBox(self.groupBoxPhase)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.radioButtonPlane = QtWidgets.QRadioButton(self.groupBox_4)
        self.radioButtonPlane.setChecked(True)
        self.radioButtonPlane.setObjectName("radioButtonPlane")
        self.gridLayout_2.addWidget(self.radioButtonPlane, 0, 0, 1, 1)
        self.radioButtonMirror = QtWidgets.QRadioButton(self.groupBox_4)
        self.radioButtonMirror.setObjectName("radioButtonMirror")
        self.gridLayout_2.addWidget(self.radioButtonMirror, 1, 0, 1, 1)
        self.doubleSpinBoxMirrorDistance = QtWidgets.QDoubleSpinBox(self.groupBox_4)
        self.doubleSpinBoxMirrorDistance.setSingleStep(0.01)
        self.doubleSpinBoxMirrorDistance.setProperty("value", 2.0)
        self.doubleSpinBoxMirrorDistance.setObjectName("doubleSpinBoxMirrorDistance")
        self.gridLayout_2.addWidget(self.doubleSpinBoxMirrorDistance, 1, 1, 1, 1)
        self.radioButtonFromFile = QtWidgets.QRadioButton(self.groupBox_4)
        self.radioButtonFromFile.setObjectName("radioButtonFromFile")
        self.gridLayout_2.addWidget(self.radioButtonFromFile, 2, 0, 1, 1)
        self.gridLayout_3.addWidget(self.groupBox_4, 8, 0, 1, 3)
        self.label_wstep = QtWidgets.QLabel(self.groupBoxPhase)
        self.label_wstep.setObjectName("label_wstep")
        self.gridLayout_3.addWidget(self.label_wstep, 6, 0, 1, 1)
        self.doubleSpinBox_dwave = QtWidgets.QDoubleSpinBox(self.groupBoxPhase)
        self.doubleSpinBox_dwave.setDecimals(3)
        self.doubleSpinBox_dwave.setMinimum(0.001)
        self.doubleSpinBox_dwave.setMaximum(1.0)
        self.doubleSpinBox_dwave.setSingleStep(0.001)
        self.doubleSpinBox_dwave.setProperty("value", 0.005)
        self.doubleSpinBox_dwave.setObjectName("doubleSpinBox_dwave")
        self.gridLayout_3.addWidget(self.doubleSpinBox_dwave, 6, 1, 1, 1)
        self.horizontalLayout.addWidget(self.groupBoxPhase)
        self.tabWidget_viewer = QtWidgets.QTabWidget(Form)
        self.tabWidget_viewer.setObjectName("tabWidget_viewer")
        self.tab_image = QtWidgets.QWidget()
        self.tab_image.setObjectName("tab_image")
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(self.tab_image)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(10, 540, 388, 524))
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.pushButton_snap = QtWidgets.QPushButton(self.horizontalLayoutWidget_3)
        self.pushButton_snap.setObjectName("pushButton_snap")
        self.horizontalLayout_4.addWidget(self.pushButton_snap)
        self.pushButton = QtWidgets.QPushButton(self.horizontalLayoutWidget_3)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout_4.addWidget(self.pushButton)
        self.labelDisplay = QtWidgets.QLabel(self.tab_image)
        self.labelDisplay.setEnabled(True)
        self.labelDisplay.setGeometry(QtCore.QRect(10, 10, 512, 512))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.labelDisplay.sizePolicy().hasHeightForWidth())
        self.labelDisplay.setSizePolicy(sizePolicy)
        self.labelDisplay.setMinimumSize(QtCore.QSize(512, 512))
        self.labelDisplay.setMaximumSize(QtCore.QSize(512, 512))
        self.labelDisplay.setFrameShape(QtWidgets.QFrame.Box)
        self.labelDisplay.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.labelDisplay.setText("")
        self.labelDisplay.setObjectName("labelDisplay")
        self.tabWidget_viewer.addTab(self.tab_image, "")
        self.tab_pupil = QtWidgets.QWidget()
        self.tab_pupil.setObjectName("tab_pupil")
        self.gridLayoutWidget = QtWidgets.QWidget(self.tab_pupil)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 10, 521, 521))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.mplwidget_ampli = MatplotlibWidget(self.gridLayoutWidget)
        self.mplwidget_ampli.setObjectName("mplwidget_ampli")
        self.gridLayout.addWidget(self.mplwidget_ampli, 0, 1, 1, 1)
        self.mplwidget_sync = MatplotlibWidget(self.gridLayoutWidget)
        self.mplwidget_sync.setObjectName("mplwidget_sync")
        self.gridLayout.addWidget(self.mplwidget_sync, 1, 0, 1, 1)
        self.mplwidget_seg = MatplotlibWidget(self.gridLayoutWidget)
        self.mplwidget_seg.setObjectName("mplwidget_seg")
        self.gridLayout.addWidget(self.mplwidget_seg, 1, 1, 1, 1)
        self.mplwidget_phase = MatplotlibWidget(self.gridLayoutWidget)
        self.mplwidget_phase.setObjectName("mplwidget_phase")
        self.gridLayout.addWidget(self.mplwidget_phase, 0, 0, 1, 1)
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.tab_pupil)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 580, 521, 31))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButton_ampli = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_ampli.sizePolicy().hasHeightForWidth())
        self.pushButton_ampli.setSizePolicy(sizePolicy)
        self.pushButton_ampli.setObjectName("pushButton_ampli")
        self.horizontalLayout_2.addWidget(self.pushButton_ampli)
        self.pushButton_zfit = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_zfit.sizePolicy().hasHeightForWidth())
        self.pushButton_zfit.setSizePolicy(sizePolicy)
        self.pushButton_zfit.setObjectName("pushButton_zfit")
        self.horizontalLayout_2.addWidget(self.pushButton_zfit)
        self.pushButton_modulate = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_modulate.sizePolicy().hasHeightForWidth())
        self.pushButton_modulate.setSizePolicy(sizePolicy)
        self.pushButton_modulate.setObjectName("pushButton_modulate")
        self.horizontalLayout_2.addWidget(self.pushButton_modulate)
        self.pushButton_save = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_save.sizePolicy().hasHeightForWidth())
        self.pushButton_save.setSizePolicy(sizePolicy)
        self.pushButton_save.setObjectName("pushButton_save")
        self.horizontalLayout_2.addWidget(self.pushButton_save)
        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(self.tab_pupil)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(10, 540, 521, 31))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_pupilOperation = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_pupilOperation.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_pupilOperation.setObjectName("horizontalLayout_pupilOperation")
        self.pushButton_flr = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_flr.sizePolicy().hasHeightForWidth())
        self.pushButton_flr.setSizePolicy(sizePolicy)
        self.pushButton_flr.setObjectName("pushButton_flr")
        self.horizontalLayout_pupilOperation.addWidget(self.pushButton_flr)
        self.pushButton_fud = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_fud.sizePolicy().hasHeightForWidth())
        self.pushButton_fud.setSizePolicy(sizePolicy)
        self.pushButton_fud.setObjectName("pushButton_fud")
        self.horizontalLayout_pupilOperation.addWidget(self.pushButton_fud)
        self.pushButton_rot90 = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_rot90.sizePolicy().hasHeightForWidth())
        self.pushButton_rot90.setSizePolicy(sizePolicy)
        self.pushButton_rot90.setObjectName("pushButton_rot90")
        self.horizontalLayout_pupilOperation.addWidget(self.pushButton_rot90)
        self.pushButton_inv = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_inv.sizePolicy().hasHeightForWidth())
        self.pushButton_inv.setSizePolicy(sizePolicy)
        self.pushButton_inv.setObjectName("pushButton_inv")
        self.horizontalLayout_pupilOperation.addWidget(self.pushButton_inv)
        self.tabWidget_viewer.addTab(self.tab_pupil, "")
        self.horizontalLayout.addWidget(self.tabWidget_viewer)
        self.groupBoxModulations = QtWidgets.QGroupBox(Form)
        self.groupBoxModulations.setCheckable(True)
        self.groupBoxModulations.setObjectName("groupBoxModulations")
        self.verticalLayoutModulations = QtWidgets.QVBoxLayout(self.groupBoxModulations)
        self.verticalLayoutModulations.setObjectName("verticalLayoutModulations")
        spacerItem2 = QtWidgets.QSpacerItem(20, 317, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayoutModulations.addItem(spacerItem2)
        self.pushButton_ADM = QtWidgets.QPushButton(self.groupBoxModulations)
        self.pushButton_ADM.setObjectName("pushButton_ADM")
        self.verticalLayoutModulations.addWidget(self.pushButton_ADM)
        self.horizontalLayout.addWidget(self.groupBoxModulations)

        self.retranslateUi(Form)
        self.tabWidget_viewer.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Adaptive Optics"))
        self.groupBox.setTitle(_translate("Form", "PSF"))
        self.checkBoxCenterLateral.setText(_translate("Form", "Center PSF laterally"))
        self.label_3.setText(_translate("Form", "# Slices:"))
        self.label_4.setText(_translate("Form", "# Frames/slice:"))
        self.pushButton_loadconfig.setText(_translate("Form", "Load configuration"))
        self.label_CY.setText(_translate("Form", "Center Y:"))
        self.label_2.setText(_translate("Form", "Range [um]:"))
        self.label_fname.setText(_translate("Form", "Image destination"))
        self.pushButton_saveconfig.setText(_translate("Form", "Save configuration"))
        self.label_mask.setText(_translate("Form", "Mask radius:"))
        self.pushButton_Acquire.setText(_translate("Form", "Acquire PSF"))
        self.lineEdit_cX.setText(_translate("Form", "128"))
        self.label_CX.setText(_translate("Form", "Center X:"))
        self.lineEdit_cY.setText(_translate("Form", "128"))
        self.label.setText(_translate("Form", "Configuration path"))
        self.checkBox_deskew.setText(_translate("Form", "Deskew"))
        self.groupBoxPhase.setTitle(_translate("Form", "Phase retrieval"))
        self.checkBox_ignore4.setText(_translate("Form", "Ignore 1st 4?"))
        self.pushButton_zernFitUnwrapped.setText(_translate("Form", "Fit Unwrpd"))
        self.pushButton_unwrap.setText(_translate("Form", "Unwrap"))
        self.label_13.setText(_translate("Form", "# Wavelengths:"))
        self.checkBox_resetAmp.setText(_translate("Form", "Reset amplitude?"))
        self.checkBox_invertPF.setText(_translate("Form", "Invert?"))
        self.checkBox_symmeterize.setText(_translate("Form", "Make symmetric?"))
        self.label_5.setText(_translate("Form", "Pixel size [um]:"))
        self.label_14.setText(_translate("Form", "Num modes to fit to:"))
        self.label_10.setText(_translate("Form", "# Iterations:"))
        self.pushButtonPhase.setText(_translate("Form", "Retrieve phase"))
        self.pushButton_modUnwrapped.setText(_translate("Form", "Mod. Unwrpd"))
        self.checkBoxNeglectDefocus.setText(_translate("Form", "Neglect defocus"))
        self.label_9.setText(_translate("Form", "Focal length:"))
        self.label_7.setText(_translate("Form", "Refractive index:"))
        self.label_8.setText(_translate("Form", "NA:"))
        self.label_6.setText(_translate("Form", "Wavelength [um]:"))
        self.groupBox_4.setTitle(_translate("Form", "Initial guess"))
        self.radioButtonPlane.setText(_translate("Form", "Plane wave"))
        self.radioButtonMirror.setText(_translate("Form", "Mirror"))
        self.radioButtonFromFile.setText(_translate("Form", "From file"))
        self.label_wstep.setText(_translate("Form", "Wavelength step"))
        self.pushButton_snap.setText(_translate("Form", "Snap shot"))
        self.pushButton.setText(_translate("Form", "XZ-view"))
        self.tabWidget_viewer.setTabText(self.tabWidget_viewer.indexOf(self.tab_image), _translate("Form", "Images"))
        self.pushButton_ampli.setText(_translate("Form", "Amplitude"))
        self.pushButton_zfit.setText(_translate("Form", "Fit"))
        self.pushButton_modulate.setText(_translate("Form", "Modulate"))
        self.pushButton_save.setText(_translate("Form", "Save"))
        self.pushButton_flr.setText(_translate("Form", "Flip L-R"))
        self.pushButton_fud.setText(_translate("Form", "Flip U-D"))
        self.pushButton_rot90.setText(_translate("Form", "Rotate 90"))
        self.pushButton_inv.setText(_translate("Form", "Invert"))
        self.tabWidget_viewer.setTabText(self.tabWidget_viewer.indexOf(self.tab_pupil), _translate("Form", "Pupils"))
        self.groupBoxModulations.setTitle(_translate("Form", "Modulations"))
        self.pushButton_ADM.setText(_translate("Form", "Apply!"))

#from matplotlibwidget import MatplotlibWidget
