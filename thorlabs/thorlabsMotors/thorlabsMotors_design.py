# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'thorlabsMotors_design.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(604, 379)
        self.gridLayoutWidget = QtWidgets.QWidget(Form)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(60, 59, 283, 131))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_position = QtWidgets.QLabel(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_position.sizePolicy().hasHeightForWidth())
        self.label_position.setSizePolicy(sizePolicy)
        self.label_position.setObjectName("label_position")
        self.gridLayout.addWidget(self.label_position, 0, 1, 1, 1)
        self.lineEdit_stepsize = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_stepsize.setObjectName("lineEdit_stepsize")
        self.gridLayout.addWidget(self.lineEdit_stepsize, 1, 2, 1, 1)
        self.lineEdit_absolute = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_absolute.setObjectName("lineEdit_absolute")
        self.gridLayout.addWidget(self.lineEdit_absolute, 1, 1, 1, 1)
        self.label_stepsize = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_stepsize.setObjectName("label_stepsize")
        self.gridLayout.addWidget(self.label_stepsize, 0, 2, 1, 1)
        self.pushButton_stepsize = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_stepsize.sizePolicy().hasHeightForWidth())
        self.pushButton_stepsize.setSizePolicy(sizePolicy)
        self.pushButton_stepsize.setObjectName("pushButton_stepsize")
        self.gridLayout.addWidget(self.pushButton_stepsize, 2, 2, 1, 1)
        self.pushButton_absolute = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_absolute.sizePolicy().hasHeightForWidth())
        self.pushButton_absolute.setSizePolicy(sizePolicy)
        self.pushButton_absolute.setObjectName("pushButton_absolute")
        self.gridLayout.addWidget(self.pushButton_absolute, 2, 1, 1, 1)
        self.lineEdit_relative = QtWidgets.QLineEdit(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_relative.sizePolicy().hasHeightForWidth())
        self.lineEdit_relative.setSizePolicy(sizePolicy)
        self.lineEdit_relative.setObjectName("lineEdit_relative")
        self.gridLayout.addWidget(self.lineEdit_relative, 1, 3, 1, 1)
        self.label_relative = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_relative.setObjectName("label_relative")
        self.gridLayout.addWidget(self.label_relative, 0, 3, 1, 1)
        self.pushButton_relative = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_relative.sizePolicy().hasHeightForWidth())
        self.pushButton_relative.setSizePolicy(sizePolicy)
        self.pushButton_relative.setObjectName("pushButton_relative")
        self.gridLayout.addWidget(self.pushButton_relative, 2, 3, 1, 1)
        self.verticalLayoutWidget = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(370, 60, 91, 201))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.pushButton_up = QtWidgets.QPushButton(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_up.sizePolicy().hasHeightForWidth())
        self.pushButton_up.setSizePolicy(sizePolicy)
        self.pushButton_up.setObjectName("pushButton_up")
        self.verticalLayout.addWidget(self.pushButton_up)
        self.pushButton_home = QtWidgets.QPushButton(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_home.sizePolicy().hasHeightForWidth())
        self.pushButton_home.setSizePolicy(sizePolicy)
        self.pushButton_home.setObjectName("pushButton_home")
        self.verticalLayout.addWidget(self.pushButton_home)
        self.pushButton_down = QtWidgets.QPushButton(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_down.sizePolicy().hasHeightForWidth())
        self.pushButton_down.setSizePolicy(sizePolicy)
        self.pushButton_down.setObjectName("pushButton_down")
        self.verticalLayout.addWidget(self.pushButton_down)
        self.horizontalLayoutWidget = QtWidgets.QWidget(Form)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(60, 200, 281, 80))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_BL = QtWidgets.QLabel(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_BL.sizePolicy().hasHeightForWidth())
        self.label_BL.setSizePolicy(sizePolicy)
        self.label_BL.setObjectName("label_BL")
        self.verticalLayout_3.addWidget(self.label_BL)
        self.lineEdit_BL = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit_BL.setObjectName("lineEdit_BL")
        self.verticalLayout_3.addWidget(self.lineEdit_BL)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.pushButton_BL = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_BL.sizePolicy().hasHeightForWidth())
        self.pushButton_BL.setSizePolicy(sizePolicy)
        self.pushButton_BL.setObjectName("pushButton_BL")
        self.horizontalLayout.addWidget(self.pushButton_BL)
        self.lcd_step = QtWidgets.QLCDNumber(Form)
        self.lcd_step.setGeometry(QtCore.QRect(60, 290, 131, 61))
        self.lcd_step.setSmallDecimalPoint(True)
        self.lcd_step.setProperty("value", 0.1)
        self.lcd_step.setProperty("intValue", 0)
        self.lcd_step.setObjectName("lcd_step")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Dialog"))
        self.label_position.setText(_translate("Form", "Absolute"))
        self.label_stepsize.setText(_translate("Form", "Stepsize"))
        self.pushButton_stepsize.setText(_translate("Form", "Stepsize"))
        self.pushButton_absolute.setText(_translate("Form", "Position"))
        self.label_relative.setText(_translate("Form", "Relative"))
        self.pushButton_relative.setText(_translate("Form", "Relative"))
        self.pushButton_up.setText(_translate("Form", "up"))
        self.pushButton_home.setText(_translate("Form", "Home"))
        self.pushButton_down.setText(_translate("Form", "down"))
        self.label_BL.setText(_translate("Form", "Backlash_begin"))
        self.pushButton_BL.setText(_translate("Form", "BL_correct"))

