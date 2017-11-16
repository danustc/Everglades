'''
This is the control program package designed for the Open-top SPIM, Everglades in Huang lab at UCSF. The program is evolved from the InControl Package designed by Dr. Ryan McGorty and updated/maintained by Dan Xie.

Requirements:   Python 3.6
                PyQt5

Last update: 11/07/2017
'''
import sys
import os
from PyQt5 import QtGui
import inLib

class Everglades(object):
    '''
    Add more instructions here.
    '''

    def __init__(self):
        print('Everglades: starting up ...')
        if len(sys.argv) > 1:
            settings_filename = sys.argv[1]
        else:
            print('Loading default settings.')
            settings_filename = 'settings_default.yaml'
        self._settings_ = inLib.load_settings(settings_filename)



