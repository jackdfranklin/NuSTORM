import numpy as np
import pandas as pd
from scipy import interpolate

import sys
import os
sys.path.append(os.path.abspath("../src"))

import flux

import Machado

# Inherits get_flux and get_units functions
# L is given in meters and converted to cm in the __init__ function
class LED_Flux(flux.Flux):

    '''Handler class which contains the raw flux data and an interpolated differential flux'''

    def __init__(self, flavour, flux_file):
        '''Initialises object with energy midpoints and fluxes from given flux file, then calculates
        and interpolates the differential flux'''
        

        self.flavour = flavour
        #Load flux file into numpy array
        raw_data = np.loadtxt(flux_file, skiprows=1)

        #Store values of energy midpoints and integrated flux
        self.energies = raw_data[:,0]
        self.int_flux = raw_data[:,1]*1e4


        #Dictionary to keep track of units
        self.units = {
            "Energy":"GeV",
            "Integrated Flux":"cm^-2 POT^-1",
            "Differential Flux":"cm^-2 GeV^-1 POT^-1"}

    def get_flux(self):

        return self.int_flux
    
       
    def get_units(self):
        '''Prints out the stored values and their units'''

        print(self.units)


