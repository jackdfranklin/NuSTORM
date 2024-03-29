import numpy as np
import pandas as pd
from scipy import interpolate
import Machado

# Pass flavour as nubar_mu or nu_e into the flux class
class Flux:

    '''Handler class which contains the raw flux data'''

    def __init__(self, flavour, flux_file):
        '''Initialises object with energy midpoints and fluxes from given flux file'''

        self.flavour = flavour
        #Load flux file into dataframe
        raw_data = np.loadtxt(flux_file, skiprows=1)

        #Store values of energy midpoints and integrated flux
        self.energies = raw_data[:,0]
        self.int_flux = raw_data[:,1]*1e4

        self.deltaE = self.energies[1]-self.energies[0]

        #Dictionary to keep track of units
        self.units = {
            "Energy":"GeV",
            "Integrated Flux":"cm^-2 POT^-1"}

    def get_flux(self):

        return self.int_flux

    def get_units(self):
        '''Prints out the stored values and their units'''

        print(self.units)




# Inherits get_flux and get_units functions
# L is given in meters and converted to cm in the __init__ function
class LED_Flux(Flux):

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









# self.deltaE = 0.06
        #Calculate differential flux from integrated flux
        # diff_flux = self.int_flux/self.deltaE

        
        

        #Interpolate differential flux to overcome binning differences
        # self.interp_diff_flux = interpolate.UnivariateSpline(self.energies, diff_flux)

    # def diff_flux(self, nu_energy):
    #     '''Wrapper function for the interpolated differential flux, energy should be in GeV'''

    #     return self.interp_diff_flux(nu_energy)
