import numpy as np
import pandas as pd
from scipy import interpolate

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

