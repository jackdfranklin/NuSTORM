import numpy as np

from scipy import interpolate

class Flux:

    '''Handler class which contains the raw flux data and an interpolated differential flux'''

    def __init__(self, flux_file):
        '''Initialises object with energy midpoints and fluxes from given flux file, then calculates
        and interpolates the differential flux'''

        #Load flux file into numpy array
        raw_data = np.loadtxt(flux_file, skiprows=1)

        #Store values of energy midpoints and integrated flux
        self.energies = raw_data[:,0]
        self.int_flux = raw_data[:,1]

        self.deltaE = 0.06
        #Calculate differential flux from integrated flux
        self.diff_flux = self.int_flux/self.deltaE

        #Dictionary to keep track of units
        self.units = dict{
            "Energy":"GeV",
            "Integrated Flux":"cm^-2 POT^-1",
            "Differential Flux":"cm^-2 GeV^-1 POT^-1"}

        #Interpolate differential flux to overcome binning differences
        self.interp_diff_flux = interpolate.UnivariateSpline(self.energies, self.diff_flux)

    def diff_flux(self, nu_energy);
        '''Wrapper function for the interpolated differential flux, energy should be in GeV'''

        return self.interp_diff_flux(nu_energy)

    def get_units(self):
        '''Prints out the stored values and their units'''

        print(self.units)
