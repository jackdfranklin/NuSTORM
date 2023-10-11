import numpy as np

import flux
import cross_sections

class Experiment:

    '''Class which contains information about the experimental setup. This includes fluxes and detector properties'''

    def __init__(self):
        '''Returns the total cross section (CC + NC) at the given energy

        Parameters
        ----------
        nu_energy: float or numpy.ndarray
                The neutrino energy or energies that the cross section will be evaluated at
        
        Returns
        -------
        cross_section: float or numpy.ndarray
        '''
        return None

   def get_events(self):
        return 0 

    def get_fluxes(self):
        return self.flux_dict

     
class ProtoDuneLike(Experiment):
    
    def __init__(self, length):

        self.length = length #Length of detector volume in cm
        self.area = (50*100)**2 #Area of face of detector volume in cm^2
        self.volume = self.length*self.area #in cm^3

        self.flux_dict = { #Store flux instances in a dictionary. Naming convention is nu/nubar for neutrino/antineutrino, followed by the flavour (e,mu,tau)
            "nubar_mu": flux.Flux("../resources/fluxes/E6spectraMuSig557Numu.txt"),
            "nu_e": flux.Flux("../resources/fluxes/E6spectraMuSig557Nue.txt")
        }

        self.cross_sec_dict = { #Store a cross section object for each flavour of neutrino
            "nubar_mu": cross_sections.Cross_Section("../resources/cross_sections/numu_Ar_xsec.txt"), #REPLACE WITH NUBAR_MU CROSS SECTIONS WHEN AVAILABLE
            "nu_e": cross_sections.Cross_Section("../resources/cross_sections/nue_Ar_xsec.txt"),
        } 

    LAr_density = 1400 * 10**-6 #kg cm^-3
    self.target_mass = self.volume*LAr_density
    self.N_Ar = 1000/40 * 6.022*10**23 * 1000 #Number of Argon atoms per kg of Argon
    self.N_e = 18*self.N_Ar #Number of electrons per kg of Argon
    self.N_p = self.N_e #Number of protons per kg of Argon
    self.N_N = 22*self.N_Ar #Number of neutrons per kg of Argon

     


    
