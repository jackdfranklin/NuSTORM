import numpy as np
import scipy.integrate as integrate

import flux
import cross_section

class Experiment:

    '''Class which contains information about the experimental setup. This includes fluxes and detector properties'''

    def get_events(self):
        return 0 

    def get_fluxes(self):
        return self.flux_dict

     
class ProtoDuneLike(Experiment):
    
    def __init__(self, length):

        self.length = length #Length of detector volume in cm
        self.area = (50*50)**2 #Area of face of detector volume in cm^2
        self.volume = self.length*self.area #in cm^3

        self.flux_dict = { #Store flux instances in a dictionary. Naming convention is nu/nubar for neutrino/antineutrino, followed by the flavour (e,mu,tau)
            "nubar_mu": flux.Flux("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Numu.txt"),
            "nu_e": flux.Flux("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Nue.txt")
        }

        self.cross_sec_dict = { #Store a cross section object for each flavour of neutrino
            "nubar_mu_Ar": cross_section.Cross_Section("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/cross_sections/numu_Ar_xsec.txt"), #REPLACE WITH NUBAR_MU CROSS SECTIONS WHEN AVAILABLE
            #"nubar_mu_e":
            "nu_e_Ar": cross_section.Cross_Section("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/cross_sections/nue_Ar_xsec.txt"),
            #"nu_e_e":
        } 

        LAr_density = 1400 * 10**-6 #kg cm^-3
        self.target_mass = self.volume*LAr_density
        self.N_Ar = 1000/40 * 6.022*10**23 * self.target_mass#Number of Argon atoms in the detector
        self.N_e = 18*self.N_Ar #Number of electrons 
        self.N_p = self.N_e #Number of protons 
        self.N_N = 22*self.N_Ar #Number of neutrons 

    def get_events(self, POT):

        Emin = 0.1
        Emax = 20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), Emin, Emax)[0]*POT*self.N_Ar
        return events


if __name__ == "__main__":

    protoDune = ProtoDuneLike(5000)
    
    N_events = protoDune.get_events(10**7)
    print(str(N_events) + " expected!")
