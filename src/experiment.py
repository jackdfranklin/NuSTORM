import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

from flux import Flux
import cross_section

class ProtoDuneLike():
    
    def __init__(self, flavour, flux_file): #Pass flavor and flux file as inputs

        self.flux = Flux(flavour, flux_file)
        #self.length = length #Length of detector volume in cm
        #self.area = (50*50)**2 #Area of face of detector volume in cm^2
        self.volume = 50*50*100*10**6 #cm^3

        directory = "../resources/cross_sections/"

        self.cross_sec_dict = { #Store a cross section object for each flavour of neutrino
            "nubar_mu": cross_section.Cross_Section(directory + "nu_mu_bar_Ar40.csv"),
            "nu_e": cross_section.Cross_Section(directory + "nu_e_Ar40.csv")
        } 
        

        LAr_density = 1.38 #g cm^-3, density at 124kPa
        self.target_mass = 100*10**6 #self.volume*LAr_density
        self.N_Ar = self.target_mass/39.948 * 6.022*10**23  #Number of Argon atoms in the detector mass/m_r *Avogadro
        self.N_e = 18*self.N_Ar #Number of electrons 
        self.N_p = self.N_e #Number of protons 
        self.N_n = 22*self.N_Ar #Number of neutrons
        self.N_targets =  self.N_p + self.N_n #number of nucleons

    #Changed N_Ar to N_targets 
    #Events for standard model
    def get_tot_Ar_events(self, POT):

        return self.get_CC_Ar_events(POT) + self.get_NC_Ar_events(POT);

    def get_CC_Ar_events(self, POT):

        cross_section = self.cross_sec_dict[self.flux.flavour]
        events = self.N_targets*POT*self.flux.get_flux()*cross_section.CC_cross_section(self.flux.energies)
        return events

    def get_NC_Ar_events(self, POT):

        cross_section = self.cross_sec_dict[self.flux.flavour]
        events = self.N_targets*POT*self.flux.get_flux()*cross_section.NC_cross_section(self.flux.energies)
        return events

