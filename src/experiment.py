import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

import flux
import cross_section
#import Machado

L = 2920 #meters

class Experiment:

    '''Class which contains information about the experimental setup. This includes fluxes and detector properties'''

    def get_events(self):
        return 0 

    def get_fluxes(self):
        return self.flux_dict

     
class ProtoDuneLike(Experiment):
    
    def __init__(self, neutrino_flux):

        #self.length = length #Length of detector volume in cm
        #self.area = (50*50)**2 #Area of face of detector volume in cm^2
        #self.volume = self.length*self.area #in cm^3
        #self.volume = 7.2*6.1*7*10**6 #cm^3
        self.volume = 50*50*100*10**6 #cm^3

        self.flux = neutrino_flux

        directory = "../resources/cross_sections/"

        self.cross_sec_dict = { #Store a cross section object for each flavour of neutrino
            "nubar_mu_Ar": cross_section.Cross_Section(directory + "nu_mu_bar_Ar40.csv"),
            "nu_e_Ar": cross_section.Cross_Section(directory + "nu_e_Ar40.csv")
        } 

        LAr_density = 1.38 #g cm^-3, density at 124kPa
        self.target_mass = 100*10**6 #self.volume*LAr_density
        self.N_Ar = self.target_mass/39.948 * 6.022*10**23  #Number of Argon atoms in the detector mass/m_r *Avogadro
        self.N_e = 18*self.N_Ar #Number of electrons 
        self.N_p = self.N_e #Number of protons 
        self.N_n = 22*self.N_Ar #Number of neutrons
        self.N_targets =  self.N_p + self.N_n

    
    #Events for standard model
    def get_tot_Ar_events(self, POT):

        cross_section = self.cross_sec_dict[self.flux.flavour]
        events = self.N_Ar*POT*self.flux.get_flux()*(cross_section.CC_cross_section(self.flux.energies)+cross_section.NC_cross_section(self.flux.energies))
        return events

    def get_CC_Ar_events(self, POT):

        cross_section = self.cross_sec_dict[self.flux.flavour]
        events = self.N_Ar*POT*self.flux.get_flux()*cross_section.CC_cross_section(self.flux.energies)
        return events

    def get_NC_Ar_events(self, POT):

        cross_section = self.cross_sec_dict[self.flux.flavour]
        events = self.N_Ar*POT*self.flux.get_flux()*cross_section.NC_cross_section(self.flux.energies)
        return events


    #Events for LED model for surviving neurinos
    #These should be contained within the flux class ideally
    def get_events_LED_nubar_mu(self, POT, a_LED):

        Emin = 0.3
        Emax = 5.5#20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'mu', 'mu'), Emin, Emax)[0]*POT*self.N_targets
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'mu', 'mu'), Emin, Emax)[0]*POT*self.N_targets

        return events_NH, events_IH
    
    def get_events_LED_nu_e(self, POT, a_LED):

        Emin = 0.3
        Emax = 5.5#20
        flux = self.flux_dict["nu_e"]   
        cross_section = self.cross_sec_dict["nu_e_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'e', 'e'), Emin, Emax)[0]*POT*self.N_targets
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'e', 'e'), Emin, Emax)[0]*POT*self.N_targets

        return events_NH, events_IH
    


    

    def get_event_bins_nubar_mu(self, POT, a_LED, num_bins):
        Emin = 0.1
        Emax = 20

        # Create an array of energy bins
        energy_bins = np.linspace(Emin, Emax, num_bins + 1)

        # Initialize arrays to store the number of events in each bin for NH and IH
        events_NH = np.zeros(num_bins)
        events_IH = np.zeros(num_bins)
        events_standard = np.zeros(num_bins)
        bin_midpoints = np.zeros(num_bins)

        flux = self.flux_dict["nubar_mu"]
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]

        for i in range(num_bins):
            # Define the integration limits for each bin
            bin_min = energy_bins[i]
            bin_max = energy_bins[i + 1]
            bin_midpoints[i] = bin_min + (bin_max-bin_min)/2

            # Integrate over the current energy bin for NH and IH
            events_NH[i] = integrate.quad(lambda E: flux.diff_flux(E) * cross_section.total_cross_section(E) * Machado.probability_NH(a_LED, E * 10**9, 'mu', 'mu'), bin_min, bin_max)[0] * POT * self.N_Ar
            events_IH[i] = integrate.quad(lambda E: flux.diff_flux(E) * cross_section.total_cross_section(E) * Machado.probability_IH(a_LED, E * 10**9, 'mu', 'mu'), bin_min, bin_max)[0] * POT * self.N_Ar
            events_standard[i] = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), bin_min, bin_max)[0]*POT*self.N_Ar

        return bin_midpoints, events_NH, events_IH, events_standard
    

    def get_event_bins_nu_e(self, POT, a_LED, num_bins):
        Emin = 0.1
        Emax = 20

        # Create an array of energy bins
        energy_bins = np.linspace(Emin, Emax, num_bins + 1)

        # Initialize arrays to store the number of events in each bin for NH and IH
        events_NH = np.zeros(num_bins)
        events_IH = np.zeros(num_bins)
        events_standard = np.zeros(num_bins)
        bin_midpoints = np.zeros(num_bins)

        flux = self.flux_dict["nu_e"]
        cross_section = self.cross_sec_dict["nu_e_Ar"]

        for i in range(num_bins):
            # Define the integration limits for each bin
            bin_min = energy_bins[i]
            bin_max = energy_bins[i + 1]
            bin_midpoints[i] = bin_min + (bin_max-bin_min)/2


            # Integrate over the current energy bin for NH and IH
            events_NH[i] = integrate.quad(lambda E: flux.diff_flux(E) * cross_section.total_cross_section(E) * Machado.probability_NH(a_LED, E * 10**9, 'e', 'e'), bin_min, bin_max)[0] * POT * self.N_Ar
            events_IH[i] = integrate.quad(lambda E: flux.diff_flux(E) * cross_section.total_cross_section(E) * Machado.probability_IH(a_LED, E * 10**9, 'e', 'e'), bin_min, bin_max)[0] * POT * self.N_Ar
            events_standard[i] = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), bin_min, bin_max)[0]*POT*self.N_Ar

        return bin_midpoints, events_NH, events_IH, events_standard
    

    



#----------------------------------------------------------------------------------------------------------------------------------------------------- 
    


if __name__ == "__main__":


    protoDune = ProtoDuneLike()
    # print(protoDune.target_mass)
    print('number of Ar atoms:' + str(protoDune.N_Ar))

    POT =  1e21 #1e30 #2e18   #10e7 
    N_points = Machado.N_points

    conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
    # R_LED_ = np.linspace(1e-3, 5, N_points)
    # R_LED = R_LED_/conversion_factor_μm_to_eV
    R_LED_ = np.logspace(np.log10(1e-3), np.log10(5), N_points)
    R_LED = R_LED_/conversion_factor_μm_to_eV

      #5000
    
    N_events_nubar_mu = protoDune.get_events_nubar_mu(POT)
    N_events_nu_e = protoDune.get_events_nu_e(POT)

    file_standard = open('src/N_events_standard.txt', 'w')
   
    print('Number of events ν_μ_bar:' + str(N_events_nubar_mu))
    print('Number of events ν_e:' + str(N_events_nu_e))

    file_standard.write(str(N_events_nubar_mu)+ '-' + str(N_events_nu_e))

    # N_events_LED = protoDune.get_events_LED_nubar_mu(POT, 5/1.9732705*10**(-1))



    file_nubar_mu = open('src/N_events_nubar_mu.txt', 'w')
    file_nu_e = open('src/N_events_nu_e.txt', 'w')



    N_events_LED_nubar_mu_NH = np.empty(N_points)
    N_events_LED_nubar_mu_IH = np.empty(N_points)
    N_events_standard_nubar_mu = np.empty(N_points)


    for i in range(N_points):
        N_events_LED_nubar_mu_NH[i] = protoDune.get_events_LED_nubar_mu(POT, R_LED[i])[0]
        N_events_LED_nubar_mu_IH[i] = protoDune.get_events_LED_nubar_mu(POT, R_LED[i])[1]
        N_events_standard_nubar_mu[i] = N_events_nubar_mu
        entry = str(R_LED_[i]) + '-' + str(N_events_LED_nubar_mu_NH[i])+ '-' + str(N_events_LED_nubar_mu_IH[i]) +'\n'
        file_nubar_mu.write(entry)


    
    
    N_events_LED_nu_e_NH = np.empty(N_points)
    N_events_LED_nu_e_IH = np.empty(N_points)
    N_events_standard_nu_e = np.empty(N_points)


    for i in range(N_points):
        N_events_LED_nu_e_NH[i] = protoDune.get_events_LED_nu_e(POT, R_LED[i])[0]
        N_events_LED_nu_e_IH[i] = protoDune.get_events_LED_nu_e(POT, R_LED[i])[1]
        N_events_standard_nu_e[i] = N_events_nu_e

        entry = str(R_LED_[i]) + '-' + str(N_events_LED_nu_e_NH[i])+ '-' + str(N_events_LED_nu_e_IH[i]) +'\n'
        file_nu_e.write(entry)

    
    
    def chi_squared(N_exp, N_theo):
        χ_squared = np.empty(N_points)

        for i in range(N_points):
            χ_squared[i] = (N_exp - N_theo[i])**2 /N_exp

        return χ_squared
    




    #Nu_bar mu figure of number of events as a function of LED length
    plt.figure(1)
    plt.plot(R_LED_, N_events_LED_nubar_mu_NH, label = 'NH', color = 'blue')
    plt.plot(R_LED_, N_events_LED_nubar_mu_IH, label = 'IH', color = 'orange' )
    plt.plot(R_LED_, N_events_standard_nubar_mu, label= 'SM', color = 'green')
    plt.xlabel(r'$R_{LED}(μm)$')
    plt.xscale('log')
    plt.legend()
    plt.ylabel(r'$Number~of~ \bar \nu_{\mu}~events$')
    plt.show()


    #Nu e figure of number of events as a function of LED length
    plt.figure(2)
    plt.plot(R_LED_, N_events_LED_nu_e_NH, label = 'NH', color = 'blue')
    plt.plot(R_LED_, N_events_LED_nu_e_IH, label = 'IH', color = 'orange' )
    plt.plot(R_LED_, N_events_standard_nu_e, label = 'SM', color = 'green')
    plt.xlabel(r'$R_{LED}(μm)$')
    plt.ylabel(r'$Number~of~ \nu_{e}~events$')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.show()





    def Plot(R_LED, N_events_NH, N_events_IH, N_events_standard, flavor):
        
        χ_squared_NH = chi_squared(N_events_standard, N_events_NH)
        χ_squared_IH = chi_squared(N_events_standard, N_events_IH)

        plt.figure()
        plt.plot(R_LED, χ_squared_NH, label = 'NH', color = 'blue')
        plt.plot(R_LED, χ_squared_IH, label = 'IH', color = 'orange')
        plt.legend()

        plt.xlabel(r'$R_{LED}(μm)$')  #l_0 = 1μm
        plt.ylabel(r'$\chi^2~difference~for~$'+ flavor + " events")
        plt.xscale('log')
        #plt.title(title)
        plt.show()
        return

    Plot(R_LED, N_events_LED_nubar_mu_NH, N_events_LED_nubar_mu_IH, N_events_nubar_mu, r'$\bar \nu_{\mu} $')
    Plot(R_LED, N_events_LED_nu_e_NH, N_events_LED_nu_e_IH, N_events_nu_e, r'$ \nu_e$')
        



    
    

    




    
    

