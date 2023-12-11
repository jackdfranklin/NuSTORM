import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

import flux
import cross_section
import Machado

class Experiment:

    '''Class which contains information about the experimental setup. This includes fluxes and detector properties'''

    def get_events(self):
        return 0 

    def get_fluxes(self):
        return self.flux_dict

     
class ProtoDuneLike(Experiment):
    
    def __init__(self):

        #self.length = length #Length of detector volume in cm
        #self.area = (50*50)**2 #Area of face of detector volume in cm^2
        #self.volume = self.length*self.area #in cm^3
        self.volume = 7.2*6.1*7*10**6 #cm^3

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

        LAr_density = 1.38 #g cm^-3, density at 124kPa
        self.target_mass = self.volume*LAr_density
        self.N_Ar = self.target_mass/40 * 6.022*10**23  #Number of Argon atoms in the detector mass/m_r *Avogadro
        self.N_e = 18*self.N_Ar #Number of electrons 
        self.N_p = self.N_e #Number of protons 
        self.N_N = 22*self.N_Ar #Number of neutrons 

    
    #Events for standard model
    def get_events_nubar_mu(self, POT):

        Emin = 0.1
        Emax = 20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), Emin, Emax)[0]*POT*self.N_Ar
        return events
    
    def get_events_nu_e(self, POT):

        Emin = 0.1
        Emax = 20
        flux = self.flux_dict["nu_e"]   
        cross_section = self.cross_sec_dict["nu_e_Ar"]
        events = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), Emin, Emax)[0]*POT*self.N_Ar
        return events


    #Events for LED model for surviving neurinos
    def get_events_LED_nubar_mu(self, POT, a_LED):

        Emin = 0.1
        Emax = 20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'mu', 'mu'), Emin, Emax)[0]*POT*self.N_Ar
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'mu', 'mu'), Emin, Emax)[0]*POT*self.N_Ar

        return events_NH, events_IH
    
    def get_events_LED_nu_e(self, POT, a_LED):

        Emin = 0.1
        Emax = 20
        flux = self.flux_dict["nu_e"]   
        cross_section = self.cross_sec_dict["nu_e_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'e', 'e'), Emin, Emax)[0]*POT*self.N_Ar
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'e', 'e'), Emin, Emax)[0]*POT*self.N_Ar

        return events_NH, events_IH
    

    #Events for LED model for oscillation neutrinos
    # def get_events_LED_nubar_mu_to_nubar_e(self, POT, a_LED):

    #     Emin = 0.1
    #     Emax = 20
    #     flux = self.flux_dict["nubar_mu"]   
    #     cross_section = self.cross_sec_dict["nu_e_Ar"]
    #     events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'mu', 'e'), Emin, Emax)[0]*POT*self.N_Ar
    #     events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'mu', 'e'), Emin, Emax)[0]*POT*self.N_Ar

    #     return events_NH, events_IH
    
    # def get_events_LED_nu_e_to_nubar_mu(self, POT, a_LED):

    #     Emin = 0.1
    #     Emax = 20
    #     flux = self.flux_dict["nu_e"]   
    #     cross_section = self.cross_sec_dict["nubar_mu_Ar"]
    #     events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'e', 'mu'), Emin, Emax)[0]*POT*self.N_Ar
    #     events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'e', 'mu'), Emin, Emax)[0]*POT*self.N_Ar

    #     return events_NH, events_IH








#----------------------------------------------------------------------------------------------------------------------------------------------------- 
    


if __name__ == "__main__":


    protoDune = ProtoDuneLike()
    # print(protoDune.target_mass)
    print('number of Ar atoms:' + str(protoDune.N_Ar))

    POT =  10e7 #10e20
    N_points = Machado.N_points

    conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
    # R_LED_ = np.linspace(1e-3, 5, N_points)
    # R_LED = R_LED_/conversion_factor_μm_to_eV
    R_LED_ = np.logspace(np.log10(1e-3/conversion_factor_μm_to_eV), np.log10(5/conversion_factor_μm_to_eV), N_points)
    R_LED = 10**R_LED_

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

    
    




#Nu_bar mu figure of number of events as a function of LED length
plt.figure(1)
plt.plot(R_LED_, N_events_LED_nubar_mu_NH, label = 'NH', color = 'blue')
plt.plot(R_LED_, N_events_LED_nubar_mu_IH, label = 'IH', color = 'orange' )
plt.plot(R_LED_, N_events_standard_nubar_mu, label= 'exp', color = 'green')
plt.xlabel(r'$Log(R_{LED})$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('log(N)')
plt.title(r'$ ν_μ \rightarrow ν_μ$')
plt.show()


#Nu e figure of number of events as a function of LED length
plt.figure(2)
plt.plot(R_LED_, N_events_LED_nu_e_NH, label = 'NH', color = 'blue')
plt.plot(R_LED_, N_events_LED_nu_e_IH, label = 'IH', color = 'orange' )
plt.plot(R_LED_, N_events_standard_nu_e, label = 'exp', color = 'green')
plt.xlabel(r'$Log(R_{LED})$')
plt.ylabel('log(N)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.title(r'$ν_e \rightarrow ν_e$')
plt.show()





    
    

    




    
    

