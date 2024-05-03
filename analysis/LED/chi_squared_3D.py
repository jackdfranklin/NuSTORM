import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.animation as animation
import networkx as nx
from matplotlib.animation import FuncAnimation, PillowWriter 
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

import sys
import os
sys.path.append(os.path.abspath("../src"))

import flux
import cross_section
import Machado_chi_squared_L_LED_L as Machado


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class Experiment:

    '''Class which contains information about the experimental setup. This includes fluxes and detector properties'''

    def get_events(self):
        return 0 

    def get_fluxes(self):
        return self.flux_dict

     
class ProtoDuneLike(Experiment):
    
    def __init__(self, L):

        #self.length = length #Length of detector volume in cm
        #self.area = (50*50)**2 #Area of face of detector volume in cm^2
        #self.volume = self.length*self.area #in cm^3
        self.volume = 7.2*6.1*7*10**6 #cm^3

        self.flux_dict = { #Store flux instances in a dictionary. Naming convention is nu/nubar for neutrino/antineutrino, followed by the flavour (e,mu,tau)
            "nubar_mu": flux.Flux("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Numu.txt", L),
            "nu_e": flux.Flux("/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Nue.txt", L)
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
        self.N_targets =  self.N_p + self.N_N

    
    #Events for standard model
    def get_events_nubar_mu(self, POT):

        Emin = 0.3
        Emax = 5.5 #20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), Emin, Emax)[0]*POT*self.N_targets
        return events
    
    def get_events_nu_e(self, POT):

        Emin = 0.3
        Emax = 5.5#20
        flux = self.flux_dict["nu_e"]   
        cross_section = self.cross_sec_dict["nu_e_Ar"]
        events = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E), Emin, Emax)[0]*POT*self.N_targets
        return events


    #Events for LED model for surviving neurinos
    def get_events_LED_nubar_mu(self, POT, a_LED, L):

        Emin = 0.3
        Emax = 5.5#20
        flux = self.flux_dict["nubar_mu"]   
        cross_section = self.cross_sec_dict["nubar_mu_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'mu', 'mu', L), Emin, Emax)[0]*POT*self.N_targets
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'mu', 'mu', L), Emin, Emax)[0]*POT*self.N_targets

        return events_NH, events_IH
    
    def get_events_LED_nu_e(self, POT, a_LED, L):

        Emin = 0.3
        Emax = 5.5#20
        flux = self.flux_dict["nu_e"]   
        cross_section = self.cross_sec_dict["nu_e_Ar"]
        events_NH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_NH(a_LED, E*10**9, 'e', 'e', L), Emin, Emax)[0]*POT*self.N_targets
        events_IH = integrate.quad(lambda E: flux.diff_flux(E)*cross_section.total_cross_section(E)*Machado.probability_IH(a_LED, E*10**9, 'e', 'e', L), Emin, Emax)[0]*POT*self.N_targets

        return events_NH, events_IH
    



#----------------------------------------------------------------------------------------------------------------------------------------------------- 
    


if __name__ == "__main__":

    N_points_LED = Machado.N_points #20
    N_points_L = Machado.N_points
    L_array = np.linspace(100, 3000, N_points_L) #keep in meters. Gets converted to cm in flux and to eV-1 in Machado_chi_squared.
    POT =  1e21 
   
    conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
    R_LED_ = np.logspace(np.log10(1e-3), np.log10(5), N_points_LED)
    R_LED = R_LED_/conversion_factor_μm_to_eV
    #print(R_LED_)
#-----------------------------------------------------------------------------------------
    def chi_squared(N_exp, N_theo):
        χ_squared = np.empty([N_points_LED, N_points_L])

        for i in range(N_points_LED):
            for j in range(N_points_L):
                χ_squared[i, j] = (N_exp[i, j] - N_theo[i, j])**2 /(N_exp[i,j])

        χ_squared = χ_squared.flatten()
        

        return χ_squared
    
#-------------------------------------------------------------------------------------------
    N_events_LED_nubar_mu_NH = np.empty([N_points_LED, N_points_L])
    N_events_LED_nubar_mu_IH = np.empty([N_points_LED, N_points_L])
    N_events_standard_nubar_mu = np.empty([N_points_LED, N_points_L])

    N_events_LED_nu_e_NH = np.empty([N_points_LED, N_points_L])
    N_events_LED_nu_e_IH = np.empty([N_points_LED, N_points_L])
    N_events_standard_nu_e = np.empty([N_points_LED, N_points_L])


    for j in range(N_points_L):   
        protoDune = ProtoDuneLike(L_array[j])
        #Standard model number of events
        N_events_nubar_mu = protoDune.get_events_nubar_mu(POT)
        N_events_nu_e = protoDune.get_events_nu_e(POT)

        for i in range(N_points_LED):  
            events = protoDune.get_events_LED_nubar_mu(POT, R_LED[i], L_array[j])
            N_events_LED_nubar_mu_NH[i, j] = events[0]
            N_events_LED_nubar_mu_IH[i, j] = events[1]
            N_events_standard_nubar_mu[i, j] = N_events_nubar_mu
            

        for i in range(N_points_LED):   
            events = protoDune.get_events_LED_nu_e(POT, R_LED[i], L_array[j])
            N_events_LED_nu_e_NH[i, j] = events[0]
            N_events_LED_nu_e_IH[i, j] = events[1]
            N_events_standard_nu_e[i, j] = N_events_nu_e

    
    
    L_array_repeated = np.tile(L_array, N_points_L) #repeats the entitre array N times
    R_LED_repeated = np.repeat(R_LED_, N_points_LED) #repeat each element after themselves

    
        
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
    def Plot_3D_chi(L_array_repeated, R_LED_repeated, N_events_LED, N_events_standard , name):
        χ_squared = chi_squared(N_events_LED, N_events_standard)
        sns.set_style("whitegrid")  

        fig = plt.figure(1)
        ax = plt.axes(projection='3d')
        trsrf = ax.plot_trisurf(L_array_repeated, R_LED_repeated, χ_squared, 
                        cmap='viridis_r', edgecolor='none', label = 'flavor')

        fig.colorbar(trsrf, ax = ax)
        
        ax.set_xlabel(r'Detector~distance~(m)')
        ax.set_ylabel(r'$R_{LED}(μm)$')
        ax.set_zlabel(r'Number~of~events')
        #ax.set_yscale('log')
        
        plt.savefig('/Users/athinavogiatzi/Documents/level 4/project/Graphs/3D graphs/3d_plot_chi-squared' + name + '.png')
        plt.show()

        return
    


    def scatter_plot_chi(L_array_repeated, R_LED_repeated, N_events_LED, N_events_standard,  name):
        χ_squared = chi_squared(N_events_LED, N_events_standard)
        data_dict = {'chi_squared': χ_squared,
             'Detector distance(m)': L_array_repeated,
             'R_LED': R_LED_repeated}
        
        data = pd.DataFrame(data_dict)
        fig, axes = plt.subplots(1, 1)
        sns.scatterplot(data=data, x='R_LED', y='chi_squared', hue='Detector distance(m)',  palette='viridis')
        axes.set_xscale('log')
        axes.grid(False)
        plt.savefig('/Users/athinavogiatzi/Documents/level 4/project/Graphs/3D graphs/scatter_plot_chi' + name + '.png')
        plt.show()

        return

    


    def plot_3D_events(L_array, R_LED_, N_events_LED, name):
        sns.set_style("whitegrid")
        N_events = N_events_LED.flatten()
        fig = plt.figure(1)
        ax = plt.axes(projection='3d')
        trsrf = ax.plot_trisurf(L_array, R_LED_, N_events, 
                        cmap='viridis_r', edgecolor='none', label = 'flavor')
        ax.set_xlabel(r'Detector~distance~(m)')
        ax.set_ylabel(r'$R_{LED}(μm)$')
        ax.set_zlabel(r'Number~of~events')
        plt.savefig('/Users/athinavogiatzi/Documents/level 4/project/Graphs/3D graphs/3d_plot_events' + name + '.png')   
        plt.show()
        return
    



    def scatter_plot_events(L_array, R_LED_, N_events_LED, name):
        N_events = N_events_LED.flatten()
        data_dict = {'N events': N_events,
             'Detector distance(m)': L_array,
             'R_LED': R_LED_}
        
        data = pd.DataFrame(data_dict)
        fig, axes = plt.subplots(1, 1)

        sns.scatterplot(data=data, x='R_LED', y='N events', hue='Detector distance(m)',  palette='viridis')
        axes.set_xscale('log')
        axes.grid(False)
        plt.savefig('/Users/athinavogiatzi/Documents/level 4/project/Graphs/3D graphs/scatter_plots_events' + name + '.png')
        plt.show()

        return


    #nu_bar_mu normal hierarchy plots
    Plot_3D_chi(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_NH, N_events_standard_nubar_mu, 'nu_mu_bar_NH')
    scatter_plot_chi(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_NH, N_events_standard_nubar_mu,  'nu_mu_bar_NH')
    plot_3D_events(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_NH,   'nu_mu_bar_NH')
    scatter_plot_events(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_NH,  'nu_mu_bar_NH')



    #nu_bar_mu inverted hierarchy plots
    Plot_3D_chi(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_IH, N_events_standard_nubar_mu, 'nu_mu_bar_IH')
    scatter_plot_chi(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_IH, N_events_standard_nubar_mu,  'nu_mu_bar_IH')
    plot_3D_events(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_IH,   'nu_mu_bar_IH')
    scatter_plot_events(L_array_repeated, R_LED_repeated, N_events_LED_nubar_mu_IH,  'nu_mu_bar_IH')
   

   #nu_e normal hierarchy plots
    Plot_3D_chi(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_NH, N_events_standard_nu_e, 'nu_e_NH')
    scatter_plot_chi(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_NH, N_events_standard_nu_e,  'nu_e_NH')
    plot_3D_events(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_NH,   'nu_e_NH')
    scatter_plot_events(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_NH,  'nu_e_NH')



    #nu_e inverted hierarchy plots
    Plot_3D_chi(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_IH, N_events_standard_nu_e, 'nu_e_IH')
    scatter_plot_chi(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_IH, N_events_standard_nu_e,  'nu_e_IH')
    plot_3D_events(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_IH,   'nu_e_IH')
    scatter_plot_events(L_array_repeated, R_LED_repeated, N_events_LED_nu_e_IH,  'nu_me_IH')
   





    



    
    

    




    
    
# fig.colorbar(trsrf, ax = ax)

        # #plt.legend( loc = 'upper left')#data, 'data points')

        # ax.view_init( azim = 163, elev = 24)

        # #data.add_legend()
        # ax.set_xlabel(r'Detector~distance~(m)')
        # ax.set_ylabel(r'$R_{LED}(μm)$')
        # ax.set_zlabel(r'Number~of~events')
        # #ax.set_yscale('log')
        # #plt.yscale('log')


        # def update(i, fig, ax):
        #     ax.view_init(elev= 5., azim=i)
        #     return fig, ax
        
        # # Writer = animation.writers['ffmpeg']
        # # writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        
        # anim = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), repeat=True, fargs=(fig, ax))
        
        # anim.save('/Users/athinavogiatzi/Documents/level 4/project/N_events_' + name + '.gif' , dpi=80, writer='ffmpeg', fps=24)
