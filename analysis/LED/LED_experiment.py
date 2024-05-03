#Length in m to eV^-1
conversion_factor_m_to_eV = 1.9732705*10**(-7)
L_1km = 1000/conversion_factor_m_to_eV
L_500m = 500/conversion_factor_m_to_eV
L_200m = 200/conversion_factor_m_to_eV

flux_dict = { #Store flux instances in a dictionary. Naming convention is nu/nubar for neutrino/antineutrino, followed by the flavour (e,mu,tau)
            "nubar_mu_1000": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D1000Numu1148.txt",
            "nu_e_1000": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D1000Nue1148.txt",

            "nubar_mu_500": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D500Numu1147.txt",
            "nu_e_500": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D500Nue1147.txt",

            "nubar_mu_200": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D200Numu1146.txt",
            "nu_e_200": "/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6D200Nue1146.txt"
        }

    #Events for LED model for surviving neurinos
    #These should be contained within the flux class ideally
    def get_events_LED_nubar_mu(self, POT, a_LED, L):
        cross_section = self.cross_sec_dict[self.flavour]
        events_NH = 0 
        events_IH = 0
        #Summing over fluxe*cross_secty
        for j in range(len(self.energies)):    
            events_NH += self.int_flux[j]*Machado.probability_NH(a_LED, self.energies[j]*1e9, 'mu', 'mu', L)*(cross_section.total_cross_section(self.energies)[j])
            events_IH += self.int_flux[j]*Machado.probability_IH(a_LED, self.energies[j]*1e9, 'mu', 'mu', L)*(cross_section.total_cross_section(self.energies)[j])
        
        events_NH = self.N_targets*POT*events_NH
        events_IH = self.N_targets*POT*events_IH
        
        return events_NH, events_IH
    
    
    def get_events_LED_nu_e(self, POT, a_LED, L):
        cross_section = self.cross_sec_dict[self.flavour]
        events_NH = 0 
        events_IH = 0
        #Summing over fluxe*cross_secty
        for j in range(len(self.energies)):    
            events_NH += self.int_flux[j]*Machado.probability_NH(a_LED, self.energies[j]*1e9, 'e', 'e', L)*(cross_section.total_cross_section(self.energies)[j])
            events_IH += self.int_flux[j]*Machado.probability_IH(a_LED, self.energies[j]*1e9, 'e', 'e', L)*(cross_section.total_cross_section(self.energies)[j])
        
        events_NH = self.N_targets*POT*events_NH
        events_IH = self.N_targets*POT*events_IH
        
        return events_NH, events_IH
    


#----------------------------------------------------------------------------------------------------------------------------------------------------- 
    
def chi_squared(N_exp, N_theo):
        χ_squared = np.empty(N_points)

        for i in range(N_points):
            χ_squared[i] = (N_exp - N_theo[i])**2 /N_exp

        return χ_squared

def Plot_N_events(R_LED, N_events_NH, N_events_IH, N_events_standard, flavor, distance ):
        
        plt.figure()
        plt.plot(R_LED, N_events_NH, label = 'NH', color = 'blue')
        plt.plot(R_LED, N_events_IH, label = 'IH', color = 'orange' )
        plt.plot(R_LED, N_events_standard, label= 'SM', color = 'green')
        plt.xlabel(r'$R_{LED}(μm)$')
        plt.xscale('log')
        plt.legend()
        plt.ylabel(r'$Number~of~$' + flavor + r"$events$" + distance)
        plt.show()

        return


def Plot_chi_squared(R_LED, N_events_NH, N_events_IH, N_events_standard, flavor, distance):
       
    χ_squared_NH = chi_squared(N_events_standard, N_events_NH)
    χ_squared_IH = chi_squared(N_events_standard, N_events_IH)

    plt.figure()
    plt.plot(R_LED, χ_squared_NH, label = 'NH', color = 'blue')
    plt.plot(R_LED, χ_squared_IH, label = 'IH', color = 'orange')
    plt.legend()

    plt.xlabel(r'$R_{LED}(μm)$')  #l_0 = 1μm
    plt.ylabel(r'$\chi^2~difference~for~$'+ flavor + r"$events$" + distance)
    plt.xscale('log')
        
    plt.show()
    return
#----------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    #Constants
    POT =  1e21 #1e30 #2e18   #10e7 
    N_points = 1000

    conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
    # R_LED_ = np.linspace(1e-3, 5, N_points)
    # R_LED = R_LED_/conversion_factor_μm_to_eV
    R_LED_ = np.logspace(np.log10(1e-3), np.log10(5), N_points) #Plot with this
    R_LED = R_LED_/conversion_factor_μm_to_eV #Entry to get_events functions



#Creating ProtoDune class instances
    
#For 1000m
    #nubar_mu
    protoDune_SM_nubar_mu_1km = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_1000"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nubar_mu_1km = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_1000"]) #flavour, flux_file are the variables of LED_Flux

    #nu_e
    protoDune_SM_nu_e_1km = ProtoDuneLike("nu_e", flux_dict["nu_e_1000"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nu_e_1km = ProtoDuneLike("nu_e", flux_dict["nu_e_1000"]) #flavour, flux_file are the variables of LED_Flux

    
#For 500m
    #nubar_mu
    protoDune_SM_nubar_mu_500m = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_500"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nubar_mu_500m = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_500"]) #flavour, flux_file are the variables of LED_Flux

    #nu_e
    protoDune_SM_nu_e_500m = ProtoDuneLike("nu_e", flux_dict["nu_e_500"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nu_e_500m = ProtoDuneLike("nu_e", flux_dict["nu_e_500"]) #flavour, flux_file are the variables of LED_Flux


#For 200m
    #nubar_mu
    protoDune_SM_nubar_mu_200m = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_200"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nubar_mu_200m = ProtoDuneLike("nubar_mu", flux_dict["nubar_mu_200"]) #flavour, flux_file are the variables of LED_Flux

    #nu_e
    protoDune_SM_nu_e_200m = ProtoDuneLike("nu_e", flux_dict["nu_e_200"] ) #flavour, flux_file are the variables of Flux
    protoDune_LED_nu_e_200m = ProtoDuneLike("nu_e", flux_dict["nu_e_200"]) #flavour, flux_file are the variables of LED_Flux

    
    
    

   

#calculating number of events
#1km
    #Standard Model events
    N_events_SM_nubar_mu_1km = protoDune_SM_nubar_mu_1km.get_tot_Ar_events(POT)
    N_events_SM_nu_e_1km = protoDune_SM_nu_e_1km.get_tot_Ar_events(POT)

    #
    N_events_LED_nubar_mu_NH_1km = np.empty(N_points)
    N_events_LED_nubar_mu_IH_1km = np.empty(N_points)
    N_events_standard_nubar_mu_1km = np.empty(N_points)

    N_events_LED_nu_e_NH_1km = np.empty(N_points)
    N_events_LED_nu_e_IH_1km = np.empty(N_points)
    N_events_standard_nu_e_1km = np.empty(N_points)



    
    for i in range(N_points):
        #nubar_mu
        #print(R_LED[i])
        #print(L_1km)
        protoDune_LED_nubar_mu = protoDune_LED_nubar_mu_1km.get_events_LED_nubar_mu(POT, R_LED[i], L_1km)
        N_events_LED_nubar_mu_NH_1km[i] = protoDune_LED_nubar_mu[0]
        N_events_LED_nubar_mu_IH_1km[i] = protoDune_LED_nubar_mu[1]
        N_events_standard_nubar_mu_1km[i] = N_events_SM_nubar_mu_1km
        
        #nu_e
        N_events_LED_nu_e = protoDune_LED_nu_e_1km.get_events_LED_nu_e(POT, R_LED[i], L_1km)
        N_events_LED_nu_e_NH_1km[i] = N_events_LED_nu_e[0]
        N_events_LED_nu_e_IH_1km[i] = N_events_LED_nu_e[1]
        N_events_standard_nu_e_1km[i] = N_events_SM_nu_e_1km


#500m
    N_events_SM_nubar_mu_500m = protoDune_SM_nubar_mu_500m.get_tot_Ar_events(POT)
    N_events_SM_nu_e_500m = protoDune_SM_nu_e_500m.get_tot_Ar_events(POT)
    N_events_LED_nubar_mu_NH_500m = np.empty(N_points)
    N_events_LED_nubar_mu_IH_500m = np.empty(N_points)
    N_events_standard_nubar_mu_500m = np.empty(N_points)

    N_events_LED_nu_e_NH_500m = np.empty(N_points)
    N_events_LED_nu_e_IH_500m = np.empty(N_points)
    N_events_standard_nu_e_500m = np.empty(N_points)



    
    for i in range(N_points):
        #nubar_mu
        protoDune_LED_nubar_mu = protoDune_LED_nubar_mu_500m.get_events_LED_nubar_mu(POT, R_LED[i], L_500m)
        N_events_LED_nubar_mu_NH_500m[i] = protoDune_LED_nubar_mu[0]
        N_events_LED_nubar_mu_IH_500m[i] = protoDune_LED_nubar_mu[1]
        N_events_standard_nubar_mu_500m[i] = N_events_SM_nubar_mu_500m
        
        #nu_e
        N_events_LED_nu_e = protoDune_LED_nu_e_500m.get_events_LED_nu_e(POT, R_LED[i], L_500m)
        N_events_LED_nu_e_NH_500m[i] = N_events_LED_nu_e[0]
        N_events_LED_nu_e_IH_500m[i] = N_events_LED_nu_e[1]
        N_events_standard_nu_e_500m[i] = N_events_SM_nu_e_500m


#200m
    N_events_SM_nubar_mu_200m = protoDune_SM_nubar_mu_200m.get_tot_Ar_events(POT)
    N_events_SM_nu_e_200m = protoDune_SM_nu_e_200m.get_tot_Ar_events(POT)
    N_events_LED_nubar_mu_NH_200m = np.empty(N_points)
    N_events_LED_nubar_mu_IH_200m = np.empty(N_points)
    N_events_standard_nubar_mu_200m = np.empty(N_points)

    N_events_LED_nu_e_NH_200m = np.empty(N_points)
    N_events_LED_nu_e_IH_200m = np.empty(N_points)
    N_events_standard_nu_e_200m = np.empty(N_points)



    
    for i in range(N_points):
        #nubar_mu
        protoDune_LED_nubar_mu = protoDune_LED_nubar_mu_200m.get_events_LED_nubar_mu(POT, R_LED[i], L_200m)
        N_events_LED_nubar_mu_NH_200m[i] = protoDune_LED_nubar_mu[0]
        N_events_LED_nubar_mu_IH_200m[i] = protoDune_LED_nubar_mu[1]
        N_events_standard_nubar_mu_200m[i] = N_events_SM_nubar_mu_200m
        
        #nu_e
        N_events_LED_nu_e = protoDune_LED_nu_e_200m.get_events_LED_nu_e(POT, R_LED[i], L_200m)
        N_events_LED_nu_e_NH_200m[i] = N_events_LED_nu_e[0]
        N_events_LED_nu_e_IH_200m[i] = N_events_LED_nu_e[1]
        N_events_standard_nu_e_200m[i] = N_events_SM_nu_e_200m


    #1km 

    #nubar_mu  
    #Plot_N_events(R_LED_, N_events_LED_nubar_mu_NH_1km, N_events_LED_nubar_mu_IH_1km, N_events_standard_nubar_mu_1km, r'$\bar \nu_{\mu}~$', r'1km')
    Plot_chi_squared(R_LED_, N_events_LED_nubar_mu_NH_1km, N_events_LED_nubar_mu_IH_1km, N_events_SM_nubar_mu_1km, r'$\bar \nu_{\mu}~$', r'1km')
    #nu_e
    # Plot_N_events(R_LED_, N_events_LED_nu_e_NH_1km, N_events_LED_nu_e_IH_1km, N_events_standard_nu_e_1km, r'$ \nu_e~$', r'1km')
    Plot_chi_squared(R_LED_, N_events_LED_nu_e_NH_1km, N_events_LED_nu_e_IH_1km, N_events_SM_nu_e_1km, r'$ \nu_e~$', r'1km')


    #500

    #nubar_mu  
    # Plot_N_events(R_LED_, N_events_LED_nubar_mu_NH_500m, N_events_LED_nubar_mu_IH_500m, N_events_standard_nubar_mu_500m, r'$\bar \nu_{\mu}~$', r'500m')
    Plot_chi_squared(R_LED_, N_events_LED_nubar_mu_NH_500m, N_events_LED_nubar_mu_IH_500m, N_events_SM_nubar_mu_500m, r'$\bar \nu_{\mu}~$', r'500m')

    #nu_e
    #Plot_N_events(R_LED_, N_events_LED_nu_e_NH_500m, N_events_LED_nu_e_IH_500m, N_events_standard_nu_e_500m, r'$ \nu_e~$', r'500km')
    Plot_chi_squared(R_LED_, N_events_LED_nu_e_NH_500m, N_events_LED_nu_e_IH_500m, N_events_SM_nu_e_500m, r'$ \nu_e~$', r'500km')

    #200

    #nubar_mu  
    #Plot_N_events(R_LED_, N_events_LED_nubar_mu_NH_200m, N_events_LED_nubar_mu_IH_200m, N_events_standard_nubar_mu_200m, r'$\bar \nu_{\mu}~$', r'200m')
    Plot_chi_squared(R_LED_, N_events_LED_nubar_mu_NH_200m, N_events_LED_nubar_mu_IH_200m, N_events_SM_nubar_mu_200m, r'$\bar \nu_{\mu}~$', r'200m')

    #nu_e
    #Plot_N_events(R_LED_, N_events_LED_nu_e_NH_200m, N_events_LED_nu_e_IH_200m, N_events_standard_nu_e_200m, r'$ \nu_e~$', r'200km')
    Plot_chi_squared(R_LED_, N_events_LED_nu_e_NH_200m, N_events_LED_nu_e_IH_200m, N_events_SM_nu_e_200m, r'$ \nu_e~$', r'200km')

    
    




    # #Nu_bar mu figure of number of events as a function of LED length
    # plt.figure(1)
    # plt.plot(R_LED_, N_events_LED_nubar_mu_NH, label = 'NH', color = 'blue')
    # plt.plot(R_LED_, N_events_LED_nubar_mu_IH, label = 'IH', color = 'orange' )
    # plt.plot(R_LED_, N_events_standard_nubar_mu, label= 'SM', color = 'green')
    # plt.xlabel(r'$R_{LED}(μm)$')
    # plt.xscale('log')
    # plt.legend()
    # plt.ylabel(r'$Number~of~ \bar \nu_{\mu}~events$')
    # plt.show()


    # #Nu e figure of number of events as a function of LED length
    # plt.figure(2)
    # plt.plot(R_LED_, N_events_LED_nu_e_NH, label = 'NH', color = 'blue')
    # plt.plot(R_LED_, N_events_LED_nu_e_IH, label = 'IH', color = 'orange' )
    # plt.plot(R_LED_, N_events_standard_nu_e, label = 'SM', color = 'green')
    # plt.xlabel(r'$R_{LED}(μm)$')
    # plt.ylabel(r'$Number~of~ \nu_{e}~events$')
    # plt.xscale('log')
    # #plt.yscale('log')
    # plt.legend()
    # plt.show()


    

    # Plot(R_LED, N_events_LED_nubar_mu_NH, N_events_LED_nubar_mu_IH, N_events_nubar_mu, r'$\bar \nu_{\mu} $')
    # Plot(R_LED, N_events_LED_nu_e_NH, N_events_LED_nu_e_IH, N_events_nu_e, r'$ \nu_e$')
        



    
    

    




    
    

