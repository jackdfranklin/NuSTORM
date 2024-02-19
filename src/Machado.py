import numpy as np
import cmath
import matplotlib.pyplot as plt




#--------------------------------------------------------------------------------------------------------------------------
#CONSTANTS

#Normal hierarchy m3>m2>m1=m0
#Inverted hierarchy  m2>m1>m3=m0

m0 = 0 #lightest mass5*10**(-2)#
N = 4
N_points = 1000

#mass differences in eVs
delta_21 = 7.53*10**(-5)
delta_31 = 2.51*10**(-3)

#circular radius of extra dimensions in μm conveted to eV^-1
#a_LED = 0.5/(1.9732705*10**(-1))
conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
R_LED_ = np.logspace(np.log10(1e-3), np.log10(5), N_points)
R_LED = R_LED_/conversion_factor_μm_to_eV
# R_LED_ = np.linspace(1e-3, 5, N_points)
# R_LED = R_LED_/conversion_factor_μm_to_eV

#Length in m to eV^-1
conversion_factor_m_to_eV = 1.9732705*10**(-7)
L_1km = 1000/conversion_factor_m_to_eV
L_500m = 500/conversion_factor_m_to_eV
L_200m = 200/conversion_factor_m_to_eV





#----------------------------------------------------------------------------------------------------------------------


#Determining masses for NH and IH

NH_dict = {
    "m1": m0,
    "m2": np.sqrt(delta_21),
    "m3": np.sqrt(delta_31)
}

IH_dict = {
    "m3": m0,
    "m1": np.sqrt(delta_31),
    "m2": np.sqrt(delta_31+delta_21)   
}



#---------------------------------------------------------------------------------------------------------------------
#PMNS matrix
#Dirac phase taken to be π
delta = -1
x = 0.0178


sin_squared_theta_12 = 0.32  
sin_squared_theta_23 = 0.5  
cos_squared_theta_12 = 1-sin_squared_theta_12
cos_squared_theta_23 = 1-sin_squared_theta_23

s12 = np.sqrt(sin_squared_theta_12)
s23 = np.sqrt(sin_squared_theta_23)
s13 = np.sqrt(x)             

c12 = np.sqrt(cos_squared_theta_12)
c23 = np.sqrt(cos_squared_theta_23)
c13 = np.sqrt(1-x)   


PMNS = np.empty((3,3))



PMNS[0,0] = c12*c13
PMNS[0,1] = s12*c13
PMNS[0,2] = s13*delta

PMNS[1,0] = -s12*c23-c12*s23*s13*delta
PMNS[1,1] = c12*c23 -s12*s23*s13*delta
PMNS[1,2] = s23*c13

PMNS[2,0] = s12*s23 -c12*c23*s13*delta
PMNS[2,1] = -c12*s23 -s12*c23*s13*delta
PMNS[2,2] = c23*c13

PMNS_dict = {
    "e1": PMNS[0,0],
    "e2": PMNS[0,1],
    "e3": PMNS[0,2],

    "mu1": PMNS[1,0],
    "mu2": PMNS[1,1],
    "mu3": PMNS[1,2],

    "tau1": PMNS[2,0],
    "tau2": PMNS[2,1],
    "tau3": PMNS[2,2]}
    





#------------------------------------------------------------------------------------------------------------------
#FUNCTIONS FOR EXTRA DIMENSIONS

def ξ_i(m_i, a):
    return np.sqrt(2)*m_i*a


#Calculate ξs depending on the hirarchy
def calculate_ξ(dictionary, a):

    ξ = np.empty((3))

    for i in range(1,4):
        mi = "m{}".format(i)
        ξ[i-1] = ξ_i(dictionary[mi],a)

    
    return ξ


def η(N, ξ_i):
    return (N+1/2)*(ξ_i**2)



def Amplitude(alpha, beta, L, E, a, eigenvalue_array, eigenvector_matrix):
    total_sum = 0
    N_max = 5  
    
    for i in range(1, 4):
        for J in range(1, 4):
            for k in range(1, 4):
                for N in range(0, N_max):
                    eigenvector = eigenvector_matrix[:, N * 3:(N + 1) * 3]
                    eigenvalue = eigenvalue_array[N * 3 + J - 1]
                    
                    total_sum += (
                        PMNS_dict["{}{}".format(alpha, i)]
                        *np.conj(PMNS_dict["{}{}".format(beta, k)])
                        * np.conj(eigenvector[i - 1, J - 1])
                        * eigenvector[k - 1, J - 1]
                        * cmath.exp(complex(0, (eigenvalue * L) / (2 * E * (a ** 2))  ))  
                    )
    
    return total_sum


def Probability(amplitude):
    absolute = abs(amplitude)
    return absolute**2


def Hamiltonian(ξ_order, N):
    H = np.empty((15,15))

    H[0] = [η(N, ξ_order[0] ), 0, 0, ξ_order[0], 0, 0, 2*ξ_order[0], 0, 0, 3*ξ_order[0], 0, 0, 4*ξ_order[0], 0, 0]
    H[1] = [0, η(N, ξ_order[1] ), 0, 0, ξ_order[1], 0, 0,2*ξ_order[1], 0,  0,3*ξ_order[1], 0,  0,4*ξ_order[1], 0]
    H[2] = [0, 0, η(N, ξ_order[2] ), 0, 0,ξ_order[2],  0, 0, 2*ξ_order[2], 0, 0, 3*ξ_order[2], 0, 0, 4*ξ_order[2]]
    H[3] = [ξ_order[0], 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0]
    H[4] = [0, ξ_order[1], 0, 0, 1, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0]
    H[5] = [0, 0, ξ_order[2], 0, 0, 1, 0, 0, 0,  0, 0, 0,  0, 0, 0]
    H[6] = [2*ξ_order[0], 0, 0, 0, 0, 0, 4, 0, 0,  0, 0, 0,  0, 0, 0]
    H[7] = [0, 2*ξ_order[1], 0, 0, 0, 0, 0, 4, 0,  0, 0, 0,  0, 0, 0]
    H[8] = [0, 0, 2*ξ_order[2], 0, 0, 0, 0, 0, 4,  0, 0, 0,  0, 0, 0]
    H[9] = [3*ξ_order[0], 0, 0, 0, 0, 0, 0, 0, 0,  9, 0, 0,  0, 0, 0]
    H[10] = [0, 3*ξ_order[1], 0, 0, 0, 0, 0, 0, 0,  0, 9, 0,  0, 0, 0]
    H[11] = [0, 0, 3*ξ_order[2], 0, 0, 0, 0, 0, 0,  0, 0, 9,  0, 0, 0]
    H[12] = [4*ξ_order[0], 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  16, 0, 0]
    H[13] = [0, 4*ξ_order[1], 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 16, 0]
    H[14] = [0, 0, 4*ξ_order[2], 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 16]

    return H


def find_eigen(mass_dictionary, a):

    ξ = calculate_ξ(mass_dictionary, a)
    H = Hamiltonian(ξ, N)
    eigenvalues, eigenvector_matrix  = np.linalg.eig(H)
    eigenvector_matrix = eigenvector_matrix[0:3, :]
    return eigenvalues, eigenvector_matrix




#----------------------------------------------------------------------------------------------------------------------------------
#Calculating the standard probabilities without extra dimensions

def P_e(E,L):

    m1 =  m0
    m2 = np.sqrt(delta_21)
    m3 = np.sqrt(delta_31)

    Delta_21 = (m2**2-m1**2)*L/(4*E)
    Delta_31 = (m3**2-m1**2)*L/(4*E)
    Delta_32 = (m3**2-m2**2)*L/(4*E)

    U_e1 = abs(PMNS_dict["{}".format("e")+"{}".format(1)])**2
    U_e2 = abs(PMNS_dict["{}".format("e")+"{}".format(2)])**2
    U_e3 = abs(PMNS_dict["{}".format("e")+"{}".format(3)])**2

    P = 1 -4*U_e1*U_e2*(np.sin(Delta_21)**2) -4*U_e1*U_e3*(np.sin(Delta_31)**2) -4*U_e2*U_e3*(np.sin(Delta_32)**2)

    return P


def P_mu(E,L):

    m1 =  m0
    m2 = np.sqrt(delta_21)
    m3 = np.sqrt(delta_31)

    Delta_21 = (m2**2-m1**2)*L/(4*E)
    Delta_31 = (m3**2-m1**2)*L/(4*E)
    Delta_32 = (m3**2-m2**2)*L/(4*E)

    U_mu1 = abs(PMNS_dict["{}".format("mu")+"{}".format(1)])**2
    U_mu2 = abs(PMNS_dict["{}".format("mu")+"{}".format(2)])**2
    U_mu3 = abs(PMNS_dict["{}".format("mu")+"{}".format(3)])**2

    P = 1 -4*U_mu1*U_mu2*(np.sin(Delta_21)**2) -4*U_mu1*U_mu3*(np.sin(Delta_31)**2) -4*U_mu2*U_mu3*(np.sin(Delta_32)**2)

    return P


def Standard_prob(alpha, Energy, length):
    
    Probabilities = np.empty(N_points)

    if alpha == "mu":
        for i in range(N_points):
            Probabilities[i] = P_mu(Energy, length)

    elif alpha == 'e':
        for i in range(N_points):
            Probabilities[i] = P_e(Energy, length)
    
    return Probabilities

#-----------------------------------------------------------------------------------------------------------------

def Plot(a, length, Energy, alpha, beta, label):
    

    probabilities_NH = np.empty(N_points)
    probabilities_IH = np.empty(N_points)
    probabilities_standard = Standard_prob(alpha, Energy ,length)

    #NH
    for n in range(N_points):
        eigenvalues_NH, eigenvectors_NH = find_eigen(NH_dict, a[n])
        probabilities_NH[n] = Probability(Amplitude(alpha, beta, length, Energy, a[n], eigenvalues_NH, eigenvectors_NH))  
        

    #ΙH
    for n in range(N_points):
        eigenvalues_IH, eigenvectors_IH = find_eigen(IH_dict, a[n])
        probabilities_IH[n] = Probability(Amplitude(alpha, beta, length,Energy, a[n], eigenvalues_IH, eigenvectors_IH))

    #x = a*conversion_factor_μm_to_eV

    plt.figure()
    plt.plot(R_LED_, probabilities_NH, label = 'NH', color = 'blue')
    plt.plot(R_LED_, probabilities_IH, label = 'IH', color = 'orange' )
    plt.plot(R_LED_, probabilities_standard, label = 'SM', color= 'g')
    plt.xlabel(r'$R_{LED}(μm)$')  #l_0 = 1μm
    plt.xscale('log')
    plt.ylabel(label)
    plt.legend(title = '{}'.format('Energy at Max Flux= '+'{}'.format(Energy/1e9)+ 'GeV'))
    #plt.title(r'$P({} \rightarrow {})$'.format(neutrino_type1, neutrino_type2))
    plt.show()
    return 

#--------------------------------------------------------------------------------------------------------------------------------------
#E = 2e9




#--------------------------------------------------------------------------------------------------------------------------------------------
file1 = open('/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Numu.txt', 'r')
file2 = open('/Users/athinavogiatzi/Documents/level 4/project/nustorm/resources/fluxes/E6spectraMuSig557Nue.txt', 'r')

data1 = np.loadtxt(file1, skiprows=1)
data2 = np.loadtxt(file2, skiprows=1)

energies1 = data1[:,0]
flux1 = data1[:,1]

energies2 = data2[:,0]
flux2 = data2[:,1]

def find_max_flux_energy(energies, flux):
    
    max_flux_index = flux.argmax()  # Index of the maximum flux value
    
    energy_at_max_flux = energies[max_flux_index]  # Energy corresponding to the maximum flux
    return energy_at_max_flux

max_flux_energy_numu = find_max_flux_energy(energies1, flux1)
max_flux_energy_nue = find_max_flux_energy(energies2, flux2)

print(max_flux_energy_numu)
print(max_flux_energy_nue)

#--------------------------------------------------------------------------------------------------------------------------------------------------
#From GeV to eV
max_flux_energy_numu = max_flux_energy_numu*10**9
max_flux_energy_nue = max_flux_energy_nue*10**9

#Give L: 200, 500, 1000
# Plot(R_LED, L, max_flux_energy_nue, 'e', 'e', r'$P(ν_e \rightarrow ν_e)$' )
# Plot(R_LED, L, max_flux_energy_numu, 'mu', 'mu', r'$P(\bar ν_μ \rightarrow \bar ν_μ)$' )


# Plot(R_LED, L, max_flux_energy_nue, 'e', 'mu', r'$ν_e$', r'$ν_μ$' )
# Plot(R_LED, L, max_flux_energy_numu, 'mu', 'e', r'$ν_μ$', r'$ν_e$' )


#--------------------------------------------------------------------------------------------------------------------------------------------------------

def probability_NH(a,  Energy, alpha, beta, L):
    eigenvalues_NH, eigenvectors_NH = find_eigen(NH_dict, a)
    probability_NH = Probability(Amplitude(alpha, beta, L, Energy, a, eigenvalues_NH, eigenvectors_NH)) 
    return probability_NH 
        

def probability_IH(a,  Energy, alpha, beta, L):
    eigenvalues_IH, eigenvectors_IH = find_eigen(IH_dict, a)
    probability_IH = Probability(Amplitude(alpha, beta, L, Energy, a, eigenvalues_IH, eigenvectors_IH)) 
    return probability_IH 

#--------------------------------------------------------------------------------------------------------------------------------------------------------------       
