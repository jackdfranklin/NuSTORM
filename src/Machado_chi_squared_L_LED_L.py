import numpy as np
import cmath
import matplotlib.pyplot as plt


#--------------------------------------------------------------------------------------------------------------------------
#CONSTANTS

#Normal hierarchy m3>m2>m1=m0
#Inverted hierarchy  m2>m1>m3=m0

m0 = 0 #lightest mass5*10**(-2)#
N = 4
N_points = 30

#mass differences in eVs
delta_21 = 7.53*10**(-5)
delta_31 = 2.51*10**(-3)

#circular radius of extra dimensions in μm conveted to eV^-1
#a_LED = 0.5/(1.9732705*10**(-1))
conversion_factor_μm_to_eV = 1.9732705*10**(-1)  
R_LED_ = np.logspace(np.log10(1e-3), np.log10(5), N_points)
R_LED = R_LED_/conversion_factor_μm_to_eV


#L is given in meters and converted into eV-1 in the probability functions
#Length in m to eV^-1
conversion_factor_m_to_eV = 1.9732705*10**(-7)


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


#--------------------------------------------------------------------------------------------------------------------------------------------------------
#The units of L are changed from meters to eV-1 here
def probability_NH(a,  Energy, alpha, beta, L):
    L = L/conversion_factor_m_to_eV
    eigenvalues_NH, eigenvectors_NH = find_eigen(NH_dict, a)
    probability_NH = Probability(Amplitude(alpha, beta, L, Energy, a, eigenvalues_NH, eigenvectors_NH)) 
    return probability_NH 
        

def probability_IH(a,  Energy, alpha, beta, L):
    L = L/conversion_factor_m_to_eV
    eigenvalues_IH, eigenvectors_IH = find_eigen(IH_dict, a)
    probability_IH = Probability(Amplitude(alpha, beta, L, Energy, a, eigenvalues_IH, eigenvectors_IH)) 
    return probability_IH 

#--------------------------------------------------------------------------------------------------------------------------------------------------------------       
