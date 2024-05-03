import numpy as np
import cmath
from scipy.optimize import newton_krylov
import matplotlib.pyplot as plt


#--------------------------------------------------------------------------------------------------------------------------
#CONSTANTS

#Normal hierarchy m3>m2>m1=m0
#Inverted hierarchy  m2>m1>m3=m0

m0 = 0 #lightest mass
N = 3

#mass differences in eVs
delta_21 = np.sqrt(7.59*10**(-5))
delta_31 = np.sqrt(2.46*10**(-3))

#circular radius of extra dimensions in μm conveted to eV^-1
a_LED = 0.5/(1.97*10**(-1))
a_0 = 0    #for no extra dimensions

#Length in km to eV^-1
conversion_factor = 1.9732705*10**(-10)
L_1 = 735/conversion_factor
L_2 = 180/conversion_factor
L_3 = 1/conversion_factor

#----------------------------------------------------------------------------------------------------------------------


#Determining masses for NH and IH

NH_dict = {
    "m1": m0,
    "m2": delta_21,
    "m3": delta_31
}

IH_dict = {
    "m3": m0,
    "m1": delta_31,
    "m2": delta_31+delta_21
}


#---------------------------------------------------------------------------------------------------------------------
#PMNS matrix
#Dirac phase taken to be π
delta = 1

#x = 0.9822
x = 0.0178

sin_squared_theta_12 = 0.32
sin_squared_theta_23 = 0.5
sin_squared_2theta_13 = 0.07

cos_squared_theta_12 = 1-sin_squared_theta_12
cos_squared_theta_23 = 1-sin_squared_theta_23
#cos_squared_theta_13 = 1-sin_squared_theta_13

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
    
#-----------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#FUNCTIONS

def ξ_i(m_i, a):
    return np.sqrt(2)*m_i*a


#Calculate ξs depending on the hirarchy
def calculate_ξ(dictionary, a):

    ξ = np.empty((3))

    for i in range(1,3):
        mi = "m{}".format(i)
        ξ[i-1] = ξ_i(dictionary[mi],a)

    #print(ξ)
    return ξ


def η(N, ξ_i):
    return (N+1/2)*(ξ_i**2)








#V_ij is zero in vaccum, because V_a = δ_eαV_CC +V_NC and the electron and neutron number density are zero. 
#T_ij= [-λ^2 + ...]δ_ij +V_ij. Only diagonal elements survive due to kronecker delta.
#When T matrix is diagonal find the eigenvalue for a given ξi numerically. 

def find_eigenvalue_diagonal_matrix(ξ_i, trials):
    
    def F(λ):
        T_ij = -(λ**2) + np.pi*0.5*λ * (ξ_i**2) / (np.tan(np.pi*λ))
        return T_ij**3
    
    lamda1 = newton_krylov(F, trials[0])
    lamda2 = None
    lamda3 = None
    
    for i in range(1, len(trials)):
        lamda = newton_krylov(F, trials[i])
        
        # Check if lamda is different from lamda1 in at least the first 3 decimal places
        if abs(lamda - lamda1) >= 0.1 and lamda2 is None:
            lamda2 = lamda
        # Check if lamda is different from lamda1 and lamda2 in at least the first 3 decimal places
        elif abs(lamda - lamda1) >= 0.1 and abs(lamda - lamda2) >= 0.1 and lamda3 is None:
            lamda3 = lamda
        
        # Exit the loop if all three unique lambdas have been found
        if lamda3 is not None:
            break
    
    return lamda1, lamda2, lamda3

#trials= np.linspace(0.1, 5, 1000)
#print(find_eigenvalue_diagonal_matrix(1, trials))




#--------------------------------------------------------------------------------------------------------
#Do not use!!!
#Calculating eigenvector elements
#This is an approximation. 
def W_i(ξ_i, N):

    if N==0:
        w = np.sqrt(1-(np.pi**2)*(ξ_i**2)/6)

    else:
        w = ξ_i/N

    return w
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––--    


W_eigenvectors = np.empty((3,3))












def W_populate(matrix, ξs, Ns):

    for ξ in range(len(ξs)):
        for n in range(len(Ns)):
            matrix[ξ,n] = W_i(ξs[ξ], Ns[n])

    return matrix

#W_eigenvectors = W_populate(W_eigenvectors, [0 ,    0.03080179, 1   ], [0,1,2])
#W_eigenvectors[2,0] =0 
#print(W_eigenvectors)


#Solving for W_i^(0N) for a specific eigenvalue
def W_i_trial(λ, ξ1, ξ2, ξ3):
    
    def F_test(λ, ξ):
        T_ij = -(λ**2) + np.pi*0.5*λ * (ξ**2) / (np.tan(np.pi*λ))
        return T_ij
    sum = F_test(λ, ξ1)*F_test(λ, ξ2)*F_test(λ, ξ3)

    return sum


#print(W_i_trial(-0.11574075642450818 ,0 ,    0.03080179, 1))








#Amplitude as a function of energy. 
#α,β = e,μ,τ
#

def Amplitude(alpha,beta,L,E,a, eigenvalue_array, eigenvector_matrix):
    total_sum = 0


    for i in range(1,3):
        for J in range(1,3):
            for k in range(1,3):
                for N in  range(0,2):  #Will make it 5 later
                    #print(total_sum)

                    if N==0:
                        eigenvector = eigenvector_matrix[:, 0:3]
                        total_sum+= PMNS_dict["{}".format(alpha)+"{}".format(i)]*np.conj(PMNS_dict["{}".format(beta)+"{}".format(k)])*np.conj(eigenvector[i-1, J-1])*eigenvector[k-1, J-1]*cmath.exp(complex(0,(eigenvalue_array[J-1, N])*L/(2*E*(a**2))))
                        #print(eigenvector[i-1, J-1])
                    elif N==1:
                        eigenvector = eigenvector_matrix[:, 3:6]
                        total_sum+= PMNS_dict["{}".format(alpha)+"{}".format(i)]*np.conj(PMNS_dict["{}".format(beta)+"{}".format(k)])*np.conj(eigenvector[i-1, J-1])*eigenvector[k-1, J-1]*cmath.exp(complex(0,(eigenvalue_array[J-1, N])*L/(2*E*(a**2))))
                        #print(eigenvector[i-1, J-1])
                    elif N==2:
                        eigenvector = eigenvector[:, 6:9]
                        total_sum+= PMNS_dict["{}".format(alpha)+"{}".format(i)]*np.conj(PMNS_dict["{}".format(beta)+"{}".format(k)])*np.conj(eigenvector[i-1, J-1])*eigenvector[k-1, J-1]*cmath.exp(complex(0,(eigenvalue_array[J-1, N])*L/(2*E*(a**2))))
                        #print(eigenvector[i-1, J-1])    
    return total_sum




# Amplitude("e","e",10, 3, NH_dict, a_LED) #works
def Probability(amplitude):
    absolute = abs(amplitude)
    return absolute**2

#print(Probability(x)) #works



def Plot(Energy_start, Energy_stop, number_of_points, length, alpha, beta, eigenvalues_NH, eigenvectors_NH, eigenvalues_IH, eigenvectors_IH, legend):
    
    E = np.linspace(Energy_start,Energy_stop, number_of_points )
    probabilities_NH = np.empty(number_of_points)
    probabilities_IH = np.empty(number_of_points)

    #NH
    for n in range(number_of_points):
        probabilities_NH[n] = Probability(Amplitude(alpha, beta, length, E[n], a_LED, eigenvalues_NH, eigenvectors_NH))  #a, eigenvalue_array, eigenvector_matrix)
        

    #ΙH
    for n in range(number_of_points):
        probabilities_IH[n] = Probability(Amplitude(alpha, beta, length,E[n], a_LED, eigenvalues_IH, eigenvectors_IH))

    

    plt.figure()
    plt.scatter(E, probabilities_NH, label = '{}'.format('NH'+ '{}'.format(legend)))
    plt.scatter(E, probabilities_IH, label = '{}'.format('IH'+ '{}'.format(legend)))
    plt.legend()
    plt.show()
    return 



#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#ξ_NH = calculate_ξ(NH_dict, a_LED)
#print(ξ_NH)
ξ_nh =  [0 ,    0.03080179, 1        ]

NH_eigenvalues_ξ1 = find_eigenvalue_diagonal_matrix(ξ_nh[0], np.linspace(-0.2, 3, 10000))
NH_eigenvalues_ξ2=find_eigenvalue_diagonal_matrix(ξ_nh[1], np.linspace(-0.02, -1.9, 70))
NH_eigenvalues_ξ3 = find_eigenvalue_diagonal_matrix(ξ_nh[2], np.linspace(0.1, 5, 1000) )

#print(NH_eigenvalues_ξ2)
#print(NH_eigenvalues_ξ3)
#ξ_IH = calculate_ξ(IH_dict, a_LED)
#print(ξ_IH)

ξ_ih = [0.17535678, 0.20615856, 0    ]
IH_eigenvalues_ξ1 = find_eigenvalue_diagonal_matrix(ξ_ih[0], np.linspace(0.001, 3, 10000))
IH_eigenvalues_ξ2 = find_eigenvalue_diagonal_matrix(ξ_ih[1], np.linspace(-0.02, -3, 7000))    
IH_eigenvalues_ξ3 = find_eigenvalue_diagonal_matrix(ξ_ih[2], np.linspace(-0.3, 0.1, 1000) )





NH_eigenvalue_array = np.empty((3,3))   #Row is n eigenvalues for a particular ξ   #Row is J index, column is n index

NH_eigenvalue_array[0,0] = NH_eigenvalues_ξ1[0]
NH_eigenvalue_array[0,1] = NH_eigenvalues_ξ1[1]
NH_eigenvalue_array[0,2] = NH_eigenvalues_ξ1[2]


NH_eigenvalue_array[1,0] = NH_eigenvalues_ξ2[0]
NH_eigenvalue_array[1,1] = NH_eigenvalues_ξ2[1]
NH_eigenvalue_array[1,2] = NH_eigenvalues_ξ2[2]


NH_eigenvalue_array[2,0] = NH_eigenvalues_ξ3[0]
NH_eigenvalue_array[2,1] = NH_eigenvalues_ξ3[1]
NH_eigenvalue_array[2,2] = NH_eigenvalues_ξ3[2]


#print(NH_eigenvalue_array)


IH_eigenvalue_array = np.empty((3,3))   #Row is n eigenvalues for a particular ξ   #Row is J index, column is n index

IH_eigenvalue_array[0,0] = IH_eigenvalues_ξ1[0]
IH_eigenvalue_array[0,1] = IH_eigenvalues_ξ1[1]
IH_eigenvalue_array[0,2] = IH_eigenvalues_ξ1[2]


IH_eigenvalue_array[1,0] = IH_eigenvalues_ξ2[0]
IH_eigenvalue_array[1,1] = IH_eigenvalues_ξ2[1]
IH_eigenvalue_array[1,2] = IH_eigenvalues_ξ2[2]


IH_eigenvalue_array[2,0] = IH_eigenvalues_ξ3[0]
IH_eigenvalue_array[2,1] = IH_eigenvalues_ξ3[1]
IH_eigenvalue_array[2,2] = IH_eigenvalues_ξ3[2]

#print(IH_eigenvalue_array)


def Hamiltonian(ξ_order, N):
    H = np.empty((9,9))

    H[0] = [η(N, ξ_order[0] ), 0, 0, ξ_order[0], 0, 0, 2*ξ_order[0], 0, 0]
    H[1] = [0, η(N, ξ_order[1] ), 0, 0, ξ_order[1], 0, 0,2*ξ_order[1], 0]
    H[2] = [0, 0, η(N, ξ_order[2] ), 0, 0,ξ_order[2],  0, 0, 2*ξ_order[2]]
    H[3] = [ξ_order[0], 0, 0, 1, 0, 0, 0, 0, 0]
    H[4] = [0, ξ_order[1], 0, 0, 1, 0, 0, 0, 0]
    H[5] = [0, 0, ξ_order[2], 0, 0, 1, 0, 0, 0]
    H[6] = [2*ξ_order[0], 0, 0, 0, 0, 0, 4, 0, 0]
    H[7] = [0, 2*ξ_order[1], 0, 0, 0, 0, 0, 4, 0]
    H[8] = [0, 0, 2*ξ_order[2], 0, 0, 0, 0, 0, 4]

    return H

H_NH = Hamiltonian(ξ_nh, N)
NH_eigenvalues_matrix, NH_eigenvector_matrix  = np.linalg.eig(H_NH)
#NH_eigenvalues_matrix = np.square(NH_eigenvalues_matrix)




H_IH = Hamiltonian(ξ_ih, N)
IH_eigenvalues_matrix,IH_eigenvector_matrix = np.linalg.eig(H_IH)
#IH_eigenvalues_matrix = np.square(IH_eigenvalues_matrix)

#Taking the M=0 elements and splitting into N=1,2,3

NH_eigenvector_matrix = NH_eigenvector_matrix[0:3, :]
IH_eigenvector_matrix = IH_eigenvector_matrix[0:3, :]

#print(NH_eigenvector_matrix)
#print(IH_eigenvector_matrix)



def normalize(W, ξ, eigenvalue_array):
    def F(ξ_l, λ_i):
        return 1 +(ξ_l**2)*((np.pi**2)/(4*np.tan(np.pi*λ_i)**2) -np.pi/(4*λ_i*np.tan(np.pi*λ_i)) +np.pi**2/4)

    W[:, 0] = W[:, 0]/np.sqrt(F(ξ[0], eigenvalue_array[0]))
    W[:, 1] = W[:, 1]/np.sqrt(F(ξ[1], eigenvalue_array[1]))
    W[:, 2] = W[:, 2]/np.sqrt(F(ξ[2], eigenvalue_array[2]))

    return W

def normalize2(W, ξ, eigenvalue_array):
    def F(ξ_l, λ_i):
        return 1 +(ξ_l**2)*((np.pi**2)/(4*np.tan(np.pi*λ_i)**2) -np.pi/(4*λ_i*np.tan(np.pi*λ_i)) +np.pi**2/4)

    n = F(ξ[0], eigenvalue_array[0]) +F(ξ[1], eigenvalue_array[1]) + F(ξ[2], eigenvalue_array[2])

    W[:, 0] = W[:, 0]/np.sqrt(n)
    W[:, 1] = W[:, 1]/np.sqrt(n)
    W[:, 2] = W[:, 2]/np.sqrt(n)

    return W



NH_eigenvector_matrix[:, 0:3] = normalize(NH_eigenvector_matrix[:, 0:3], ξ_nh, NH_eigenvalue_array[:, 0])
NH_eigenvector_matrix[:, 3:6] = normalize(NH_eigenvector_matrix[:, 3:6], ξ_nh, NH_eigenvalue_array[:, 1])
NH_eigenvector_matrix[:, 6:9] = normalize(NH_eigenvector_matrix[:, 6:9], ξ_nh, NH_eigenvalue_array[:, 2])

#print(NH_eigenvector_matrix)

IH_eigenvector_matrix[:, 0:3] = normalize(IH_eigenvector_matrix[:, 0:3], ξ_ih, IH_eigenvalue_array[:, 0])
IH_eigenvector_matrix[:, 3:6] = normalize(IH_eigenvector_matrix[:, 3:6], ξ_ih, IH_eigenvalue_array[:, 1])
IH_eigenvector_matrix[:, 6:9] = normalize(IH_eigenvector_matrix[:, 6:9], ξ_ih, IH_eigenvalue_array[:, 2])





#--------------------------------------------------------------------------------------------------------------------------------------------
#1km Plot #E_stop = 9*10**(6) , 500
Plot(1.5*10**(6),9*10**(6) , 500, L_3, 'e', 'e', NH_eigenvalue_array, NH_eigenvector_matrix,IH_eigenvalue_array, IH_eigenvector_matrix, 'L3' )
Plot(1.5*10**(6),9*10**(6) , 500, L_2, 'e', 'e', NH_eigenvalue_array, NH_eigenvector_matrix,IH_eigenvalue_array, IH_eigenvector_matrix, 'L2' )
Plot(1*10**(9),5*10**(9) , 500, L_1, 'mu', 'mu', NH_eigenvalue_array, NH_eigenvector_matrix,IH_eigenvalue_array, IH_eigenvector_matrix, 'L1' )


#print( NH_eigenvector_matrix[:, 3:6])





#-----------------------------------------------------------------------------------------------------------------------

def P_e(E,L):

    m1 =  m0
    m2 = delta_21
    m3 = delta_31

    Delta_21 = (m2**2-m1**2)*L/(4*E)
    Delta_31 = (m3**2-m1**2)*L/(4*E)
    Delta_32 = (m3**2-m2**2)*L/(4*E)

    U_e1 = abs(PMNS_dict["{}".format("e")+"{}".format(1)])**2
    U_e2 = abs(PMNS_dict["{}".format("e")+"{}".format(2)])**2
    U_e3 = abs(PMNS_dict["{}".format("e")+"{}".format(3)])**2

    P = 1 -4*U_e1*U_e2*(np.sin(Delta_21)**2) -4*U_e1*U_e3*(np.sin(Delta_31)**2) -4*U_e2*U_e3*(np.sin(Delta_32)**2)

    return P


def Plot_probability_e(Energy_stop, number_of_points, length, legend):
    
    E = np.linspace(1,Energy_stop, number_of_points )
    Probabilities = np.empty(number_of_points)

    for i in range(number_of_points):
        Probabilities[i] = P_e(E[i], length)

    plt.figure()
    plt.scatter(E, Probabilities, label = '{}'.format('standard'+ '{}'.format(legend)))
    plt.legend()
    plt.show()
    return 



#Plot_probability_e(9*10**(6) , 500, L_2, "L2")
#Plot_probability_e(9*10**(6) , 500, L_3, "L3")



def P_mu(E,L):

    m1 =  m0
    m2 = delta_21
    m3 = delta_31

    Delta_21 = (m2**2-m1**2)*L/(4*E)
    Delta_31 = (m3**2-m1**2)*L/(4*E)
    Delta_32 = (m3**2-m2**2)*L/(4*E)

    U_mu1 = abs(PMNS_dict["{}".format("mu")+"{}".format(1)])**2
    U_mu2 = abs(PMNS_dict["{}".format("mu")+"{}".format(2)])**2
    U_mu3 = abs(PMNS_dict["{}".format("mu")+"{}".format(3)])**2

    P = 1 -4*U_mu1*U_mu2*(np.sin(Delta_21)**2) -4*U_mu1*U_mu3*(np.sin(Delta_31)**2) -4*U_mu2*U_mu3*(np.sin(Delta_32)**2)

    return P


def Plot_probability_mu(Energy_stop, number_of_points, length, legend):
    
    E = np.linspace(1,Energy_stop, number_of_points )
    Probabilities = np.empty(number_of_points)

    for i in range(number_of_points):
        Probabilities[i] = P_mu(E[i], length)

    plt.figure()
    plt.scatter(E, Probabilities, label = '{}'.format('standard'+ '{}'.format(legend)))
    plt.legend()
    plt.show()
    return 



#Plot_probability_mu(9*10**(9) , 500, L_1, 'L1')
