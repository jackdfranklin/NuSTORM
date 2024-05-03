import numpy as np
import cmath
from scipy.optimize import newton_krylov
import matplotlib.pyplot as plt


#--------------------------------------------------------------------------------------------------------------------------
#CONSTANTS

#Normal hierarchy m3>m2>m1=m0
#Inverted hierarchy  m2>m1>m3=m0

m0 = 5e-2 #lightest mass5*10**(-2)#
N = 4

#mass differences in eVs
delta_21 = 7.53*10**(-5)#7.59*10**(-5)
delta_31 = 2.51*10**(-3)#2.46*10**(-3)

#circular radius of extra dimensions in μm conveted to eV^-1
a_LED = 50/(1.97*10**(-1))  #0.38#


#Length in km to eV^-1
conversion_factor = 1.9732705*10**(-10)
L_1 = 735/conversion_factor
L_2 = 180/conversion_factor
L_3 = 1/conversion_factor

L_4 = 13/conversion_factor
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
    "m2": np.sqrt(delta_31+delta_21)   #np.sqrt(delta_31) + np.sqrt(+delta_21)
}


#---------------------------------------------------------------------------------------------------------------------
#PMNS matrix
#Dirac phase taken to be π
delta = 0

#x = 0.9822
x = 0.0178
#x=2.18*10**(-2)


sin_squared_theta_12 = 0.32  #0.304#
sin_squared_theta_23 = 0.5  #0.514
#sin_squared_2theta_13 = #0.07

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





        


def find_eigenvalues(ξ_i, trials):
    
    def F(λ):
        T_ij = λ - np.pi*0.5*(ξ_i**2) / (np.tan(np.pi*λ))
        return T_ij
    
    lamda1 = None
    lamda2 = None
    lamda3 = None
    lamda4 = None
    lamda5 = None
    
    for i in range(len(trials)):
        try:
            lamda = newton_krylov(F, trials[i])
            
            # Check if lamda is different from lamda1 in at least the first 3 decimal places
            if lamda1 is None:
                lamda1 = lamda
            # Check if lamda is different from lamda1 and lamda2 in at least the first 3 decimal places
            elif abs(lamda - lamda1) >= 0.1 and lamda2 is None:
                lamda2 = lamda
            # Check if lamda is different from lamda1, lamda2, and lamda3 in at least the first 3 decimal places
            elif abs(lamda - lamda1) >= 0.1 and abs(lamda - lamda2) >= 0.1 and lamda3 is None:
                lamda3 = lamda
            # Check if lamda is different from lamda1, lamda2, lamda3, and lamda4 in at least the first 3 decimal places
            elif abs(lamda - lamda1) >= 0.1 and abs(lamda - lamda2) >= 0.1 and abs(lamda - lamda3) >= 0.1 and lamda4 is None:
                lamda4 = lamda
            # Check if lamda is different from lamda1, lamda2, lamda3, lamda4, and lamda5 in at least the first 3 decimal places
            elif abs(lamda - lamda1) >= 0.1 and abs(lamda - lamda2) >= 0.1 and abs(lamda - lamda3) >= 0.1 and abs(lamda - lamda4) >= 0.1 and lamda5 is None:
                lamda5 = lamda
            
            # Break out of the loop if all five unique lambdas have been found
            if lamda5 is not None:
                break
        
        except Exception as e:
            # Handle exception (e.g., newton_krylov failed for the current trial)
            #print(f"Failed for trial {i+1}: {e}")
            continue
    
    return lamda1, lamda2, lamda3, lamda4, lamda5


#Amplitude as a function of energy. 
#α,β = e,μ,τ
   
def S_elements(ξ, λ):

    x = 1 +(np.pi**2)*0.5*ξ**2 +2*(λ/ξ)**2
    return 2/x

def S_array(ξs, λ):
    S = np.empty((3,5))
    for i in range(1, 3):    
        for N in range(4): 
            S[i-1,N] = S_elements(ξs[i-1], λ[i-1,N])
    return S



def Amplitude(alpha, beta, L, E, a,  eigenvalues, S_squared):
    total_sum = 0
    N_max = 4  # Adjust N_max as needed
    
    for i in range(1, 3):    
        for N in range(N_max):      
            total_sum += (
                PMNS_dict["{}{}".format(alpha, i)]
                *np.conj(PMNS_dict["{}{}".format(beta, i)])
                * S_squared[i-1, N]
                * cmath.exp(complex(0, ((eigenvalues[i-1, N]**2) *L ) / (2 * E * (a ** 2))  ))  
                )
    
    return total_sum





def Probability(amplitude):
    absolute = abs(amplitude)
    return absolute**2





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
    plt.plot(E, probabilities_NH, label = '{}'.format('NH'+ '{}'.format(legend)))
    plt.plot(E, probabilities_IH, label = '{}'.format('IH'+ '{}'.format(legend)))
    plt.legend()
    # plt.ylim(0, 1)
    plt.show()
    return 



#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
ξ_NH = calculate_ξ(NH_dict, a_LED)
ξ_IH = calculate_ξ(IH_dict, a_LED)


# ξ_NH = [0.   ,    0.031147, 1.      ]
# ξ_IH = [0.17982731, 0.18250478, 1.        ]

print(ξ_NH)
print(ξ_IH)

#NH_eigenvalues_ξ1 = find_eigenvalues(ξ_NH[0], np.linspace(0.0001,10, 10000))
NH_eigenvalues_ξ1 = np.zeros(5) #We know this is just zero
NH_eigenvalues_ξ2 = find_eigenvalues(ξ_NH[1], np.linspace(-0.022006697, 10, 10000))
NH_eigenvalues_ξ3 = find_eigenvalues(ξ_NH[2], np.linspace(0.417339091, 5, 1000) )


IH_eigenvalues_ξ1 = find_eigenvalues(ξ_IH[0], np.linspace(0.123872303, 5, 1000))
IH_eigenvalues_ξ2 = find_eigenvalues(ξ_IH[1], np.linspace(0.1256196299, 5, 1000))    
IH_eigenvalues_ξ3 = find_eigenvalues(ξ_IH[2], np.linspace(0.417, 5, 1000) )


print(NH_eigenvalues_ξ1)
print(NH_eigenvalues_ξ2)
print(NH_eigenvalues_ξ3)

print(IH_eigenvalues_ξ1)
print(IH_eigenvalues_ξ2)
print(IH_eigenvalues_ξ3)


NH_eigenvalues = np.empty((3,5))
IH_eigenvalues = np.empty((3,5))


NH_eigenvalues[0] = NH_eigenvalues_ξ1
NH_eigenvalues[1] = np.array([arr.item() for arr in NH_eigenvalues_ξ2])
NH_eigenvalues[2] = np.array([arr.item() for arr in NH_eigenvalues_ξ3])

IH_eigenvalues[0] = np.array([arr.item() for arr in IH_eigenvalues_ξ1])  
IH_eigenvalues[1] = np.array([arr.item() for arr in IH_eigenvalues_ξ2])
IH_eigenvalues[2] = np.array([arr.item() for arr in IH_eigenvalues_ξ3])



S_array_NH = S_array(ξ_NH, NH_eigenvalues)
S_array_IH = S_array(ξ_IH, IH_eigenvalues)







#Plot(1.5*10**(6),9*10**(6) , 5000, L_3, 'e', 'e', NH_eigenvalues, S_array_NH, IH_eigenvalues, S_array_IH, 'L3' )
#Plot(1.5*10**(6),9*10**(6) , 5000, L_2, 'e', 'e',NH_eigenvalues,S_array_NH, IH_eigenvalues, S_array_IH,'L2' )
Plot(1*10**(9),5*10**(9) , 5000, L_4, 'mu', 'mu', NH_eigenvalues, S_array_NH, IH_eigenvalues, S_array_IH, 'L4' )












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



Plot_probability_mu(9*10**(9) , 500, L_4, 'L4')





#-------------------------------------------------------------------------------------------------------
# ms = [m0, np.sqrt(delta_21),np.sqrt(delta_31)]

# for i in range(len(ms)):
     
#     print(np.sqrt( 2/(1+np.pi**2*(ms[i]*a_LED)**2 +NH_eigenvalues_matrix[i]**2/((ms[i]*a_LED)**2))  )  )


# print(NH_eigenvector_matrix[:, 0:3])


# for i in range(len(ms)):
     
#     print(np.sqrt( 2/(1+np.pi**2*(ms[i]*a_LED)**2 +NH_eigenvalues_matrix[3+i]**2/((ms[i]*a_LED)**2))  )  )


# print(NH_eigenvector_matrix[:, 3:6])



# for i in range(len(ms)):
     
#     print(np.sqrt( 2/(1+np.pi**2*(ms[i]*a_LED)**2 +NH_eigenvalues_matrix[6+i]**2/((ms[i]*a_LED)**2))  )  )


# print(NH_eigenvector_matrix[:, 6:9])






# print(NH_dict["m1"]*a_LED)
# print(NH_eigenvalues_matrix[0])

# print(NH_dict["m2"]*a_LED)
# print(NH_eigenvalues_matrix[1])

# print(NH_dict["m3"]*a_LED)
# print(NH_eigenvalues_matrix[2])


# print(delta_21)
# print((NH_eigenvalues_matrix[1]**2 - NH_eigenvalues_matrix[0]**2)/(a_LED**2))
# print((IH_eigenvalues_matrix[1]**2 - IH_eigenvalues_matrix[0]**2)/(a_LED**2))

# print(delta_31)
# print((NH_eigenvalues_matrix[2]**2 - NH_eigenvalues_matrix[0]**2)/(a_LED**2))
# print((IH_eigenvalues_matrix[2]**2 - IH_eigenvalues_matrix[0]**2)/(a_LED**2))

# print('end')

# def S_in(mi, lamda):
#     x = 1 + (np.pi**2)*((mi*a_LED)**2) + lamda**2/((mi*a_LED)**2)
#     return 2/x

# print(S_in(NH_dict["m1"], NH_eigenvalues_matrix[0]))
# print(S_in(NH_dict["m2"], NH_eigenvalues_matrix[1]))
# print(S_in(NH_dict["m3"], NH_eigenvalues_matrix[2]))


# print(S_in(NH_dict["m1"], NH_eigenvalues_matrix[3]))
# print(S_in(NH_dict["m2"], NH_eigenvalues_matrix[4]))
# print(S_in(NH_dict["m3"], NH_eigenvalues_matrix[5]))


# print(S_in(NH_dict["m1"], NH_eigenvalues_matrix[6]))
# print(S_in(NH_dict["m2"], NH_eigenvalues_matrix[7]))
# print(S_in(NH_dict["m3"], NH_eigenvalues_matrix[8]))