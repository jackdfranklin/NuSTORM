import numpy as np
import cmath
from scipy.optimize import newton_krylov
import matplotlib.pyplot as plt


#--------------------------------------------------------------------------------------------------------------------------
#CONSTANTS

#Normal hierarchy m3>m2>m1=m0
#Inverted hierarchy  m2>m1>m3=m0

m0 = 0 #lightest mass5*10**(-2)#
N = 4

#mass differences in eVs
delta_21 = 7.53*10**(-5)#7.59*10**(-5)
delta_31 = 2.51*10**(-3)#2.46*10**(-3)

#circular radius of extra dimensions in μm conveted to eV^-1
a_LED = 0.5/(1.97*10**(-1))  #0.38#


#Length in km to eV^-1
conversion_factor = 1.9732705*10**(-10)
L_1 = 735/conversion_factor
L_2 = 180/conversion_factor
L_3 = 1/conversion_factor

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
delta = 1

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





#Amplitude as a function of energy. 
#α,β = e,μ,τ
   


def Amplitude(alpha, beta, L, E, a, eigenvalue_array, eigenvector_matrix):
    total_sum = 0
    N_max = 4  # Adjust N_max as needed
    
    for i in range(1, 3):
        for J in range(1, 3):
            for k in range(1, 3):
                for N in range(N_max):
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

#print(ξ_NH)
#print(ξ_IH)









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





H_NH = Hamiltonian(ξ_NH, N)
NH_eigenvalues_matrix, NH_eigenvector_matrix  = np.linalg.eig(H_NH)


H_IH = Hamiltonian(ξ_IH, N)
IH_eigenvalues_matrix,IH_eigenvector_matrix = np.linalg.eig(H_IH)



# print(NH_eigenvector_matrix)
# print(IH_eigenvector_matrix)







def check(Hamiltonian, eigenvector_matrix, eigenvalues):
    for i in range(0, 15):
        lhs = np.dot(Hamiltonian, eigenvector_matrix[:, i])
        print(eigenvalues[i])
# Calculate the right-hand side of the equation
        rhs = eigenvalues[i] * eigenvector_matrix[:, i]

# Check if both sides are close within a specified tolerance
        are_equal = np.allclose(lhs, rhs)

# Print the result
        print(are_equal)

    return

#print(check(H_NH, NH_eigenvector_matrix, NH_eigenvalues_matrix))
#print(check(H_IH, IH_eigenvector_matrix, IH_eigenvalues_matrix))

print(NH_eigenvector_matrix)





#Taking the M=0 elements and splitting into N=1,2,3

NH_eigenvector_matrix = NH_eigenvector_matrix[0:3, :] 
IH_eigenvector_matrix = IH_eigenvector_matrix[0:3, :]



    



#print(NH_eigenvector_matrix)
# print(IH_eigenvector_matrix)



Plot(1.5*10**(6),9*10**(6) , 1000, L_3, 'e', 'e', NH_eigenvalues_matrix, NH_eigenvector_matrix,IH_eigenvalues_matrix, IH_eigenvector_matrix, 'L3' )
Plot(1.5*10**(6),9*10**(6) , 1000, L_2, 'e', 'e', NH_eigenvalues_matrix, NH_eigenvector_matrix,IH_eigenvalues_matrix, IH_eigenvector_matrix, 'L2' )
Plot(1*10**(9),5*10**(9) , 1000, L_1, 'mu', 'mu', NH_eigenvalues_matrix, NH_eigenvector_matrix,IH_eigenvalues_matrix, IH_eigenvector_matrix, 'L1' )












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



print(NH_dict)
print(IH_dict)


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