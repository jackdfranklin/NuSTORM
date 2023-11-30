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
delta_21 = 7.59*10**(-5)#7.59*10**(-5)
delta_31 = 2.46*10**(-3)#2.46*10**(-3)

#circular radius of extra dimensions in μm conveted to eV^-1
a_LED = 0.5/(1.97*10**(-1))  #0.38#-1


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
    "m2": np.sqrt(delta_31+delta_21)   
}


#---------------------------------------------------------------------------------------------------------------------
#PMNS matrix
#Dirac phase taken to be π
delta = -1

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

    for i in range(1,4):
        mi = "m{}".format(i)
        ξ[i-1] = ξ_i(dictionary[mi],a)
        
    return ξ


def η(N, ξ_i):
    return (N+1/2)*(ξ_i**2)

#Has been multiplied by a^2
def M_dagger_M(ξ, N):
    L = np.empty((5, 5))
    L[0] = [(N+1/2)*(ξ**2), ξ, 2*ξ, 3*ξ, 4*ξ]
    L[1] = [ξ, 1, 0, 0, 0]
    L[2] = [2*ξ, 0, 4, 0, 0]
    L[3] = [3*ξ, 0, 0, 9, 0]
    L[4] = [4*ξ, 0, 0, 0, 16]

    return L/(a_LED**2)  


#Amplitude as a function of energy. 
#α,β = e,μ,τ
   


def Amplitude(alpha, beta, L, E,a, eigenvalues, eigenvectors):
    total_sum = 0
    N_max = 5  # Adjust N_max 
    
    for i in range(1, 4):
       
        for N in range(N_max):      
            total_sum += (
                PMNS_dict["{}{}".format(alpha, i)]
                *np.conj(PMNS_dict["{}{}".format(beta, i)])
                #*np.conj(eigenvectors[i-1,N])
                * eigenvectors[i-1, N]**2
                * cmath.exp(complex(0, (eigenvalues[i-1, N]*L ) / (2 * E)  )  )
                )
    
    return total_sum





def Probability(amplitude):
    absolute = abs(amplitude)
    return absolute**2



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


def Standard_prob(alpha, Energy_start, Energy_stop, number_of_points, length):
    E = np.linspace(Energy_start,Energy_stop, number_of_points )
    Probabilities = np.empty(number_of_points)

    if alpha == "mu":
        for i in range(number_of_points):
            Probabilities[i] = P_mu(E[i], length)

    elif alpha == 'e':
        for i in range(number_of_points):
            Probabilities[i] = P_e(E[i], length)

    return Probabilities




# print(P_e(1.6e6, L_3))
# print(P_e(8e6, L_3))

# print(P_e(1e5, L_3))
# print(P_e(1e5, L_3))

#exit()
#print(Probability(Amplitude('e', 'e', L_3, 1e6, a_LED, eigenvalues_NH, eigenvectors_NH))  #a, eigenvalue_array, eigenvector_matrix)


















def Plot(Energy_start, Energy_stop, number_of_points, length, alpha, beta, eigenvalues_NH, eigenvectors_NH, eigenvalues_IH, eigenvectors_IH, legend, x_label):
    
    #E = np.logspace(np.log10(Energy_start),np.log10(Energy_stop), number_of_points )
    E = np.linspace(Energy_start, Energy_stop, number_of_points)
    probabilities_NH = np.empty(number_of_points)
    probabilities_IH = np.empty(number_of_points)
    #probabilities_standard = Standard_prob(alpha, Energy_start, Energy_stop, number_of_points, length)
    
    #NH
    for n in range(number_of_points):
        probabilities_NH[n] = Probability(Amplitude(alpha, beta, length, E[n], a_LED, eigenvalues_NH, eigenvectors_NH))  #a, eigenvalue_array, eigenvector_matrix)
        

    #ΙH
    for n in range(number_of_points):
        probabilities_IH[n] = Probability(Amplitude(alpha, beta, length,E[n], a_LED, eigenvalues_IH, eigenvectors_IH))

    


    plt.figure()
    plt.plot(E, probabilities_NH, label = '{}'.format('NH'+ '{}'.format(legend)))
    plt.plot(E, probabilities_IH, label = '{}'.format('IH'+ '{}'.format(legend)))
    #plt.plot(E, probabilities_standard, label = '{}'.format('standard'+ '{}'.format(legend)))
    plt.ylabel("Probability")
    plt.xlabel(x_label)
    #plt.xscale('log')
    plt.legend()
    #plt.ylim(0.7, 1)
    plt.show()
    return 



def check(Hamiltonian, eigenvector_matrix, eigenvalues):
    for i in range(0, 5):
        lhs = np.dot(Hamiltonian, eigenvector_matrix[:, i])
        print(eigenvalues[i])
# Calculate the right-hand side of the equation
        rhs = eigenvalues[i] * eigenvector_matrix[:, i]

# Check if both sides are close within a specified tolerance
        are_equal = np.allclose(lhs, rhs)

# Print the result
        print(are_equal)

    return



#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
ξ_NH = calculate_ξ(NH_dict, a_LED)
ξ_IH = calculate_ξ(IH_dict, a_LED)


print(ξ_NH)
print(ξ_IH)


M1_M1_dagger_NH = M_dagger_M(ξ_NH[0], N)
M2_M2_dagger_NH = M_dagger_M(ξ_NH[1], N)
M3_M3_dagger_NH = M_dagger_M(ξ_NH[2], N)

M1_M1_dagger_IH = M_dagger_M(ξ_IH[0], N)
M2_M2_dagger_IH = M_dagger_M(ξ_IH[1], N)
M3_M3_dagger_IH = M_dagger_M(ξ_IH[2], N)

eigenval_M1_NH, eigenvector_M1_NH = np.linalg.eig(M1_M1_dagger_NH)
eigenval_M2_NH, eigenvector_M2_NH = np.linalg.eig(M2_M2_dagger_NH)
eigenval_M3_NH, eigenvector_M3_NH = np.linalg.eig(M3_M3_dagger_NH)

eigenval_M1_IH, eigenvector_M1_IH = np.linalg.eig(M1_M1_dagger_IH)
eigenval_M2_IH, eigenvector_M2_IH = np.linalg.eig(M2_M2_dagger_IH)
eigenval_M3_IH, eigenvector_M3_IH = np.linalg.eig(M3_M3_dagger_IH)

# print(np.matmul(np.matmul(eigenvector_M1_NH, M1_M1_dagger_NH), eigenvector_M1_NH.transpose()))
# print(np.matmul(np.matmul(eigenvector_M2_NH, M2_M2_dagger_NH), eigenvector_M2_NH.transpose()))
# print(np.matmul(np.matmul(eigenvector_M3_NH, M3_M3_dagger_NH), eigenvector_M3_NH.transpose()))



# print(np.matmul(np.matmul(eigenvector_M2_NH.transpose(), M2_M2_dagger_NH), eigenvector_M2_NH))
# print(np.matmul(np.matmul(eigenvector_M3_NH.transpose(), M3_M3_dagger_NH), eigenvector_M3_NH))

#print(check(M3_M3_dagger_NH,eigenvector_M3_NH,eigenval_M3_NH ))


# eigenvector_M1_NH = eigenvector_M1_NH.transpose()
# eigenvector_M2_NH = eigenvector_M2_NH.transpose()
# eigenvector_M3_NH = eigenvector_M3_NH.transpose()

# eigenvector_M1_IH = eigenvector_M1_IH.transpose()
# eigenvector_M2_IH = eigenvector_M2_IH.transpose()
# eigenvector_M3_IH = eigenvector_M3_IH.transpose()




eigenvector_M1_NH = eigenvector_M1_NH[0, :]
eigenvector_M2_NH = eigenvector_M2_NH[0, :]
eigenvector_M3_NH = eigenvector_M3_NH[0, :]

eigenvector_M1_IH = eigenvector_M1_IH[0, :]
eigenvector_M2_IH = eigenvector_M2_IH[0, :]
eigenvector_M3_IH = eigenvector_M3_IH[0, :]


eigenvector_NH = np.empty((3, 5))
eigenvector_NH[0] = eigenvector_M1_NH
eigenvector_NH[1] = eigenvector_M2_NH
eigenvector_NH[2] = eigenvector_M3_NH

eigenvector_IH = np.empty((3, 5))
eigenvector_IH[0] = eigenvector_M1_IH
eigenvector_IH[1] = eigenvector_M2_IH
eigenvector_IH[2] = eigenvector_M3_IH


eigenvalues_NH = np.empty((3, 5))
eigenvalues_NH[0] = eigenval_M1_NH
eigenvalues_NH[1] = eigenval_M2_NH
eigenvalues_NH[2] = eigenval_M3_NH


eigenvalues_IH = np.empty((3, 5))
eigenvalues_IH[0] = eigenval_M1_IH
eigenvalues_IH[1] = eigenval_M2_IH
eigenvalues_IH[2] = eigenval_M3_IH



# print(eigenvector_NH)

# print(eigenvector_IH)






Energy_start1 = 1e9
Energy_stop1 = 5e9
Energy_start2 = 1e6
Energy_stop2 = 9e6
Energy_start3 = 1e6
Energy_stop3 = 9e6


file_standard_Figure1 = open('standard_probability_to_plot_Figure1.txt', 'w')
file_standard_Figure2 = open('standard_probability_to_plot_Figure2.txt', 'w')
file_standard_Figure3 = open('standard_probability_to_plot_Figure3.txt', 'w')

def probability_mu(file, Energy_start, Energy_stop, number_of_points, length):
    
    E = np.linspace(Energy_start,Energy_stop, number_of_points )
    Probabilities = np.empty(number_of_points)

    for i in range(number_of_points):
        Probabilities[i] = P_mu(E[i], length)
        entry = str(E[i]) + '-' + str(Probabilities[i])+'\n'
        file.write(entry)
    return 


def probability_e(file, Energy_start, Energy_stop, number_of_points, length):
    
    E = np.linspace(Energy_start,Energy_stop, number_of_points )
    Probabilities = np.empty(number_of_points)

    for i in range(number_of_points):
        Probabilities[i] = P_e(E[i], length)
        entry = str(E[i]) + '-' + str(Probabilities[i])+'\n'
        file.write(entry)
    return 

probability_mu(file_standard_Figure1, Energy_start1, Energy_stop1, 5000, L_1)
probability_e(file_standard_Figure2,Energy_start2, Energy_stop2,5000,  L_2 )
probability_e(file_standard_Figure3, Energy_start3, Energy_stop3,5000, L_3)

file_standard_Figure1.close()
file_standard_Figure2.close()
file_standard_Figure3.close()


NH_file_Figure3 = open('NH_probability_to_plot_Figure3.txt', 'w')
IH_file_Figure3 = open('IH_probability_to_plot_Figure3.txt', 'w')
NH_file_Figure2 = open('NH_probability_to_plot_Figure2.txt', 'w')
IH_file_Figure2 = open('IH_probability_to_plot_Figure2.txt', 'w')
NH_file_Figure1 = open('NH_probability_to_plot_Figure1.txt', 'w')
IH_file_Figure1 = open('IH_probability_to_plot_Figure1.txt', 'w')



def Write(file1, file2, Energy_start, Energy_stop, number_of_points, length, alpha, beta, eigenvalues_NH, eigenvectors_NH, eigenvalues_IH, eigenvectors_IH):
    


    E = np.linspace(Energy_start, Energy_stop, number_of_points)
    probabilities_NH = np.empty(number_of_points)
    probabilities_IH = np.empty(number_of_points)
    
    #NH
    for n in range(number_of_points):
        probabilities_NH[n] = Probability(Amplitude(alpha, beta, length, E[n], a_LED, eigenvalues_NH, eigenvectors_NH))  #a, eigenvalue_array, eigenvector_matrix)
        entry = str(E[n]) + '-' + str(probabilities_NH[n])+'\n'
        file1.write(entry)
    #ΙH
    for n in range(number_of_points):
        probabilities_IH[n] = Probability(Amplitude(alpha, beta, length,E[n], a_LED, eigenvalues_IH, eigenvectors_IH))
        entry = str(E[n]) + '-' + str(probabilities_IH[n])+'\n'
        file2.write(entry)
        #print(probabilities_IH[n])
    
    
    return 

Write(NH_file_Figure3, IH_file_Figure3, Energy_start3, Energy_stop3,  5000, L_3, 'e', 'e', eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH )
Write(NH_file_Figure2, IH_file_Figure2,Energy_start2, Energy_stop2,  5000, L_2, 'e', 'e',eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH)
Write(NH_file_Figure1, IH_file_Figure1,Energy_start1, Energy_stop1, 5000, L_1, 'mu', 'mu', eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH )




NH_file_Figure3.close()
IH_file_Figure3.close()
NH_file_Figure2.close()
IH_file_Figure2.close()
NH_file_Figure1.close()
IH_file_Figure1.close()








    



#print(NH_eigenvector_matrix)
# print(IH_eigenvector_matrix)



Plot(Energy_start3, Energy_stop3 , 5000, L_3, 'e', 'e', eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH, 'L3', "E(MeV)" )
Plot(Energy_start2, Energy_stop2 , 5000, L_2, 'e', 'e',eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH ,'L2', "E(MeV)" )
Plot(Energy_start1, Energy_stop1 , 5000, L_1, 'mu', 'mu', eigenvalues_NH, eigenvector_NH,eigenvalues_IH, eigenvector_IH, 'L1', "E(GeV)" )












#-----------------------------------------------------------------------------------------------------------------------





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





