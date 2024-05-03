import numpy as np
import matplotlib.pyplot as plt
import Machado 


N_points = Machado.N_points


def Sort_values(file_name):
    R_LED_values = []
    N_events_NH_values = []
    N_events_IH_values = []

    # Read the file
    with open(file_name, 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Split the line into R_LED, N_events_LED_nu_e_NH, and N_events_LED_nu_e_IH using '-' as a delimiter
            #R_LED, N_events_NH, N_events_IH = map(float, line.strip().split('-'))
            R_LED, N_events_NH, N_events_IH = map(
            float,
            map(lambda x: float(x) if ('e' in x or 'E' in x) else x, line.strip().split('-'))
            )
            # Append the values to their respective arrays
            R_LED_values.append(R_LED)
            N_events_NH_values.append(N_events_NH)
            N_events_IH_values.append(N_events_IH)

    # Sort the arrays based on R_LED values
    sorted_indices = sorted(range(len(R_LED_values)), key=lambda k: R_LED_values[k])
    sorted_R_LED_values = [R_LED_values[i] for i in sorted_indices]
    sorted_N_events_NH_values = [N_events_NH_values[i] for i in sorted_indices]
    sorted_N_events_IH_values = [N_events_IH_values[i] for i in sorted_indices]

    return sorted_R_LED_values, sorted_N_events_NH_values, sorted_N_events_IH_values



R_LED, N_events_nubar_mu_NH, N_events_nubar_mu_IH = Sort_values('src/N_events_nubar_mu.txt')
R_LED, N_events_nu_e_NH, N_events_nu_e_IH = Sort_values('src/N_events_nu_e.txt')

with open('src/N_events_standard.txt', 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Split the line into N_events_nubar_mu and N_events_nu_e using '-' as a delimiter
            N_events_nubar_mu, N_events_nu_e = map(float, line.strip().split('-'))

#print(N_events_nubar_mu)


def chi_squared(N_exp, N_theo):
    χ_squared = np.empty(N_points)

    for i in range(N_points):
        χ_squared[i] = (N_exp - N_theo[i])**2 /N_exp

    return χ_squared


def Plot(R_LED, N_events_NH, N_events_IH, N_events_standard, title):
     
    χ_squared_NH = chi_squared(N_events_standard, N_events_NH)
    χ_squared_IH = chi_squared(N_events_standard, N_events_IH)

    plt.figure()
    plt.plot(R_LED, χ_squared_NH, label = 'NH')
    plt.plot(R_LED, χ_squared_IH, label = 'IH')
    plt.legend()

    plt.xlabel(r'$R_{LED}(μm)$')  #l_0 = 1μm
    plt.ylabel(r'$\chi^2$')
    plt.xscale('log')
    plt.title(title)
    plt.show()
    return

Plot(R_LED, N_events_nubar_mu_NH, N_events_nubar_mu_IH, N_events_nubar_mu, r'$\bar ν_μ \rightarrow \bar ν_μ$')
Plot(R_LED, N_events_nu_e_NH, N_events_nu_e_IH, N_events_nu_e, r'$ ν_e \rightarrow  ν_e$')
     
