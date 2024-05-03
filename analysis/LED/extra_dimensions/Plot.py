import matplotlib.pyplot as plt

def Sort_values(file_name):
    energy_values = []
    probability_values = []

    # Read the file
    with open(file_name, 'r') as file:
    # Iterate through each line in the file
        for line in file:
        # Split the line into energy and probability using the '-' as a delimiter
            energy, probability = map(float, line.strip().split('-'))
        
        # Append the values to their respective arrays
            energy_values.append(energy)
            probability_values.append(probability)

    # Sort the arrays based on energy values
    sorted_indices = sorted(range(len(energy_values)), key=lambda k: energy_values[k])
    sorted_energy_values = [energy_values[i] for i in sorted_indices]
    sorted_probability_values = [probability_values[i] for i in sorted_indices]

    return sorted_energy_values, sorted_probability_values

#print(Sort_values('NH_probability_to_plot.txt'))

Figure1_NH = Sort_values('NH_probability_to_plot_Figure1.txt')
Figure2_NH = Sort_values('NH_probability_to_plot_Figure2.txt')
Figure3_NH = Sort_values('NH_probability_to_plot_Figure3.txt')

Figure1_IH = Sort_values('IH_probability_to_plot_Figure1.txt')
Figure2_IH = Sort_values('IH_probability_to_plot_Figure2.txt')
Figure3_IH = Sort_values('IH_probability_to_plot_Figure3.txt')

Figure1_standard = Sort_values('standard_probability_to_plot_Figure1.txt')
Figure2_standard = Sort_values('standard_probability_to_plot_Figure2.txt')
Figure3_standard = Sort_values('standard_probability_to_plot_Figure3.txt')

plt.figure(1)
plt.plot(Figure1_NH[0], Figure1_NH[1], label= 'NH_L1')
plt.plot(Figure1_IH[0], Figure1_IH[1], label= 'IH_L1')
plt.plot(Figure1_standard[0], Figure1_standard[1],  label= 'standard_L1')
plt.legend()
plt.ylabel("Probability")
plt.xlabel('E(GeV)')
plt.show()


plt.figure(2)
plt.plot(Figure2_NH[0], Figure2_NH[1], label= 'NH_L2')
plt.plot(Figure2_IH[0], Figure2_IH[1], label= 'IH_L2')
plt.plot(Figure2_standard[0], Figure2_standard[1],  label= 'standard_L2')
plt.legend()
plt.ylabel("Probability")
plt.xlabel('E(MeV)')
plt.show()


plt.figure(3)
plt.plot(Figure3_NH[0], Figure3_NH[1], label= 'NH_L3')
plt.plot(Figure3_IH[0], Figure3_IH[1], label= 'IH_L3')
plt.plot(Figure3_standard[0], Figure3_standard[1],  label= 'standard_L3')
plt.legend()
plt.ylabel("Probability")
plt.xlabel('E(MeV)')
plt.show()