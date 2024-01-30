import numpy as np
import matplotlib.pyplot as plt
from experiment import Experiment, ProtoDuneLike


protoDune = ProtoDuneLike()
POT =  10e7
conversion_factor_μm_to_eV = 1.9732705*10**(-1)
a = 0.6/conversion_factor_μm_to_eV

bin_midpoints_nubar_mu, events_NH_nubar_mu, events_IH_nubar_mu, events_standard_nubar_mu = protoDune.get_event_bins_nubar_mu(POT, a, 20)
bin_midpoints_nu_e, events_NH_nu_e, events_IH_nu_e, events_standard_nu_e = protoDune.get_event_bins_nu_e(POT, a, 20)

plt.figure(1)
plt.scatter(bin_midpoints_nubar_mu, events_NH_nubar_mu, label = 'NH', color = 'blue',  s = 30, marker = 'd')
plt.scatter(bin_midpoints_nubar_mu, events_IH_nubar_mu, label = 'IH', color = 'orange',  s = 30, marker = 'x')
plt.scatter(bin_midpoints_nubar_mu, events_standard_nubar_mu, label = 'SM', color = 'green',  s = 30)
plt.ylabel('Number of events')
plt.xlabel('Energy bin midpoint(GeV)')
plt.legend()
plt.title(r'$\bar ν_μ \rightarrow \bar ν_μ$')
plt.show()

plt.figure(2)
plt.scatter(bin_midpoints_nu_e, events_NH_nu_e, label = 'NH', color = 'blue',   s = 30, marker = 'd')
plt.scatter(bin_midpoints_nu_e, events_IH_nu_e, label = 'IH', color = 'orange',  s = 30, marker = 'x')
plt.scatter(bin_midpoints_nu_e, events_standard_nu_e, label = 'SM', color = 'green',  s = 30)
plt.ylabel('Number of events')
plt.xlabel('Energy bin midpoint(GeV)')
plt.legend()
plt.title(r'$ ν_e \rightarrow ν_e$')
plt.show()