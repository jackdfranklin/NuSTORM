import numpy as np

import sys
import os

#Add src directory to the path so python can find the files
sys.path.append(os.path.abspath("../src"))
import flux
import experiment



def events_test(flavour, flux_file, results):

    pd = experiment.ProtoDuneLike(flavour , flux_file) 

    tot_events = np.sum(pd.get_tot_Ar_events(10**15)) 

    if (abs(tot_events - results[0]) < 1e-6):

        print("Total events test passed! N = {:.3f}".format(tot_events))
    
    else:

        print("Total events test failed... expected N = 4.6401753, got N = {:.7f}".format(tot_events))

    CC_events = np.sum(pd.get_CC_Ar_events(10**15)) 

    if (abs(CC_events - results[1]) < 1e-6):

        print("CC events test passed! N = {:.3f}".format(CC_events))

    else:

        print("CC events test failed... expected N = 3.2050076, got N = {:.7f}".format(CC_events))

    NC_events = np.sum(pd.get_NC_Ar_events(10**15)) 

    if (abs(NC_events - results[2]) < 1e-6):

        print("NC events test passed! N = {:.3f}".format(NC_events))

    else:

        print("NC events test failed... expected N = 1.4351677, got N = {:.7f}".format(NC_events))


if __name__ == '__main__':

    events_test("nubar_mu", "../resources/fluxes/E5spectraMuSig558Numu.txt", [5.4527597, 3.7944726, 1.6582871])

    events_test("nu_e", "../resources/fluxes/E5spectraMuSig558Nue.txt", [12.7733192, 9.6435353, 3.1297839])
