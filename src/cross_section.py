import numpy as np
import pandas as pd
from scipy import interpolate

class Cross_Section:

    def __init__(self, cross_section_file):

        #Read table in as dataframe
        df = pd.read_csv(cross_section_file, header=0, index_col=0)

        #Store values of energy midpoints and integrated cross section(s)
        self.energies = df["bins"].to_numpy()
        self._CC_cross_section = df["tot_cc"].to_numpy()

        self._NC_cross_section = df["tot_nc"].to_numpy()

        #Dictionary to keep track of units
        self.units = {
            "Energy":"GeV",
            "Cross Section":"cm^2 nucleus^-1"
        } 

    def get_cross_section(self, energies, cross_sections):
        '''Aligns cross section array to energy bins, averages if bin width smaller than that of energies'''

        aligned_cross_sections = cross_sections[(self.energies>=energies[0])&(self.energies<=energies[-1])] #Only look at bins within valid range
        
        bin_width_ratio = int(round((energies[1]-energies[0])/(self.energies[1]-self.energies[0]))) #Assumes bin width of cross sections <= bin width of fluxes
        return np.add.reduceat(aligned_cross_sections, np.arange(0, aligned_cross_sections.size, bin_width_ratio))/bin_width_ratio #Averages over valid bin width

    def total_cross_section(self, energies):
        '''Returns the total cross section (CC + NC) array

        Returns
        -------
        cross_section: numpy.ndarray
        '''
        return self.get_cross_section(energies, self._CC_cross_section + self._NC_cross_section)

    def CC_cross_section(self, energies):
        '''Returns the CC cross section
        
        Returns
        -------
        cross_section: numpy.ndarray
        '''

        return self.get_cross_section(energies, self._CC_cross_section)


    def NC_cross_section(self, energies):
        '''Returns the NC cross section
        
        Returns
        -------
        cross_section: numpy.ndarray
        '''

        return self.get_cross_section(energies, self._NC_cross_section)

class Electron_Scattering(Cross_Section):

    sin2_thetaW = 0.2315
    sin4_thetaW = sin2_thetaW**2
    m_e = 0.511 * 10**-3 #GeV, PDG21
    m_e_2 = m_e**2

    const_dict = {
        "12": (1+2*sin2_thetaW)**2 + (4/3)*sin4_thetaW,
        "-12": (1/3)*(1+2*sin2_thetaW)**2 + 4*sin4_thetaW,
        "13": (1-2*sin2_thetaW)**2 + (4/3)*sin4_thetaW,
        "-13": (1/3)*(1-2*sin2_thetaW)**2 + 4*sin4_thetaW,
        "14": (1-2*sin2_thetaW)**2 + (4/3)*sin4_thetaW,
        "-14": (1/3)*(1-2*sin2_thetaW)**2 + 4*sin4_thetaW,
    }
    
    lepton_mass_dict = {

        "12": m_e,
        "-12": m_e,
        "13": 0.106, #Gev, PDG21
        "-13": 0.106, #GeV, PDG21
        "14": 1.76, #GeV, PDG14
        "-14": 1.76, #GeV, PDG14

    }

    def __init__(self, PDG_ID):

        
        self.total_cross_section_const = self.const_dict[PDG_ID]
        self.lepton_mass = self.lepton_mass_dict[PDG_ID]

