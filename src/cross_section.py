import numpy as np
import pandas as pd
from scipy import interpolate

class Cross_Section:

    def __init__(self, cross_section_file):

        #Read table in as dataframe
        df = pd.read_csv(cross_section_file, header=0, index_col=0)

        #Store values of energy midpoints and integrated cross section(s)
        self.energies = df["bins"].to_numpy()
        self.CC_cross_section = df["tot_cc"].to_numpy()
        self.NC_cross_section = df["tot_nc"].to_numpy()

        #Dictionary to keep track of units
        self.units = {
            "Energy":"GeV",
            "Cross Section":"cm^2 nucleus^-1"
        } 

    def total_cross_section(self):
        '''Returns the total cross section (CC + NC) array

        Returns
        -------
        cross_section: numpy.ndarray
        '''
        return self.CC_cross_section + self.NC_cross_section()

    def CC_cross_section(self):
        '''Returns the CC cross section
        
        Returns
        -------
        cross_section: numpy.ndarray
        '''

        return self.CC_cross_section


    def NC_cross_section(self):
        '''Returns the NC cross section
        
        Returns
        -------
        cross_section: numpy.ndarray
        '''

        return self.NC_cross_section

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

