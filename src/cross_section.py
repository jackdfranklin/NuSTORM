import numpy as np
from scipy import interpolate

class Cross_Section:

    def __init__(self, cross_section_file):

        raw_data = np.loadtxt(cross_section_file, skiprows=3)

        #Store values of energy midpoints and integrated cross section(s)
        self.energies = raw_data[:,0]
        self.CC_cross_section = raw_data[:,1]
        self.NC_cross_section = raw_data[:,2]

        #Interpolate cross sections
        self.interp_CC = interpolate.UnivariateSpline(self.energies, self.CC_cross_section)
        self.interp_NC = interpolate.UnivariateSpline(self.energies, self.NC_cross_section)

        #Dictionary to keep track of units
        self.units = {
            "Energy":"GeV",
            "Cross Section":"cm^2 nucleus^-1"
        } 

    def total_cross_section(self, nu_energy):
        '''Returns the total cross section (CC + NC) at the given energy

        Parameters
        ----------
        nu_energy: float or numpy.ndarray
                The neutrino energy or energies that the cross section will be evaluated at
        
        Returns
        -------
        cross_section: float or numpy.ndarray
        '''
        return self.interp_CC(nu_energy) + self.interp_NC(nu_energy)

    def CC_cross_section(self, nu_energy):
        '''Returns the CC cross section at the given energy
        Parameters
        ----------
        nu_energy: float or numpy.ndarray
                The neutrino energy or energies that the cross section will be evaluated at
        
        Returns
        -------
        cross_section: float or numpy.ndarray
        '''

        return self.interp_CC(nu_energy)


    def NC_cross_section(self, nu_energy):
        '''Returns the NC cross section at the given energy
        Parameters
        ----------
        nu_energy: float or numpy.ndarray
                The energy or energies that the cross section will be evaluated at
        
        Returns
        -------
        cross_section: float or numpy.ndarray
        '''

        return self.interp_NC(nu_energy)

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

