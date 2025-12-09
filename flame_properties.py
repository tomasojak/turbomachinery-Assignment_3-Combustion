# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:23:39 2022

@author: luuk + Rafael + Adam

premixed flame properties
"""

#%% IMPORT PACKAGES:
import sys
import os
import glob
import cantera as ct
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle as pkl
from progressbar import progressbar
script_dir = os.path.dirname(__file__) # The script directory
script_name = os.path.basename(__file__) # The script name


#%% MATPLOTLIB PARAMETERS
# Close all figures
plt.close('all')

# Set font parameters
plt.rcParams.update({
    "text.usetex": True,
    #"text.latex.preamble": r'\usepackage{sfmath}',
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    "font.sans-serif":"Computer Modern Sans Serif",
    "font.size": 11})

# Set the axis and grid to be always below the points
plt.rcParams['axes.axisbelow'] = True

# Set colormap
jet = mpl.colormaps['jet']
viridis = mpl.colormaps['viridis']
tab10 = mpl.colormaps['tab10']
tab20 = mpl.colormaps['tab20']
colors = tab20.colors


#%% FUNCTIONS

# This function creates flame objects for every mixture provided. 
# It returns a list of the flame objects as well as a library containing their properties.
def create_flame_library(script_dir, filename, mixtures, chemical_mechanism = 'gri30.yaml'):
    
    sep = os.path.sep # Gets the OS path separator.
    
    file_dir = script_dir + sep + 'premixed_flame_libraries' # Defines the file directory and create it if it does not exist.
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
    
    file = file_dir + sep + filename + '.pkl' # Creates the full file path based on the directory and name.
    if os.path.isfile(file): # Adds a number to the file path if it is already used by another file.
        file_nr = 0
        while True:
            file_nr = file_nr + 1
            new_file = file.split('.pkl')[0] + '_' + str(file_nr) + '.pkl'
            if os.path.isfile(new_file):
                continue
            else:
                file = new_file
                break
    
    flames = []
    flames_library = []
    for i, mixture in enumerate(progressbar(mixtures)): 
        H2_percentage = mixture['H2_percentage'] # Sets the hydrogen volume percentage [-].
        phi = mixture['phi'] # Sets the mixture's equivalence ratio [-].
        T_u = mixture['T_u'] # Sets the mixture's unburned temperature [K].
        p_u = mixture['p_u'] # Sets the mixture's unburned pressure [Pa].
        
        flame = PremixedFlame(phi, H2_percentage, T_u, p_u, chemical_mechanism) # Creates a flame object (non-pickable).
        flames.append(flame)

        # Extracts the key properties to a dictionary which is pickable and append them to a list. 
        flame_properties = {'index': i + 1, 'phi':flame.phi, 'H2%':flame.H2_percentage, 'T_u':flame.T_u, 'p_u':flame.p_u, 'rho_u':flame.rho_u, 'mu_u':flame.mu_u,
                            'nu_u':flame.nu_u, 'alpha_u':flame.alpha_u, 'X_i_air':flame.X_i_air, 'X_fuel':flame.X_i_fuel, 'X_u':flame.X_u, 'Y_u':flame.Y_u,
                            'D_binary':flame.D_binary, 'D_mix':flame.D_mix, 'PHI': flame.PHI, 'Le_eff_bin':flame.Le_eff_bin, 'Le_eff_mix':flame.Le_eff_mix, 'Le_o':flame.Le_o, 'Delta':flame.Delta,
                            'LHV_mixture':flame.LHV_mixture, 'LHV_fuel':flame.LHV_fuel, 's_l0':flame.s_l0, 'T_ad':flame.T_ad
                    }
        flames_library.append(flame_properties)
        
        print(f'Flame calculated: phi = {round(phi, 2)} [-], H2% = {H2_percentage} [-], T_u = {T_u} [K], p_u = {p_u} [Pa].')
        
    with open(file, 'wb') as outp: # Writes the list containing the flame properties to a file (wb = write binary).
        pkl.dump(flames_library, outp, pkl.HIGHEST_PROTOCOL)
        
    return flames, flames_library

###############################################################################

# This function loads the existing flame library.
def load_flame_library(script_dir, filename):
    
    sep = os.path.sep # Gets the OS path separator.

    file_dir = script_dir + sep + 'premixed_flame_libraries' # Defines the file directory.
    
    # Find all files corresponding to the spcified filename
    matching_files = glob.glob(file_dir + sep + filename + '.pkl')
    
    # Read all the files and collect all flame properties into a single flames library
    loaded_flames_library = []
    for file in matching_files:
        
        with open(file, 'rb') as f:
            loaded_flames_library = loaded_flames_library + pkl.load(f)
        
    return loaded_flames_library

###############################################################################
# REVIEW THIS FUNCTION IF USED
# This function returns the lower heating value (LHV) and higher heating value (HHV) for a given fuel.
# https://cantera.org/dev/userguide/heating-value.html
def heating_value(fuel):
    
    T_u = 293.15 # Sets the unburned temperature of the fuel [K].
    p_u = ct.one_atm # Sets the unburned pressure of the fuel [Pa].
    
    gas = ct.Solution('gri30.yaml')
    # gas = ct.Solution('sandiego20161214_mechCK.yaml')
    gas.TP = T_u, p_u
    gas.set_equivalence_ratio(1.0, fuel, 'O2')
    
    h1 = gas.enthalpy_mass
    Y_fuel = gas[fuel].Y[0]
    
    water = ct.Water()
    water.TQ = T_u, 0 # Sets the liquid water state (vapor fraction x = 0).
    h_liquid = water.h
    water.TQ = T_u, 1 # Sets the gaseous water state (vapor fraction x = 1).
    h_gas = water.h

    x_products = {'CO2': gas.elemental_mole_fraction('C'),
                  'H2O': 0.5 * gas.elemental_mole_fraction('H'),
                  'N2': 0.5 * gas.elemental_mole_fraction('N')} # The complete combustion products. 
    
    gas.TPX = None, None, x_products
    
    Y_H2O = gas['H2O'].Y[0]
    h2 = gas.enthalpy_mass
    
    LHV = -(h2 - h1)/Y_fuel/1e6
    HHV = -(h2 - h1 + (h_liquid - h_gas)*Y_H2O)/Y_fuel/1e6
    
    return LHV, HHV


#%% CLASSES

class UnburntMixture:
    
    def __init__(self, phi, H2_percentage, T_u, p_u = ct.one_atm, chemical_mechanism = 'gri30.yaml'):
        
        ################# SET THE PARAMETERS ################
        
        # Air composition
        X_O2_air = 21/100 # The mole fraction of O2 in air [-].
        X_N2_air = 78/100 # Fraction of N2 in air [-].
        X_Ar_air = 1/100 # Fraction of Ar in air [-].
        X_i_air = {'N2':X_N2_air, 'O2':X_O2_air, 'AR':X_Ar_air} # The air composition expressed in mole fractions [-].
        
        # DNG composition
        X_CH4_DNG = 81.87/100 # The mole fraction of CH4 in the DNG (obtained from the Bronkhorst datasheet) [-].
        X_C2H6_DNG = 3.73/100 # The mole fraction of C2H6 in the DNG (obtained from the Bronkhorst datasheet) [-].
        X_N2_DNG = 14.4/100 # The mole fraction of N2 in the DNG (obtained from the Bronkhorst datasheet) [-].
        X_rest_DNG = 0/100 # The mole fraction of other species present in the DNG [-]. 
        
        # Fuel composition
        X_H2_fuel = H2_percentage/100 # The mole fraction of H2 in the fuel mixture [-].
        X_DNG_fuel = 1 - X_H2_fuel # The mole fraction of DNG in the fuel mixture [-].
        X_CH4_fuel = X_DNG_fuel*X_CH4_DNG # The mole fraction of CH4 in the fuel mixture [-].
        X_C2H6_fuel = X_DNG_fuel*X_C2H6_DNG # The mole fraction of C2H6 in the fuel mixture [-].
        X_N2_fuel = X_DNG_fuel*X_N2_DNG # The mole fraction of N2 in the fuel mixture [-].
        X_rest_fuel = X_DNG_fuel*X_rest_DNG # The mole fraction of other species present in the fuel [-].
        X_i_fuel = {'H2':X_H2_fuel, 'CH4':X_CH4_fuel, 'C2H6':X_C2H6_fuel, 'N2':X_N2_fuel} # The fuel composition expressed in mole fractions [-].
        
        # Checks if the total mole/volume fraction sums up to 1.
        check_air = sum(X_i_air.values()) # Sums up the mole/volume fractions of the air species .
        check_fuel = sum(X_i_fuel.values()) # Sums up the mole/volume fractions of the air species.
        if check_air == 1.0 and round(check_fuel,3) == 1.0:
            pass
        else:
            sys.exit('The fuel or air composition is incorrect!' + ' Sum of air volume fractions: ' + str(check_air) + ';' + ' Sum of fuel volume fractions: ' + str(check_fuel))
        
        # Define the unburnt mixture
        mixture = ct.Solution(chemical_mechanism) # Object representing the gaseous mixture.
        mixture.transport_model = 'multicomponent' # Sets the transport models.
        mixture.set_equivalence_ratio(phi, X_i_fuel, X_i_air) # Computes the final amounts of fuel and air based on the equivalence ratio.
        mixture.TP = T_u, p_u # Sets the state of the gaseous mixture.
        
        ################# COMPUTE THE UNBURNT MIXTURE PROPERTIES ################
                
        # Properties of the unburnt mixture
        X_u = mixture.mole_fraction_dict() # Returns the mole fraction for each species in the unburnt mixture.
        Y_u = mixture.mass_fraction_dict() # Returns the mass fraction for each species in the unburnt mixture.
        h_u = mixture.enthalpy_mass # The specific enthalpy of the unurned mixture [J/kg].
        cp_u = mixture.cp_mass # The specific heat capacity at constant pressure of the unburned mixture [J/kg/K].
        cv_u = mixture.cv_mass # The specific heat capacity at constant volume of the unburned mixture [J/kg/K].
        rho_u = mixture.density_mass # The density of the unburned mixture [kg/m^3].
        mu_u = mixture.viscosity # The dynamic visocity of the unburned mixture [Pa.s].
        nu_u = mu_u/rho_u # The kinematic viscosity of the unburned mixture [m^(2)/s].
        lambda_u = mixture.thermal_conductivity # The thermal conductivity of the unburned mixture [W/m/K].
        alpha_u = lambda_u/(rho_u*cp_u) # The thermal diffusivity of the unburned mixture [m^2/s].
        
        # Get the diffusion coefficients of the mixture
        D_binary = mixture.binary_diff_coeffs # The binary diffusion coefficients [m^2/s].
        D_mix = mixture.mix_diff_coeffs # The mixture-averaged diffusion coefficients [m^2/s].
        
        # Extract the relevant individual coefficients (species indices are those of GRI30)
        D_H2_N2, D_CH4_N2, D_C2H6_N2, D_N2_N2 = D_binary[0][47], D_binary[13][47], D_binary[26][47], D_binary[47][47]
        D_H2_mix, D_CH4_mix, D_C2H6_mix = D_mix[0], D_mix[13], D_mix[26]
        D_O2_mix, D_N2_mix, D_Ar_mix = D_mix[3], D_mix[47], D_mix[48]
        
        # Compute the excess-to-deficient reactant ratio
        if phi <= 1:
            PHI = 1/phi
        else:
            PHI = phi
            
        ############### COMPUTE THE LOWER AND HIGHER HEATING VALUES ###############
        # REVIEW IF USED
        
        Y_fuel = 0 # Initializes a variable for storing the total mass fraction of the fuel [-].
        LHV_fuel = 0 # Initializes a variable for storing the LHV of the fuel [MJ/kg].
        LHV_mixture = 0 # Initializes a variable for storing the LHV of the mixture [MJ/kg].
        
        for fuel_key in X_i_fuel.keys(): # Calculates the total mass fraction of the fuel.
            if Y_u.get(fuel_key) is not None and fuel_key in ['H2', 'CH4', 'C2H6']:
                Y_fuel = Y_fuel + Y_u[fuel_key]
        
        for fuel_key in X_i_fuel.keys(): # Calculates the LHV and HHV of the fuel.
            if Y_u.get(fuel_key) is not None and fuel_key in ['H2', 'CH4', 'C2H6']:
                LHV, HHV = heating_value(fuel_key)
                LHV_fuel = LHV_fuel + (LHV*Y_u[fuel_key])/Y_fuel
                LHV_mixture = LHV_mixture + (LHV*Y_u[fuel_key])
        
        ################# COMPUTE THE LEWIS NUMBER AND MASS DIFFUSIVITY RATIO #################
        
        # Compute the effective Lewis number (laminar flames), Luuk's approach (simple estimate)
        Le_eff_bin = alpha_u/(1 - X_N2_fuel) * (X_H2_fuel/D_H2_N2 + X_CH4_fuel/D_CH4_N2 + X_C2H6_fuel/D_C2H6_N2)
        
        # Compute the fuel and oxidiser Lewis numbers (laminar flames), Lapalme's approach
        Le_f = alpha_u/(1 - X_N2_fuel) * (X_H2_fuel/D_H2_mix + X_CH4_fuel/D_CH4_mix + X_C2H6_fuel/D_C2H6_mix)
        Le_o = alpha_u/D_O2_mix
        
        # Compute the fuel to oxidiser mass diffusivity ratio
        D_f_mix = 1/(1 - X_N2_fuel) * (X_H2_fuel*D_H2_mix + X_CH4_fuel*D_CH4_mix + X_C2H6_fuel*D_C2H6_mix)
        D_o_mix = D_O2_mix
        Delta = D_f_mix/D_o_mix
        
        ############### OBJECT ATTRIBUTES ###############
        
        self.mixture = mixture
        self.X_i_air = X_i_air
        self.X_i_fuel = X_i_fuel
        self.X_u = X_u
        self.Y_u = Y_u
        self.h_u = h_u
        self.cp_u = cp_u
        self.cv_u = cv_u
        self.rho_u = rho_u
        self.mu_u = mu_u
        self.nu_u = nu_u
        self.lambda_u = lambda_u
        self.alpha_u = alpha_u
        self.D_binary = D_binary
        self.D_mix = D_mix
        self.PHI = PHI
        self.LHV_fuel = LHV_fuel
        self.LHV_mixture = LHV_mixture
        self.Le_eff_bin = Le_eff_bin
        self.Le_f = Le_f
        self.Le_o = Le_o
        self.Delta = Delta
        
        

class PremixedFlame:
    
    def __init__(self, phi, H2_percentage, T_u, p_u = ct.one_atm, chemical_mechanism = 'gri30.yaml'):
        
        # Ideal gas constant
        R = ct.gas_constant*1e-3 # [J/mol/K]
        
        # Create the mixture (in an unburnt state)
        unburnt_mixture = UnburntMixture(phi, H2_percentage, T_u, p_u, chemical_mechanism)
        mixture = unburnt_mixture.mixture # Cantera mixture object (in unburnt state)
        
################# COMBUST THE MIXTURE USING A LAMINAR UNSTRETCHED FREE FLAME (1D MODEL) ################
    
        width = 0.03 # Sets the 1-D domain size [m].
        flame = ct.FreeFlame(mixture, width = width) # Creates an object for the freely-propagating premixed flames.
        flame.set_refine_criteria(ratio = 3, slope = 0.06, curve = 0.12) # Sets the domain refining criteria.
        flame.solve(loglevel = 1, auto = True)
        
        s_l0 = flame.velocity[0] # The unstretched laminar flame speed [m/s].
        T_ad = mixture.T # The adiabatic flame temperature [K].

######### COMPUTE THE BURNT MIXTURE PROPERTIES ###############################
        
        h_b = mixture.enthalpy_mass # The specific enthalpy of the burned mixture [J/kg].
        cp_b = mixture.cp_mass # The specific heat capacity at constant pressure of the burned mixture [J/kg/K].
        cv_b = mixture.cv_mass # The specific heat capacity at constant volume of the burned mixture [J/kg/K].
        rho_b = mixture.density_mass # The density of the burned mixture [kg/m^3].
        mu_b = mixture.viscosity # The dynamic visocity of the burned mixture [Pa.s].
        nu_b = mu_b/rho_b # The kinematic viscosity of the burned mixture [m^(2)/s].
        lambda_b = mixture.thermal_conductivity # The thermal conductivity of the burned mixture [W/m/K].
        alpha_b = lambda_b/(rho_b*cp_b) # The thermal diffusivity of the burned mixture [m^2/s].
        #Ze = E_a*(T_ad - T_u)/(R*T_ad**2) # Zeldovich number [-]

################# COMPUTE THE LEWIS NUMBER AND MASS DIFFUSIVITY RATIO OF THE UNBURNT MIXTURE #################
        
        # Compute the effective Lewis number (laminar flames), Lapalme's approach
        Le_eff_mix = unburnt_mixture.Le_f # Good estimate for mixtures where phi < 0.8
        
        # Correction for near-stochiometry mixtures not implemented yet (not able to clearly define activation energy)
        # if phi < 0.8:
        #     Le_eff_mix = Le_f
        # elif (0.8 <= phi) & (phi <= 1):
        #     A_1 = 1 + Ze*(PHI - 1)
        #     Le_eff_mix = 1 + ((Le_o - 1) + (Le_f - 1)*A_1)/(1 + A_1)

################# OBJECT ATTRIBUTES ####################

        # Experimental parameters
        self.phi = phi # The equivalence ratio [-].
        self.H2_percentage = H2_percentage # The hydrogen percentage by volume (in comparison to the DNG) [-].
        self.T_u = T_u # The unburned temperature of the fuel-air mixture [-].
        self.p_u = p_u # The unburned pressure of the fuel-air mixture [-].
        
        # Properties of the unburnt mixture
        self.X_i_air = unburnt_mixture.X_i_air
        self.X_i_fuel = unburnt_mixture.X_i_fuel
        self.X_u = unburnt_mixture.X_u
        self.Y_u = unburnt_mixture.Y_u
        self.LHV_fuel = unburnt_mixture.LHV_fuel
        self.LHV_mixture = unburnt_mixture.LHV_mixture     
        self.h_u = unburnt_mixture.h_u
        self.cp_u = unburnt_mixture.cp_u
        self.cv_u = unburnt_mixture.cv_u
        self.rho_u = unburnt_mixture.rho_u
        self.mu_u = unburnt_mixture.mu_u
        self.nu_u = unburnt_mixture.nu_u
        self.lambda_u = unburnt_mixture.lambda_u
        self.alpha_u = unburnt_mixture.alpha_u
        self.D_binary = unburnt_mixture.D_binary
        self.D_mix = unburnt_mixture.D_mix
        self.PHI = unburnt_mixture.PHI
        self.Le_eff_bin = unburnt_mixture.Le_eff_bin
        self.Le_o = unburnt_mixture.Le_o
        self.Delta = unburnt_mixture.Delta     
        
        # Flame properties
        self.flame = flame
        self.T_ad = T_ad
        self.s_l0 = s_l0
        self.Le_eff_mix = Le_eff_mix
        
        # Burnt mixture properties
        self.h_b = h_b
        self.cp_b = cp_b
        self.cv_b = cv_b
        self.rho_b = rho_b
        self.mu_b = mu_b
        self.nu_b = nu_b
        self.lambda_b = lambda_b
        self.alpha_b = alpha_b

#%% CREATE FLAME LIBRARY

if __name__ == "__main__":

##################################################################################
#            Creates a new flame library (COMMENT OUT IF NOT IN USE):            # 
##################################################################################
        
    print(f"Running Cantera Version: {ct.__version__}")
    
    mixtures = []
    phis = np.arange(0.2, 1.05, 0.05)
    for phi in phis:
        mixture = {'phi': phi, 'H2_percentage': 100, 'T_u': 423.15, 'p_u':ct.one_atm} # 293.15/423.15
        mixtures.append(mixture)
    
    flames, flames_lib = create_flame_library(script_dir, 'phi_020_100_H2_100_Tu_150_pu_101325_GRI30_rafael', mixtures)
    
    
#%% LOAD FLAME LIBRARIES FOR LE AND DELTA
    
    flames_lib = load_flame_library(script_dir, 'phi_*_H2_000_100_Tu_*_pu_101325_GRI30_rafael')
    
    phis = []
    H2s = []
    T_us = []
    Le_effs = []
    Le_eff_bins = []
    Le_os = []
    Deltas = []
    for flame_properties in flames_lib:
        phis.append(flame_properties['phi'])
        H2s.append(flame_properties['H2%'])
        T_us.append(flame_properties['T_u'] - 273.15) # convert to celsius
        Le_effs.append(flame_properties['Le_eff_mix'])
        Le_eff_bins.append(flame_properties['Le_eff_bin'])
        Le_os.append(flame_properties['Le_o'])
        Deltas.append(flame_properties['Delta'])
    
    phis = np.array(phis)
    H2s = np.array(H2s)
    T_us = np.array(T_us)
    Le_effs = np.array(Le_effs)
    Le_eff_bins = np.array(Le_eff_bins)
    Le_os = np.array(Le_os)
    Deltas = np.array(Deltas)
        
#%% PLOTS FOR LE AND DELTA
    
    # Create empty lists of figures, axes and figure names
    figs = []
    axs = []
    fig_names = []
    
    
    # Plot effective Lewis number versus fuel hydrogen content
    fig_name = 'Le_vs_X_H2'
    fig, ax = plt.subplots()
    phi_unique = np.unique(phis)
    for i, phi in enumerate(phi_unique):
        x_1 = H2s[(phis == phi) & (T_us == 20)]
        y_1 = Le_effs[(phis == phi) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'20, {}'.format(phi), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
    for i, phi in enumerate(phi_unique):    
        x_2 = H2s[(phis == phi) & (T_us == 150)]
        y_2 = Le_effs[(phis == phi) & (T_us == 150)]
        ax.scatter(x_2, y_2, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'150, {}'.format(phi), zorder=2)
        ax.plot(x_2, y_2, c=colors[2*i+1], ls='--', zorder=1)
        
    # Set the limits, labels and legend
    ax.set_aspect(140)
    xlabel = r'$X_{\mathrm{H}_2, f}$ [\%]'
    ylabel = r'$\mathrm{Le_{eff}}$ [--]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$T_u$ [$^{\circ}$C], $\phi$ [--]', loc='best', ncol=2, handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot oxidiser Lewis number versus fuel hydrogen content
    fig_name = 'Le_o_vs_X_H2'
    fig, ax = plt.subplots()
    phi_unique = np.unique(phis)
    for i, phi in enumerate(phi_unique):
        x_1 = H2s[(phis == phi) & (T_us == 20)]
        y_1 = Le_os[(phis == phi) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'20, {}'.format(phi), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
    for i, phi in enumerate(phi_unique):    
        x_2 = H2s[(phis == phi) & (T_us == 150)]
        y_2 = Le_os[(phis == phi) & (T_us == 150)]
        ax.scatter(x_2, y_2, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'150, {}'.format(phi), zorder=2)
        ax.plot(x_2, y_2, c=colors[2*i+1], ls='--', zorder=1)
        
    # Set the limits, labels and legend
    ax.set_aspect(200)
    xlabel = r'$X_{\mathrm{H}_2, f}$ [\%]'
    ylabel = r'$\mathrm{Le}_{o}$ [--]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$T_u$ [$^{\circ}$C], $\phi$ [--]', loc='best', ncol=2, handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot diffusivity ratio versus fuel hydrogen content
    fig_name = 'Delta_vs_X_H2'
    fig, ax = plt.subplots()
    phi_unique = np.unique(phis)
    for i, phi in enumerate(phi_unique):
        x_1 = H2s[(phis == phi) & (T_us == 20)]
        y_1 = Deltas[(phis == phi) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'20, {}'.format(phi), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
    for i, phi in enumerate(phi_unique):    
        x_2 = H2s[(phis == phi) & (T_us == 150)]
        y_2 = Deltas[(phis == phi) & (T_us == 150)]
        ax.scatter(x_2, y_2, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'150, {}'.format(phi), zorder=2)
        ax.plot(x_2, y_2, c=colors[2*i+1], ls='--', zorder=1)
        
    # Set the limits, labels and legend
    ax.set_aspect(33)
    xlabel = r'$X_{\mathrm{H}_2, f}$ [\%]'
    ylabel = r'$\Delta$ [--]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$T_u$ [$^{\circ}$C], $\phi$ [--]', loc='best', ncol=2, handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    #%% LOAD FLAME LIBRARIES FOR NU
        
    flames_lib = load_flame_library(script_dir, 'phi_*_100_H2_*_pu_101325_GRI30_rafael')
    
    phis = []
    H2s = []
    T_us = []
    nus = []
    for flame_properties in flames_lib:
        phis.append(flame_properties['phi'])
        H2s.append(flame_properties['H2%'])
        T_us.append(flame_properties['T_u'] - 273.15) # convert to celsius
        nus.append(flame_properties['nu_u'])
    
    phis = np.array(phis)
    H2s = np.array(H2s)
    T_us = np.array(T_us)
    nus = np.array(nus)
    
    # Compute kinematic viscosity ratio between heated and ambient temperature mixtures
    phis_20 = phis[T_us == 20]
    H2s_20 = H2s[T_us == 20]
    nu_ratios = nus[T_us == 150]/nus[T_us == 20]
    
            
    #%% PLOTS FOR NU
        
    # Create empty lists of figures, axes and figure names
    figs = []
    axs = []
    fig_names = []
    
    
    # Plot kinematic viscosity versus equivalence ratio for ambient temperature mixtures
    fig_name = 'nu_vs_phi_T_u_20'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):
        x_1 = phis[(H2s == H2) & (T_us == 20)]
        y_1 = nus[(H2s == H2) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
    
    # Set the limits, labels and legend
    ax.set_aspect(1.42e5)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$\nu_u$ [m$^2$/s]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', ncol=1,
              handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot kinematic viscosity versus equivalence ratio for heated mixtures
    fig_name = 'nu_vs_phi_T_u_150'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):
        x_1 = phis[(H2s == H2) & (T_us == 150)]
        y_1 = nus[(H2s == H2) & (T_us == 150)]
        ax.scatter(x_1, y_1, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], ls='--', zorder=1)
   
    # Set the limits, labels and legend
    ax.set_aspect(7.6e4)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$\nu_u$ [m$^2$/s]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', ncol=1,
              handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot kinematic viscosity ratio (heated to ambient) versus equivalence ratio
    fig_name = 'nu_ratio_vs_phi'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):
        x_1 = phis_20[H2s_20 == H2]
        y_1 = nu_ratios[H2s_20 == H2]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
    
    # Set the limits, labels and legend
    ax.set_aspect(85)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$\nu_{u,h}/\nu_{u,0}$ [--]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', ncol=1,
              handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
            
#%% LOAD FLAME LIBRARIES FOR T_AD AND S_L0
    
    flames_lib = load_flame_library(script_dir, 'phi_*_100_H2_*_pu_101325_refined_SGN_2')
    
    phis = []
    H2s = []
    T_us = []
    T_ads = []
    s_l0s = []
    for flame_properties in flames_lib:
        phis.append(flame_properties['phi'])
        H2s.append(flame_properties['H2%'])
        T_us.append(flame_properties['T_u'] - 273.15) # convert to celsius
        T_ads.append(flame_properties['T_ad'])
        s_l0s.append(flame_properties['S_L0'])
    
    phis = np.array(phis)
    H2s = np.array(H2s)
    T_us = np.array(T_us)
    T_ads = np.array(T_ads)
    s_l0s = np.array(s_l0s)
    

#%% PLOTS FOR T_AD AND S_L0
    
    # Create empty lists of figures, axes
    figs = []
    axs = []
    fig_names = []
    
    
    # Plot adiabatic temperature versus equivalence ratio for ambient temperature mixtures
    fig_name = 'T_ad_vs_phi_T_u_20'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):
        x_1 = phis[(H2s == H2) & (T_us == 20)]
        y_1 = T_ads[(H2s == H2) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
        
    # Set the limits, labels and legend
    ax.set_ylim([850, 2500])
    ax.set_aspect(4.85e-4)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$T_{\mathrm{ad}}$ [K]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot adiabatic temperature versus equivalence ratio for heated mixtures
    fig_name = 'T_ad_vs_phi_T_u_150'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):   
        x_1 = phis[(H2s == H2) & (T_us == 150)]
        y_1 = T_ads[(H2s == H2) & (T_us == 150)]
        ax.scatter(x_1, y_1, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], ls='--', zorder=1)
        
    # Set the limits, labels and legend
    ax.set_ylim([850, 2500])
    ax.set_aspect(4.85e-4)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$T_{\mathrm{ad}}$ [K]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot unstretched laminar flame speed versus equivalence ratio for ambient temperature mixtures
    fig_name = 's_l0_vs_phi_T_u_20'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):
        x_1 = phis[(H2s == H2) & (T_us == 20)]
        y_1 = s_l0s[(H2s == H2) & (T_us == 20)]
        ax.scatter(x_1, y_1, marker='o', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], zorder=1)
         
    # Set the limits, labels and legend
    ax.set_ylim([-0.2, 4.3])
    ax.set_aspect(0.2)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$s_{l, 0}$ [m/s]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='upper left', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot unstretched laminar flame speed versus equivalence ratio for heated mixtures
    fig_name = 's_l0_vs_phi_T_u_150'
    fig, ax = plt.subplots()
    H2_unique = np.unique(H2s)
    for i, H2 in enumerate(H2_unique):   
        x_1 = phis[(H2s == H2) & (T_us == 150)]
        y_1 = s_l0s[(H2s == H2) & (T_us == 150)]
        ax.scatter(x_1, y_1, marker='^', s=None, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=r'{}'.format(H2), zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], ls='--', zorder=1)
        
    # Set the limits, labels and legend
    ax.set_ylim([-0.2, 4.3])
    ax.set_aspect(0.2)
    xlabel = r'$\phi$ [--]'
    ylabel = r'$s_{l, 0}$ [m/s]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=r'$X_{\mathrm{H}_2, f}$ [\%]', loc='best', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)

        
#%% SAVE THE FIGURES
    
    # Specify the figures directory
    sep = os.path.sep # Gets the OS path separator.
    sub_dir = script_name.removesuffix('.py')
    fig_dir = 'U:' + sep + 'HELIOS' + sep + 'rpichler' + sep + 'Research' + sep + 'Figures' + sep + sub_dir
    # Create the figures directory if it does not exist
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
        
    # Apply tight_layout to each figure and save it to the specified folder
    for i, fig in enumerate(figs):
        fig_nr = i + 1
        # Specify which figures to save
        if 1 <= fig_nr <= len(figs):
            fig.tight_layout()
            fig_name_suffix = '.svg'
            fig_path = fig_dir + sep + fig_names[i] + fig_name_suffix
            fig.savefig(fig_path, format='svg', dpi=500, bbox_inches='tight')

