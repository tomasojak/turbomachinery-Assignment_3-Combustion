# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:23:39 2022

@author: luuk + Rafael + Adam

premixed flame properties
"""

#%% PACKAGES
import sys
import os
import numpy as np
import cantera as ct
import matplotlib as mpl
from matplotlib import pyplot as plt
script_dir = os.path.dirname(__file__) # The script directory
script_name = os.path.basename(__file__) # The script name


#%% MATPLOTLIB PARAMETERS
# Close all figures
plt.close('all')

# Set plot parameters
plt.rcParams.update({
    'axes.axisbelow': True, # Set the axis and grid to be always below the points
    "font.size": 11})

# Set colormap
tab20 = mpl.colormaps['tab20']
colors = tab20.colors


#%% CLASSES

class UnburntMixture:
    
    def __init__(self, phi, H2_percentage, T_u = 293.15, p_u = ct.one_atm, chemical_mechanism = 'gri30.yaml'):
        
        ################# SET THE PARAMETERS ################
        
        # Air composition
        X_O2_air = 21/100 # The mole fraction of O2 in air [-].
        X_N2_air = 79/100 # Fraction of N2 in air [-].
        X_i_air = {'N2':X_N2_air, 'O2':X_O2_air} # The air composition expressed in mole fractions [-].
        
        # Natural gas (NG) composition
        X_CH4_NG = 1 # The mole fraction of CH4 in the NG [-].
       
        # Fuel composition
        X_H2_fuel = H2_percentage/100 # The mole fraction of H2 in the fuel mixture [-].
        X_NG_fuel = 1 - X_H2_fuel # The mole fraction of NG in the fuel mixture [-].
        X_CH4_fuel = X_NG_fuel*X_CH4_NG # The mole fraction of CH4 in the fuel mixture [-].
        X_i_fuel = {'H2':X_H2_fuel, 'CH4':X_CH4_fuel} # The fuel composition expressed in mole fractions [-].
        
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
        D_H2_N2, D_CH4_N2 = D_binary[0][47], D_binary[13][47]
        D_H2_mix, D_CH4_mix = D_mix[0], D_mix[13]
        D_O2_mix = D_mix[3]
        
        # Compute the excess-to-deficient reactant ratio
        if phi <= 1:
            PHI = 1/phi
        else:
            PHI = phi
            
        
        ################# COMPUTE THE LEWIS NUMBER AND MASS DIFFUSIVITY RATIO #################
        
        # Compute the effective Lewis number (laminar flames), Luuk's approach (simple estimate)
        Le_eff_bin = alpha_u * (X_H2_fuel/D_H2_N2 + X_CH4_fuel/D_CH4_N2)
        
        # Compute the fuel and oxidiser Lewis numbers (laminar flames), Lapalme's approach
        Le_f = alpha_u * (X_H2_fuel/D_H2_mix + X_CH4_fuel/D_CH4_mix)
        Le_o = alpha_u/D_O2_mix
        
        # Compute the fuel to oxidiser mass diffusivity ratio
        D_f_mix = (X_H2_fuel*D_H2_mix + X_CH4_fuel*D_CH4_mix)
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
        self.Le_eff_bin = Le_eff_bin
        self.Le_f = Le_f
        self.Le_o = Le_o
        self.Delta = Delta
        
        

class PremixedFlame:
    
    def __init__(self, phi, H2_percentage, T_u = 293.15, p_u = ct.one_atm, chemical_mechanism = 'gri30.yaml'):
        
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

################# COMPUTE THE LEWIS NUMBER AND MASS DIFFUSIVITY RATIO OF THE UNBURNT MIXTURE #################
        
        # Compute the effective Lewis number (laminar flames), Lapalme's approach
        Le_eff_mix = unburnt_mixture.Le_f # Good estimate for mixtures where phi < 0.8

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

#%% 1-D FLAME SIMULATIONS

if __name__ == "__main__":

##################################################################################
#            Creates a new flame library (COMMENT OUT IF NOT IN USE):            # 
##################################################################################
        
    print(f"Running Cantera Version: {ct.__version__}")
    
    # Simulation parameters
    X_H2s = np.array([0, 100], dtype=int) # fuel hydrogen fraction [%]
    phis = np.array([0.8, 1, 1.3]) # equivalence ratios [-]
    
    # Run the simulations and collect T_ad and s_l0 in arrays
    T_ads = []
    s_l0s = []
    for X_H2 in X_H2s:
        for phi in phis:
            flame = PremixedFlame(phi, X_H2)
            T_ad = flame.T_ad # Adiabatic temperature [K]
            s_l0 = flame.s_l0 # Laminar flame speed [m/s]
            T_ads.append(T_ad)
            s_l0s.append(s_l0)
    T_ads = np.array(T_ads)
    s_l0s = np.array(s_l0s)
    

#%% PLOTS
    
    # Create empty lists of figures, axes
    figs = []
    axs = []
    fig_names = []
    
    
    # Plot adiabatic temperature versus equivalence ratio
    fig_name = 'T_ad_vs_phi'
    fig, ax = plt.subplots()
    for i, X_H2 in enumerate(X_H2s):
        x_1 = phis
        y_1 = T_ads[3*i:3*i+3]
        ax.scatter(x_1, y_1, marker='o', s=50, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=f'{X_H2}', zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], lw=2, zorder=1)
        
    # Set the labels and legend
    xlabel = '$\phi$ [-]'
    ylabel = '$T_{\mathrm{ad}}$ [K]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title='$X_{\mathrm{H}_2, f}$ [%]', loc='best', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)
    
    
    # Plot unstretched laminar flame speed versus equivalence ratio
    fig_name = 's_l0_vs_phi'
    fig, ax = plt.subplots()
    for i, X_H2 in enumerate(X_H2s):
        x_1 = phis
        y_1 = s_l0s[3*i:3*i+3]
        ax.scatter(x_1, y_1, marker='o', s=50, color=colors[2*i+1], edgecolors=colors[2*i], linewidths=None, label=f'{X_H2}', zorder=2)
        ax.plot(x_1, y_1, c=colors[2*i+1], lw=2, zorder=1)
         
    # Set the labels and legend
    xlabel = '$\phi$ [-]'
    ylabel = '$s_{l, 0}$ [m/s]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title='$X_{\mathrm{H}_2, f}$ [%]', loc='best', handletextpad = 0.1, columnspacing=0.1, framealpha = 0)
    fig.tight_layout()
    figs.append(fig)
    axs.append(ax)
    fig_names.append(fig_name)

        
#%% SAVE THE FIGURES
    
    # Specify the figures directory
    sep = os.path.sep # Gets the OS path separator.
    sub_dir = script_name.removesuffix('.py')
    fig_dir = script_dir + sep + sub_dir
    # Create the figures directory if it does not exist
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
        
    # Apply tight_layout to each figure and save it to the specified folder
    for i, fig in enumerate(figs):
        fig_nr = i + 1
        # Specify which figures to save
        if 1 <= fig_nr <= len(figs):
            fig.tight_layout()
            fig_name_suffix = '.pdf'
            fig_path = fig_dir + sep + fig_names[i] + fig_name_suffix
            fig.savefig(fig_path, format='pdf', dpi=500, bbox_inches='tight')

