# Combustion Assignement

## 2025-2026 Turbomachinery Course

## Preliminary preparations

Before starting the assignement, you need to install a Python distribution on your machine. You can do so by following the instructions given [here](./installation_instructions.md). If you do not want or are not able to install a Python distribution on your PC, you may run the code directly by using an online IDE, such as [Jupyter](https://jupyter.org/), [Online Python](https://www.online-python.com/) or [PyTogether](https://pytogether.org/). 

## Introduction

A key component of gas turbines is the combustor. In the combustor, the fuel and oxidiser (air) will be mixed and combusted to increase the enthalpy of the working fluid. Therefore, understanding the physics which govern combustion is essential for gas turbine design. This assignement aims to enhance your knowledge and understanding of combustion through the use of a 1-D flame simulation. The base code is provided with the questions. It contains a 1-D premixed flame simulation which computes the main flame parameters, such as the adiabatic flame temperature $T_{\mathrm{ad}}$ and the unstretched laminar flame speed $s_{l,0}$. 

In order to run the simulation, you will have to specify the following unburnt mixture parameters: the equivalence ratio $\phi$ and hydrogen volume fraction in the fuel $X_{\mathrm{H_2},f}$. The unburnt mixture temperature and pressure are set to default values $T_u=20$ °C and $p_u=1$ atm. The equivalence ratio is an important parameter in premixed combustion, which reflects how close a certain mixture is to stoichiometry. It is defined as

$$
\phi = \frac{n_f/n_o}{(n_f/n_o)_\mathrm{st}} = \frac{m_f/m_o}{(m_f/m_o)_\mathrm{st}},
$$

where $n$ and $m$ represent the quantity (mole) and mass, respectively, and the subscripts $f$, $o$ and $\mathrm{st}$ stand for fuel, oxidiser and stoichiometry, respectively. Thus, $\phi=1$ means that the mixture is at stoichiometric conditions, while $\phi<1$ and $\phi>1$ indicate fuel-lean and fuel-rich mixtures, respectively.

Throughout this assignment, it will be assumed, for simplicity, that natural gas is composed of pure methane ($\mathrm{CH_4}$). Hence, $X_{\mathrm{H_2},f}=0$ % means that the fuel is pure methane, while $X_{\mathrm{H_2},f}=100$ % means that it is pure hydrogen.

Air will always be used as an oxidiser, and it will be assumed that its composition (in volume fraction) is 79 % nitrogen and 21 % oxygen.

## Question 1

Depending on your student number, use the following equivalence ratios:

- Even number: $\phi_1=0.8$, $\phi_2=1$ and $\phi_3=1.15$

- Odd number: $\phi_1=0.9$, $\phi_2=1$ and $\phi_3=1.3$

Run the 1-D flame simulation for all three equivalence ratios, both using pure natural gas and pure hydrogen as a fuel (6 cases). 

Plot the adiabatic temperatures obtained for both fuels against the equivalence ratio. What do you observe? Comment on the effect of the equivalence ratio and the fuel composition. What does this mean in the context of a gas turbine combustor? What are  possible challenges caused by replacing natural gas by hydrogen in the fuel?

Plot the laminar flame speeds obtained for both fuels against the equivalence ratio. What do you observe? Comment on the effect of the equivalence ratio and the fuel composition. What does this mean in the context of a gas turbine combustor? What are possible challenges caused by replacing natural gas by hydrogen in the fuel?

## Question 2

Real flames encoutered in engineering applications are rarely 1-dimensional. The laminar flame speed obtained in question 1 is therefore a rather theoretical model, which does not account for much of the actual physics encountered in real flames. 

Hence, let us now consider the case of a Bunsen burner. Assume that the Bunsen burner consists of a simple pipe of inner diameter $d$, with the flame anchored at the burner rim. Also assume that fuel and oxidiser are homogeneously mixed upstream of the burner inlet, and that the flow inside the burner is turbulent and fully developed.

The Karlovitz number is a non-dimensional number which relates the laminar flame timescale to the Kolmogorov timescale of turbulent flow. It is defined as

$$
\mathrm{Ka} = \frac{\tau_l}{\tau_\eta} = \left( \frac{u_\eta}{s_{l,0}} \right)^2,
$$

where $u_\eta$ is the Kolmogorov velocity.

## Authors

- Luuk Altenburg
- Rafael Pichler

## License

This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Project status

Under development and in use.