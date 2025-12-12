# Installation instructions

## Installation

The code in this repository has been written using Python. Therefore, to be able to use the code, you need to install a Python distribution. 

### Installing a Python distribution

To use Python on your computer, you will need to create a Python environment into which you will install your desired Python distribution (i.e. version), as well as any packages necessary to run the code. The easiest way to do this is by installing an environment and package manager such as Miniforge. The steps below summarise how to install a Python distribution using Miniforge.

1. [Download](https://conda-forge.org/download/) and install Miniforge. It is recommended to install Miniforge only for your local user account, unless several people are expected to collaborate on the same project using the same computer. When installing Miniforge, make sure to **add it to your PATH environment variable and register it as your default Python 3.12**.
2. Open a Miniforge command prompt. You should be able to find it using the search bar.
3. Create a new environment with Python 3.13. To do so, run the following command: `mamba create --name py313 python=3.13`. Choose an environment name which makes sense to you. Here `py313` is used as an example.
4. Close the Minforge command prompt and reopen it.
5. Activate the new environment by running `mamba activate py313`. If it does not work, run `conda activate py313`.
6. Install the necessary packages into the environment by running `mamba install package-name`. For this code, you will need `numpy`, `cantera` and `matplotlib`.
7. Install an IDE (Integrated Development Environment). You can install it just like any other package. Here, Spyder is recommended (`mamba install spyder`).
8. Open Spyder or any other IDE you installed, by entering the IDE name (`spyder`) in the Miniforge command prompt. 
9. From your IDE, you can open and run the code similarly as in MATLAB for example.

To re-open your IDE after installation is complete, simply open a Miniforge command prompt, activate your Python environment and enter the IDE name. The active environment is indicated between parentheses in the Miniforge command prompt.

### Packages

The packages necessary to run this code are listed below.

- NumPy: `numpy`
- MatPlotLib: `matplotlib`
- Cantera: `cantera`

## Resources

- [Mamba documentation](https://mamba.readthedocs.io/en/latest/index.html)
- [Conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) (most mamba commands are the same as conda ones and can just be run by replacing `conda` with `mamba`)
- [Python documentation](https://www.python.org/doc/)
- [Python f-string cheatsheet](https://fstring.help/cheat/)
- [NumPy documentation](https://numpy.org/doc/stable/)
- [MatPlotLib documentation](https://matplotlib.org/stable/index.html)
- [Cantera documentation](https://www.cantera.org/)