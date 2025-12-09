# HELIOS


## Description
Repository for the code linked to the experimental part of the HELIOS project; the investigation of hydrogen's fundamental combustion characteristics (flame speed, flashback, flow interaction, blowoff). This code is used for the purposes of setup design, data post-processing, modeling and the like.
The subdirectories and the files they contains are described below.

- `control_panel_verification.py`: Short script that verifies the volumetric flow rates computed by the LabView panel using cantera.

### Data visualisation
This directory contains all scripts used for data visualisation, such as plotting fields and profiles
- `uffb_ref_data`: Subdirectory containing LDA data (stored in `.mat` files) collected by Mark Tummers, necessary to validate the unreactive flow experiments on the Bunsen burner.
- `flashback_maps.py`: This module produces flashback maps according to the control panel data as well as data entered manually. It imports the `piv_data_postprocessing` package.
- `flashback_maps_triple.py`: This module produces flashback maps showing all 3 transition points, not only the downstream-bend transition. It imports the `piv_data_postprocessing` package.
- `flow_fields_bunsen.py`: This module allows for the visualisation of flow fields resulting from PIV experiments on the Bunsen burner. It imports the `piv_data_postprocessing` package.
- `flow_fields_flamesheet.py`: This module allows for the visualisation of flow fields resulting from PIV experiments on the Flamesheet burner. It imports the `piv_data_postprocessing` package.
- `temperature_graphs.py`: This module generates different kinds of plots allowing to visualise instataneous and steady-state temperarture measurements. It imports the `piv_data_postprocessing` package.
- `unreactive_flow_fields_bunsen.py`: This module allows for the visualisation of unreactive flow fields resulting from PIV experiments on the Bunsen burner. It imports the `piv_data_postprocessing` package.
- `unreactive_flow_fields_flamesheet.py`: This module allows for the visualisation of unreactive flow fields resulting from PIV experiments on the Flamesheet burner. It imports the `piv_data_postprocessing` package.

### Documentation
This directory contains all documentation files other than the readme file.
- `git_SSH_permission_denied_error_fix.md`: Describes how to solve the "Permission denied (public key)" error when trying to set up an SSH key for Git.
- `common_conda_commands.md`: Lists common Conda commands with a short description of what they do.

### Flame front detection
This package (directory) contains a series of modules that detect the flame front based on raw images in the `.tif` format.
- `__init__.py`: This module defines the `flame_front_detection` subdirectory as a package and takes care of all the relative imports such that the whole package can be imported from any script outside this directory.
- `aux_functions.py`: Contains auxiliary functions necessary to run `get_contour_procedures.py`.
- `contour.py`: Classes `Frame` and `ContourData` are defined is this module. It imports `contour_properties.py` and `get_contour_procedures.py`.
- `contour_properties.py`: Contains functions which are used to calculate the contour's properties. These are necessary in order to run `contour.py`.
- `flame.py`: Class `Flame` is defined in this module, which basically serves to create `flame` objects corresponding to real flames with their key properties. This module imports `premixed_flame_properties.py` and `contour.py`.
- `flames_info.py`: Base module in which the information corresponding to each flame analysed is entered, such that a corresponding object can be created, and the flame front can be detected using the algorithm. It also allows reading and saving flame objects from and to pickle files, such that the algorithm must only be run once for a given flame. This module imports `flame.py` and `contour_properties.py`. 
- `get_contour_procedures.py`: Contains the functions which are used to detect the flame front. The flame front is detected by running a bilateral (space and intensity) gaussian filter over the raw image, before the filtered image is converted into a binary map thanks to a histogram method based on pixel intensity. In this method, the image is divided into two regions at the minimum located between the two histogram peaks. This module imports `aux_functions.py`.
- `premixed_flame_properties.py`: Class `PremixedFlame` is defined is this module, which evaluates the transport, thermodynamic and chemical properties of the flame object that is being passed to it.

### PIV data post-processing
This package (directory) contains a series of modules which allow the PIV data to be post-processed, when it is imported from another script. Post-processing includes extracting the time-dependent velocity fields from the data files (usually `.csv`), and computing their statistics.
- `__init__.py`: This module defines the `piv_data_postprocessing` subdirectory as a package and takes care of all the relative imports such that the whole package can be imported from any script outside this directory.
- `burner.py`: Class `Burner` is defined in this module. This class ajusts the data according to the burner parameters (such as defining the reference length and shifting the coordinate system). 
- `cp_data.py`: Class `ControlPanelData` is defined in this module. This class reads the control data file corresponding to one experimental run and retrieves the relevant data. This module imports `velocity_field.py`.
- `piv_data.py`: Class `PIVData` is defined in this module. This class reads all PIV data files corresponding to one experiment by calling `VelocityFields`, and builds 3D arrays for the velocity components and magnitude. It also retrieves the x and y coordinate axes, and computes statistical variables over the time-dependent fields. This module imports `burner.py`, `cp_data.py` and `velocity_field.py`.
- `tc_data.py`: Class `ThermocoupleData` is defined in this module. This class reads the thermocouple data file corresponding to one experimental run and retrieves the relevant data.  This module imports `velocity_field.py`.
- `velocity_field.py`: Class `VelocityField` is defined in this module. This class reads one velocity field obtained through PIV from the data file it has been saved to, and extracts its components and magnitude to (numpy) arrays.

### Setup design
This directory contains all modules linked to the design of the model FlameSheet burner setup.
- `flowswitch.py`: Performs all calculations linked to the installation of the flowswitch upstream of the preheater. It most importantly returns the acceptable pipe diameters on which the flowswitch must be mountain for it to function correctly with the preheater.
- `preheater.py`: Performs all calculations linked to the design of the preheating system. It most importantly returns the required heating power in order to attin a given burner entry temperature.

## Installation
The code in this repository has been written using Python and its versions are being controlled using Git. Therefore, to be able to use and develop the code, you need both to install a Python distribution and Git. 

### Installing a Python distribution
To use Python on your computer, you will need to create a Python environment into which you will install your desired Python distribution (i.e. version), as well as any packages necessary to run the code. The easiest way to do this is by installing an environment and package manager such as Conda. There are two versions of Conda: Miniconda and Anaconda. Miniconda is only the environment and package manager, while Anaconda already includes a certain number of common packages. If you are not using that many packages, such as in this code, it is preferred to install Miniconda and then only install the few necessary packages. This will keep the instalation as simple and as small as possible. The steps below summarise how to install a Python distribution using Miniconda.
1. Install Miniconda by following the steps on: <https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>. It is recommended to install Miniconda only for your local user account, unless several people are expected to collaborate on the same project using the same computer.
2. Verify that conda has been installed correctly and check that it is up to date by referring to: <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-conda.html>.
3. Append (or add) any necessary channels to your conda configuration by referring to: <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html>. For this project, the `cantera` and `conda-forge` channels are necessary in addition to the `defaults`. It is recommended that you prioritise the `defaults` followed by the `cantera` and then the `conda-forge`, as the `conda-forge` usually has the newest releases which are still unstable.
4. Check which Python distribution you would like to install by referring to: <https://www.python.org/downloads/>.
5. Create an environment with the desired Python distribution by referring to: <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>.
6. Install any necessary packages into the environment by referring to: <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html>. Most IDEs, such as Spyder can simply be installed like any other package. Available packages are listed on: <https://anaconda.org/>.
 
### Installing and using Git
1. Download and install Git from: <https://git-scm.com/downloads>. You can additionnally refer to: <https://docs.gitlab.com/ee/topics/git/how_to_install_git/index.html>
2. Create an SSH key and add it to GitLab by referring to: <https://docs.gitlab.com/ee/user/ssh.html>. This is not required for Git to function, but it will save you from logging into GitLab each time you want to push to GitLab. If you get the `Permission denied (public key)` error when trying to set up the SSH key, refer to [this](./documentation/git_SSH_permission_denied_error_fix.md).
3. Clone the project into a local repository (folder) from which you want to work from by referring to: <https://docs.gitlab.com/ee/tutorials/make_first_git_commit/>.
4. Create your own developement branch and name it 'NetID-dev' by further referring to <https://docs.gitlab.com/ee/tutorials/make_first_git_commit/>.
5. Make changes, commit them and push them onto GitLab by further referring to <https://docs.gitlab.com/ee/tutorials/make_first_git_commit/>.

### Channels
The channels necessary to install the packages listed further are the following. They are listed in the recommended order of descending priority.
- `defaults`
- `cantera`
- `conda-forge`

### Packages
The packages necessary to run this code are listed below.
- NumPy: `numpy`
- SciPy: `scipy`
- MatPlotLib: `matplotlib`
- Pandas: `pandas`
- Cantera: `cantera`
- OpenCV: `opencv`
- Progress Bar 2: `progressbar2`
- imutils: `imutils`
- tqdm: `tqdm`
- scikit-learn: `scikit-learn`

## Resources
- Conda documentation: <https://docs.conda.io/projects/conda/en/latest/index.html>
- Common Conda commands: [here](./documentation/common_conda_commands.md)
- Anaconda package repository: <https://anaconda.org/>
- Git documentation: <https://git-scm.com/doc>
- Git everyday commands: <https://git-scm.com/docs/giteveryday>
- GitLab documentation: <https://docs.gitlab.com/>
- Python documentation: <https://www.python.org/doc/>
- Python f-string cheatsheet: <https://fstring.help/cheat/>
- NumPy documentation: <https://numpy.org/doc/>
- SciPy documentation: <https://docs.scipy.org/doc//scipy/index.html>
- MatPlotLib documentation: <https://matplotlib.org/stable/index.html>
- Pandas documentation: <https://pandas.pydata.org/docs/>
- Cantera documentation: <https://cantera.org/documentation/>

## Troubleshooting

### Conda
- Weird error occurs after updating Conda or some package(s): Most likely this is caused by a CONDA_TRASH file that Conda has created in your Python environment while updating. To solve the issue simply navigate to your Python environment (myenv) in Windows explorer (usually the path will be `C:\Users\rpichler\AppData\Local\miniconda3\envs\myenv`) and delete the CONDA_TRASH file (finishes by `.conda_trash`).

## Authors
- Luuk Altenburg
- Rafael Pichler

## License
This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Project status
Under development and in use.