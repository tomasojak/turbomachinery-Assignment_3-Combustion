# Installation instructions

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