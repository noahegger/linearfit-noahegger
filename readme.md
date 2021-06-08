# Program Goal

This project contains multiple components. The first part of the project consists of finding the best fit parameters from the linear semi-empirical mass formula (henceforth SEMF). In order to do this, one feeds the experimental data `EXPERIMENT_AMES2016.dat ` into the program. From there, the program constructs the relevant matrices and vectors and through the use of LU decomposition, finds the best fit parameters for the SEMF as well as the corresponding uncertainties in binding energy. The program then prints to the screen the best fit parameters and corresponding uncertainties, as well as writes a file called `results_basic.dat ` containing the experimental values, experimental uncertainties, theoretical values, and theoretical uncertainties. The second part of the project utilizes the theoretical model (SEMF with found best fit parameters) in order to find the most stable isotope; that is, the value of N for which the binding energy per nucleon is the lowest (the valley of stability). Additionally, the program finds for every value of Z between 1 and 118 the position of the neutron drip-line; that is, the largest value of N for which the separation energy is positive. Moreover, the program contained in the `plots_analysis.ipynb ` file plots the deviation of the theoretical binding energy from the experimental binding energy. Within this program, there exists a discussion of the pattern seen from this difference. In addition, there are two more plots that show the experimental uncertainty as a function of Z and the theoretical uncertainties as a function of Z. There exists a discussion comparing these two plots as well.

Finally, the program makes an attempt at improving the SEMF via addition of another parameter. At the end of the `plots_analysis.ipynb ` file, the reduced chi-squared for the data set is computed.

# Directions for Usage

Navigate to the relevant directory through terminal. Make sure all `.f90` files are in the same directory as well as the `EXPERIMENT_AME2016.ipynb ` file. Type "make" into the terminal as this will compile the files and create an executable called `nuclear_energies`. In terminal type "./nuclear_energies" to run the program. The program will ask for a data file. Type "EXPERIMENT_AME2016.dat". The program will then display the best fit parameters for the SEMF as well as the additional ansatz parameter as well as the corresponding uncertainties. The program will write into a file called results_basic.dat' containing the experimental values, experimental uncertainties, theoretical values, and theoretical uncertainties. Additionally, the program will write to a file called `results_advanced.dat ` containing the positions of valley of stability and the neutron drip-line. You can visualize all this data by running the jupyter notebook `plots_analysis.ipynb` file. Make sure the `plots_analysis.ipynb ` file is contained within the same directory as all the `.f90 ` files and the `EXPERIMENT_AME2016.dat ` file.

# Contents

This program contains the files `linear_algebra.f90`, `read_write.f90`, `types.f90`, `main.f90`, `nuclear_model.f90`, as well as `plots_analysis.ipynb`, this `readme.md`, the `EXPERIMENT_AME2016.dat `, and the `makefile`. 

`linear_algebra.f90` This module computes the inverse of the alpha matrix constructed within the nuclear model module, then solves the system of equations to obtain the parameters needed for the binding energy calculations.

`nuclear_model.f90` Constructs the relevant matrices and vectors to be used in the `linear_algebra.f90` module. Also, constructs the parameters associated with each term in the SEMF. Finally, also locates the positions of the most stable isotope and constructs the neutron driplines.

`read_write.f90` takes inputs from user and writes data to results files. 

`types.f90` contains types of inputs to be used (real double precision, pi ~ 3.14159..., etc)

`main.f90` calls to read/write module to execute the program. 

`plots_analysis.ipynb` can be opened with jupyter notebook. Contains plots comparing deviation in theoretical and experimental binding energies as a function of N and Z. Also contains plots of experimental and theoretical uncertainty as a function of Z. Also constructs a plot indicating the position of the neutron drip line and most stable isotopes (how Z depends on N). Moreover, a reduced chi-squared for the theoretical model (SEMF) is calculated. 

`makefile` allows user to type "make" into terminal to compile the program.

`readme.md` this file. Contains instructions for the use of this program. 

`EXPERIMENT_AME2016.dat ` Contains experimental data. Most importantly, the associated binding energies for various isotopes up to Z = 118.





