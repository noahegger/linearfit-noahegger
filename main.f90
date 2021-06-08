! Program: nuclear_energies
! By: Noah Egger
!------------------------------------------------------------------------------
! This program reads experimental data provided and calculates the binding
! energies for isotopes given proton number up to 118. The program then computes
! the most stable isotope for a given Z as well as the neutorn dripline for the
! range of isotopes.
!
!------------------------------------------------------------------------------
program nuclear_energies
use types
use read_write, only : read_exp_data, write_predictions, write_nuclear_lines
use nuclear_model, only : find_best_parameters !, most_stable_n
implicit none

integer, allocatable :: n_protons(:), n_neutrons(:)
real(dp), allocatable :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)

! Reads the input data
call read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
! Finds best parameters for binding energy expression
call find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
! Writes predictions of binding energy as well as the associated error
call write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)

! Writes data corresponding to position of stable isotopes and 
! position of the neutron dripline

call write_nuclear_lines(c_parameters, n_protons, n_neutrons)
!call most_stable_n(c_parameters, n_stable)

end program nuclear_energies