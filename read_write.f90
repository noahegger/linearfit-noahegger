!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module reads the input file containing the experimental data.
!! Additionally, the module writes the calculated BEs and corresponding
!! errors to a file 'results.dat'. Finally, the module writes the neutron
!! drip line and isotope data to a file 'advanced_results.dat'.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_exp_data
!! write_predictions
!! write_nuclear_lines
!!----------------------------------------------------------------------
module read_write

use types
use nuclear_model, only : semi_empirical_mass, semi_empirical_error, most_stable_n, neutron_drip_line

implicit none

private
public :: read_exp_data, write_predictions, write_nuclear_lines

contains

!-----------------------------------------------------------------------
!! Subroutine: read_exp_data
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine reads the input file specified by user, 
!! 'EXPERIMENT_AME2016.dat'.
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
subroutine read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
    implicit none
    integer, intent(out), allocatable :: n_protons(:), n_neutrons(:)
    real(dp), intent(out), allocatable :: exp_values(:), uncertainties(:)

    character(len=128) :: filename, string_trash
    logical :: file_exists
    integer :: file_unit, number_data_points, integer_trash1, integer_trash2, i
    real(dp) :: real_trash


    ! use the print statement to display a message explaining what the program
    ! does

    print *, 'This program will read information given by a file.'
    print *, 'please provide the file name with the experimental data'
    print *, 'The program will take the experimental data and perform'
    print *, 'a semi-empirical calculation of the BEs for a number of'
    print *, 'isotopes. Additionally, it will calculate the neutron'
    print *, 'drip lines and most stable isotopes.'
    print *, 'Please provide the file name with the experimental data:'
    read(*, '(a)') filename 

    ! when trying to allocate an array that was passed as an argument, it's
    ! always good idea to deallocate them if for whatever reason they're
    ! already allocated
    
    if(allocated(n_protons)) deallocate(n_protons)
    if(allocated(n_neutrons)) deallocate(n_neutrons)
    if(allocated(exp_values)) deallocate(exp_values)
    if(allocated(uncertainties)) deallocate(uncertainties)

    ! when trying to open a file provided by the user it's good practice to
    ! check if the file exists in the current directory
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
        ! Open the file and read the data. (don't forget to close the file)
        open(unit = 1, file=filename)
        read(1,*) number_data_points

        allocate(n_protons(1:number_data_points))
        allocate(n_neutrons(1:number_data_points))
        allocate(exp_values(1:number_data_points))
        allocate(uncertainties(1:number_data_points))
    do i=1,2
        read(1,*)
    enddo
    ! Loop for reading each row
        do i=1,number_data_points
    ! Data in colums we do not need is stored in a variable and never used. 
        read(1,*) integer_trash1, string_trash, integer_trash2, n_neutrons(i), & 
        n_protons(i), exp_values(i), real_trash, uncertainties(i)
        enddo 
        close(1)
    else
        ! If file does not exist, print error and abort program. 
        print*, 'The file named ', trim(filename),' could not be found.'
        stop
    endif
end subroutine read_exp_data

!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine records in a file 6 columns containing the number of
!! protons, the number of neutrons, the experimental binding energy,
!! the experimental error, the theoretical binding energy, and the 
!! theoretical error.
!!----------------------------------------------------------------------
!! Input:
!!
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! covariance       real        Array containing the elements of the covariance matrix
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    character(len=*), parameter :: file_name = 'results.dat'
    integer :: unit1, i, length
    real(dp) :: theoretical_error, theoretical_BE
    length = size(n_neutrons)

    open(newunit = unit1, file = file_name)
    write(unit1,28) 'Protons', 'Neutrons', 'Experimental BE', 'Experimental Error', 'Theoretical BE', 'Theoretical Error'
    28 format(5x, a, 5x, a, 5x, a, 14x, a, 12x, a, 19x, a)
    do i = 1, length

        theoretical_BE = semi_empirical_mass(c_parameters, n_protons(i), n_neutrons(i))
        theoretical_error = semi_empirical_error(covariance, n_protons(i), n_neutrons(i))
        write(unit1,*) n_protons(i), n_neutrons(i), exp_values(i), uncertainties(i), theoretical_BE, theoretical_error

    enddo
    close(unit1)

    print *, 'Theoretical binding energies and their errors were written in ', file_name

end subroutine write_predictions

!--------------------------------------------------------------------------------
! Your advanced subroutine goes here.
! Remember to document your new subroutine

! In order to work with the jupyter notebook. Write down the results in a file named 
! 'results_advanced.dat'.
!
! The file should have 3 columns. The number of protons, the position (value
! of N) for the stable isotopes, and the position (value of N) for the neutron
! drip-line
!
! Don't forget to put a header in the file.

!-----------------------------------------------------------------------
!! Subroutine: write_nuclear_lines
!-----------------------------------------------------------------------
!! By: Noah Egger 
!!
!! This subroutine writes to a file names 'results_advanced.dat'. The file
!! contains 3 columns: the number of protons, the position (value of N) for
!! stable isotopes, and the position (value of N) for the neutron drip line.
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_nuclear_lines(c_parameters, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    character(len=*), parameter :: file_name = 'results_advanced.dat'
    integer :: unit2, i, p_max
    integer, allocatable :: n_stable(:), n_drip(:), z_values(:)
    
    allocate(n_stable(1:118))
    allocate(n_drip(1:118))
    allocate(z_values(1:118))


    ! Retrieve array containing neutron number associated with stable isotopes. 

    call most_stable_n(c_parameters, n_stable)
    p_max = size(n_stable)

    ! Call subroutine to calculate largest neutron count, per z value, for
    ! which the separation energy is positive. 

    call neutron_drip_line(c_parameters, p_max, n_drip)
    
    ! Do loop to fill array representing all possible proton numbers in the isotopes. 

    do i=1,p_max
        z_values(i) = i
    enddo

    ! Open file to write to. 

    open(newunit=unit2,file=file_name)

    write(unit2,2) 'Proton Number: Z', 'Neutron Number: N', 'Drip Line: N'
    2 format(4x,a,4x,a,6x,a)

    ! Loop to write Z, N, and drip line.

    do i=1,p_max      
    write(unit2,*) z_values(i), n_stable(i), n_drip(i)

    enddo
    close(unit2)
    print *, 'Stable isotopes and neutron dripline were written to ', file_name
end subroutine write_nuclear_lines
    
end module read_write