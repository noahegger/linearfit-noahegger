!-----------------------------------------------------------------------
!Module: nuclear_model
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module addresses the relevant physics. The linear terms within
!! the semi-empirical mass formula (SEMF henceforth) are calculated in
!! individual functions. The alpha matrix is constructed using the product
!! of neighboring linear terms divided by the square of the uncertainty for
!! a given iteration, summed over all iterations "k". Additionally, we
!! construct the beta vector which consists of the parameters indexed
!! by "i", multiplied by the experimental BE indexed by "k" and divided by 
!! the uncertainty, also indexed by "k", and summed over all k for each i
!! within the beta vector. The matrix and vector constructed here are solved
!! within the linear_algebra module in order to construct the parameters
!! associated with each term in the SEMF. Finally, we also locate
!! the positions of the most stable isotope and construct the neutron driplines.
!!
!!
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! find_best_parameters
!! construct_alpha_beta
!! calculate_linear_termns
!! print_best_parameters
!! most_stable_n
!! neutron_drip_line
!!----------------------------------------------------------------------
!! Included functions:
!!
!! volume_term
!! surface_term
!! asymmetry_term
!! coulomb_term
!! pairing_term
!! my_extra_term
!! semi_empirical_mass
!! semi_empirical_error
!! delta
!! find_largest_n
!!----------------------------------------------------------------------
module nuclear_model
use types
use linear_algebra, only : solve_linear_system
implicit none

private

public :: find_best_parameters, semi_empirical_mass, semi_empirical_error, most_stable_n, neutron_drip_line
contains


!-----------------------------------------------------------------------
!! Subroutine: find_best_parameters
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Describe what the subroutine does
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! c_parameters     real        Array containing the semi-empirical mass formula parameters
!! covariance       real        Array containing the covariance matrix of the parameters
!-----------------------------------------------------------------------
subroutine find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out), allocatable ::  c_parameters(:), covariance(:,:)

    ! Number of terms in binding energy calculation (6 to do EC).

    integer, parameter :: n_parameters = 6

    ! Initialize alpha matrix and beta vector. Must have size n parameters.

    real(dp) :: alpha(1:n_parameters,1:n_parameters), beta(1:n_parameters)

    ! c_parameters and covariance were passed as arguments and need to be
    ! allocated. Deallocate them if they're allocated and then allocate them
    ! with the correct size

    if (allocated(c_parameters)) deallocate(c_parameters)
    if (allocated(covariance)) deallocate(covariance)

    ! Allocate arrays with correct size. Covariance is used in solve_linear_system.

    allocate(c_parameters(n_parameters))
    allocate(covariance(n_parameters, n_parameters))
    
    call construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)

    ! The subroutine below (defined in the linear_algebra module) should solve
    ! the matrix equation in the README and return the c_parameters and the
    ! inverse of alpha (the covariance matrix)

    call solve_linear_system(alpha,beta,c_parameters,covariance)

    ! Now just print the parameters (with it's uncertainties) to screen

    call print_best_parameters(c_parameters,covariance)
end subroutine find_best_parameters

!-----------------------------------------------------------------------
!! Subroutine: construct_alpha_beta
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine constructs the alpha matrix and beta vector in order
!! to form the linear system to be solved.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! alpha            real        Array containing the alpha matrix
!! beta             real        Array containing the beta vector
!-----------------------------------------------------------------------
subroutine construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out) :: alpha(:,:), beta(:)

    integer :: n_data, n_parameters, alpha_shape(1:2), i, j, k, beta_size
    real(dp):: linear_terms(1:size(beta))

    ! Check if the alpha array is a square matrix
    ! Also check that beta has the same number of elements as alpha has rows (or columns)

    alpha_shape = shape(alpha)
    beta_size = size(beta)

    if (alpha_shape(2) /= beta_size) then
        print*, 'The size of beta does not equal the column size of alpha.'
        stop
    endif

    if (alpha_shape(1) /= alpha_shape(2)) then
        print*, 'The alpha matrix is not square.'
        stop
    endif

    ! Number of data points in file.

    n_data = size(uncertainties)

    ! Number of terms (parameters in SEMF)

    n_parameters = alpha_shape(1)

    alpha = 0._dp
    beta = 0._dp

    ! Constructing alpha matrix.

    ! "i" will represent rows
    ! "j" will represent columns
    ! "k" will sum through all data points for any given i,j
    ! The subroutine below will return the f_\alpha(Z_i,N_i) terms defined in
    ! the README file

    ! The i'th and j'th terms of f_\alpha(Z_i,N_i) are multiplied and divided by 
    ! the experimental uncertainty for a specific k, all summed over k to form
    ! the alpha matrix

    do k = 1, n_data
        do i= 1, n_parameters
            do j = 1, n_parameters
                call calculate_linear_terms(n_protons(k), n_neutrons(k), linear_terms)
                alpha(i,j) = alpha(i,j) + (linear_terms(i)*linear_terms(j))/((uncertainties(k))**2)
            enddo
        enddo
    enddo

    ! Constructing beta vector. 
    ! The i'th term of f_\alpha(Z_i,N_i) is multiplied by the experimental binding
    ! energies indexed with k, then divided by experimental uncertainty for a 
    ! a specific k, summed over k to form the beta vector
!  

    do i = 1, n_parameters
        call calculate_linear_terms(n_protons(i), n_neutrons(i), linear_terms)
        do k = 1, n_data
            call calculate_linear_terms(n_protons(k), n_neutrons(k), linear_terms)
            beta(i) = beta(i) + linear_terms(i)*exp_values(k)/(uncertainties(k)**2)
        enddo
    enddo

            

end subroutine construct_alpha_beta

!-----------------------------------------------------------------------
!! Subroutine: calculate_linear_terms
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine calls the individual terms within the SEMF, each 
!! calculated separately in separate routines. Here, "terms" refers to
!! the vlaues in front of the individual parameters c_vol, c_coul, etc.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! N                integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! linear_terms        real        Array containing the linear terms in the semi-empirical mass formula
!-----------------------------------------------------------------------
subroutine calculate_linear_terms(Z, N, linear_terms)
    implicit none
    integer, intent(in) :: Z, N
    real(dp), intent(out) :: linear_terms(:)

    ! We could write down all the formulas for each term here. However, in
    ! order to keep the code readable and easy to understand  we'll  separate
    ! them into different functions

    linear_terms(1) = volume_term(Z,N)
    linear_terms(2) = surface_term(Z,N)
    linear_terms(3) = asymmetry_term(Z,N)
    linear_terms(4) = coulomb_term(Z,N)
    linear_terms(5) = pairing_term(Z,N)
    linear_terms(6) = my_extra_term(Z,N)

end subroutine calculate_linear_terms

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Calculates the volume term in the SEMF to be used in the linear_terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function volume_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real(Z+N, kind=dp)

end function volume_term

!-----------------------------------------------------------------------
!! function: surface_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Calculates the surface term in the SEMF to be used in the linear_terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        surface term
!-----------------------------------------------------------------------
real(dp) function surface_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real(Z+N, kind=dp)**(2._dp/3._dp)

end function surface_term

!-----------------------------------------------------------------------
!! function: asymmetry_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Calculates the asymmetry term in the SEMF to be used in the linear_terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        asymmetry term
!-----------------------------------------------------------------------
real(dp) function asymmetry_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    
    r = (real(N-Z, kind=dp)**2._dp)/real(N+Z,kind=dp)

end function asymmetry_term

!-----------------------------------------------------------------------
!! function: coulomb_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Calculates the coulomb term in the SEMF to be used in the linear_terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        coulomb term
!-----------------------------------------------------------------------
real(dp) function coulomb_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    
    r = real(Z*(Z - 1),kind=dp)/(real(Z+N,kind=dp)**(1._dp/3._dp))

end function coulomb_term

!-----------------------------------------------------------------------
!! function: pairing_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!! 
!! Calculates the pairing term in the SEMF to be used in the linear_terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        pairing term
!-----------------------------------------------------------------------
real(dp) function pairing_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    
    r = ((real(Z+N, kind=dp)**(-3._dp/4._dp)))*delta(Z, N)

end function pairing_term

!-----------------------------------------------------------------------
!! function: my_extra_term
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Calculates the additional ansatz term in the SEMF to be used in the 
!! linear_terms array. At this point, a total shot in the dark.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function my_extra_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    !r = real((Z**2), kind=dp)/(real((Z+N),kind=dp))
    r = real((N - Z)**2, kind=dp)/real((Z+N)**(4._dp/3._dp),kind=dp)
    

end function my_extra_term

!-----------------------------------------------------------------------
!! function: delta
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Returns -1, 0, or 1 depending on whether if Z and N are both even,
!! A is odd, or Z and N are both odd, respectively.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        result of 1, -1, or 0
!-----------------------------------------------------------------------
real(dp) function delta(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    ! If number of protons and neutron are both even, return a 1 for delta. 
    ! If number of protons and neutrons are both odd, return a -1 for delta.
    ! If A = N+Z is odd, return a 0 for delta.

    if ( modulo(Z,2) == 0 .and. modulo(N,2) == 0 ) then
        r = 1
    else if (modulo(Z,2) == 1 .and. modulo(N,2) == 1) then
        r = -1
    else
        r = 0
    endif


end function delta

!-----------------------------------------------------------------------
!! Subroutine: print_best_parameters
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine prints the best fit parameters within the SEMF and the 
!! associated uncertainty.
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! covariance       real        Array containing covariance matrix
!-----------------------------------------------------------------------
subroutine print_best_parameters(c_parameters, covariance)
    implicit none
    real(dp), intent(in) :: c_parameters(:), covariance(:,:)

    ! This is an example of how to define and use formats

    ! How can you use the error formula in the README to calculate the error
    ! bar in each parameter?
    
    print *, ' Best fit values:             value                 uncertainty'
    print 1, ' Volume parameter:   ', c_parameters(1),         sqrt(covariance(1,1))
    print 1, ' Surface parameter:  ', c_parameters(2),         sqrt(covariance(2,2))
    print 1, ' Asymmetry parameter:', c_parameters(3),         sqrt(covariance(3,3))
    print 1, ' Coulomb parameter:  ', c_parameters(4),         sqrt(covariance(4,4))
    print 1, ' Pairing term:       ', c_parameters(5),         sqrt(covariance(5,5))
    print 1, ' Extra term:         ', c_parameters(5),         sqrt(covariance(6,6))

1 format(a,f15.8,e28.16)
end subroutine print_best_parameters



!-----------------------------------------------------------------------
!! function: semi_empirical_mass
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function computes the binding energy using the SEMF given the
!! number of protons and neutrons for a particular isotope. In addition,
!! the function utilizes linear_terms which supplies the coefficients
!! which were computed in separate functions. c represents the c_vol,
!! c_surf, etc.
!!----------------------------------------------------------------------
!! Input:
!!
!! c    real        Array containing the parameters of the semi-empirical mass formula
!! Z    integer     number of protons in an isotope
!! N    integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r    real        Binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_mass(c, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: c(:)
    integer, intent(in) :: Z, N
    real(dp), allocatable :: linear_terms(:)
    integer :: c_size, i

    ! Initialize variable with size equal to linear terms size. 
    c_size = size(c)

    ! Allocate array with length equal to number of terms in SEMF. 
    allocate(linear_terms(1:c_size))

    ! Calculates coefficients for a given proton and neutron number. 

    call calculate_linear_terms(Z, N, linear_terms(:)) 

    ! Compute binding energy. 
    r = 0
    do i = 1, c_size 
        r = r + c(i)*linear_terms(i)
    enddo

end function semi_empirical_mass

!-----------------------------------------------------------------------
!! function: semi_empirical_error
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Computes the error in the binding energy calculation.
!!----------------------------------------------------------------------
!! Input:
!!
!! covariance   real        2D array containing the parameters' covariance matrix
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        statistical uncertainty in the binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_error(covariance, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: covariance(:,:)
    integer, intent(in) :: Z, N
    real(dp), allocatable :: linear_terms(:), n_parameters(:)
    integer :: i, j, n_terms
    real(dp) :: s

    ! Allocate array holding row and column number of alpha and covariance matrix. 

    allocate(n_parameters(1:2))
    n_parameters = shape(covariance) 

    ! Row number of covariance will be equal to number of terms in binding energy. 
    n_terms = n_parameters(1)

    ! Allocate array to hold linear terms for given proton and neutron number.
    allocate(linear_terms(1:n_terms))

    ! Fill array of linear terms for given proton and neutron number. 
    call calculate_linear_terms(Z, N, linear_terms(:))

    ! Begin do loop to calculate error for each binding energy.
    s = 0
    do i=1,n_terms
    ! First term in linear terms array will be the i'th. "i" is also row number in 
    ! covariance matrix. 
        do j=1,n_terms
    ! Second term in linear terms array will be the j'th. "j" is also column number
    ! in covariance matrix. 
            s = s + linear_terms(i)*linear_terms(j)*covariance(i,j)
        enddo
    enddo
    ! Take the square root to find the error. 
    r = sqrt(s)

end function semi_empirical_error

!--------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!! Subroutine: most_stable_n
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine finds the neutron number that produces
!! the lowest binding energy per nucleon for a given proton number. 
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters    real      1D array consisting of parameters which multiply each term in linear_term to get BE.
!-----------------------------------------------------------------------
!! Output:
!!
!! n_stable       integer    1D array containing neutron number associated with lowest BE/A for a particular proton number
!-----------------------------------------------------------------------
subroutine most_stable_n(c_parameters, n_stable)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(out), allocatable :: n_stable(:)
    integer :: n, z, index, m, n_max, p_max
    real(dp) :: b, BE_theory
    allocate(n_stable(1:118))

! p_max is largest number of protons seen in input file

    p_max = size(n_stable)

z=1
index = 1

! Find the largest number of neutrons for an given proton number
! by calling the function find_largest_n .

n_max = find_largest_n(c_parameters, p_max )

! Start loop to find neutron number corresponding to lowest binding energy per nucleon. 
! Outer loop will loop over proton number, starting at one and ending at the highest
! number of protons achieved, p_max. This is 118 for EXPERIMENT_AME2016.dat input file.

do z = 1, p_max

! Variable "b" to store binding energy per nucleon, starting with first binding 
! energy per nucleon value (corresponding to protons number = neutrons number = 1).

    b = semi_empirical_mass(c_parameters, z, z)/real(z+z,kind=dp)

! Loop over all neutron numbers for given z value defined in outer loop.

    do n = 1, n_max

! Calculate binding energy per nucleon (BE/A) 

    BE_theory = semi_empirical_mass(c_parameters, z, n)/real(n+z,kind=dp)
        if (BE_theory < b) then
! If this second BE/A is smaller than the outside loop's BE/A, we store the second BE/A
! as the new "b" value. 
            b = BE_theory

! Store the current neutron number (that achieved this low BE/A) as "m".
            m = n   
        endif

! This is repeated until all neutron numbers have been checked for a given proton number.

    enddo

! Store the neutron number associated with lowest BE/A into array n_stable
! The index of the array starts at one and increases by one everytime we enter a value.  

    n_stable(index) = m
    index = index + 1    
enddo
end subroutine most_stable_n

!-----------------------------------------------------------------------
!! Subroutine: neutron_drip_line
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine finds the largest neutron number for which the
!! neutron separation energy is positive for a given proton number. 
!! Neutron separation energy is defined as BE(Z,N-1) - BE(Z,N).
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters   real        1D array containing the parameters
!! p_max          integer     largest number of protons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! n_dripline    integer      1D array containing the neutron number associated with lowest separation energy for given proton number
!------------------------------------------------------------------------
subroutine neutron_drip_line(c_parameters, p_max, n_dripline)
    implicit none
    integer, intent(in) :: p_max
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(out), allocatable :: n_dripline(:)
    integer :: z, n, m, n_max
    real(dp) :: term_1, term_2, s

! Allocate array to hold the largest neutron numbers. 

    allocate(n_dripline(1:p_max))

! Start loop to retrieve this largest neutron number for positive separation energy.

do z = 1, p_max
    n = 1

! Calculates first term 
    term_1 =  semi_empirical_mass(c_parameters, z, n)
! Calculates second term  
    term_2 = semi_empirical_mass(c_parameters, z, n+1)
    do while (term_1 - term_2 > 0)

! While the difference is positive, perform these actions, then check to see if 
! they are still positive for the next pair of neutron numbers in this 
! specific isotope. 

        term_1 =  semi_empirical_mass(c_parameters, z, n)
        m = n + 1
        term_2 = semi_empirical_mass(c_parameters, z, m)
        n = n + 1 
    enddo

! Once this loop exits, the last value of "m" we stopped on will be one less than the 
! largest neutron number achieved for this particular proton number.
! We increase "z" by one and repeat the process, until we reach z = p_max.

    n_dripline(z) = m + 1

enddo

end subroutine neutron_drip_line

!-----------------------------------------------------------------------
!! function: find_largest_n
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function uses the largest proton number from the input file and
!! impliments a while loop to determine the largest neutron number for 
!! which the difference in BE of two neighboring isotopes becomes 
!! negative. This number is then used to determine how large the loops over 
!! "n" should be in the dripine and stable isotope subroutines.
!!  
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters   real        1D array containing the parameters
!! p_max          integer     largest number of protons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        neutron number for largest isotope of 118 proton element
!-----------------------------------------------------------------------
integer function find_largest_n(c_parameters, p_max) result(r)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(in) :: p_max
    integer :: p, n, m
    real(dp) :: term_2, term_1

! P is the proton number for the isotope. 
    p = p_max
    n = 1
! Calculate term_1 
    term_1 =  semi_empirical_mass(c_parameters, p, n)
! Calculate term_2 
    m = n + 1
    term_2 = semi_empirical_mass(c_parameters, p, m)
    do while (term_1 - term_2 > 0)
! While the difference is positive, perfore these actions, then check to see if 
! they are still positive for the next pair of neutron numbers in this 
! specific isotope. 
        term_1 =  semi_empirical_mass(c_parameters, p, n)
        m = n + 1
        term_2 = semi_empirical_mass(c_parameters, p, m)
        n = n + 1
    enddo

! Once this loop exits, the last value of "n" we stopped on will be the largest 
! neutron number achieved. Now we know how high our loop over neutron number 
! needs to be for each isotope.
    r = m

end function find_largest_n

end module nuclear_model