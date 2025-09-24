!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! firm_price_wrapper.f90
!
! Call the code to solve the example dynamic price adjustment
! model for the Structural Summer School.
!
! Stephen Terry
!
! This Version: August 2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program firm_price_wrapper
use params
use base_lib
use omp_lib
use solve
implicit none

double precision :: x(nparam),fx

x(1) = 0.0 !dummy input

!call the wrapper function
fx = fGMM(x)

end program firm_price_wrapper
