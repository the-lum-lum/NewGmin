! styblinski_tang
! Reference for the data.f90 file. 
!
! Copy this into the working directory by running GMIN -c
!
! The values in comments are the defaults used.
module user_variables
  ! Set of global variables controlled by the user
  !
  use global_variables
  ! Array parameters, these must be declared.
  integer, parameter :: N = 100
  integer, parameter :: M = 15

  ! Functions to override 
  logical, parameter :: run_init = .false.
  logical, parameter :: take_step_override = .false.

  ! Declare model globals here.
  ! This replaces the commons file.

  double precision, dimension(N) :: X, G
  contains
    subroutine set_defaults()
      ! Assign values in the model here.
      ! Replaces the data file in old GMIN


      ! Can also override variables related to L-BFGS or basin hopping 
      ! These are defined in global_variables

      !###################!
      ! L-BFGS Parameters !
      !###################!

      integer :: max_iterations           = 5000
      double precision :: convergence_rms = 1d-10
      double precision :: max_step_size   = 0.25d0
      double precision :: H0init          = 1d-1
      ! Use the relative (or absoulte) energy to test the max_erise 
      ! during the L-BFGS step.
      logical :: relative_energy_check    = .true.
      double precision :: dE_max          = 1d-1
      ! Print convergence information to GMIN_OUT
      logical :: debug                    = .true.

      !##########################!
      ! Basin Hopping Parameters !
      !##########################!
      
      ! Number of basin hopping steps to take
      integer :: max_runs                = 5
      ! Store unique energy values in RAM (minima are on disk)
      double precision, allocatable :: e_store(:)
      ! Relative error for minima to be considered unique when writing 
      ! to disk
      double precision :: max_similarity = 1d-4
      double precision :: temperature    = 100d0
      ! Output results to binary files, if not, text is used.
      logical :: binary_io               = .false.
      ! Save only the lowest i minima. 0 Saves all 
      integer :: minima_save_limit       = 0
    end subroutine set_defaults
end module 
