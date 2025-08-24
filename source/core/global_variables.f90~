module global_variables
  ! This contains a set of sensible global defaults and global variables needed
  ! by the minimiser.
  !
  ! These are not designed to be overwritten by the user, instead add the 
  ! values to the list in the users data.f90 and they will be used instead.
  !
  ! If such a file doesn't exist, call GMIN -c to create one.
  ! 
  implicit none 

  integer :: N, M
  ! Global Variables, should not be changed by the user.
  double precision :: E, step_size, rms
  ! Current iteration of basin hopping. 
  integer :: basin_hopping_iter 

  !JRP 30/06/17 A dirty way to force new droplet outputs into the lowest file
  double precision :: VOL_OUT, AREA_LV

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
  ! Store unique energy values in RAM (minima on disk)
  double precision, allocatable :: e_store(:)
  ! Relative error for minima to be considered unique when writing 
  ! to disk
  double precision :: max_similarity = 1d-4
  double precision :: temperature    = 100d0
  ! Output results to binary files, if not, text is used.
  logical :: binary_io               = .false.
  ! Save only the lowest i minima. 0 Saves all 
  integer :: minima_save_limit       = 0
  
  double precision, allocatable, dimension(:) :: X, G
  contains
    subroutine allocate_arrays()
      ! Allocate the arrays after the size has been set by user.
      ! To save memory, X and G are compiled in and so have to be set 
      ! in the users data.f90 file.
      !
      call set_defaults()
      allocate(e_store(max_runs))
      allocate(X(N))
      allocate(G(N))
    
    end subroutine
end module
