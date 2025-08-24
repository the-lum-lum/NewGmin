subroutine perturb_coordinates_default(amplitude)
  use user_variables
  ! Default perturbation of coordinates.
  !
  ! This uses a uniform random perturbation across all coordinates
  ! This may be overridden by the user by updating ...
  !
  double precision, intent(in) :: amplitude
  double precision :: shift(N)

  call random_number(shift)
  X = X + shift*amplitude*2d0 - amplitude
end subroutine perturb_coordinates_default
