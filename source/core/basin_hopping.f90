module basin_hopping
  use bfgs_logic
  use input_output
  implicit none
  ! Note that we do have e_store, however, this will not necessarily 
  ! correspond to the last accepted energy.
  double precision, private :: EOld
  double precision, private :: EPrev = 1d10
  ! High starting value ensures that the first step is accepted, as 
  ! we would expect.
  contains
    subroutine basin_hop()
      ! Run the basin hopping algorithm, adding unique minima to storage.
      !
      do while (basin_hopping_iter .lt. max_runs)
        call basin_hop_step()
      end do
    end subroutine basin_hop

    subroutine basin_hop_step()
      ! Take a step of the basin hopping run
      !
      ! Create an backup of the coordinates in case the step is rejected.
      ! TODO: This energy check makes no sense
      double precision, dimension(N) :: XOld
      XOld = X
      EOld = E


      ! Perturb coordinates 
      ! Call from external module
      if (perturb_coordinates_override) then
        call perturb_coordinates()
      else
        call perturb_coordinates_default(4d0)
      end if

      ! Minimise the function
      ! This sets X and E

      call quench()

      ! Accept or reject step according to criterion
      if (is_step_accept()) then

        ! Keep X and E 
        ! Store if unique
        call store_minima()
        EPrev = E

      else 

        ! Revert back to the old coordinates so they may be perturbed again
        ! Todo: do we add these to file?
        X = XOld
        E = EOld

      end if
      basin_hopping_iter = basin_hopping_iter + 1
    end subroutine basin_hop_step

    subroutine perturb_coordinates_default(amplitude)
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

    logical function is_step_accept()
      ! Accept or reject the basin hopping step, based on the Metropolis criterion.
      !
      ! If the energy of new minima is less than the last it is accepted
      ! automatically, if not, then there is a probability of success based
      ! off of the energy difference.
      !
      double precision :: rand, accept_prob, dprand

      if (E .lt. EPrev) then
        is_step_accept = .true.
        write(50, '(A)') "Accepted: E < EOld"
        return
      else 
        rand = DPRAND()
        accept_prob = exp( (EPrev - E) / temperature ) 
        if (rand < accept_prob) then
          is_step_accept = .true.
          write(50, '(A, G16.10)') "Accepted: on odds of", accept_prob
        else
          is_step_accept = .false.
          write(50, '(A, G16.10)') "Rejected: on odds of", accept_prob
        end if
      end if
    end function is_step_accept
end module basin_hopping
