module bfgs_logic
  ! Contains the logic used to drive the L-BFGS algorithm.
  ! This should not need to be modified.
  !
  ! The 'heavy lifting' is done the subroutines in lbfgs_maths,
  ! it is those steps that are open to parallelisation, or the calculation
  ! of the gradient, which is handled in the potential file.
  !
  ! These steps have been modified from the pele code to run in pure fortran,
  ! lbfgs_maths has been used directly.
  use potential
  use input_output
  implicit none

  ! Result Values
  integer :: n_steps, n_fev, iteration, point, CONV_COUNT, n_grad
  double precision :: H0

  ! Internal Values
  ! Some of these may be replaced with the working array W
  double precision, dimension(:), allocatable :: W
  double precision, dimension(:), private, allocatable :: dXold, dGold, stp, H0vec
  double precision, dimension(:), private, allocatable :: X0, G0
  double precision, private :: E0

  contains
    subroutine quench()
      ! Run the minimiser.
      ! Keep iterating until we reach preset limit or reach
      ! convergence.
      !
      ! If this is the first quench.
      if (.not. allocated(dXold)) then
        allocate(dXold(N))
        allocate(dGold(N))
        allocate(stp(N))
        allocate(H0vec(N))
        allocate(X0(N))
        allocate(G0(N))
        allocate(W(0:N*(2*M+1) + 2*M-1))
      end if

      call calc_energy_gradient()
      rms = norm2(G)/sqrt(dble(N))

      n_fev = 1

      iteration = 0
      point = 0
      H0 = H0init
      H0vec = H0init

      do while ((iteration .lt. max_iterations) .and. (.not. is_stop_criteron()))
        call one_iteration()
      end do
      write(convergence_io,'(A,I6)') "Quench Complete", basin_hopping_iter
      write(*,*) "Quench Complete", basin_hopping_iter
      write(convergence_io,'(A)') " " 
      write(convergence_io,'(A)') " " 
    end subroutine quench

    subroutine one_iteration()
      ! Take a single step in the L-BFGS algorithm
      !

      call get_step(stp)

      ! We have to do this somewhat differently
      ! (or at least it seems) as the pele code doesn't modify the arrays in
      ! place.
      E0 = E
      X0 = X
      G0 = G
      call adjust_step_size()

      dXold = X - X0
      dGold = G - G0

      !rms = norm2(G)/sqrt(dble(N))
      rms = norm2(G)/sqrt(dble(N))
      if (debug) then
        write(convergence_io,'(A,ES15.8,A,ES15.8,A,I6,A,ES10.4)') &
          'Energy   ',E,'   RMS   ', rms, '   after step', iteration , '  of size  ', step_size
      end if
    end subroutine one_iteration

    subroutine get_step(stp)
      ! This is overriding the python code
      ! Wrapper to calculations.
      ! This may all get moved to the external file
      !
      ! See if the initial step has been taken
      ! There is fortran code of this already
      double precision, intent(out), dimension(N) :: stp
      integer :: ISPT, IYPT, NPT

      ! This all deals with making sure the array is the
      ! is allocated correctly.
      ! However this is copied from python...
      if (iteration .ne. 0) then
        ISPT = N + 2 * M
        IYPT = ISPT + N * M

        !NPT  = N * ((point + M - 1) MOD M)
        NPT  = N * MOD(point + M - 1, M)
        ! This seems pretty inefficient
        ! TODO: Maybe write to W directly?
        W(ISPT + NPT : ISPT + NPT + N - 1) = dXold
        W(IYPT + NPT : IYPT + NPT + N - 1) = dGold
      end if

      call MYLBFGS_UPDATESTEP(iteration, N, M, G, W, H0vec, point, stp)

      H0 = H0vec(1)
      iteration = iteration + 1
      point = MOD(iteration, M)
    end subroutine get_step

    subroutine adjust_step_size()
      ! Make sure that the step is sensible, and take it.
      !
      ! This is known as a backtracking line search
      !
      ! variables
      ! =internal=
      ! f: step scale
      !
      ! =global=
      ! X, G, E
      !
      double precision :: f, step_size_l
      integer :: n_decrease

      ! Assign copies
      X0 = X
      G0 = G
      E0 = E
      f = 1d0

      ! Check that we aren't going uphill
      if (dot_product(G, stp) .gt. 0) then
        write(*,*) "Returned uphill step. Reversing"
        stp = -stp
      end if

      step_size_l = norm2(stp)

      if (step_size_l * f .gt. max_step_size) then
        f = max_step_size / step_size_l
      end if


      ! decrease step size until it is accepted
      n_decrease = 0
      do while (n_decrease < 10)
        X = X0 + f * stp
        call calc_energy_gradient()
        n_fev = n_fev + 1

        ! If the energy rise is too great then reduce step size.
        if (is_accept_step(E, E0)) then
          ! Accepted, all is good
          !write(*,*) "Step Accepted"
          exit
        else
          ! Rejected, panic.
          f = f/10d0
          n_decrease = n_decrease + 1
          write(50, '(A)') "Energy increased, finding new step size."
          if (n_decrease .gt. 10) then
            write(*,*) "FAILED"
            write(*,*) "Too many failures in adjust_step_size"
            exit
          end if
          E = E0
          G = G0
          X = X0
        end if
      end do
      ! The step_size global is just for debugging.
      step_size = step_size_l
    end subroutine adjust_step_size

    logical function is_stop_criteron()
	double precision :: e0
      ! Test to see if convergence has been reached.
      ! At the moment this is done by means of a simple RMS check
      !
      ! N.B. stop_signal == .true. implies convergence.
      if (rms .lt. convergence_rms) then
		CONV_COUNT=0
		!JRP 05/08/16: Convergence is now defined as the gradient of each node being < convergence_rms
		do n_grad=1,N
			if (G(n_grad) > convergence_rms) then
				CONV_COUNT = CONV_COUNT+1
			endif
		enddo
		IF (CONV_COUNT == 0) THEN
			is_stop_criteron = .true.
        		!write(50, '(A)') "Convergence Reached"
			PRINT *, E
			!PRINT *, E, VOLUME, VOLCONST
			!Print *, E, E_STRETCH, E_BEND
			!open(57,file='energies.out')
			!write(57,*) E, E_STRETCH, E_BEND
			!close(57)
			!CALL COMPUTE_ENERGY(E0)
			!print *, 'Energy = ', E0
			!open(57,file='g2_e.out')
			!write(57,*) E0, E, A_SV, E_LV
			!CLOSE(57)
		ELSE
			is_stop_criteron = .false.
		ENDIF
      else
        is_stop_criteron = .false.
      end if
    end function is_stop_criteron

    logical function is_accept_step(E_new, E_old)
      ! Do we accept the new step. By default this is dE < dE_max?
      ! We may change this later on.
      double precision, intent(inout) :: E_old, E_new
      double precision :: dE

      ! Relative or absolute energy check
      if (relative_energy_check .eqv. .true.) then
        if (E_old .eq. 0d0) then
          E_old = 1d-40
        end if
        dE = (E_new - E_old)/abs(E_old)
      else
        dE = E_new - E_old
      end if
      is_accept_step = (dE .lt. dE_max)
    end function
end module
