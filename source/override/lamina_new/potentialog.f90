!MODULE PFWetting
MODULE potential
  use global_variables
  implicit none
  DOUBLE PRECISION, ALLOCATABLE :: PHI(:,:,:), COORDSG(:), GRAD(:), FINITE_DIFF(:)
  DOUBLE PRECISION AREA, VOLUME, DROPVOL
  DOUBLE PRECISION, ALLOCATABLE :: WetEnergy(:), ConAngle(:)
  DOUBLE PRECISION, ALLOCATABLE :: DGRADSQPHI(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE ::  AREAPHI(:)
  INTEGER :: DROPNODES, n_calls
  LOGICAL :: CONSTVOL=.FALSE.
  LOGICAL :: RANGEVOL=.FALSE.
  LOGICAL :: TEST_COLLAPSE=.FALSE.
  LOGICAL :: QUARTERSYM=.FALSE.
  LOGICAL :: E_CONSTRAINT=.FALSE.

  ! Functions to override
  logical, parameter :: run_init = .true.
  logical, parameter :: perturb_coordinates_override = .true.

  ! Declare model globals here.
  ! This replaces the commons file.

  DOUBLE PRECISION ::  PFCA1, PFCA2, KX, KY, MAG_POT

  INTEGER :: NPOSTX, NPOSTY, WIDTHX, WIDTHY, TOPEXX, TOPEXY, HEIGHT1, HEIGHT2, LIPX, LIPY, LIPZ
  INTEGER, ALLOCATABLE :: SURFLIST(:,:), NODESTATE(:,:,:), NORMLIST(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: WEIGHTLIST(:), WETCONTACT(:), SURFWEIGHT(:), CORNERLIST(:,:)
  INTEGER :: GRIDX, GRIDY, GRIDZ, SURFNODES, MINIMISE_ENERGY
  DOUBLE PRECISION CM1, CM2, CM3, CM4, CM5, VOLCONST, PFepsilon, PFGRIDSIZE, PRESSURE, VMIN, VMAX, VWIDTH, VAV
  DOUBLE PRECISION ::  E_0, E_TARGET, E_TEST, DELTA_PHI
  DOUBLE PRECISION, ALLOCATABLE :: WEIGHT(:)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS,COORDSO
  INTEGER MCSTEPS(3), ISTEP
CONTAINS
  !Compute the energy, energy derivatives and gradient (if called for)
  subroutine calc_energy_gradient()
    implicit none
    ! Wrapper to phase field model
    call PFGRAD2_implement(X, G, E, .true.)
  end subroutine

  subroutine init()
    ! Wrapped for code that needs to be called
    write(*,*) "InitWetting"
    call INIPFWetting()
    write(*,*) "InitWetting end"
  end subroutine

  subroutine perturb_coordinates()
    ! Wrapper
    IMPLICIT NONE
    INTEGER J1, J2, J3, Cur, Cur2, CAPFILL
    INTEGER JP, SHIFT, RANX, RANY, XSTART, XSTOP, YSTART, YSTOP
    DOUBLE PRECISION PI, DPRAND, COORDS0(GRIDX*GRIDY*GRIDZ,1), &
      COORDS2(GRIDX*GRIDY*GRIDZ), COORDS3(GRIDX*GRIDY*GRIDZ,1), &
      COORDS4(GRIDX*GRIDY*GRIDZ,1), RANZ
    !DPRAND is a double-precision random number

 OPEN(30, FILE = 'coords')
  READ(30, *) COORDS0
  CLOSE(30)

!RETURN
n_calls= 0.0
DO J1=1,GRIDX
	DO J2=1,GRIDY
		DO J3=1,GRIDZ
			Cur  = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			IF (COORDS0(CUR,1) ==0.0) THEN
				COORDS0(CUR,1) = 1.0
!				COORDS0(CUR,1) = COORDS0( 0+0+J3 , 1 )
			ENDIF
		ENDDO
	ENDDO
ENDDO

X=COORDS0(:,1)
RETURN
    ISTEP = basin_hopping_iter


 IF (ISTEP == 10) THEN
	RANZ=GRIDZ+1
	CAPFILL=-1
  ELSEIF (ISTEP == 20) THEN
	RANZ = 0
	CAPFILL=-1
  ELSEIF (ISTEP ==30) THEN
	RANZ = 0
	CAPFILL=1
  ELSEIF (ISTEP == 40) THEN
	RANZ = HEIGHT1-LIPZ+0
	CAPFILL=-1
  ELSEIF (ISTEP == 50) THEN
	RANZ = HEIGHT1-LIPZ+0
	CAPFILL=1
  ELSEIF (ISTEP >=0 ) THEN
	RANZ = HEIGHT1+HEIGHT2
	CAPFILL=-1
  ELSEIF (ISTEP == 70) THEN
	RANZ = HEIGHT1+HEIGHT2
	CAPFILL=1
  ELSE
	RANZ = INT(DPRAND()*(HEIGHT1+HEIGHT2+6))
        CAPFILL=1
  ENDIF


  XSTART = INT(1 + ((GRIDX/NPOSTX)-WIDTHX)/2)
  YSTART = INT(1 + ((GRIDY/NPOSTY)-WIDTHY)/2)
  XSTOP = XSTART+WIDTHX
  YSTOP = YSTART+WIDTHY

  PRINT *, RANZ, CAPFILL
  DO J1 = 1, GRIDX
	DO J2 = 1, GRIDY
		DO J3 = 1, GRIDZ
			Cur  = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			IF (J3 > RANZ) THEN
				COORDSG(Cur) = 1
			ELSEIF (J3 == RANZ) THEN
				 COORDSG(Cur) = 0.5
			ELSE
				COORDSG(Cur) = -1
			ENDIF
			
			IF ((J1 >= XSTART-TOPEXX + LIPX) .AND. (J1 <= XSTOP + TOPEXX -LIPX) .AND. (J2 >= YSTART-TOPEXY + LIPY) &
		            .AND. (J2 <= YSTOP + TOPEXY -LIPY) .AND. (J3 <= HEIGHT1) .AND. (J3 >= HEIGHT1-2*LIPZ-8)) THEN
				COORDSG(CUR) = CAPFILL !*DPRAND()
			ENDIF
			
			!IF ( (J2 > YSTART-TOPEXY+LIPY) .AND. (J3 <= HEIGHT1+HEIGHT2-3) ) THEN
			!	COORDSG(CUR,JP) = CAPFILL*DPRAND()
			!ENDIF

			!IF ( (J2>=GRIDY/2) .AND. (J3 <= HEIGHT1+HEIGHT2) ) THEN
			!	COORDSG(CUR,JP) = -1.0
			!ENDIF
		ENDDO
	ENDDO
  ENDDO
X = COORDSG



  end subroutine perturb_coordinates

!Compute the energy, energy derivatives and gradient (if called for)
Subroutine PFGRAD2_IMPLEMENT(COORDS2,V,E,GTEST)
  
  IMPLICIT NONE
  LOGICAL GTEST
  INTEGER J1, J2, J3, Cur, S1
  DOUBLE PRECISION COORDS2(N), E, ENERGY_GRADIENT, MODE2_GRADIENT
  DOUBLE PRECISION :: V(N)
  

  E = 0.0
  GRAD(:) = 0.0
  V(:) = 0.0
  MINIMISE_ENERGY=0
 n_calls = n_calls+1


  !Loop over the grid to fill the 3d array PHI with the phi values input as the 1D array COORDSG
  DO J1 = 1, GRIDX
     DO J2 = 1, GRIDY
        DO J3 = 1, GRIDZ
           Cur = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3 !Assign each point in 3D to an element in the 1D array PHI
           PHI(J1,J2,J3) = COORDS2(Cur) !Fill PHI with the phi values at each point input into the function
	   IF (NODESTATE(J1,J2,J3) ==-1) THEN
		PHI(J1,J2,J3) = 0
		COORDS2(CUR) = 0
	   ENDIF

	   IF (QUARTERSYM) THEN
		IF ((J1 > GRIDX/2) .OR. (J2 > GRIDY/2) ) THEN
			PHI(J1,J2,J3) = 0.0
			COORDS2(CUR) = 0.0
		ENDIF
	   ENDIF


        ENDDO
     ENDDO
  ENDDO
  
  COORDSG = COORDS2

MINIMISE_ENERGY=4
!For minimise energy:
! 0 - perform a grad^2 minimisation
! 1 - perform a 'normal' energy minimisation
! 2 - compute the energy of an input state. If this is not converged in the grad^2 minimisation scheme, grad^2 minimisation will occur
! 3 - perform a grad^2 minimisation with an approximate derivative of grad^2
! 4 - test code. Compare G with a finite-difference of E



 IF (MINIMISE_ENERGY==1) THEN
	E=0
 	CALL COMPUTE_ENERGY(E)
 	!print *, 'E=', E

  	!This computes the values of (del Energy)

  	!(del Energy) needs to be stored for use in COMPUTE_DERIVATIVES
  	!but here the quantity to be minimise 'E' = |del Energy|^2

  	IF (GTEST) THEN
 	   	CALL COMPUTE_GRAD()
  	ENDIF

	!OPEN(26, FILE='grad.out')
	!write(26,'(F20.10)') GRAD*CM1
	!CLOSE(26)

	!STOP
	
	!IF (n_calls==1500) then
		!open(13,file='energies.out')
		!write(13,*) E*CM1
	!	open(12, file='coords.out')
	!	write(12,'(f20.10)') coordsg
	!	close(12)
	!	STOP
	!Else
		!stop
	!endif


	IF (CONSTVOL) THEN
		E = CM2 *(VOLUME-VOLCONST)**2
	ENDIF

	IF (RANGEVOL) THEN
		E = E + CM2*( 2 + tanh(1/VWIDTH*(VOLUME-VMAX)) - tanh(1/VWIDTH*(VOLUME-VMin)) )*(VOLUME-VAV)**2
	ENDIF


        IF (TEST_COLLAPSE) THEN
		CALL LABEL_COLLAPSE()
	ENDIF

  	V = GRAD * CM1
  	E = E * CM1

ELSEIF (MINIMISE_ENERGY==0) THEN

  IF (RANGEVOL) THEN
	! (This is just to get the total volume for the system)
	VOLUME = 0.0
	DO J1=1,GRIDX
		DO J2=1,GRIDY
			DO J3=1,GRIDZ
				CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				IF (NODESTATE(J1,J2,J3) >=0) THEN	
					VOLUME = VOLUME + (COORDSG(CUR) + 1.0)/2.0*WEIGHT(CUR)
				ENDIF
			ENDDO
		ENDDO
	ENDDO
  ENDIF



  CALL COMPUTE_GRAD()
   
  E=0
  DO J1 = 1,GRIDX*GRIDY*GRIDZ
	E = E + GRAD(J1)**2
  ENDDO
 
   IF (RANGEVOL) THEN
	E = E + CM2*( 2 + tanh(1/VWIDTH*(VOLUME-VMAX)) - tanh(1/VWIDTH*(VOLUME-VMin)) )*(VOLUME-VAV)**2
   ENDIF

  IF (E_CONSTRAINT) THEN
	E_0=0.0
	CALL COMPUTE_ENERGY(E_0)
	E = E + CM3*(E_0 - E_TARGET)**2
  ENDIF

  !This computes the values of each d/d(phi_i) |del Energy|^2
  IF (GTEST) THEN
     CALL COMPUTE_DERIVATIVES(V)

     IF (E_CONSTRAINT) THEN
	V = V+2*CM3*(E_0-E_TARGET)*GRAD
     ENDIF

  ENDIF

	!IF (n_calls==2275) then
		!open(13,file='energies.out')
		!write(13,*) E*CM1
	!	open(12, file='coords.out')
	!	write(12,'(f20.10)') coordsg
	!	close(12)
	!	STOP
	!Else
		!stop
	!endif

  V = V * CM1
  E = E * CM1

	!IF (n_calls==3000) then
		!open(13,file='energies.out')
		!write(13,*) E*CM1
	!	open(12, file='coords.out')
	!	write(12,'(f20.10)') coordsg
	!	close(12)
	!	print *, 'stopped at n=', n_calls
	!	STOP
	!Else
		!stop
	!endif



ELSEIF (MINIMISE_ENERGY == 2) THEN

 CALL COMPUTE_ENERGY(E)
 print *, 'E=', E

  CALL COMPUTE_GRAD()
   
  E=0
  DO J1 = 1,GRIDX*GRIDY*GRIDZ
	E = E + GRAD(J1)**2
  ENDDO
 

  !This computes the values of each d/d(phi_i) |del Energy|^2
  IF (GTEST) THEN
     CALL COMPUTE_DERIVATIVES(V)
  ENDIF


  V = V * CM1
  E = E * CM1

ELSEIF (MINIMISE_ENERGY == 3) THEN

  CALL COMPUTE_GRAD()
   
  E=0
  DO J1 = 1,GRIDX*GRIDY*GRIDZ
	E = E + GRAD(J1)**2
  ENDDO
 !PRINT*, E

  !This computes the values of each d/d(phi_i) |del Energy|^2
  IF (GTEST) THEN
     CALL ALTERNATIVE_COMPUTE_DERIVATIVES(V)
  ENDIF


  V = V * CM1
  E = E * CM1

ELSEIF (MINIMISE_ENERGY == 4) THEN

   IF (ALLOCATED(FINITE_DIFF)) THEN
   ELSE 
   ALLOCATE(FINITE_DIFF(GRIDX*GRIDY*GRIDZ))
   ENDIF
   FINITE_DIFF=0.0

   !CALL COMPUTE_GRAD()
   !OPEN(10, FILE='grad.out')
   !write(10,'(1F20.10)') grad*CM1
   !close(10)

   CALL COMPUTE_ENERGY(E)
   E_TEST=E
   DELTA_PHI=1.0D-8
   DO J1=1,GRIDX
	DO J2=1,GRIDY
		DO J3=1,GRIDZ
			CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3	
			COORDSG(CUR)=COORDSG(CUR) + DELTA_PHI
			CALL COMPUTE_ENERGY(E)
			FINITE_DIFF(CUR) = (E-E_TEST)/DELTA_PHI * CM1
			!print *, FINITE_DIFF(CUR)
			COORDSG(CUR)=COORDSG(CUR) - DELTA_PHI
		ENDDO
	ENDDO
   ENDDO

   !OPEN(11, FILE='finite.out')
   !write(11,'(1F20.10)') finite_diff
   !close(11)

   
   V=FINITE_DIFF
   E = E_TEST*CM1
	
				

ENDIF

END SUBROUTINE PFGRAD2_IMPLEMENT




SUBROUTINE COMPUTE_ENERGY(E)

  
  IMPLICIT NONE
  INTEGER J1, J2, J3, J1P, J1M, J2M, J2P, CUR, CUR1P, CUR1M, CUR2P, CUR2M, CUR3P, CUR3M, norm(3)
  DOUBLE PRECISION :: A_ELEMENT, grtemp, grtemp_p, grtemp_m, E
  
  VOLUME=0.0
  E=0.0

  DO J1 = 1, GRIDX
      DO J2 = 1, GRIDY
          DO J3 = 1, GRIDZ

                CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1P = (J1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1M = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR2P = (J1-1)*GRIDY*GRIDZ + (J2)*GRIDZ + J3
		CUR2M = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
		CUR3P = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3+1
		CUR3M = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3-1

		
		IF (QUARTERSYM) THEN
			
			IF (J1 == 1) THEN
				!CUR1M = (2-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				CUR1M=CUR
			ENDIF
			IF (J1 == GRIDX/2) THEN
				CUR1P = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3

			ENDIF
			IF (J2 == 1) THEN
				!CUR2M =  (2-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
				CUR2M = CUR
			ENDIF
			IF (J2 == GRIDY/2) THEN
				CUR2P = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
			ENDIF
			IF (J3 == GRIDZ) THEN
				CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
			ENDIF


		ELSE

			IF (J1 == 1) THEN
				CUR1M = (GRIDX-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				!CUR1M=CUR
			ENDIF
			IF (J1 == GRIDX) THEN
				CUR1P = (J2-1)*GRIDZ + J3
				!CUR1P=CUR
			ENDIF
			IF (J2 == 1) THEN
				CUR2M =  (J1-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
				!CUR2M=CUR
			ENDIF
			IF (J2 == GRIDY) THEN
				CUR2P = (J1-1)*GRIDY*GRIDZ + J3
				!CUR2P=CUR
			ENDIF
			IF (J3 == GRIDZ) THEN
				CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
				!CUR3M = CUR
			ENDIF

		ENDIF

		IF (NODESTATE(J1,J2,J3) >= 0) THEN
		IF (CONSTVOL) THEN
			VOLUME=VOLUME+COORDSG(CUR)*WEIGHT(CUR)
		ELSE
			VOLUME=VOLUME+((COORDSG(CUR)+1.0)/2.0)*WEIGHT(CUR)
		ENDIF
		ENDIF

		IF ( NODESTATE(J1,J2,J3) == 0 ) THEN

			!Bulk contribution
			E = E + ( 0.25*COORDSG(CUR)**4 - 0.5*COORDSG(CUR)**2+0.25)*WEIGHT(CUR)/PFEPSILON
			!Pressure contribution
			E = E + (COORDSG(CUR)+1)*WEIGHT(CUR)*PRESSURE/2.0
			!Potential contribution
			!E = E + MAG_POT*((COORDSG(CUR)+1.0)/2.0*WEIGHT(CUR)*J3)**2
			!Spatial gradient contribution from element i
			IF (GRIDX /= 1) THEN
			E = E + ( (COORDSG(CUR1P)- COORDSG(CUR))**2 + (COORDSG(CUR1M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)
			ENDIF
			E = E + ( (COORDSG(CUR2P)- COORDSG(CUR))**2 + (COORDSG(CUR2M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)
			E = E + ( (COORDSG(CUR3P)- COORDSG(CUR))**2 + (COORDSG(CUR3M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)

		ELSEIF( NODESTATE(J1,J2,J3) > 0 ) THEN		
			
			!Bulk contribution
			E = E + ( 0.25*COORDSG(CUR)**4 - 0.5*COORDSG(CUR)**2+0.25)*WEIGHT(CUR)/PFEPSILON

			!Pressure contribution
			E = E + (COORDSG(CUR)+1)*WEIGHT(CUR)*PRESSURE/2.0
			!Potential contribution
			!E = E + MAG_POT*((COORDSG(CUR)+1.0)/2.0*WEIGHT(CUR)*J3)**2
			!Surface contribution
			E = E + WETENERGY(CUR)*SURFWEIGHT(NODESTATE(J1,J2,J3))*( (-1.0/6.0)*COORDSG(CUR)**3 + (1.0/2.0)*COORDSG(CUR) + (1.0/3.0) )
			norm = NORMLIST(:,NODESTATE(J1, J2, J3))

			IF (GRIDX /= 1) THEN
			IF (norm(1) == 0) THEN
				E = E + ( (COORDSG(CUR1P)- COORDSG(CUR))**2 + (COORDSG(CUR1M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)
			ELSEIF (norm(1) == 1) THEN
				E = E + ( (COORDSG(CUR1P)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)				
			ELSEIF (norm(1) == -1) THEN
				E = E + ( (COORDSG(CUR1M)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF
			ENDIF

			IF (norm(2) == 0) THEN
				E = E + ( (COORDSG(CUR2P)- COORDSG(CUR))**2 + (COORDSG(CUR2M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)
			ELSEIF (norm(2) == 1) THEN
				E = E + ( (COORDSG(CUR2P)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)				
			ELSEIF (norm(2) == -1) THEN
				E = E + ( (COORDSG(CUR2M)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF

			IF (norm(3) == 0) THEN
				E = E + ( (COORDSG(CUR3P)- COORDSG(CUR))**2 + (COORDSG(CUR3M) - COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(4.0*PFGRIDSIZE**2)
			ELSEIF (norm(3) == 1) THEN
				E = E + ( (COORDSG(CUR3P)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)				
			ELSEIF (norm(3) == -1) THEN
				E = E + ( (COORDSG(CUR3M)- COORDSG(CUR))**2 )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF

		ENDIF
		

	    ENDDO
	ENDDO
    ENDDO

 OPEN(24,FILE='volume.out')
 write(24,*) VOLUME
 close(24)

END SUBROUTINE COMPUTE_ENERGY




! Compute del E (the gradient of the energy with respect to variations in each individual phi_i)
SUBROUTINE COMPUTE_GRAD()

  
  IMPLICIT NONE
  INTEGER J1, J2, J3, J1P, J1M, J2M, J2P, CUR, CUR1P, CUR1M, CUR2P, CUR2M, CUR3P, CUR3M, norm(3), J3P, J3M, JN
  DOUBLE PRECISION :: A_ELEMENT, grtemp, grtemp_p, grtemp_m

  

  DO J1 = 1, GRIDX
      DO J2 = 1, GRIDY
          DO J3 = 1, GRIDZ

                CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1P = (J1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1M = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR2P = (J1-1)*GRIDY*GRIDZ + (J2)*GRIDZ + J3
		CUR2M = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
		CUR3P = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3+1
		CUR3M = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3-1

		J1P=J1+1
		J1M=J1-1
		J2P=J2+1
		J2M=J2-1
		J3P=J3+1
		J3M=J3-1



		IF (QUARTERSYM) THEN
			
			IF (J1 == 1) THEN
				!CUR1M = (2-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				!J1M=2
				CUR1M=CUR
				J1M=J1
			ENDIF
			IF (J1 == GRIDX/2) THEN
				CUR1P = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				J1P=J1-1
			ENDIF
			IF (J2 == 1) THEN
				CUR2M =  (2-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
				J2M=2
				CUR2M=CUR
				J2M=J2
			ENDIF
			IF (J2 == GRIDY/2) THEN
				CUR2P = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
				J2P=J2-1
			ENDIF
			IF (J3 == GRIDZ) THEN
				CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
				J3P=GRIDZ
			ENDIF


		ELSE

			IF (J1 == 1) THEN
				CUR1M = (GRIDX-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
				J1M = GRIDX
				!CUR1M=CUR
				!J1M=J1
			ENDIF
			IF (J1 == GRIDX) THEN
				CUR1P = (J2-1)*GRIDZ + J3
				J1P=1
				!CUR1P=CUR
				!J1P=J1
			ENDIF
			IF (J2 == 1) THEN
				CUR2M =  (J1-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
				J2M=GRIDY
				!CUR2M=CUR
				!J2M=J2
			ENDIF
			IF (J2 == GRIDY) THEN
				CUR2P = (J1-1)*GRIDY*GRIDZ + J3
				J2P=1
				!CUR2P=CUR
				!J2P=J2
			ENDIF
			IF (J3 == GRIDZ) THEN
				CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
				!CUR3M = CUR
				J3P=GRIDZ
			ENDIF

		ENDIF

		grtemp = 0
		grtemp_p = 0
		grtemp_m = 0
		IF ( NODESTATE(J1,J2,J3) == 0 ) THEN

			!Bulk contribution
			grtemp = grtemp + ( COORDSG(CUR)**3 - COORDSG(CUR) )*WEIGHT(CUR)/PFEPSILON
			!Pressure contribution
			grtemp = grtemp + WEIGHT(CUR)*PRESSURE/2.0
    			IF (CONSTVOL) THEN
    				grtemp = grtemp + 2.0D0 * CM2 * (Volume-VOLCONST) * WEIGHT(Cur)
   	 		ENDIF
		

!    			IF (RANGEVOL) THEN
!    				grtemp = grtemp + CM2*WEIGHT(CUR)/(2.0*VWIDTH) * &
!& 					( 1/(cosh(1/VWIDTH*(VMAX-VOLUME)))**2 - 1/(cosh(1/VWIDTH*(VMIN-VOLUME)))**2 )
!   	 		ENDIF


			!Spatial gradient contribution from element i
			IF (GRIDX /= 1) THEN
			grtemp = grtemp - ( COORDSG(CUR1P) + COORDSG(CUR1M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF
			grtemp = grtemp - ( COORDSG(CUR2P) + COORDSG(CUR2M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			IF (J3 /= GRIDZ) THEN
			grtemp = grtemp - ( COORDSG(CUR3P) + COORDSG(CUR3M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF

			!Spatial gradient contribution from elements i+1
			IF (GRIDX /= 1) THEN
			IF (nodestate(J1P,J2,J3) == 0) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1P,J2,J3) > 0) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(PFGRIDSIZE**2)
			ENDIF
			ENDIF
			IF (nodestate(J1,J2P,J3) == 0) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2P,J3) > 0) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(PFGRIDSIZE**2)
			ENDIF
			IF ((nodestate(J1,J2,J3P) == 0) .AND. (J3P /= GRIDZ)) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)
			ELSEIF ((nodestate(J1,J2,J3P) > 0) .AND. (J3P /= GRIDZ)) THEN
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(PFGRIDSIZE**2)
			ENDIF

			!Spatial gradient contribution from elements i-1
			IF (GRIDX /= 1) THEN
			IF (nodestate(J1M,J2,J3) == 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1M,J2,J3) > 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(PFGRIDSIZE**2)
			ENDIF
			ENDIF
			IF (nodestate(J1,J2M,J3) == 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2M,J3) > 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(PFGRIDSIZE**2)
			ENDIF
			IF (nodestate(J1,J2,J3M) == 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2,J3M) > 0) THEN
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(PFGRIDSIZE**2)
			ENDIF

		ELSEIF( NODESTATE(J1,J2,J3) > 0 ) THEN		

			!Bulk contribution
			grtemp = grtemp + ( COORDSG(CUR)**3 - COORDSG(CUR) )*WEIGHT(CUR)/PFEPSILON

			!Pressure contribution
			grtemp = grtemp + WEIGHT(CUR)*PRESSURE/2.0
			!Surface contribution
			grtemp = grtemp + WETENERGY(CUR)*SURFWEIGHT(NODESTATE(J1,J2,J3))*( (-1.0/2.0)*COORDSG(CUR)**2 + (1.0/2.0) )
			!External potential contribution
			!grtemp = grtemp + MAG_POT*(COORDSG(CUR)+1.0) * (J3*PFGRIDSIZE)**2 * WEIGHT(CUR)**2 / 2.0
    			IF (CONSTVOL) THEN
    				grtemp = grtemp + 2.0D0 * CM2 * (Volume-VOLCONST) * WEIGHT(Cur)
   	 		ENDIF
!    			IF (RANGEVOL) THEN
!    				grtemp = grtemp + CM2*WEIGHT(CUR)/(2.0*VWIDTH) * &
!& 					( 1/(cosh(1/VWIDTH*(VMAX-VOLUME)))**2 - 1/(cosh(1/VWIDTH*(VMIN-VOLUME)))**2 )
!   	 		ENDIF


			norm = NORMLIST(:,NODESTATE(J1, J2, J3))
			JN = NODESTATE(J1,J2,J3)

			IF (GRIDX /= 1) THEN
			IF (norm(1) == 0) THEN
				grtemp = grtemp - ( COORDSG(CUR1P) + COORDSG(CUR1M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
			ELSEIF (norm(1) == 1) THEN
				grtemp = grtemp - ( COORDSG(CUR1P) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)				
			ELSEIF (norm(1) == -1) THEN
				grtemp = grtemp - ( COORDSG(CUR1M) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
			ENDIF
			ENDIF

			IF (norm(2) == 0) THEN
				grtemp = grtemp - ( COORDSG(CUR2P) + COORDSG(CUR2M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
			ELSEIF (norm(2) == 1) THEN
				grtemp = grtemp - ( COORDSG(CUR2P) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)				
			ELSEIF (norm(2) == -1) THEN
				grtemp = grtemp - ( COORDSG(CUR2M) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
			ENDIF

			IF (norm(3) == 0) THEN
				grtemp = grtemp - ( COORDSG(CUR3P) + COORDSG(CUR3M) - 2*COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
			ELSEIF (norm(3) == 1) THEN
				grtemp = grtemp - ( COORDSG(CUR3P) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p = grtemp_p + ( COORDSG(CUR) - COORDSG(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)				
			ELSEIF (norm(3) == -1) THEN
				grtemp = grtemp - ( COORDSG(CUR3M) - COORDSG(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m = grtemp_m + ( COORDSG(CUR) - COORDSG(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
			ENDIF

		ENDIF
		
		GRAD(CUR) = grtemp + grtemp_p + grtemp_m

    			IF ((RANGEVOL) .AND. (NODESTATE(J1,J2,J3) >=0) .AND. (MINIMISE_ENERGY == 1) ) THEN
    				GRAD(CUR) = GRAD(CUR) + CM2*WEIGHT(CUR)/(2.0*VWIDTH) * &
& 					( 1/(cosh(1/VWIDTH*(VMAX-VOLUME)))**2 - 1/(cosh(1/VWIDTH*(VMIN-VOLUME)))**2 )*(VOLUME-VAV)**2
    				GRAD(CUR) = GRAD(CUR) + CM2*WEIGHT(CUR) * &
& 					( 2 + tanh(1/VWIDTH*(VOLUME-VMAX)) - tanh(1/VWIDTH*(VOLUME-VMin)) )*(VOLUME-VAV)
   	 		ENDIF


	    ENDDO
	ENDDO
    ENDDO


END SUBROUTINE COMPUTE_GRAD


!Compute the derivative of |del Energy|^2 with respect to each phi_i
SUBROUTINE COMPUTE_DERIVATIVES(V)

  
  IMPLICIT NONE
  INTEGER :: J1, J2, J3, CUR, CUR1P, CUR1M, CUR2P, CUR2M, CUR3P, CUR3M, J1P, J1M, J2P, J2M, J3P, J3M, JN
  INTEGER :: J1PP, J1MM, J2PP, J2MM, J3PP, J3MM
  DOUBLE PRECISION :: V(N), grtemp, grtemp_p1, grtemp_p2, grtemp_p3, grtemp_m1, grtemp_m2, grtemp_m3, grrang, multiplier 
!  DOUBLE PRECISION, ALLOCATABLE :: HESSG(:,:)
  INTEGER ::  norm(3)

!  ALLOCATE(HESSG(GRIDX*GRIDY*GRIDZ,GRIDX*GRIDY*GRIDZ))
!  HESSG(:,:) = 0.0	

  DO J1 = 1, GRIDX
     DO J2 = 1, GRIDY
        DO J3 = 1, GRIDZ
		!PRINT *, J1,J2,J3,NODESTATE(J1,J2,J3)
                CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1P = (J1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1M = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR2P = (J1-1)*GRIDY*GRIDZ + (J2)*GRIDZ + J3
		CUR2M = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
		CUR3P = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3+1
		CUR3M = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3-1

		J1P=J1+1
		J1M=J1-1
		J2P=J2+1
		J2M=J2-1
		J3P=J3+1
		J3M=J3-1

		J1PP=J1+2
		J1MM=J1-2
		J2PP=J2+2
		J2MM=J2-2
		J3PP=J3+2
		J3MM=J3-2

		IF (J1 == 1) THEN
			CUR1M = (GRIDX-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			J1M=GRIDX
			J1MM=GRIDX-1
		ENDIF
		IF (J1 == GRIDX) THEN
			CUR1P = (J2-1)*GRIDZ + J3
			J1P=1
			J1PP=2
		ENDIF
		IF (J2 == 1) THEN
			CUR2M =  (J1-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
			J2M=GRIDY
			J2MM=GRIDY-1
		ENDIF
		IF (J2 == GRIDY) THEN
			CUR2P = (J1-1)*GRIDY*GRIDZ + J3
			J2P=1
			J2PP=2
		ENDIF
		IF (J3 == GRIDZ) THEN
			CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
			!CUR3M = CUR
		ENDIF

		IF (J1 ==2) THEN
			J1MM = GRIDX
		ENDIF
		IF (J1 == GRIDX-1) THEN
			J1PP = 1
		ENDIF
		IF (J2 == 2) THEN
			J2MM = GRIDY
		ENDIF
		IF (J2 == GRIDY-1) THEN
			J2PP = 1
		ENDIF
		
		grtemp = 0
	
		grtemp_p1 = 0
		grtemp_p2 = 0
		grtemp_p3 = 0

		grtemp_m1 = 0
		grtemp_m2 = 0
		grtemp_m3 = 0

		IF ( NODESTATE(J1,J2,J3) == 0 ) THEN

			! Calculate d^2E/phi_i^2
			! Bulk term
			grtemp = grtemp + ( 3*COORDSG(CUR)**2 - 1.0 )*WEIGHT(CUR)/PFEPSILON

			!External potential contribution
			!grtemp = grtemp + MAG_POT* (J3*PFGRIDSIZE)**2 * WEIGHT(CUR)**2 / 2.0



			!Gradient terms arising from phi_i gradient
			IF (GRIDX /= 1) THEN
			grtemp = grtemp + 1 * PFEPSILON*WEIGHT(CUR)/PFGRIDSIZE**2
			ENDIF
			grtemp = grtemp + 2 * PFEPSILON*WEIGHT(CUR)/PFGRIDSIZE**2

			IF (GRIDX /= 1) THEN
			IF (nodestate(J1P,J2,J3) >= 0) THEN
				IF (nodestate(J1P,J2,J3) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)
					grtemp_p1 = grtemp_p1 - PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)
					IF (nodestate(J1PP,J2,J3) == 0) THEN
						grtemp_p1 = grtemp_p1 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1PP,J2,J3) > 0) THEN
						grtemp_p1 = grtemp_p1 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
				ELSEIF (nodestate(J1P,J2,J3) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR1P)/(PFGRIDSIZE**2)
					grtemp_p1 = grtemp_p1 - PFEPSILON*WEIGHT(CUR1P)/(PFGRIDSIZE**2)
					grtemp_p1 = grtemp_p1 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF
			ENDIF
			ENDIF

			IF (nodestate(J1,J2P,J3) >= 0) THEN
				IF (nodestate(J1,J2P,J3) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)
					grtemp_p2 = grtemp_p2 - PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)
					IF (nodestate(J1,J2PP,J3) == 0) THEN
						grtemp_p2 = grtemp_p2 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1,J2PP,J3) > 0) THEN
						grtemp_p2 = grtemp_p2 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
				ELSEIF (nodestate(J1,J2P,J3) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR2P)/(PFGRIDSIZE**2)
					grtemp_p2 = grtemp_p2 - PFEPSILON*WEIGHT(CUR2P)/(PFGRIDSIZE**2)
					grtemp_p2 = grtemp_p2 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF			
			ENDIF

			IF (J3P <= GRIDZ) THEN
			IF (nodestate(J1,J2,J3P) >= 0) THEN			
				IF (nodestate(J1,J2,J3P) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)
					grtemp_p3 = grtemp_p3 - PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)
					IF (J3PP <= GRIDZ) THEN
					IF (nodestate(J1,J2,J3PP) == 0) THEN
						grtemp_p3 = grtemp_p3 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1,J2,J3PP) > 0) THEN
						grtemp_p3 = grtemp_p3 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
					ENDIF
				ELSEIF (nodestate(J1,J2,J3P) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR3P)/(PFGRIDSIZE**2)
					grtemp_p3 = grtemp_p3 - PFEPSILON*WEIGHT(CUR3P)/(PFGRIDSIZE**2)
					grtemp_p3 = grtemp_p3 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF
			ENDIF
			ENDIF

			IF (GRIDX /= 1) THEN
			IF (nodestate(J1M,J2,J3) >= 0) THEN
				IF (nodestate(J1M,J2,J3) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)
					grtemp_m1 = grtemp_m1 - PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)
					IF (nodestate(J1MM,J2,J3) == 0) THEN
						grtemp_m1 = grtemp_m1 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1MM,J2,J3) > 0) THEN
						grtemp_m1 = grtemp_m1 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
				ELSEIF (nodestate(J1M,J2,J3) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR1M)/(PFGRIDSIZE**2)
					grtemp_m1 = grtemp_m1 - PFEPSILON*WEIGHT(CUR1M)/(PFGRIDSIZE**2)
					grtemp_m1 = grtemp_m1 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF
			ENDIF
			ENDIF

			IF (nodestate(J1,J2M,J3) >= 0) THEN
				IF (nodestate(J1,J2M,J3) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)
					grtemp_m2 = grtemp_m2 - PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)
					IF (nodestate(J1,J2MM,J3) == 0) THEN
						grtemp_m2 = grtemp_m2 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1,J2MM,J3) > 0) THEN
						grtemp_m2 = grtemp_m2 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
				ELSEIF (nodestate(J1,J2M,J3) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR2M)/(PFGRIDSIZE**2)
					grtemp_m2 = grtemp_m2 - PFEPSILON*WEIGHT(CUR2M)/(PFGRIDSIZE**2)
					grtemp_m2 = grtemp_m2 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF
			ENDIF

			IF (nodestate(J1,J2,J3M) >= 0) THEN
				IF (nodestate(J1,J2,J3M) == 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)
					grtemp_m3 = grtemp_m3 - PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)
					IF (nodestate(J1,J2,J3MM) == 0) THEN
						grtemp_m3 = grtemp_m3 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
					ELSEIF (nodestate(J1,J2,J3MM) > 0) THEN
						grtemp_m3 = grtemp_m3 - PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
					ENDIF
				ELSEIF (nodestate(J1,J2,J3M) > 0) THEN
					grtemp = grtemp + PFEPSILON*WEIGHT(CUR3M)/(PFGRIDSIZE**2)
					grtemp_m3 = grtemp_m3 - PFEPSILON*WEIGHT(CUR3M)/(PFGRIDSIZE**2)	
					grtemp_m3 = grtemp_m3 - PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				ENDIF
			ENDIF



		ELSEIF( NODESTATE(J1,J2,J3) > 0 ) THEN		

			! Bulk term
			grtemp = grtemp + ( 3*COORDSG(CUR)**2 - 1.0 )*WEIGHT(CUR)/PFEPSILON
			! Arising from phi_i gradient
			!IF (GRIDX /= 1) THEN
			!	grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
			!ELSE
			!	grtemp = grtemp + 2*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
			!ENDIF
			
			!External potential contribution
			!grtemp = grtemp + MAG_POT* (J3*PFGRIDSIZE)**2 * WEIGHT(CUR)**2 / 2.0

			!Surface contribution
			grtemp = grtemp + WETENERGY(CUR)*SURFWEIGHT(NODESTATE(J1,J2,J3))*(-1.0*COORDSG(CUR))

			norm = NORMLIST(:,NODESTATE(J1, J2, J3))
			JN = NODESTATE(J1,J2,J3)
			
			IF (GRIDX /= 1) THEN
			IF (norm(1) == 0) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p1 = grtemp_p1 + PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)
				grtemp_p1 = grtemp_p1 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,NODESTATE(J1P,J2,J3))

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m1 = grtemp_m1 + PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
				grtemp_m1 = grtemp_m1 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,NODESTATE(J1M,J2,J3))

			ELSEIF (norm(1) == 1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p1 = grtemp_p1 + PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)
				grtemp_p1 = grtemp_p1 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				grtemp_m1=0.0

			ELSEIF (norm(1) == -1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m1 = grtemp_m1 + PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
				grtemp_m1 = grtemp_m1 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				grtemp_p1=0.0

			ENDIF
			ENDIF


			
			IF (norm(2) == 0) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p2 = grtemp_p2 + PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)
				grtemp_p2 = grtemp_p2 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,NODESTATE(J1,J2P,J3))

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m2 = grtemp_m2 + PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
				grtemp_m2 = grtemp_m2 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,NODESTATE(J1,J2M,J3))

			ELSEIF (norm(2) == 1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p2 = grtemp_p2 + PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)
				grtemp_p2 = grtemp_p2 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				grtemp_m2=0.0

			ELSEIF (norm(2) == -1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m2 = grtemp_m2 + PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
				grtemp_m2 = grtemp_m2 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				grtemp_p2=0.0
			ENDIF



			IF (norm(3) == 0) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp +  PFEPSILON*WEIGHT(CUR)/PFGRIDSIZE**2
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p3 = grtemp_p3 + PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)
				grtemp_p3 = grtemp_p3 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,NODESTATE(J1,J2,J3P))

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m3 = grtemp_m3 + PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
				grtemp_m3 = grtemp_m3 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,NODESTATE(J1,J2,J3M))

			ELSEIF (norm(3) == 1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)

				!Gradient terms arising from phi_i+1 gradient
				grtemp_p3 = grtemp_p3 + PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)
				grtemp_p3 = grtemp_p3 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				!grtemp_m3=0.0


			ELSEIF (norm(3) == -1) THEN
				!Gradient terms arising from phi_i gradient
				grtemp = grtemp + 1*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp = grtemp + PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)

				!Gradient terms arising from phi_i-1 gradient
				grtemp_m3 = grtemp_m3 + PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
				grtemp_m3 = grtemp_m3 + PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)*2.0

				!Solid node contribution is supressed
				grtemp_p3=0.0

			ENDIF

		ENDIF



!		HESSG(CUR, CUR) = HESSG(CUR, CUR) + grtemp

!		HESSG(CUR1P,CUR) = HESSG(CUR1P,CUR) + grtemp_p1
!		HESSG(CUR,CUR1P) = HESSG(CUR,CUR1P) + grtemp_p1

!		HESSG(CUR2P,CUR) = HESSG(CUR2P,CUR) + grtemp_p2
!		HESSG(CUR,CUR2P) = HESSG(CUR,CUR2P) + grtemp_p2
		
!		HESSG(CUR3P,CUR) = HESSG(CUR3P,CUR) + grtemp_p3
!		HESSG(CUR,CUR3P) = HESSG(CUR,CUR3P) + grtemp_p3


!		HESSG(CUR1M,CUR) = HESSG(CUR1M,CUR) + grtemp_m1
!		HESSG(CUR,CUR1M) = HESSG(CUR,CUR1M) + grtemp_m1

!		HESSG(CUR2M,CUR) = HESSG(CUR2M,CUR) + grtemp_m2
!		HESSG(CUR,CUR2M) = HESSG(CUR,CUR2M) + grtemp_m2
		
!		HESSG(CUR3M,CUR) = HESSG(CUR3M,CUR) + grtemp_m3
	!	HESSG(CUR,CUR3M) = HESSG(CUR,CUR3M) + grtemp_m3


		IF (NODESTATE(J1,J2,J3) >= 0) THEN
		grtemp = 2*GRAD(CUR)*grtemp

		grtemp_p1 = 2*GRAD(CUR1P)*grtemp_p1
		grtemp_p2 = 2*GRAD(CUR2P)*grtemp_p2
		grtemp_p3 = 2*GRAD(CUR3P)*grtemp_p3

		grtemp_m1 = 2*GRAD(CUR1M)*grtemp_m1
		grtemp_m2 = 2*GRAD(CUR2M)*grtemp_m2
		grtemp_m3 = 2*GRAD(CUR3M)*grtemp_m3
	
	  	V(CUR) = grtemp + grtemp_p1 + grtemp_p2 + grtemp_p3 + grtemp_m1 + grtemp_m2 + grtemp_m3

		ENDIF

    			IF ((RANGEVOL) .AND. (NODESTATE(J1,J2,J3) >=0) ) THEN
    				V(CUR) = V(CUR) + CM2*WEIGHT(CUR)/(2.0*VWIDTH) * &
& 					( 1/(cosh((VMAX-VOLUME)/VWIDTH))**2 - 1/(cosh((VMIN-VOLUME)/VWIDTH))**2 )*(VOLUME-VAV)**2
    				V(CUR) = V(CUR) + CM2*WEIGHT(CUR) * &
& 					( 2 + tanh((VOLUME-VMAX)/VWIDTH) - tanh((VOLUME-VMin)/VWIDTH) )*(VOLUME-VAV)
   	 		ENDIF



        ENDDO
     ENDDO
  ENDDO

IF (RANGEVOL) THEN
	multiplier = 0.0
	!Start by making the multiplier common to all the V(cur)
	multiplier = CM2/(2.0*VWIDTH**2)*( tanh(1/VWIDTH*(VMAX-VOLUME))/( COSH(1/VWIDTH*(VMAX-VOLUME)) )**2 &
	&			   - tanh(1/VWIDTH*(VMIN-VOLUME))/( COSH(1/VWIDTH*(VMIN-VOLUME)) )**2 )
	multiplier = multiplier*dot_product(WEIGHT,GRAD)

	!Now add the volume constraint onto the gradient V.
	DO J1=1, GRIDX
		DO J2=1,GRIDY
			DO J3=1,GRIDZ
               			CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			
				grrang = multiplier*WEIGHT(CUR)				
				!V(CUR) = V(CUR) + 2*GRAD(CUR)*grrang
		
			ENDDO
		ENDDO
	ENDDO					
ENDIF


!  HessG = HESSG * CM1

 ! OPEN(96, FILE='hess.out')
 ! WRITE(96,'(144F20.10)') HESSG
 ! CLOSE(96)
 ! STOP

END SUBROUTINE COMPUTE_DERIVATIVES






!---------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE ALTERNATIVE_COMPUTE_DERIVATIVES(V)

  
  IMPLICIT NONE
  INTEGER :: J1, J2, J3, CUR, CUR1P, CUR1M, CUR2P, CUR2M, CUR3P, CUR3M, J1P, J1M, J2P, J2M, J3P, J3M, JN
  INTEGER :: J1PP, J1MM, J2PP, J2MM, J3PP, J3MM, JK
  DOUBLE PRECISION :: V(N), grtemp, grtemp_p1, grtemp_p2, grtemp_p3, grtemp_m1, grtemp_m2, grtemp_m3, DELTA_PHI
  DOUBLE PRECISION ::  GRAD_PLUS(7), GRAD_MINUS(7), GSQP, GSQM, COORDS_SAVED(GRIDX*GRIDY*GRIDZ)
!  DOUBLE PRECISION, ALLOCATABLE :: HESSG(:,:)
  INTEGER ::  norm(3)

 DELTA_PHI=1.0D-8
 V(:) = 0.0
  COORDS_SAVED=COORDSG
  DO J1=1,GRIDX
	DO J2=1,GRIDY
		DO J3=1,GRIDZ
			
                	CUR = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3

			!1) Get the grad+
			GRAD_PLUS=DIFF_GRAD(CUR,COORDSG,DELTA_PHI, J1, J2, J3)
			COORDSG=COORDS_SAVED

			!2) Get the grad-
			GRAD_MINUS=DIFF_GRAD(CUR,COORDSG,-1*DELTA_PHI, J1, J2, J3)
			COORDSG=COORDS_SAVED

			!3) Compute the derivative
			DO JK=1,7
				V(CUR) = V(CUR) - GRAD_PLUS(JK)**2 - GRAD_MINUS(JK)**2
			ENDDO
		ENDDO
	ENDDO
  ENDDO
V=V/(2*DELTA_PHI)


END SUBROUTINE ALTERNATIVE_COMPUTE_DERIVATIVES



FUNCTION DIFF_GRAD(CUR,COORDS_IN, DELTA_PHI, J1, J2, J3)
  
  IMPLICIT NONE
  INTEGER J1, J2, J3, J1P, J1M, J2M, J2P, CUR, CUR1P, CUR1M, CUR2P, CUR2M, CUR3P, CUR3M, norm(3), J3P, J3M, JN
  DOUBLE PRECISION :: A_ELEMENT, grtemp, COORDS_IN(GRIDX*GRIDY*GRIDZ), DIFF_GRAD(7), DELTA_PHI
  DOUBLE PRECISION :: grtemp_p1, grtemp_m1, grtemp_p2, grtemp_m2, grtemp_p3, grtemp_m3

		COORDS_IN(CUR) = COORDS_IN(CUR) + DELTA_PHI  

		CUR1P = (J1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR1M = (J1-2)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
		CUR2P = (J1-1)*GRIDY*GRIDZ + (J2)*GRIDZ + J3
		CUR2M = (J1-1)*GRIDY*GRIDZ + (J2-2)*GRIDZ + J3
		CUR3P = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3+1
		CUR3M = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3-1

		J1P=J1+1
		J1M=J1-1
		J2P=J2+1
		J2M=J2-1
		J3P=J3+1
		J3M=J3-1

		IF (J1 == 1) THEN
			CUR1M = (GRIDX-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			J1M = GRIDX
		ENDIF
		IF (J1 == GRIDX) THEN
			CUR1P = (J2-1)*GRIDZ + J3
			J1P=1
		ENDIF
		IF (J2 == 1) THEN
			CUR2M =  (J1-1)*GRIDY*GRIDZ + (GRIDY-1)*GRIDZ + J3
			J2M=GRIDY
		ENDIF
		IF (J2 == GRIDY) THEN
			CUR2P = (J1-1)*GRIDY*GRIDZ + J3
			J2P=1
		ENDIF
		IF (J3 == GRIDZ) THEN
			CUR3P = CUR !This will set the gradients at z = GRIDZ to zero
			!CUR3M = CUR
		ENDIF

		grtemp = 0

		grtemp_p1 = 0
		grtemp_p2 = 0
		grtemp_p3 = 0

		grtemp_m1 = 0
		grtemp_m2 = 0
		grtemp_m3 = 0

		IF ( NODESTATE(J1,J2,J3) == 0 ) THEN

			!Bulk contribution
			grtemp = grtemp + ( COORDS_IN(CUR)**3 - COORDS_IN(CUR) )*WEIGHT(CUR)/PFEPSILON
			!Pressure contribution
			grtemp = grtemp + WEIGHT(CUR)*PRESSURE/2.0
			!External potential contribution
			grtemp = grtemp + MAG_POT*(COORDS_IN(CUR)+1.0) * (J3*PFGRIDSIZE)**2 * WEIGHT(CUR)**2 / 2.0

			!Spatial gradient contribution from element i
			IF (GRIDX /= 1) THEN
			grtemp = grtemp - ( COORDS_IN(CUR1P) + COORDS_IN(CUR1M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF
			grtemp = grtemp - ( COORDS_IN(CUR2P) + COORDS_IN(CUR2M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			IF (J3 /= GRIDZ) THEN
			grtemp = grtemp - ( COORDS_IN(CUR3P) + COORDS_IN(CUR3M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
			ENDIF

			!Spatial gradient contribution from elements i+1
			IF (GRIDX /= 1) THEN
			IF (nodestate(J1P,J2,J3) == 0) THEN
				grtemp_p1 = grtemp_p1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1P,J2,J3) > 0) THEN
				grtemp_p1 = grtemp_p1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(PFGRIDSIZE**2)
			ENDIF
			ENDIF
			IF (nodestate(J1,J2P,J3) == 0) THEN
				grtemp_p2 = grtemp_P2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2P,J3) > 0) THEN
				grtemp_p2 = grtemp_p2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(PFGRIDSIZE**2)
			ENDIF
			IF ((nodestate(J1,J2,J3P) == 0) .AND. (J3P /= GRIDZ)) THEN
				grtemp_p3 = grtemp_p3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)
			ELSEIF ((nodestate(J1,J2,J3P) > 0) .AND. (J3P /= GRIDZ)) THEN
				grtemp_p3 = grtemp_p3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(PFGRIDSIZE**2)
			ENDIF

			!Spatial gradient contribution from elements i-1
			IF (GRIDX /= 1) THEN
			IF (nodestate(J1M,J2,J3) == 0) THEN
				grtemp_m1 = grtemp_m1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1M,J2,J3) > 0) THEN
				grtemp_m1 = grtemp_m1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(PFGRIDSIZE**2)
			ENDIF
			ENDIF
			IF (nodestate(J1,J2M,J3) == 0) THEN
				grtemp_m2 = grtemp_m2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2M,J3) > 0) THEN
				grtemp_m2 = grtemp_m2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(PFGRIDSIZE**2)
			ENDIF
			IF (nodestate(J1,J2,J3M) == 0) THEN
				grtemp_m3 = grtemp_m3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)
			ELSEIF (nodestate(J1,J2,J3M) > 0) THEN
				grtemp_m3 = grtemp_m3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(PFGRIDSIZE**2)
			ENDIF

		ELSEIF( NODESTATE(J1,J2,J3) > 0 ) THEN		

			!Bulk contribution
			grtemp = grtemp + ( COORDS_IN(CUR)**3 - COORDS_IN(CUR) )*WEIGHT(CUR)/PFEPSILON

			!Pressure contribution
			grtemp = grtemp + WEIGHT(CUR)*PRESSURE/2.0
			!Surface contribution
			grtemp = grtemp - WETENERGY(CUR)*SURFWEIGHT(NODESTATE(J1,J2,J3))
			!External potential contribution
			grtemp = grtemp + MAG_POT*(COORDS_IN(CUR)+1.0) * (J3*PFGRIDSIZE)**2 * WEIGHT(CUR)**2 / 2.0


			norm = NORMLIST(:,NODESTATE(J1, J2, J3))
			JN = NODESTATE(J1,J2,J3)

			IF (GRIDX /= 1) THEN
			IF (norm(1) == 0) THEN
				grtemp = grtemp - ( COORDS_IN(CUR1P) + COORDS_IN(CUR1M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p1 = grtemp_p1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)
				grtemp_m1 = grtemp_m1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
			ELSEIF (norm(1) == 1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR1P) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p1 = grtemp_p1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1P) )*PFEPSILON*WEIGHT(CUR1P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(1,JN)				
			ELSEIF (norm(1) == -1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR1M) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m1 = grtemp_m1 + ( COORDS_IN(CUR) - COORDS_IN(CUR1M) )*PFEPSILON*WEIGHT(CUR1M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(2,JN)
			ENDIF
			ENDIF

			IF (norm(2) == 0) THEN
				grtemp = grtemp - ( COORDS_IN(CUR2P) + COORDS_IN(CUR2M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p2 = grtemp_p2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)
				grtemp_m2 = grtemp_m2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
			ELSEIF (norm(2) == 1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR2P) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p2 = grtemp_p2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2P) )*PFEPSILON*WEIGHT(CUR2P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(3,JN)				
			ELSEIF (norm(2) == -1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR2M) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m2 = grtemp_m2 + ( COORDS_IN(CUR) - COORDS_IN(CUR2M) )*PFEPSILON*WEIGHT(CUR2M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(4,JN)
			ENDIF

			IF (norm(3) == 0) THEN
				grtemp = grtemp - ( COORDS_IN(CUR3P) + COORDS_IN(CUR3M) - 2*COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(2.0*PFGRIDSIZE**2)
				grtemp_p3 = grtemp_p3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)
				grtemp_m3 = grtemp_m3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
			ELSEIF (norm(3) == 1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR3P) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_p3 = grtemp_p3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3P) )*PFEPSILON*WEIGHT(CUR3P)/(2.0*PFGRIDSIZE**2)*CORNERLIST(5,JN)				
			ELSEIF (norm(3) == -1) THEN
				grtemp = grtemp - ( COORDS_IN(CUR3M) - COORDS_IN(CUR) )*PFEPSILON*WEIGHT(CUR)/(PFGRIDSIZE**2)
				grtemp_m3 = grtemp_m3 + ( COORDS_IN(CUR) - COORDS_IN(CUR3M) )*PFEPSILON*WEIGHT(CUR3M)/(2.0*PFGRIDSIZE**2)*CORNERLIST(6,JN)
			ENDIF

		ENDIF
		
		DIFF_GRAD(1) = grtemp
		DIFF_GRAD(2) = grtemp_p1
		DIFF_GRAD(3) = grtemp_m1
		DIFF_GRAD(4) = grtemp_p2
		DIFF_GRAD(5) = grtemp_m2
		DIFF_GRAD(6) = grtemp_p3
		DIFF_GRAD(7) = grtemp_m3


END FUNCTION DIFF_GRAD


!---------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------






  !Compute the energy of interaction between the droplet and the surface

  SUBROUTINE INIPFWetting()
  IMPLICIT NONE
  INTEGER J1, J2, J3, Cur, S1
  DOUBLE PRECISION PI
  character(len=10)       :: datechar,timechar,zonechar
  integer                 :: values(8),itime1


  CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
  itime1= values(7)*39 + values(8)
  CALL SDPRND(itime1)
  !print *, '1'
  CALL MAKESURFACE()
  !print *, '2'
  ALLOCATE(WEIGHT(N)) !Allocate memory to the allocatable WEIGHT 1d array

  WEIGHT(:) = PFGRIDSIZE**3
  PI = ATAN(1.0D0) * 4.0D0 !Define pi
  VOLCONST = 2.0D0 * VOLCONST - GRIDX * GRIDY * GRIDZ * PFGRIDSIZE**3
 

  ! Loop over each grid point to assign each the variable weight. Note that because element 1 and element GRIDX,Y,Z are the same, the weights associated with both is halved
  !print *, '3-2'
  DO J1 = 1, GRIDX
     DO J2 = 1, GRIDY
        DO J3 = 1, GRIDZ

           Cur = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
           IF (nodestate(J1,J2,J3) > 0) THEN
              WEIGHT(Cur) = weightlist(nodestate(J1,J2,J3))*PFGRIDSIZE**3

	   elseif (J3 == GRIDZ) THEN
		WEIGHT(Cur) = 0.5*PFGRIDSIZE**3
           ELSE
              WEIGHT(Cur) = PFGRIDSIZE**3

          ENDIF

        ENDDO
     ENDDO
  ENDDO

  IF (RANGEVOL) THEN
	PRINT *, 'Volume range', VMIN, VMAX
	VAV = (VMAX+VMIN)/2.0
  ENDIF


  !print *, '4'
  ALLOCATE(ConAngle(N), WetEnergy(N))
  ALLOCATE(PHI(GRIDX,GRIDY,GRIDZ), GRAD(N))
  ALLOCATE(COORDSG(N))
  !ConAngle is the contact angle
  !Wet energy is the surface energy density arising from the liquid drop interacting with the solid surface
  WetEnergy(:) = 0.0D0
  ConAngle(:) = 90.0/180.0 * ATAN(1.0D0) * 4.0D0
  !print *, '5'
  do S1 = 1, surfnodes
      Cur = (surflist(1,S1)-1)*GRIDY*GRIDZ + (surflist(2,S1)-1)*GRIDZ + surflist(3,S1)
      ConAngle(Cur) = WETCONTACT(S1)/180.0*ATAN(1.0)*4.0 !Covert contact angle in degrees to radians
      WETENERGY(CUR) = (-1.0*SQRT(2.0))*COS(CONANGLE(CUR))
  ENDDO
  END SUBROUTINE INIPFWETTING

!Take a random step, automatically accepting the change (basin hopping => no energy acceptance criteria)
SUBROUTINE TAKESTEPPFWETTING(JP)

  
  IMPLICIT NONE
  INTEGER J1, J2, J3, Cur, Cur2
  INTEGER JP, SHIFT, RANX, RANY, RANZ, CAPFILL, XSTART, XSTOP, YSTART, YSTOP
  DOUBLE PRECISION PI, DPRAND, COORDS0(GRIDX*GRIDY*GRIDZ,1), COORDSG(GRIDX*GRIDY*GRIDZ)
  DOUBLE PRECISION :: COORDS3(GRIDX*GRIDY*GRIDZ,1), COORDS4(GRIDX*GRIDY*GRIDZ,1)
  !DPRAND is a double-precision random number

  OPEN(30, FILE = 'coords')
  READ(30, *) COORDS0
  CLOSE(30)

  return

    ISTEP = basin_hopping_iter
 IF (ISTEP == 1) THEN
	RANZ=GRIDZ+1
	CAPFILL=1
  ELSEIF (ISTEP == 2) THEN
	RANZ = 0
	CAPFILL=-1
  ELSEIF (ISTEP == 3) THEN
	RANZ = 0
	CAPFILL=1
  ELSEIF (ISTEP == 4) THEN
	RANZ = HEIGHT1-LIPZ+0
	CAPFILL=-1
  ELSEIF (ISTEP == 5) THEN
	RANZ = HEIGHT1-LIPZ+0
	CAPFILL=1
  ELSEIF (ISTEP == 6) THEN
	RANZ = HEIGHT1+HEIGHT2
	CAPFILL=-1
  ELSEIF (ISTEP == 7) THEN
	RANZ = HEIGHT1+HEIGHT2
	CAPFILL=1
  ELSE
	RANZ = INT(DPRAND()*(HEIGHT1+HEIGHT2+6))
  ENDIF


  XSTART = INT(1 + ((GRIDX/NPOSTX)-WIDTHX)/2)
  YSTART = INT(1 + ((GRIDY/NPOSTY)-WIDTHY)/2)
  XSTOP = XSTART+WIDTHX
  YSTOP = YSTART+WIDTHY

  PRINT *, RANZ, CAPFILL
  DO J1 = 1, GRIDX
	DO J2 = 1, GRIDY
		DO J3 = 1, GRIDZ
			Cur  = (J1-1)*GRIDY*GRIDZ + (J2-1)*GRIDZ + J3
			IF (J3 > RANZ) THEN
				COORDSG(Cur) = 1
			ELSEIF (J3 == RANZ) THEN
				 COORDSG(Cur) = 0
			ELSE
				COORDSG(Cur) = -1
			ENDIF
			
			IF ((J1 >= XSTART-TOPEXX + LIPX) .AND. (J1 <= XSTOP + TOPEXX -LIPX) .AND. (J2 >= YSTART-TOPEXY + LIPY) &
		            .AND. (J2 <= YSTOP + TOPEXY -LIPY) .AND. (J3 <= HEIGHT1) .AND. (J3 >= HEIGHT1-2*LIPZ-8)) THEN
				COORDSG(CUR) = CAPFILL *DPRAND()
			ENDIF
			
			!IF ( (J2 > YSTART-TOPEXY+LIPY) .AND. (J3 <= HEIGHT1+HEIGHT2-3) ) THEN
			!	COORDSG(CUR,JP) = CAPFILL*DPRAND()
			!ENDIF

			!IF ( (J2>=GRIDY/2) .AND. (J3 <= HEIGHT1+HEIGHT2) ) THEN
			!	COORDSG(CUR,JP) = -1.0
			!ENDIF
		ENDDO
	ENDDO
  ENDDO
X = COORDSG

END SUBROUTINE TAKESTEPPFWETTING

  SUBROUTINE MAKESURFACE()

implicit none

!THESE WILL BE GLOBALLY DEFINED IN THE FULL PROGRAM
!integer GRIDX, GRIDY, GRIDZ
!looping indices
integer :: J1, J2, J3, S1, P1, P2, W1, W2, J_SOLID3
!integer :: WIDTHX, WIDTHY, TOPEXX, TOPEXY, HEIGHT1, HEIGHT2, NPOSTX, NPOSTY


!REAL :: PFCA1
integer:: xstart, xstop, ystart, ystop, XP, XM, YP, YM, ZP, ZM

!integer, allocatable :: surflist(:,:), normlist(:,:)!An array of all the node poisitions on the surface
!integer, allocatable :: nodestate(:,:,:) !An array of all the nodes, taking values of 1 if in the bulk, 0 otherwise
!real, allocatable :: weightlist(:), wetcontact(:) !An array of the normal vector at each node position and weight at each position
!integer surfnodes !number of surface nodes


!Initialise typical system variables
!GRIDX = 30
!GRIDY = 30
!GRIDZ = 30

!NPOSTX = 2
!NPOSTY = 2

!WIDTHX = 3
!WIDTHY = 3

!TOPEXX = 4
!TOPEXY = 4

!HEIGHT1 = 8
!HEIGHT2 = 2

!PFCA1 = 110.0

allocate(nodestate(GRIDX, GRIDY, GRIDZ))
NODESTATE(:,:,:) = 0


! Make the surface nodes layer by layer up the z axis

DO P1 = 1,NPOSTX
	DO P2 = 1, NPOSTY
		
		XSTART = INT((P1-1)*(GRIDX/NPOSTX)+1 + ((GRIDX/NPOSTX)-WIDTHX)/2)
		YSTART = INT((P2-1)*(GRIDY/NPOSTY)+1 + ((GRIDY/NPOSTY)-WIDTHY)/2)
		XSTOP = XSTART+WIDTHX
		YSTOP = YSTART+WIDTHY
		
		print *, xstart,xstop,ystart,ystop, GRIDX/2
		J3 = 1
		
		NODESTATE(:,:,J3) = 1

		DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
			DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
				IF ((J1 > XSTART) .AND. (J1 < XSTOP) .AND. (J2 > YSTART) .AND. (J2 < YSTOP)) THEN
					NODESTATE(J1,J2,J3) = 0
				ENDIF
			ENDDO
		ENDDO

		DO J3 = 2,HEIGHT1-LIPZ
			DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
				DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
					IF ((J1 >= XSTART) .AND. (J1 <= XSTOP) .AND. (J2 >= YSTART) .AND. (J2 <= YSTOP)) THEN
						NODESTATE(J1,J2,J3) = 1
					ENDIF
					IF ((J1 > XSTART) .AND. (J1 < XSTOP) .AND. (J2 > YSTART) .AND. (J2 < YSTOP)) THEN
						NODESTATE(J1,J2,J3) = 0
					ENDIF
				ENDDO
			ENDDO
		ENDDO

		J3 = HEIGHT1-LIPZ+1
		DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
			DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
				IF ((J1 >= XSTART-TOPEXX) .AND. (J1 <= XSTOP+TOPEXX) .AND. (J2 >= YSTART-TOPEXY) .AND. (J2 <= YSTOP+TOPEXY)) THEN
					NODESTATE(J1,J2,J3) = 1
				ENDIF

				IF ((J1 > XSTART-TOPEXX+LIPX) .AND. (J1 < XSTOP+TOPEXX-LIPX)&
					& .AND. (J2 > YSTART-TOPEXY+LIPY) .AND. (J2 < YSTOP+TOPEXY-LIPY)) THEN
					NODESTATE(J1,J2,J3) = 0
				ENDIF

				IF ((J1 >= XSTART) .AND. (J1 <= XSTOP) .AND. (J2 >= YSTART) .AND. (J2 <= YSTOP)) THEN
					NODESTATE(J1,J2,J3) = 1
				ENDIF

				IF ((J1 > XSTART) .AND. (J1 < XSTOP) .AND. (J2 > YSTART) .AND. (J2 < YSTOP)) THEN
					NODESTATE(J1,J2,J3) = 0
				ENDIF
			ENDDO
		ENDDO

		DO J3 = HEIGHT1-LIPZ+2, HEIGHT1
			DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
				DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
					IF ((J1 >= XSTART-TOPEXX) .AND. (J1 <= XSTOP+TOPEXX)&
						& .AND. (J2 >= YSTART-TOPEXY) .AND. (J2 <= YSTOP+TOPEXY)) THEN
						NODESTATE(J1,J2,J3) = 1
					ENDIF

					IF ((J1 > XSTART-TOPEXX) .AND. (J1 < XSTOP+TOPEXX) .AND. (J2 > YSTART-TOPEXY) .AND. (J2 < YSTOP+TOPEXY)) THEN
						NODESTATE(J1,J2,J3) = 0
					ENDIF

					IF ((J1 >= XSTART-TOPEXX+LIPX) .AND. (J1 <= XSTOP+TOPEXX-LIPX)&
						& .AND. (J2 >= YSTART-TOPEXY+LIPY) .AND. (J2 <= YSTOP+TOPEXY-LIPY)) THEN
						NODESTATE(J1,J2,J3) = 1
					ENDIF

					IF ((J1 > XSTART-TOPEXX+LIPX) .AND. (J1 < XSTOP+TOPEXX-LIPX)&
						& .AND. (J2 > YSTART-TOPEXY+LIPY) .AND. (J2 < YSTOP+TOPEXY-LIPY)) THEN
						NODESTATE(J1,J2,J3) = 0
					ENDIF

					IF ((J1 >= XSTART) .AND. (J1 <= XSTOP) .AND. (J2 >= YSTART) .AND. (J2 <= YSTOP)) THEN
						NODESTATE(J1,J2,J3) = 1
					ENDIF

					IF ((J1 > XSTART) .AND. (J1 < XSTOP) .AND. (J2 > YSTART) .AND. (J2 < YSTOP)) THEN
						NODESTATE(J1,J2,J3) = 0
					ENDIF
				ENDDO
			ENDDO
		ENDDO

		J3 = HEIGHT1+1
		DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
			DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
				IF ((J1 >= XSTART-TOPEXX) .AND. (J1 <= XSTOP+TOPEXX)&
					& .AND. (J2 >= YSTART-TOPEXY) .AND. (J2 <= YSTOP+TOPEXY)) THEN
					NODESTATE(J1,J2,J3) = 1
				ENDIF

				IF ((J1 > XSTART-TOPEXX) .AND. (J1 < XSTOP+TOPEXX) .AND. (J2 > YSTART-TOPEXY) .AND. (J2 < YSTOP+TOPEXY)) THEN
					NODESTATE(J1,J2,J3) = 0
				ENDIF

				IF ((J1 >= XSTART-TOPEXX+LIPX) .AND. (J1 <= XSTOP+TOPEXX-LIPX)&
					& .AND. (J2 >= YSTART-TOPEXY+LIPY) .AND. (J2 <= YSTOP+TOPEXY-LIPY)) THEN
					NODESTATE(J1,J2,J3) = 1
				ENDIF
				IF ((J1 > XSTART) .AND. (J1 < XSTOP) .AND. (J2 > YSTART) .AND. (J2 < YSTOP)) THEN
					NODESTATE(J1,J2,J3) = 0
				ENDIF
			ENDDO
		ENDDO
		
		DO J3 = HEIGHT1+2,HEIGHT1+HEIGHT2
			DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
				DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
					IF ((J1 >= XSTART-TOPEXX) .AND. (J1 <= XSTOP+TOPEXX) .AND. (J2 >= YSTART-TOPEXY) .AND. (J2 <= YSTOP+TOPEXY)) THEN
						NODESTATE(J1,J2,J3) = 1
					ENDIF
					IF ((J1 > XSTART-TOPEXX) .AND. (J1 < XSTOP+TOPEXX) .AND. (J2 > YSTART-TOPEXY) .AND. (J2 < YSTOP+TOPEXY)) THEN
						NODESTATE(J1,J2,J3) = 0
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		
		J3 = HEIGHT1+HEIGHT2+1
		DO J1 = (P1-1)*(GRIDX/NPOSTX)+1,P1*(GRIDX/NPOSTX)
			DO J2 = (P2-1)*(GRIDY/NPOSTY)+1,P2*(GRIDY/NPOSTY)
				IF ((J1 >= XSTART-TOPEXX) .AND. (J1 <= XSTOP+TOPEXX) .AND. (J2 >= YSTART-TOPEXY) .AND. (J2 <= YSTOP+TOPEXY)) THEN
					NODESTATE(J1,J2,J3) = 1
				ENDIF
			ENDDO
		ENDDO
	ENDDO
ENDDO


! Calculate where the solid nodes are and update nodestate accordingly
DO J1 = 1,GRIDX
	DO J2 = 1,GRIDY
		
		J3 = GRIDZ
		W1 = 0
		W2 = 0
		
		DO WHILE (W1 == 0)
			If (J3 == 1) THEN	
				W1 = 1
			ELSEIF (NODESTATE(J1,J2,J3) == 1) THEN
				J_SOLID3 = J3-1
				DO WHILE (W2 == 0)
					IF (J_SOLID3 == 1) THEN
						NODESTATE(J1,J2,J_SOLID3) = -1
						W2 = 1
						W1 = 1
					ELSEIF (NODESTATE(J1,J2,J_SOLID3) == 0) THEN
						NODESTATE(J1,J2,J_SOLID3) = -1
					ELSEIF (NODESTATE(J1,J2,J_SOLID3) == 1) THEN
						W2 = 1
						W1 = 1
					ENDIF
					J_SOLID3 = J_SOLID3-1
				ENDDO
			ENDIF
			J3 = J3-1
		ENDDO		
	ENDDO
ENDDO


!Quartersym - Only work with one quarter of the simulation system, and enforce mirror symmetry
IF (QUARTERSYM) THEN
PRINT *, 'Quartersym operational'
  DO J1 = 1, GRIDX
     DO J2 = 1, GRIDY
        DO J3 = 1, GRIDZ
		IF ( (J1 > GRIDX/2) .OR. (J2 > GRIDY/2) ) THEN
			NODESTATE(J1,J2,J3) = -1
		ENDIF
        ENDDO
     ENDDO
  ENDDO
ENDIF


!Count the surface nodes
S1 = 0
DO J1 = 1,GRIDX
	DO J2 = 1,GRIDY
		DO J3 = 1,GRIDZ
			IF (NODESTATE(J1,J2,J3) == 1) THEN
				S1 = S1+1
			ENDIF
		ENDDO
	ENDDO
ENDDO


SURFNODES = S1



allocate(surflist(3,surfnodes))
allocate(normlist(3,surfnodes))
allocate(weightlist(surfnodes))
allocate(Wetcontact(surfnodes))
allocate(surfweight(surfnodes))
allocate(cornerlist(6,surfnodes))

surflist(:,:) = 0.0
normlist(:,:) = 0.0
WETCONTACT(:) = PFCA1
weightlist(:) = 0.0
surfweight(:)=0.0
cornerlist(:,:) = 1.0

!Calculate the normals and weights of each surface node, and label each surface node with a positive integer
S1 = 0
DO J1 = 1,GRIDX
	DO J2 = 1,GRIDY
		DO J3 = 1,GRIDZ
			
			IF (NODESTATE(J1,J2,J3) == 1) THEN
				S1 = S1+1
				surflist(:,S1) = (/J1, J2, J3/)
				NODESTATE(J1,J2,J3) = S1

				XM = J1-1
				XP = J1+1
				YM = J2-1
				YP = J2+1
				ZM = J3-1
				ZP = J3+1

				IF (J1 == 1) THEN
					XM = GRIDX
				ENDIF
				IF (J1 == GRIDX) THEN
					XP = 1
				ENDIF
				IF (J2 == 1) THEN
					YM = GRIDY
				ENDIF
				IF (J2 == GRIDY) THEN
					YP = 1
				ENDIF

				
				

				!Now calculate the normals
				IF (NODESTATE(XM,J2,J3) == -1) THEN
					NORMLIST(1,S1) = 1
				ELSEIF (NODESTATE(XP,J2,J3) == -1) THEN
					NORMLIST(1,S1) = -1
				ELSEIF ((NODESTATE(XM,J2,J3) >= 0) .AND. (NODESTATE(XP,J2,J3) == 0)) THEN
					NORMLIST(1,S1) = 1
				ELSEIF ((NODESTATE(XP,J2,J3) >= 0) .AND. (NODESTATE(XM,J2,J3) == 0)) THEN
					NORMLIST(1,S1) = -1
				ENDIF

				IF (NODESTATE(J1,YM,J3) == -1) THEN
					NORMLIST(2,S1) = 1
				ELSEIF (NODESTATE(J1,YP,J3) == -1) THEN
					NORMLIST(2,S1) = -1
				ELSEIF ((NODESTATE(J1,YM,J3) >= 0) .AND. (NODESTATE(J1,YP,J3) == 0)) THEN
					NORMLIST(2,S1) = 1
				ELSEIF ((NODESTATE(J1,YP,J3) >= 0) .AND. (NODESTATE(J1,YM,J3) == 0)) THEN
					NORMLIST(2,S1) = -1
				ENDIF

				IF ((J3 == 1) .AND. (NODESTATE(XM,J2,J3) > 0) .AND. (NODESTATE(XP,J2,J3) > 0)&
					& .AND. (NODESTATE(J1,YM,J3) > 0) .AND. (NODESTATE(J1,YP,J3) > 0))  THEN
					IF (NODESTATE(XM,YM,J3) == -1) THEN
						NORMLIST(:,S1) = (/1,1,1/)
					ELSEIF (NODESTATE(XM,YP,J3) == -1) THEN
						NORMLIST(:,S1) = (/1,-1,1/)
					ELSEIF (NODESTATE(XP,YM,J3) == -1) THEN
						NORMLIST(:,S1) = (/-1,1,1/)
					ELSEIF (NODESTATE(XP,YP,J3) == -1) THEN
						NORMLIST(:,S1) = (/-1,-1,1/)
					ELSEIF (NODESTATE(XM,YM,J3) == 0) THEN
						NORMLIST(:,S1) = (/-1,-1,1/)
					ELSEIF (NODESTATE(XM,YP,J3) == 0) THEN
						NORMLIST(:,S1) = (/-1,1,1/)
					ELSEIF (NODESTATE(XP,YM,J3) == 0) THEN
						NORMLIST(:,S1) = (/1,-1,1/)
					ELSEIF (NODESTATE(XP,YP,J3) == 0) THEN
						NORMLIST(:,S1) = (/1,1,1/)
					ENDIF

				ENDIF

				IF ((J3 == HEIGHT1+1) .AND. (NODESTATE(XM,J2,J3) > 0) .AND. (NODESTATE(XP,J2,J3) > 0)&
					& .AND. (NODESTATE(J1,YM,J3) > 0) .AND. (NODESTATE(J1,YP,J3) > 0))  THEN
					IF (NODESTATE(XM,YM,J3) == -1) THEN
						NORMLIST(:,S1) = (/1,1,-1/)
					ELSEIF (NODESTATE(XM,YP,J3) == -1) THEN
						NORMLIST(:,S1) = (/1,-1,-1/)
					ELSEIF (NODESTATE(XP,YM,J3) == -1) THEN
						NORMLIST(:,S1) = (/-1,1,-1/)
					ELSEIF (NODESTATE(XP,YP,J3) == -1) THEN
						NORMLIST(:,S1) = (/-1,-1,-1/)
					ELSEIF (NODESTATE(XM,YM,J3) == 0) THEN
						NORMLIST(:,S1) = (/-1,-1,-1/)
					ELSEIF (NODESTATE(XM,YP,J3) == 0) THEN
						NORMLIST(:,S1) = (/-1,1,-1/)
					ELSEIF (NODESTATE(XP,YM,J3) == 0) THEN
						NORMLIST(:,S1) = (/1,-1,-1/)
					ELSEIF (NODESTATE(XP,YP,J3) == 0) THEN
						NORMLIST(:,S1) = (/1,1,-1/)
					ENDIF

				ENDIF

				IF (J3 == 1) THEN
					NORMLIST(3,S1) = 1
				ELSEIF (J3 == GRIDZ) THEN
					
				ELSEIF (NODESTATE(J1,J2,ZM) == -1) THEN
					NORMLIST(3,S1) = 1
				ELSEIF (NODESTATE(J1,J2,ZP) == -1) THEN
					NORMLIST(3,S1) = -1	
				ELSEIF ((NODESTATE(J1,J2,ZM) >= 0) .AND. (NODESTATE(J1,J2,ZP) == 0)) THEN
					NORMLIST(3,S1) = 1
				ELSEIF ((NODESTATE(J1,J2,ZP) >= 0) .AND. (NODESTATE(J1,J2,ZM) == 0)) THEN
					NORMLIST(3,S1) = -1
				ENDIF

				!The weights of each node
				IF (J3==1) THEN
					IF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 3) THEN
						WEIGHTLIST(S1) = 0.5
						IF ((NODESTATE(XM,YM,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XM,YP,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YM,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YP,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF

					ELSEIF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 2) THEN
						WEIGHTLIST(S1) = 0.25
					ELSEIF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 1) THEN
						WEIGHTLIST(S1) = 0.5
					ENDIF
				ELSE
					IF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 3) THEN
						WEIGHTLIST(S1) = 1.0
						IF ((NODESTATE(XM,YM,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XM,YP,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YM,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YP,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XM,YM,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XM,YP,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YM,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF
						IF ((NODESTATE(XP,YP,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.125
						ENDIF

					ELSEIF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 2) THEN

						WEIGHTLIST(S1) = 1
						IF ((NODESTATE(XM,YM,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XM,YP,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XP,YM,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XP,YP,J3) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF

						IF ((NODESTATE(XM,J2,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XM,J2,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XP,J2,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(XP,J2,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF

						IF ((NODESTATE(J1,YM,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(J1,YM,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(J1,YP,ZM) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						IF ((NODESTATE(J1,YP,ZP) == -1)) THEN
							WEIGHTLIST(S1) = WEIGHTLIST(S1)-0.25
						ENDIF
						
						If (WEIGHTLIST(S1) < 0.75) THEN
							WEIGHTLIST(S1) = 0.25
						ENDIF

					ELSEIF ((ABS(NORMLIST(1,S1))+ABS(NORMLIST(2,S1))+ABS(NORMLIST(3,S1))) == 1) THEN
						WEIGHTLIST(S1) = 0.5
					ENDIF




				ENDIF


			ENDIF


			IF (weightlist(nodestate(J1,J2,J3)) == 0.5) THEN
				SURFWEIGHT(S1) = 1.0*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.25) THEN
				SURFWEIGHT(S1) = 1.0*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.75) THEN
				SURFWEIGHT(S1) = 1.0*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.125) THEN
				SURFWEIGHT(S1) = 0.75*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.875) THEN
				SURFWEIGHT(S1) = 0.75*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.625) THEN
				SURFWEIGHT(S1) = 1.25*PFGRIDSIZE**2
			ELSEIF (weightlist(nodestate(J1,J2,J3)) == 0.375) THEN
				SURFWEIGHT(S1) = 1.25*PFGRIDSIZE**2
			ENDIF


		ENDDO
	ENDDO
ENDDO
				
!Calculate if a surface node is adjacent to a corner, and in which direction
!If the adjacent surface node is a corner, the corresponding element in cornerlist is 0
!otherwise, it is 1
DO S1=1,SURFNODES
	J1 = SURFLIST(1,S1)+1
	IF (J1 > GRIDX) THEN
		J1 = 1
	ENDIF
	J2 = SURFLIST(2,S1)
	J3 = SURFLIST(3,S1)
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(1,NODESTATE(J1,J2,J3)) ==-1) THEN !-1
			CORNERLIST(1,S1) = 2.0
		ENDIF
	ENDIF

	J1 = SURFLIST(1,S1)-1
	IF (J1 < 1) THEN
		J1 = GRIDX
	ENDIF
	J2 = SURFLIST(2,S1)
	J3 = SURFLIST(3,S1)
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(1,NODESTATE(J1,J2,J3)) ==1) THEN !1
			CORNERLIST(2,S1) = 2.0
		ENDIF
	ENDIF

	J1 = SURFLIST(1,S1)
	J2 = SURFLIST(2,S1)+1
	IF (J2 > GRIDY) THEN
		J2 = 1
	ENDIF
	J3 = SURFLIST(3,S1)
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(2,NODESTATE(J1,J2,J3)) ==-1) THEN !-1
			CORNERLIST(3,S1) = 2.0
		ENDIF
	ENDIF

	J1 = SURFLIST(1,S1)
	J2 = SURFLIST(2,S1)-1
	IF (J2 < 1) THEN
		J2 = GRIDX
	ENDIF
	J3 = SURFLIST(3,S1)
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(2,NODESTATE(J1,J2,J3)) ==1) THEN !1
			CORNERLIST(4,S1) = 2.0
		ENDIF
	ENDIF
	
	J1 = SURFLIST(1,S1)
	J2 = SURFLIST(2,S1)
	J3 = SURFLIST(3,S1)+1
	IF (J3 <= GRIDZ) THEN
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(3,NODESTATE(J1,J2,J3)) ==-1) THEN !-1
			CORNERLIST(5,S1) = 2.0
		ENDIF
	ENDIF	
	ENDIF

	J1 = SURFLIST(1,S1)
	J2 = SURFLIST(2,S1)
	J3 = SURFLIST(3,S1)-1
	IF (J3 >= 1) THEN
	IF (NODESTATE(J1,J2,J3) > 0) THEN
		IF (NORMLIST(3,NODESTATE(J1,J2,J3)) ==1) THEN !1
			CORNERLIST(6,S1) = 2.0
		ENDIF
	ENDIF	
	ENDIF
ENDDO



			
open(12, file = 'nodestateG')
write(12,'(30I5)') nodestate
close(12)

open(10, file = 'surflistG')
write(10,'(3I5)') surflist
close(10)

open(11, file = 'normlistG')
write(11,'(3I5)') normlist
close(11)

open(13, file = 'weightlistG')
write(13,'(F7.4)') weightlist
close(13)

open(18, file='surfweightlist')
write(18, '(F20.10)') SURFWEIGHT
 close(18)

open(19,file='cornerlist')
write(19,'(6f5.3)') cornerlist
close(19)

END SUBROUTINE MAKESURFACE







!A subroutine for use in the critical pressure tests. It will find by which mechanism the system collapses
SUBROUTINE LABEL_COLLAPSE()

IMPLICIT NONE
INTEGER :: J1, J2, J3, XSTART, XSTOP, YSTART, YSTOP
DOUBLE PRECISION :: PHI_BASE, PHI_PIL, PHI_CAP
 CHARACTER(len=3) :: OUTPUT

!nodestate(xstart-1,ystart-1,2) = 1
!nodestate(xstart-1,gridy/2,2) = 2
!nodestate(xstart-1,ystart-1,HEIGHT1) = 3

!nodestate(xstart-1-topexx+LIPX+2, ystart-1-topexy+LIPY+2, HEIGHT1) = 4
!nodestate(xstart-1-topexx+LIPX+2, GRIDY/2, HEIGHT1) = 5


XSTART = INT(((GRIDX/NPOSTX)-WIDTHX)/2)
YSTART = INT(((GRIDY/NPOSTY)-WIDTHY)/2)
XSTOP = XSTART+WIDTHX
YSTOP = YSTART+WIDTHY

!Phi at the centre of the posts
PHI_BASE = PHI(1,1,1)

!Phi max along the pillar edge and face
J1 = XSTART-1
J2 = YSTART-1
PHI_PIL = PHI(J1,J2,2)
DO J3=3,HEIGHT1
	IF (PHI(J1,J2,J3) > PHI_PIL) THEN
		PHI_PIL = PHI(J1,J2,J3)
	ENDIF
ENDDO

J1 = XSTART-1
J2 = GRIDY/2
DO J3=2,HEIGHT1
	IF (PHI(J1,J2,J3) > PHI_PIL) THEN
		PHI_PIL = PHI(J1,J2,J3)
	ENDIF
ENDDO

!Phi at the corners of the cap underside



IF (LIPX /= 0) THEN
	PHI_CAP = PHI(GRIDX/2, GRIDY/2-WIDTHY/2-TOPEXY/2, HEIGHT1)
	!PHI_CAP = PHI(xstart-1-topexx+LIPX+2, ystart-1-topexy+LIPY+2, HEIGHT1)
	!IF ( PHI(xstart-1-topexx+LIPX+2, GRIDY/2, HEIGHT1) > PHI_CAP) THEN
	!	PHI_CAP = PHI(xstart-1-topexx+LIPX+2, GRIDY/2, HEIGHT1)
	!ENDIF
ELSE
!If the geometry is a reentrant structure (not doubly reentrant) the cap contact state does not exist
	PHI_CAP = -1
ENDIF




!PRINT *, PHI_BASE, PHI_PIL, PHI_CAP
!PRINT *, MAX(PHI_BASE, PHI_PIL, PHI_CAP)

IF ((PHI_BASE >= PHI_PIL) .AND. (PHI_BASE >= PHI_CAP)) THEN
	OUTPUT = 'BAS'
ELSEIF ((PHI_PIL >= PHI_BASE) .AND. (PHI_PIL >= PHI_CAP)) THEN
	OUTPUT = 'PIL'
ELSEIF ((PHI_CAP >= PHI_BASE) .AND. (PHI_CAP >= PHI_PIL)) THEN
	OUTPUT = 'CAP'
ENDIF

IF (MAX(PHI_BASE, PHI_PIL, PHI_CAP) >= 0) THEN
	OPEN(89, FILE='collapse.out')
	WRITE(89,*) OUTPUT
	CLOSE(89)
	OPEN(90, FILE='collapse.coords')
	WRITE(90,'(F20.10)') COORDSG
	CLOSE(90)

	STOP
ENDIF



END SUBROUTINE LABEL_COLLAPSE 





END MODULE potential
