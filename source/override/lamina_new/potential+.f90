MODULE potential
  use global_variables
  implicit none

  ! Functions to override
  logical, parameter :: run_init = .true.
  logical, parameter :: perturb_coordinates_override = .true.

  ! Model global variables
  DOUBLE PRECISION, ALLOCATABLE :: COORDS(:,:), K_ARRAY(:), R_ARRAY(:), B_ARRAY(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: COORDS_0(:,:), ANGLES(:), GRAD(:), GRAD_SPAT(:,:)
  INTEGER, ALLOCATABLE :: BONDS(:,:), TRIANGLES(:,:)
  DOUBLE PRECISION :: F_I, F_S, L, H, K1, K2, K3, R1, R2, R3, BI, BE, LS, HS, CM1, KC, NORM_ANG, KY, WY,LOWER_RAT
  DOUBLE PRECISION :: KAREA, CAREA, AREA0
  INTEGER :: N_ROD, N_SEG, i, j, k, N_BONDS, IO



CONTAINS
!------------------------------------------------------------------------------------------
  !Wrapper subroutines
  subroutine init()
    ! Wrapped for code that needs to be called
    write(*,*) "InitWetting"
    call INIT_SYSTEM()
    write(*,*) "InitWetting end"
  end subroutine

  subroutine calc_energy_gradient()
    implicit none
    ! Wrapper to phase field model
    call IMPLEMENT_POTENTIAL(X, G, E, .true.)
  end subroutine
!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
!System initialisation
SUBROUTINE INIT_SYSTEM()

  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: LS_ARRAY(:) 
  INTEGER :: CUR, NODE1, NODE2, NODE0, S1
  

  
  ALLOCATE(COORDS(N_ROD*(N_SEG+1),2), K_ARRAY(3), R_ARRAY(3), COORDS_0(N_ROD,2), ANGLES(N), GRAD(N), GRAD_SPAT(N_ROD*(N_SEG+1),2))
  COORDS=0.0
  ANGLES=0.0
  GRAD=0.0
  COORDS_0=0.0

  !ALLOCATE(LS_ARRAY(N_SEG))
  !DO i=1,N_SEG
  !  LS_ARRAY(i)= (L / N_SEG) + 0.001*i   
  !ENDDO
  !Define the segment lengths and vertical separation 
  LS = L/(N_SEG)
  HS = H/(N_ROD-1)

  !Fill the elastic constant and equilibrium bond length arrays
  K_ARRAY(1)=K1
  K_ARRAY(2)=K2
  K_ARRAY(3)=K3
  R_ARRAY(1)=R1
  R_ARRAY(2)=R2
  R_ARRAY(3)=R3

  !Read in the coordinates of the rod starting points
  OPEN(1,FILE='initpos')
  DO i=1,N_ROD
  	READ(1,*) COORDS_0(i,:)
  ENDDO
  CLOSE(1)

  !Read in the bonds
  N_BONDS=0
  OPEN(1, FILE='bonds')
  DO
  READ(1,*, iostat=io)
  IF  (io/=0) EXIT
  	N_BONDS=N_BONDS+1
  END DO
  CLOSE(1)
  
  ALLOCATE(BONDS(N_BONDS,3))

  OPEN(1,FILE='bonds')
  DO i=1,N_BONDS
  	READ(1,*) BONDS(i,:)
  ENDDO
  CLOSE(1)
 
  !Make a list of the node triplets used in the area repulsion term
  AREA0=0.5*LS*HS !Equilibrium area of triangles (for now they are uniform)
  ALLOCATE(TRIANGLES(2*(N_ROD-1)*N_SEG,3))
  TRIANGLES=0
  S1=1
 ! DO i=1,N_ROD-1
!	DO j=1,N_SEG+1
!		IF (j<N_SEG+1) THEN
!			NODE0=(i-1)*(N_SEG+1)+j
!			NODE1=(i)*(N_SEG+1)+j+1
!			NODE2=(i)*(N_SEG+1)+j
!			TRIANGLES(S1,1)=NODE0
!			TRIANGLES(S1,2)=NODE1
!			TRIANGLES(S1,3)=NODE2	
!			S1=S1+1
!		ENDIF
!		
!		IF (j>1) THEN
!			NODE0=(i-1)*(N_SEG+1)+j
!			NODE1=(i)*(N_SEG+1)+j
!			NODE2=(i)*(N_SEG+1)+j-1
!			TRIANGLES(S1,1)=NODE0
!			TRIANGLES(S1,2)=NODE1
!			TRIANGLES(S1,3)=NODE2	
!			S1=S1+1
!		ENDIF
!		
!	ENDDO
!  ENDDO

  DO i=1,N_ROD-1
	DO j=1,N_SEG
		NODE0=(i-1)*(N_SEG+1)+j
		NODE1=(i)*(N_SEG+1)+j+1
		NODE2=(i)*(N_SEG+1)+j
		TRIANGLES(S1,1)=NODE0
		TRIANGLES(S1,2)=NODE1
		TRIANGLES(S1,3)=NODE2
		S1=S1+1	
	ENDDO
  ENDDO
			
  DO i=2,N_ROD
	DO j=2,N_SEG+1
		NODE0=(i-1)*(N_SEG+1)+j
		NODE1=(i-2)*(N_SEG+1)+j-1
		NODE2=(i-2)*(N_SEG+1)+j
		TRIANGLES(S1,1)=NODE0
		TRIANGLES(S1,2)=NODE1
		TRIANGLES(S1,3)=NODE2
		S1=S1+1	
	ENDDO
  ENDDO

!DO i=1,2*(N_ROD-1)*N_SEG
!	PRINT *, TRIANGLES(i,:)
!ENDDO
!STOP

END SUBROUTINE INIT_SYSTEM
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!System perturbation
SUBROUTINE PERTURB_COORDINATES()
  IMPLICIT NONE
  INTEGER :: i, j, CUR 

  !This is where we can make initial changes to the angles from the initial input angles X
  ANGLES=X

  DO i=1,N_ROD
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		IF (j <= N_SEG-25) THEN
			!ANGLES(CUR) = 0.00
		ENDIF
	ENDDO
  ENDDO

  X=ANGLES

OPEN(1,FILE='coords.perturb')
WRITE(1,'(F20.10)') ANGLES
CLOSE(1)



END SUBROUTINE PERTURB_COORDINATES

!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!Compute the energy, energy derivatives and gradient
Subroutine IMPLEMENT_POTENTIAL(X,G,E,GTEST) 
  IMPLICIT NONE
  LOGICAL GTEST
  DOUBLE PRECISION :: X(N), E, G(N)
  
  ANGLES=X
  CALL UNPACK_COORDINATES()

  E = 0.0
  GRAD = 0.0

  CALL COMPUTE_ENERGY(E)

  IF (GTEST) THEN
   	CALL COMPUTE_GRAD()
  ENDIF

  G = GRAD * CM1
  E = E * CM1

  !CALL TEST_GRADIENT(G)


END SUBROUTINE IMPLEMENT_POTENTIAL
!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
!Unpack the angles to form the vertex coordinates  
SUBROUTINE UNPACK_COORDINATES()
  IMPLICIT NONE

  !Coordinates of 1st and 2nd vertex of each rod
  DO i=1,N_ROD
    COORDS((i-1)*(N_SEG+1)+1,1) = COORDS_0(i,1);
    COORDS((i-1)*(N_SEG+1)+1,2) = COORDS_0(i,2); 
    
    COORDS((i-1)*(N_SEG+1)+2,1) = COORDS((i-1)*(N_SEG+1)+1,1) + LS;
    COORDS((i-1)*(N_SEG+1)+2,2) = COORDS((i-1)*(N_SEG+1)+1,2); 
  ENDDO

  !Coordinates of remaining vertices of rod 1 to (N_ROD-1)
  DO i=1,N_ROD
    DO j=3,N_SEG+1
        COORDS((i-1)*(N_SEG+1)+j,1) = COORDS((i-1)*(N_SEG+1)+j-1,1) + LS*cos(ANGLES((i-1)*(N_SEG-1) + j-2));
        COORDS((i-1)*(N_SEG+1)+j,2) = COORDS((i-1)*(N_SEG+1)+j-1,2) + LS*sin(ANGLES((i-1)*(N_SEG-1) + j-2));
    ENDDO
  ENDDO

END SUBROUTINE UNPACK_COORDINATES
!------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------
!Energy computation
SUBROUTINE COMPUTE_ENERGY(E)
  IMPLICIT NONE
  DOUBLE PRECISION :: epsilon, sigma, R_cutoff, potential
  DOUBLE PRECISION :: E, EB, ES, EY, ELI, ELS, EC, EAREA, ETEMP, COORD1(2), COORD2(2), R, NORMAL(2), TEMP
  DOUBLE PRECISION :: DIFF_0, COS_DIFF, VERTDIF, X0,XP,XM,Y0,YP,YM,AREA 
  DOUBLE PRECISION :: VNODE0(2), VNODE1(2), VNODE2(2)
  INTEGER :: CUR, NODE1, NODE2, BONDTYPE, CURP, CURM, CUR2P, CUR2M, CUR2

  E=0.0
  EB=0.0
  ES=0.0
  ELI=0.0
  ELS=0.0
  EC=0.0
  EY=0.0
  EAREA=0.0


  !1) Bending energy
  !1.1) Interior rods
  DO i=2,N_ROD-1
	DO j=2,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		EB = EB + B_ARRAY(i,j)*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2
	ENDDO
  ENDDO
  EB = 2*EB/LS

  !1.2) Exterior rod 1
  DO j=2,N_SEG-1
	CUR=j
	EB = EB + (2*BE/LS)*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2
  ENDDO
  
  DO j=2,N_SEG-1
	CUR=(N_ROD-1)*(N_SEG-1)+j
	EB = EB +(2*BE/LS)*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2
  ENDDO

  
  !2) Stretching energy
  DO i=1,N_BONDS
	NODE1=BONDS(i,1)
	NODE2=BONDS(i,2)
	BONDTYPE=BONDS(i,3)

	R=SQRT( (COORDS(NODE1,1)-COORDS(NODE2,1))**2 + (COORDS(NODE1,2)-COORDS(NODE2,2))**2 )

	!With 1/R term
	!ES = ES + K_ARRAY(BONDTYPE)*(R-R_ARRAY(BONDTYPE))**2*(1.0/R)
	!Without 1/R term
	!ES = ES + K_ARRAY(BONDTYPE)*((R-R_ARRAY(BONDTYPE))/R_ARRAY(BONDTYPE))**2
  !ENDDO
  !ES = 0.5*ES
  

  epsilon = 0.01
  sigma = R_ARRAY(BONDTYPE)/ (2.0**(1.0/6.0))
  !R_cutoff = 2.5D0 * sigma
  !Apply cutoff
  !IF (R > R_cutoff) CYCLE
  !Lennard-Jones potential energy 
  potential = (4.0 * epsilon * ((sigma/R)**12.0 - (sigma/R)**6.0) + 1.0/4.0)
  ES= potential +ES
  !PRINT *,ES  

  ENDDO
  
  !3) Compression load
  DO i=1,N_ROD
	CUR=i*(N_SEG+1)
	ELI = ELI + COORDS(CUR,1)
  ENDDO
  ELI=ELI*F_I


  !4) pinching load
  ETEMP=0.0
  DO j=2,N_SEG+1
	CUR=j
	ETEMP = ETEMP + COORDS(CUR,2)-COORDS(1,2)
  ENDDO
  ELS = ELS-ETEMP*F_S

  ETEMP=0.0
  DO j=2,N_SEG+1
  	CUR=(N_ROD-1)*(N_SEG+1)+j
  	ETEMP = ETEMP + COORDS(CUR,2)-COORDS((N_ROD-1)*(N_SEG+1)+1,2)
  ENDDO
  ELS = ELS+ETEMP*F_S


  !5) Angular soft constraints of rod ends
  !5.1) TWO ENDS SHOULD BE ALIGNED
  DO i=1,N_ROD
  CUR=(N_ROD-1)*(N_SEG-1)+j
  VERTDIF =COORDS((i-1)*(N_SEG+1)+2,2)-COORDS((i-1)*(N_SEG+1)+N_SEG+1,2)
  EC = EC + 0.5*KC*VERTDIF*VERTDIF
  ENDDO

	 
  !5.2) END ANGLES CONSTRAINED (PINNED)
  DO i=1,N_ROD
  CUR=(i-1)*(N_SEG-1)+1
  EC = EC + 0.5*KC*(1-COS(ANGLES(CUR)))*(1-COS(ANGLES(CUR)))
  ENDDO

  DO i=1,N_ROD
  CUR=(i-1)*(N_SEG-1)+N_SEG-1
  EC = EC + 0.5*KC*(1-COS(ANGLES(CUR)))*(1-COS(ANGLES(CUR)))
  ENDDO

  !6) Energy penalising a node on a rod below becoming above the same node on the rod above
  !DO i=1,N_ROD-1
	!DO j=1,N_SEG+1
	!	CUR = (i-1)*(N_SEG+1) + j
	!	CURP = (i)*(N_SEG+1) + j
	!	TEMP= 0.0
	!	TEMP = ( 1.0 - TANH((COORDS(CURP,2)-COORDS(CUR,2))/WY) )
	!	TEMP = TEMP*( 2.0 - (COORDS(CURP,2)-COORDS(CUR,2)) )
	!	TEMP = TEMP * KY * 0.5
		
	!	EY = EY+TEMP

	!ENDDO
  !ENDDO


  !6) Repulsion term that prevents the node triangles having <0 area
  DO i=1,2*(N_ROD-1)*N_SEG
	X0=COORDS(TRIANGLES(i,1),1)
	XP=COORDS(TRIANGLES(i,2),1)
	XM=COORDS(TRIANGLES(i,3),1)
	Y0=COORDS(TRIANGLES(i,1),2)
	YP=COORDS(TRIANGLES(i,2),2)
	YM=COORDS(TRIANGLES(i,3),2)

	AREA=0.5*((XP-X0)*(YM-Y0)-(XM-X0)*(YP-Y0))

	EAREA=EAREA+KAREA*EXP(-CAREA*AREA/AREA0)
  ENDDO


E = EB + ES + ELI + ELS + EC + EAREA

!PRINT *, EB,ES, ELI,ELS,EC,EAREA,E,ELI+EB+ES+ELS

END SUBROUTINE COMPUTE_ENERGY
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!Gradient computation
SUBROUTINE COMPUTE_GRAD()
  IMPLICIT NONE
  DOUBLE PRECISION :: THETAM, THETAP, TEMP, R, WEIGHTX, WEIGHTY
  DOUBLE PRECISION :: DIFF_0, DIFF_P, DIFF_M, DIFF_0_SQ, DIFF_P_SQ, DIFF_M_SQ, DIFF_0_Y, DIFF_P_Y, DIFF_M_Y
  DOUBLE PRECISION :: TEMP_P, TEMP_M, COS_DIFF, GTEMP
  DOUBLE PRECISION :: X0, XP, XM, Y0, YP, YM, AREA, VPRE
  INTEGER :: CUR, CURP, CURM,  NODE1, NODE2, BONDTYPE, NODE_J, NODE_K, END_NODE_K, CUR2, CUR2P, CUR2M, CUR2PP, CUR2MM, GSIGN

  GRAD=0.0
  GRAD_SPAT=0.0


 !1) Bending gradient
  !1.1) Internal rods, not including ends
  DO i=2,N_ROD-1
	DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*B_ARRAY(i,j)/LS)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )
	ENDDO 
  ENDDO

  !1.2) Internal rods, first segment (angle 1 = 0)
  DO i=2,N_ROD-1
	j=1
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*B_ARRAY(i,j)/LS)*( - TAN(THETAP)/(COS(THETAP))**2 )
  ENDDO

  !1.3) Internal rods, last segment (no N_SEG+1 contribution)
  DO i=2,N_ROD-1
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))

		GRAD(CUR) = GRAD(CUR) + (2*B_ARRAY(i,j)/LS)*( TAN(THETAM)/(COS(THETAM))**2 )
  ENDDO
		
  !1.4) External rods, not including ends
  i=1
	DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*BE/LS)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )
	ENDDO 
 
  i= N_ROD
	DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*BE/LS)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )
	ENDDO 

  !1.5) External rods, first segments (angle 1 = 0)
  i=1
	j=1
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*BE/LS)*( - TAN(THETAP)/(COS(THETAP))**2 )

  !1.4) External rods,  last segment (no N_SEG+1 contribution)
  i=1
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))

		GRAD(CUR) = GRAD(CUR) + (2*BE/LS)*( TAN(THETAM)/(COS(THETAM))**2 ) 
   
   i=N_ROD
	j=1
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD(CUR) = GRAD(CUR) + (2*BE/LS)*( - TAN(THETAP)/(COS(THETAP))**2 )

  !1.4) External rods,  last segment (no N_SEG+1 contribution)
  i=N_ROD
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))

		GRAD(CUR) = GRAD(CUR) + 2*BE/LS*( TAN(THETAM)/(COS(THETAM))**2 ) 



  !2) Stretching gradient
  !2.1) Start by computing the dE/dxi and dE/dyi
  DO i=1,N_BONDS
	NODE1=BONDS(i,1)
	NODE2=BONDS(i,2)
	BONDTYPE=BONDS(i,3)

	R=SQRT( (COORDS(NODE1,1)-COORDS(NODE2,1))**2 + (COORDS(NODE1,2)-COORDS(NODE2,2))**2 )
	
	! X-derivatives
	!With 1/R term
	!TEMP=K_ARRAY(BONDTYPE)*(COORDS(NODE1,1)-COORDS(NODE2,1))*(1.0-R_ARRAY(BONDTYPE)/R)*(1/R*(1-0.5*(1.0-R_ARRAY(BONDTYPE)/R)))
	!Without 1/R term
	TEMP=(K_ARRAY(BONDTYPE)/(R_ARRAY(BONDTYPE)**2))*(COORDS(NODE1,1)-COORDS(NODE2,1))*(1.0-R_ARRAY(BONDTYPE)/R)

	GRAD_SPAT(NODE1,1) = GRAD_SPAT(NODE1,1) + TEMP
	GRAD_SPAT(NODE2,1) = GRAD_SPAT(NODE2,1) - TEMP

	! Y-derivatives
	!With 1/R term
	!TEMP=(K_ARRAY(BONDTYPE)/R_ARRAY(BONDTYPE)**2)*(COORDS(NODE1,2)-COORDS(NODE2,2))*(1.0-R_ARRAY(BONDTYPE)/R)*(1/R*(1-0.5*(1.0-R_ARRAY(BONDTYPE)/R)))
	!Without 1/R term
	TEMP=(K_ARRAY(BONDTYPE)/(R_ARRAY(BONDTYPE)**2))*(COORDS(NODE1,2)-COORDS(NODE2,2))*(1.0-R_ARRAY(BONDTYPE)/R)

	GRAD_SPAT(NODE1,2) = GRAD_SPAT(NODE1,2) + TEMP
	GRAD_SPAT(NODE2,2) = GRAD_SPAT(NODE2,2) - TEMP
  ENDDO


  !2.2) Now convert these to dE/dtheta_i
  !2.2.1) Interior rods 
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		! 1D location of current angle
		CUR=(i-1)*(N_SEG-1)+j
	
		DO k=j+1,N_SEG
			! 1D location of a node which the current angle affects
			NODE_K=(i-1)*(N_SEG+1)+k+1
			
			GRAD(CUR) = GRAD(CUR) - LS*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LS*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
		ENDDO		
	ENDDO
  ENDDO


  !3) Interocular load
  !3.1)Internal rods
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		GRAD(CUR) = GRAD(CUR) - F_I*LS*SIN(ANGLES(CUR))
	ENDDO
  ENDDO


!  !4) Sclera load
  i=1
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		GRAD(CUR) = GRAD(CUR) - F_S*LS*COS(ANGLES(CUR))*(N_SEG-j)
	ENDDO
 
   i=N_ROD
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		GRAD(CUR) = GRAD(CUR) + F_S*LS*COS(ANGLES(CUR))*(N_SEG-j)
	ENDDO



!5)Angular constraint
  !5.1) End aligment
  DO i=1,N_ROD
    DO j=1,N_SEG-1
      CUR=(i-1)*(N_SEG-1)+j 
      GRAD(CUR) = GRAD(CUR) - KC*(COORDS((i-1)*(N_SEG+1)+2,2)-COORDS((i-1)*(N_SEG+1)+N_SEG+1,2))*LS*COS(ANGLES(CUR))
    ENDDO
  ENDDO

  !5.2) END ANGLES CONSTRAINED (PINNED)
  DO i=1,N_ROD
   CUR=(i-1)*(N_SEG-1)+1
   GRAD(CUR) =  GRAD(CUR)+KC*SIN(ANGLES(CUR))*(1-COS(ANGLES(CUR)))
  ENDDO

  DO i=1,N_ROD
   CUR=(i-1)*(N_SEG-1)+N_SEG-1
   GRAD(CUR) =  GRAD(CUR)+KC*SIN(ANGLES(CUR))*(1-COS(ANGLES(CUR)))
  ENDDO	

 
  !6) Repulsion term that prevents the node triangles having <0 area
  GRAD_SPAT=0.0
  DO i=1,2*(N_ROD-1)*N_SEG
	X0=COORDS(TRIANGLES(i,1),1)
	XP=COORDS(TRIANGLES(i,2),1)
	XM=COORDS(TRIANGLES(i,3),1)
	Y0=COORDS(TRIANGLES(i,1),2)
	YP=COORDS(TRIANGLES(i,2),2)
	YM=COORDS(TRIANGLES(i,3),2)

	AREA=0.5*((XP-X0)*(YM-Y0)-(XM-X0)*(YP-Y0))

	VPRE=EXP(-CAREA*AREA/AREA0)

	GRAD_SPAT(TRIANGLES(i,1),1)=GRAD_SPAT(TRIANGLES(i,1),1)+VPRE*(-YM+YP)
	GRAD_SPAT(TRIANGLES(i,2),1)=GRAD_SPAT(TRIANGLES(i,2),1)+VPRE*(YM-Y0)
	GRAD_SPAT(TRIANGLES(i,3),1)=GRAD_SPAT(TRIANGLES(i,3),1)+VPRE*(-YP+Y0)

	GRAD_SPAT(TRIANGLES(i,1),2)=GRAD_SPAT(TRIANGLES(i,1),2)-VPRE*(-XM+XP)
	GRAD_SPAT(TRIANGLES(i,2),2)=GRAD_SPAT(TRIANGLES(i,2),2)-VPRE*(XM-X0)
	GRAD_SPAT(TRIANGLES(i,3),2)=GRAD_SPAT(TRIANGLES(i,3),2)-VPRE*(-XP+X0)
  ENDDO
  GRAD_SPAT=GRAD_SPAT*(-CAREA*KAREA/AREA0*0.5)


  !6.2) Now convert these to dE/dtheta_i
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		! 1D location of current angle
		CUR=(i-1)*(N_SEG-1)+j
		DO k=j+1,N_SEG
			! 1D location of a node which the current angle affects
			NODE_K=(i-1)*(N_SEG+1)+k+1
			GRAD(CUR) = GRAD(CUR) - LS*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LS*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
		ENDDO		
	ENDDO
  ENDDO

!open(2, file='grad.out')
!write(2,'(F20.10)') GRAD
!close(2)
!stop
END SUBROUTINE COMPUTE_GRAD
!------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------
!Test the gradient vs a finite difference approximation
SUBROUTINE TEST_GRADIENT(G)
IMPLICIT NONE
INTEGER :: I
DOUBLE PRECISION :: DT, A_0, A_P, A_M, E_P, E_M, TEST_GRAD(N), GRAD_DIFF(N), G(N)

DT=0.000001

DO I=1,N
	A_0 = ANGLES(I)
	A_P = ANGLES(I) + DT
	A_M = ANGLES(I) - DT

	E=0.0

	ANGLES(I) = A_P
	CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
	E_P=E
	ANGLES(I)=A_M
	CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
	E_M=E

	TEST_GRAD(I) = (E_P-E_M)/(2*DT)*CM1	
	
	ANGLES(I) = A_0

	PRINT *, G(I), TEST_GRAD(I), G(I)-TEST_GRAD(I)
	IF (ABS(G(I)-TEST_GRAD(I))>1.0D-5) THEN
		PRINT *, 'TEST GRAD FAILED AT I=', I
	ENDIF
ENDDO

GRAD_DIFF = G-TEST_GRAD



OPEN(1,FILE='grad.diff')
WRITE(1,'(F20.10)') GRAD_DIFF
CLOSE(1)

OPEN(1,FILE='grad.orig')
WRITE(1,'(F20.10)') G
CLOSE(1)

OPEN(1,FILE='grad.test')
WRITE(1,'(F20.10)') TEST_GRAD
CLOSE(1)

STOP

END SUBROUTINE TEST_GRADIENT



END MODULE potential
