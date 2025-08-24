!lamina_new_clamped
subroutine set_defaults()
    use potential
    !#################!
    ! User Parameters !
    !#################!
  
    CM1 = 1.0D3 !Energy scale factor
   
    OPEN(1,FILE='data.in')
    READ(1,*) N_ROD
    READ(1,*) N_SEG
    READ(1,*) L
    READ(1,*) H
    READ(1,*) K3
    READ(1,*) LOWER_RAT
    READ(1,*) F_I
    READ(1,*) F_S
    READ(1,*) BE
    READ(1,*) RAND_MAG
    CLOSE(1)
  
  
  
    R3=H/(N_ROD-1.0) !LC spring equilibrium length
    
    K1=K3 !Vertical spring stiffness
    R1=H/(N_ROD-1.0) !Vertical spring equilibrium length
  
    K2=K3  !Diagonal spring stiffness
    R2=SQRT( (L/N_SEG)**2 + (H/(N_ROD-1.0))**2 )  !Diagonal spring equilibrium length
    
  ! Allocate a 2D array for bending stiffness: one row per filament, one column per bending segment
    ALLOCATE(B_ARRAY(N_ROD, N_SEG-1))
    B_ARRAY(1, :) = BE
    B_ARRAY(N_ROD, :) = BE
    
    OPEN(1, FILE='b_array.in', STATUS='OLD', ACTION='READ', IOSTAT=io_status)
    IF (io_status /= 0) THEN
        PRINT *, "Error opening b_array.in"
        STOP
    END IF
    
    DO i = 2, N_ROD-1
        READ(1, *) (B_ARRAY(i, j), j = 1, N_SEG-1)
    END DO
    CLOSE(1)
    

  
   
  
    KC=100000.0 !End constraint magnitude
    KY=10.0 !Overlap penalisation magnitude initially 10.0
    WY=0.001 !Overlap penalisation onset width
  
    KAREA=1.0D6 !Repulsion magnitude
    CAREA=40.0 !Repulsion exponential factor
  
  
    N = (N_ROD)*(N_SEG-1)
    !###################!
    ! L-BFGS Parameters !
    !###################!
  
    M = 10
    max_iterations  = 1000000 !Max number of minimisation steps
    convergence_rms = 1d-5 !RMS gradient at convergence
    max_step_size   = 0.25d0
    H0init          = 1d-1
    relative_energy_check = .true.
    dE_max          = 1d-1
    debug = .true.
  
    !##########################!
    ! Basin Hopping Parameters !
    !##########################!
  
    max_runs          = 1
    max_similarity    = 1d-4
    temperature       = 1000d4
    binary_io         = .false.
    minima_save_limit = 0
  end subroutine set_defaults
  