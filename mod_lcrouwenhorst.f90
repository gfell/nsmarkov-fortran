MODULE mod_lcrouwenhorst
    IMPLICIT NONE
    
    CONTAINS

      SUBROUTINE lcrouwenhorst(rho,sigma_eps,y_grid,trans)

! ========================================================================
!  Code for discretization method in "An Extension of Rouwenhorst Method
!  for Non-Stationary Environments" by Giulio Fella, Giovanni Gallipoli
!  and Jutong Pan
! 
!  Rouwenhurst method to approximate non-stationary AR(1) process by a discrete Markov chain
!       y(t) = rho(t)*y(t-1)+ epsilon(t),   epsilion(t)~iid N(0,sigmaeps(t))
!       with INITIAL condition y(0) = 0 (equivalently y(1)=epsilon(1) ) 
!  INPUT:  rho - Tx1 vector of serial correlation coefficients
!          sigma_eps - Tx1 vector of standard deviations of innovations
!  OUTPUT: trans  - NxNxT array of T (NxN) transition matrices  
!                   transition probabilities are arranged by row
!          y_grid - NxT matrix, each column stores the Markov state 
! 	   space for period t
! ========================================================================%
        
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), INTENT(in) :: sigma_eps(:)
    REAL(dp), INTENT(in) :: rho(:)
    REAL(dp), DIMENSION(:,:), INTENT(out) ::  y_grid
    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: trans   
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: trans0
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: sigma    ! SD of y(t) 
    REAL(dp) :: den,p,range
    INTEGER :: N        ! number of income states
    INTEGER :: LIFET    ! life time
    INTEGER :: t,i     ! index

    LIFET = SIZE(y_grid,dim=2)     ! number of periods
    N    = SIZE(y_grid,dim=1)     ! number of states
    den = real(N-1,KIND=dp) 
    
    ALLOCATE(trans0(N,N), sigma(LIFET))

! *** Step 1: construct the state space y_grid(t) in each period t.
! Evenly-spaced N-state space over [-sqrt(N-1)*sigma_y(t),sqrt(N-1)*sigma_y(t)].

    ! 1.a Compute unconditional variances of y(t)
    sigma(1) = sigma_eps(1)
    DO t=2,LIFET
        sigma(t) = SQRT(rho(t)**2*sigma(t-1)**2+sigma_eps(t)**2)
    END DO

    ! 1.b Construct state space
    DO t=1,LIFET
       range = 2*SQRT(REAL(N-1,dp))*sigma(t) ! range of grid points
       DO i = 1,N
         y_grid(i,t) = -SQRT(REAL(N-1,dp))*sigma(t) + REAL(i-1,KIND=dp)/den*range
       END DO
    END DO

!  *** Step 2: Compute the transition matrices trans(:,:,t) from
!              period (t-1) to period t
! The transition matrix for period t is defined by parameter p(t).
! p(t) = 0.5*(1+rho*sigma(t-1)/sigma(t))

! Note: trans(:,:,1) is the transition matrix from y(0)=0 to any
! gridpoint of y_grid(1) in period 1.
! Any of its rows is the (unconditional) distribution in period 1.
    
    p = 0.5d0           ! First period: p(1) = 0.5 as y(1) is white noise.

    CALL rhmat(trans(:,:,1))

    ! Step 2(b): get the transition matrices for t>1
    DO t=2,LIFET
        p = (sigma(t)+rho(t)*sigma(t-1))/(2*sigma(t))
        CALL rhmat(trans(:,:,t))  ! trans(:,:,t) is the transition matrix from t-1 to t
    END DO

    CONTAINS
    SUBROUTINE rhmat(tr)
        IMPLICIT NONE
        REAL(dp), DIMENSION(:,:), INTENT(OUT) :: tr
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: tr0,tr1
        INTEGER :: h
    
        IF (size(tr,dim=1)==2) THEN
            tr = reshape((/p,1-p,1-p,p/),(/2,2/))
        ELSE
            ALLOCATE(tr0(2,2))
            tr0 = reshape((/p,1-p,1-p,p/),(/2,2/))
            DO h = 3,size(tr,dim=1)
                ALLOCATE(tr1(h,h))
                tr1=0
                tr1(1:h-1,1:h-1)=p*tr0
                tr1(1:h-1,2:h)=tr1(1:h-1,2:h)+(1-p)*tr0
                tr1(2:h,1:h-1)=tr1(2:h,1:h-1)+(1-p)*tr0
                tr1(2:h,2:h)=tr1(2:h,2:h)+p*tr0
                tr1(2:h-1,:)=tr1(2:h-1,:)/2
                DEALLOCATE(tr0)
                ALLOCATE(tr0(h,h))
                tr0 = tr1
                DEALLOCATE(tr1)
            ENDDO
            tr = tr0
            DEALLOCATE(tr0)
        END IF
    END SUBROUTINE rhmat

END SUBROUTINE lcrouwenhorst


END MODULE
