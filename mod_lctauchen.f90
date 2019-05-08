MODULE mod_lctauchen
    IMPLICIT NONE
    
    CONTAINS

SUBROUTINE lctauchen(rho,sigma_eps,omega,y_grid,trans)

! ========================================================================
!  Code for discretization method in "An Extension of Rouwenhorst Method
!  for Non-Stationary Environments" by Giulio Fella, Giovanni Gallipoli
!  and Jutong Pan
! 
!  Tauchen method to approximate non-stationary AR(1) process by a discrete Markov chain
!       y(t) = rho(t)*y(t-1)+ epsilon(t),   epsilion(t)~iid N(0,sigmaeps(t))
!       with INITIAL condition y(0) = 0 (equivalently y(1)=epsilon(1) ) 
!  INPUT:  rho 	     - Tx1 vector of serial correlation coefficients
!          sigma_eps - Tx1 vector of standard deviations of innovations
!          omega     - Tx1 vector of scaling factors for grid range
!                      (equal 3 in Tauchen (1986) paper)
!  OUTPUT: trans     - NxNxT matrix of T (NxN) transition matrices  
!                      transition probabilities are arranged by row
!          y_grid    - an NxT matrix, each column stores the Markov state 
! 	               space for period t
! ========================================================================%

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  REAL(dp), INTENT(IN):: rho(:), sigma_eps(:), omega(:)
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: trans
  REAL(dp), DIMENSION(:,:), INTENT(OUT) :: y_grid
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: temp3d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: trans0,cdf
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: temp1d, sigma, h
  INTEGER :: N        ! number of income states
  INTEGER :: LIFET    ! life time
  INTEGER :: i,j,t      ! index

    N = SIZE(y_grid,dim=1)
    LIFET = SIZE(y_grid,dim=2)
    ALLOCATE( sigma(LIFET), h(LIFET), temp1d(N), temp3d(N,N,LIFET),&
        & cdf(N,N) )


    ! *** Step 1: construct the state space y_grid(t) in each period t,
    ! Evenly-spaced N-state space over [-omega(t)*sigma_y(t),omega(t)*sigma_y(t)].

    ! 1. Compute unconditional variances of y(t)

    ! Part 1: obtain the initial distribution
    ! initial distribution is the stationary distribution of the Markov chain for rho=0
    ! when rho=0, the stationary distribution is any row of the transition matrix
    sigma(1) = sigma_eps(1)
    DO t=2,LIFET
      sigma(t) = SQRT(rho(t)**2*sigma(t-1)**2+sigma_eps(t)**2)
    END DO
    h(1:LIFET) = 2*omega(1:LIFET)*sigma(1:LIFET)/REAL(N-1,KIND=dp)  
    ! Construct state space
    DO t=1,LIFET
        DO i = 1,N
         y_grid(i,t) = -omega(t)*sigma(t) + REAL(i-1,KIND=dp)*h(t)
        END DO
    END DO

    !  *** Step 2: Compute the transition matrices trans(:,:,t) from
    !              period (t-1) to period t

    ! Compute the transition matrix in period 1; i.e., from y(0)=0
    ! to any gridpoint of y_grid(1) in period 1.
    ! Any of its rows is the (unconditional) distribution in period 1.
    DO i=1,N
      temp1d = (y_grid(:,1)-h(1)/2)/sigma_eps(1);
      temp1d = MAX(temp1d,-37.d0); ! To avoid underflow in next line
      cdf(i,:) = cdf_normal(temp1d);
    END DO
    trans(:,1,1) = cdf(:,2); 
    trans(:,N,1) = 1-cdf(:,N);
    DO i=2,N-1
        trans(:,i,1) = cdf(:,i+1)-cdf(:,i)
    END DO
    
    ! Compute the transition matrices for t>2
    DO t=2,LIFET
        DO i=1,N
          temp3d(i,:,t) = (y_grid(:,t) - rho(t)*y_grid(i,t-1) - h(t)/2)/sigma_eps(t)
          ! The next line is to avoid an underflow when computing the normal cdf in
          ! the following line
          temp3d(i,:,t) = MAX(temp3d(i,:,t),-37.d0)
          cdf(i,:) = cdf_normal(temp3d(i,:,t))
        END DO
        trans(:,1,t) = cdf(:,2)   ! note that here trans(:,:,t) represents the trans matrix from t to t+1
        trans(:,N,t) = 1-cdf(:,N)
        DO j=2,N-1
        trans(:,j,t) = cdf(:,j+1)-cdf(:,j)
        END DO
    END DO
    
END SUBROUTINE lctauchen

FUNCTION cdf_normal(X)
  ! Returns the value of the cdf of the Standard Normal distribution at point X

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), INTENT(in), DIMENSION(:)  :: X
    REAL(dp), DIMENSION(SIZE(X)) :: cdf_normal
 
    cdf_normal =  0.5d0*( 1+ERF(x/SQRT(2.d0)) )  

  END FUNCTION cdf_normal

END MODULE mod_lctauchen
