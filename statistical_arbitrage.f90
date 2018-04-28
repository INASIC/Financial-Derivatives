!==============================================================================
! Name: Anthony Hills
! URN: 6317482
! PHY30
!
! This program analyses the long term profitability and risk of the strategy
! of a hedge fund.
! We use the Black-Scholes model to simulate the share price over time. ...
!=============================================================================
PROGRAM statistical_arbitrage
  IMPLICIT NONE

  ! Parameters
  REAL, PARAMETER :: mu = 0.16  ! annual drift
  REAL, PARAMETER :: sigma = 0.2  ! volatility
  REAL, PARAMETER :: r = 0.04  ! risk-free rate
  REAL, PARAMETER  :: k = 0.2  !

  ! Helpers
  INTEGER :: i, i_final, realisation
  REAL :: present_value, mean_v, varaince_v
  REAL :: t_final, dt, t
  REAL :: dW, N
  REAL :: share_price_today, last_price, portfolio_today
  REAL :: S, S0, dS, p

  REAL, DIMENSION(:,:), ALLOCATABLE :: port


! Plot a histogram of the random numbers generated from the Box-Muller transform
OPEN(unit=100, file='random_numbers_box_muller.dat')
i_final = 10000
DO i=1, i_final
  WRITE(100,*) rand_box_muller()
END DO
CLOSE(100)

! ! Simulate share prices for 1 year
! OPEN(unit=10, file='daily_simulated_share_prices.dat')
! OPEN(unit=20, file='portfolio.dat')
! OPEN(unit=30, file='rand.dat')
! t = 0.; t_final = 1.; dt = 1./250.
! S = 100.
! WRITE(10,*) t*250, S
! DO WHILE (t <= t_final)
!   dW = sqrt(dt) * rand_box_muller()
!   dS = 100. * (mu * dt + sigma * dW)
!   S = S + dS
!   t = t + dt
!   WRITE(10,*) t*250, S
!   WRITE(30,*) rand_box_muller()
! END DO
! CLOSE(10)
! CLOSE(20)

ALLOCATE(port(1:250, 1:3))
DO realisation=1, 1000  !
  WRITE(6,*) (realisation/1000.) * 100, '%'
  ! WRITE(6,*) realisation/1000.
  OPEN(unit=10, file='daily_simulated_share_prices.dat')
  OPEN(unit=20, file='portfolio.dat')
  OPEN(unit=30, file='rand.dat')
  t = 0.; t_final = 1.; dt = 1./250.
  S0 = 100.
  S = S0
  WRITE(10,*) t*250, S
  DO WHILE (t <= t_final)
    dW = sqrt(dt) * rand_box_muller()
    dS = S0 * (mu * dt + sigma * dW)
    S = S + dS
    t = t + dt

    IF (S >= S0 * (1. + k) * exp(r*t)) THEN
      p = S0 * k
    ELSE
      p = S * exp(-r*t) - S0
    END IF

    WRITE(20,*) p
  END DO
END DO
CLOSE(10)
CLOSE(20)







! =============================== Functions ===================================
CONTAINS
!===========================================================================
! (b) Function that generates a normally distributed random number with
! mean zero and variance one from the Box-Muller transform
!===========================================================================
REAL FUNCTION rand_box_muller()
  REAL :: x1, x2, z1, z2
  REAL :: pi = 3.1415926535897932384626433832795028841971693993751058209749445
  INTEGER :: i

  ! ----- variables for portable seed setting -----
  INTEGER :: i_seed
  INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----
  REAL :: r

  DO i=1, 2  ! Loop twice
    ! ----- Set up random seed portably -----
    CALL RANDOM_SEED(size=i_seed)
    ALLOCATE(a_seed(1:i_seed))
    CALL RANDOM_SEED(get=a_seed)
    CALL DATE_AND_TIME(values=dt_seed)
    a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    CALL RANDOM_SEED(put=a_seed)
    DEALLOCATE(a_seed)
    ! ----- Done setting up random seed -----

    IF (i == 1) CALL RANDOM_NUMBER(x1);
    IF (i == 2) CALL RANDOM_NUMBER(x2);
  END DO

  ! Generate normally distributed numbers, according to Box-Muller transform
  z1 = sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2)
  z2 = sqrt(-2.0 * log(x1)) * sin (2.0 * pi * x2)

  rand_box_muller = z1
END FUNCTION
END PROGRAM statistical_arbitrage
