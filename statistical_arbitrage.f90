!==============================================================================
! Name: Anthony Hills
! URN: 6317482
! PHY3048
!
! This program analyses the long term profitability and risk of the strategy
! of a hedge fund.
!
! The Black-Scholes model is used to simulate the share price over time, and
! the portfolio's performance is evaluated by measuring its mean, variance,
! and probability of a loss as a function of time.
!=============================================================================
PROGRAM statistical_arbitrage
  IMPLICIT NONE

  ! Parameters
  INTEGER, PARAMETER :: realisations = 1000  ! number of simulations
  REAL, PARAMETER :: mu = 0.16  ! drift
  REAL, PARAMETER :: sigma = 0.20  ! volatility
  REAL, PARAMETER :: S0 = 100.  ! share price at t=0
  REAL, PARAMETER :: dt = (1.0/250.)  ! time step
  REAL, PARAMETER :: r = 0.04  ! risk-free rate
  REAL, PARAMETER :: k = 0.2
  REAL, PARAMETER :: pi = 4.*ATAN(1.)

  ! Helpers
  REAL :: x1, x2, z1, z2, dW, dS, S, t, mean, variance, prob_loss
  INTEGER :: realisation, trading_days, i_seed, years, day, losses, i
  REAL, DIMENSION(:,:), ALLOCATABLE :: portfolio

! ----------------------------------------------------------------------------
  ! Allocate the portfolio array using the user input number of trading days
  WRITE(6,*) 'How many years would you like to simulate for?'
  READ(5,*) years
  trading_days = years / dt
  ALLOCATE(portfolio(trading_days, realisations))

  ! Open files to be written to
  OPEN(10,file='share_price_simulation.dat')
  OPEN(20,file='portfolio_value_histogram.dat')
  OPEN(30,file='probability_of_loss.dat')
  OPEN(40,file='portfolio_mean_variance.dat')
  OPEN(50,file='box_muller.dat')

! ----------------------------------------------------------------------------
  ! Simulate value of portfolio for chosen number of realisations
  portfolio(:,:) = 0.;  variance = 0.  ! Initialize variables
  DO realisation=1, realisations
     S = S0; t = 0.0  ! Reinitialize initial share price at t=0
     DO day=0, trading_days  ! Evaluate value for each day
        t = t + dt
        dW = SQRT(dt) * rand_box_muller()  ! Wiener process
        dS = S * (mu * dt + sigma * dW)
        S = S + dS  ! Update current share price
        IF (realisation == 1) THEN
           WRITE(10,*) t, S  ! Save simulated share price using 1 realisation
        END IF

        ! Define hedging strategy
        portfolio(day,realisation) = S * EXP(-r * t) - S0
        IF (S >= S0 * (1. + k) * EXP(r * t)) THEN
           portfolio(day,realisation) = S0 * k
        END IF
        WRITE(20,*) portfolio(day,realisation)  ! Save portfolio values to file
      END DO
   END DO

! ----------------------------------------------------------------------------
   ! Evaluate portfolio's mean, variance and loss probability with time
   t = 0.;
   DO day=0, trading_days
      losses = 0.; variance = 0.; mean = 0.  ! Initialize variables
      t = t + dt
      mean = SUM(portfolio(day,:)) / realisations  ! Calculte mean
      DO realisation=1, realisations
          variance = variance + (portfolio(day,realisation) - mean) ** 2
        IF (portfolio(day, realisation) < 0.) THEN
          losses = losses + 1
        END IF
        prob_loss = REAL(losses) / realisations  ! Probability of a loss
      END DO
      WRITE(30,*) prob_loss, t  ! Save probability of a loss to file
      variance = variance / (realisations-1)
      WRITE(40,*) mean, variance, t  ! Save portfolio statistics to file
  END DO

! ----------------------------------------------------------------------------
  DO i=1, 10000
    WRITE(50,*) rand_box_muller()  ! Save random numbers with mean 0 and var 1
  END DO

  CLOSE(10); CLOSE(20); CLOSE(30); CLOSE(40); CLOSE(50)  ! Finish output

  ! Display information to user
  WRITE(6,*) 'A stochastic simulation for the share price has been saved to:'
  WRITE(6,*) "                 share_price_simulation.dat               "
  WRITE(6,*) ''
  WRITE(6,*) "The portfolio value for 1000 realisations and for your chosen:"
  WRITE(6,*) "number of years has been saved to:"
  WRITE(6,*) "                portfolio_value_histogram.dat "
  WRITE(6,*) ''
  WRITE(6,*) "The probability of a loss as a function of time was saved to:"
  WRITE(6,*) "                  probability_of_loss.dat"
  WRITE(6,*) ""
  WRITE(6,*) "The mean and variance for the value of the portfolio as a "
  WRITE(6,*) "function of time has been saved to:"
  WRITE(6,*) "                 portfolio_mean_variance.dat "
  WRITE(6,*) ""

!=============================== Functions ====================================
CONTAINS
  ! Generates a normally distributed random number with mean zero and variance
  ! one from the Box-Muller transform
  REAL FUNCTION rand_box_muller()
    REAL :: x1, x2, z1, z2
    REAL, PARAMETER :: pi = 4.*ATAN(1.)
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
