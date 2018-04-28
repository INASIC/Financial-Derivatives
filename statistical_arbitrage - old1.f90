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

  ! Parameters
  REAL, PARAMETER :: mu = 0.16  ! annual drift
  REAL, PARAMETER :: sigma = 0.2  ! volatility
  REAL, PARAMETER :: r = 0.04  ! risk-free rate
  REAL, PARAMETER  :: k = 0.2  !

  ! Helpers
  INTEGER :: i, i_final, seed1, seed2, seed
  REAL :: present_value, mean_v, varaince_v
  REAL :: t_final, dt, t
  REAL :: dw, N
  REAL :: delta_share, s0, share_price_today, last_price, portfolio_today


! Plot a histogram of the random numbers generated from the Box-Muller transform
OPEN(unit=100, file='random_numbers_box_muller.dat')
i_final = 10000
DO i=1, i_final
  seed1 = i; seed2 = i_final + i
  WRITE(100,*) rand_box_muller(seed1, seed2)
END DO
CLOSE(100)

! ! Simulate daily share price over one year
! t = 0.
! t_final = 1.
! ! dt = (1./250.)
! dt = t_final / 1000
! share_price_today = 100.
! OPEN(unit=10, file='daily_simulated_share_prices.dat')
! OPEN(unit=20, file='dW.dat')
! WRITE(10,*) t*250, share_price_today
! seed1 = 1; seed2 = INT(seed1 + t_final * 250)  ! Initialize seeds
! DO WHILE (t <= t_final)
!   WRITE(6,*) t/t_final
!   seed1 = INT(seed1 + t*250.); seed2 = INT(seed2 + t*250.)  ! Different seeds
!   last_price = share_price_today
!   N = rand_box_muller(seed1, seed2)  ! Set dW to random value from Gaussian
!   ! N = rand_box_muller()
!   dw = sqrt(dt) * N
!   delta_share = 100. * (mu * dt + sigma * dw)
!   t = t + dt
!   share_price_today = last_price + delta_share  ! Evaluate price tomorrow
!   WRITE(10,*) t*250, share_price_today
!   WRITE(20,*) dw
! END DO
! CLOSE(10)
! CLOSE(20)


! Simulate daily share price over one year
t = 0.
t_final = 1.
! dt = (1./250.)
dt = t_final / 1000
share_price_today = 100.
OPEN(unit=10, file='daily_simulated_share_prices.dat')
OPEN(unit=20, file='dW.dat')
OPEN(unit=30, file='portfolio.dat')
WRITE(10,*) t*250, share_price_today
portfolio_today = share_price_today
WRITE(30,*) portfolio_today
seed1 = 1; seed2 = INT(seed1 + t_final * 250)  ! Initialize seeds
DO WHILE (t <= t_final)
  WRITE(6,*) t/t_final
  seed1 = INT(seed1 + t*250.); seed2 = INT(seed2 + t*250.)  ! Different seeds
  last_price = share_price_today
  N = rand_box_muller(seed1, seed2)  ! Set dW to random value from Gaussian
  ! N = rand_box_muller()
  dw = sqrt(dt) * N
  delta_share = 100. * (mu * dt + sigma * dw)
  t = t + dt
  share_price_today = last_price + delta_share  ! Evaluate price tomorrow
  ! IF (t <= t_star) THEN
  IF (share_price_today >= 100.*(1+k)*exp(r*t)) THEN
    portfolio_today = share_price_today*exp(-r*t) - 100.
  ELSE
    portfolio_today = 100. * k
  END IF
  WRITE(10,*) t*250, share_price_today
  WRITE(20,*) dw
  WRITE(30,*) portfolio_today
END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)




WRITE(6,*) "Answers to the coursework:"
WRITE(6,*) "(a)", present_value, " = present value of portfolio "
WRITE(6,*) ""
WRITE(6,*) "(c) Simulated daily share prices for one year have been saved to:"
WRITE(6,*) "                     daily_share_price_c.dat"
WRITE(6,*) ""
WRITE(6,*) "(d) Value of the portfolio for 1, 2, and 5 years using 1000 "
WRITE(6,*) "realisations of the stock price have been saved to:"
WRITE(6,*) "                    hist_portfolio_value.dat"
WRITE(6,*) ""
WRITE(6,*) "(e) mean =", mean_v
WRITE(6,*) "variance = ", varaince_v
WRITE(6,*) ""
WRITE(6,*) "(f) Probability of a loss as a function of time from t=0 to 20"
WRITE(6,*) "years has been saved to:"
WRITE(6,*) "                       probability_of_loss.dat"
WRITE(6,*) ""
WRITE(6,*) "(g) For a volatility of 0.6, ..."


END PROGRAM statistical_arbitrage

!================================================================
! (a) Derive the present value of the portfolio given in Eq. (2)
!================================================================
REAL FUNCTION derive_present_value(t_present)
  ! Parameters
  INTEGER, PARAMETER :: t_invest_bond = 100  ! Time sell for risk free bond
  INTEGER, PARAMETER :: t_final = 1 * 365  ! Should it not be 250?
  INTEGER, PARAMETER :: t0 = 0  ! Initial day
  ! INTEGER, PARAMETER ::
  INTEGER :: dt = 1  ! Integration step
  REAL, PARAMETER :: r = 0.04  ! Risk free rate, LIBOR
  REAL, PARAMETER :: k = 0.2  !
  INTEGER :: t_present

  ! Helpers
  INTEGER :: t
  REAL, DIMENSION(t_final-t0) :: v, share_price, discounted_value

  share_price(1) = 100  ! Share price at t = 0

  ! Portfolio discounted value as a function of time (Eq. 2)
  IF (t <= t_invest_bond) THEN
    v(t) = share_price(t) * exp(-r*t) - share_price(1)
  ELSE
    v(t) = share_price(1) * k
  END IF

  ! Equation 2:
  derive_present_value = v(t_present)

END FUNCTION

!===========================================================================
! (b) Function that generates a normally distributed random number with
! mean zero and variance one from the Box-Muller transform
!===========================================================================
! REAL FUNCTION rand_box_muller(seed1, seed2)
  REAL FUNCTION rand_box_muller(seed1, seed2)
  INTEGER, INTENT(IN) :: seed1, seed2
  ! INTEGER :: seed1=2
  REAL :: x1, x2, z1, z2
  REAL :: pi = 3.1415926535897932384626433832795028841971693993751058209749445
  ! Intrinsic uniform random number
  ! CALL SRAND(seed1); x1 = RAND()
  ! CALL SRAND(seed2); x2 = RAND()

  ! CALL SRAND(5)
  !
  ! x1 = RAND(seed1)
  ! x2 = RAND(seed2)

  CALL RANDOM_NUMBER(x1); CALL RANDOM_NUMBER(x2)


  ! Generate normally distributed numbers, according to Box-Muller transform
  z1 = sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2)
  z2 = sqrt(-2.0 * log(x1)) * sin (2.0 * pi * x2)

  ! Only return either z1 or z2
  ! RETURN z1, z2
  rand_box_muller = z1
END FUNCTION
! Perhaps visualize lots of random numbers drawn from this function for report
! Compare with the uniform distribution of x1, x2


!==============================================================================
! (b) Produce a realisation of daily share price for an annual drift, mu = 0.16
! , and volatility, sigma = 0.2. Assume 250 trading days in a year. Use this to
! simulate daily share price for one year.
!==============================================================================
!  SUBROUTINE simulate_share_prices(t_final, realisations, volatility, share_price)
!    REAL, INTENT(IN) :: t_final, volatility, drift
!    INTEGER, INTENT(IN) :: realisations
!    ! REAL,
!
!   ! Parameters
!   REAL, PARAMETER :: mu = 0.16  ! annual drift
!   REAL, PARAMETER :: sigma = 0.2  ! volatility
!   REAL, PARAMETER :: r = 0.04  ! risk-free rate
!   REAL, PARAMETER  :: k = 0.2  !
!
!   ! Simulate daily share price over one year
!   t = 0.
!   t_final = 1.
!   ! dt = (1./250.)
!   dt = t_final / 1000
!   share_price_today = 100.
!   OPEN(unit=10, file='daily_simulated_share_prices.dat')
!   OPEN(unit=20, file='dW.dat')
!   WRITE(10,*) t*250, share_price_today
!   seed1 = 1; seed2 = INT(seed1 + t_final * 250)  ! Initialize seeds
!   DO WHILE (t <= t_final)
!     WRITE(6,*) t/t_final
!     seed1 = INT(seed1 + t*250.); seed2 = INT(seed2 + t*250.)  ! Different seeds
!     last_price = share_price_today
!     N = rand_box_muller(seed1, seed2)  ! Set dW to random value from Gaussian
!     ! N = rand_box_muller()
!     dw = sqrt(dt) * N
!     delta_share = 100. * (mu * dt + sigma * dw)
!     t = t + dt
!     share_price_today = last_price + delta_share  ! Evaluate price tomorrow
!     WRITE(10,*) t*250, share_price_today
!     WRITE(20,*) dw
!   END DO
!   CLOSE(10)
!   CLOSE(20)
!
! END SUBROUTINE
