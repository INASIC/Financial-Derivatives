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
REAL, FUNCTION derive_present_value()
  ! Parameters
  INTEGER, PARAMETER :: t_invest_bond = 100  ! Time sell for risk free bond
  INTEGER, PARAMETER :: t_final = 1 * 365
  INTEGER, PARAMETER ::
  INTEGER :: t0 = 0, dt = 1  ! Final day, Integration time
  REAL, PARAMETER :: r =   ! Risk free rate
  REAL, PARAMETER :: k =   !

  ! Helpers
  INTEGER :: t
  REAL, DIMENSION(t0:t_final) :: v, share_price, discounted_value


  share_price(t0) = 100  ! Share price at t = 0
  t_invest_bond = 100

  ! Portfolio discounted value as a function of time (Eq. 2)
  IF (t <= t_invest_bond) THEN
    v(t) = share_price(t) * exp(-r*t) - share_price(t0)
  ELSE
    v(t) = share_price(t0) * k
  END IF

  ! Equation 2:
  discounted_value = v
  derive_present_value = discounted_value

END FUNCTION

!===========================================================================
! (b) Function that generates a normally distributed random number with
! mean zero and variance one from the Box-Muller transform
!===========================================================================
REAL FUNCTION rand_box_muller()
  INTEGER :: seed1 = 1, seed2 = 2
  REAL :: x1, x2, z1, z2
  REAL :: pi = 3.1415926535897932384626433832795028841971693993751058209749445
  ! Intrinsic uniform random number
  CALL SRAND(seed1); x1 = RAND()
  CALL SRAND(seed2); x2 = RAND()

  ! Generate normally distributed numbers, according to Box-Muller transform
  z1 = sqrt(-2.0 * ln(x1) * cos(2.0 * pi * x2))
  z2 = sqrt(-2.0 * ln(x1)) * sin (2.0 * pi * x2)

  ! Only return either z1 or z2
  ! RETURN z1, z2
  rand_box_muller = z1
END FUNCTION
! Perhaps visualize lots of random numbers drawn from this function for report


!==============================================================================
! (b) Produce a realisation of daily share price for an annual drift, mu = 0.16
! , and volatility, sigma = 0.2. Assume 250 trading days in a year. Use this to
! simulate daily share price for one year.
!==============================================================================
 SUBROUTINE simulated_share_prices(t0, tfinal)
  ! Parameters
  REAL, PARAMETER :: mu = 0.16  ! annual drift
  REAL, PARAMETER :: sigma = 0.2  ! volatility
  REAL, PARAMETER :: r = 0.04  ! risk-free rate
  REAL, PARAMETER  :: k = 0.2  !
  INTEGER, PARAMETER :: t_final = 250  ! trading days in a year
  INTEGER, PARAMETER :: dt = 1  ! daily intervals
  INTEGER, PARAMETER :: t0 = 0

  ! Helpers
  REAL, DIMENSION(t0:t_final) :: share_price
  REAL simulated_share_prices(t0:tfinal)  ! Set as a vector
  REAL :: dw =  ! Wiener process
  REAL :: delta_share
  INTEGER :: t,

  share_price(t0) = 100  ! Starting share price at t=0
  dw = ...  ! Wiener process

  ! Simulate daily share price over one year
  t = 0
  DO WHILE (t <= t_final)
    delta_share = (mu * dt + sigma * dw) * share_price(t)
    t = t + dt
    share_price(t) = share_price(t-1) + delta_share  ! Evaluate price tomorrow
  END DO

  ! share_price_after(1, years, realisations)

  realisations = 1000
  dt = realisations / t_final

  t_final = 1 * 250  ! 1 year
  t_final = 2 * 250  ! 2 years
  t_final = 5 * 250  ! 5 years

  ! Evaluate v from t=0 up to 20 years
  years = 20
  mean_v = SUM(v) / realisations
  the_sum = 0
  n = realisations

  DO i=0, n
    the_sum = the_sum (v(i) - mean_v) ** 2.
  ENDDO

  variance_v = 1.0/(n-1) * sum(x - mean_v) ** 2

  OPEN(unit=10, file='simulated_share_prices.dat')
    DO t=t0, t_final
      WRITE(10,*) t, share_price(t)
    END DO
  CLOSE(10)
END SUBROUTINE
