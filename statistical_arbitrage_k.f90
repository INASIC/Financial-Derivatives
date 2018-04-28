PROGRAM statistical_arbitrage
IMPLICIT NONE
!http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html


!State all variables
REAL :: r, k, x1, x2, z1, z2, pi, dW, dS, S, S0, time, mean, variance
REAL, PARAMETER :: mu = 0.16, dt = (1.0/250.0), sigma = 0.2
REAL, ALLOCATABLE :: port(:,:)

INTEGER :: realisation, realisations, tr_days, i_seed, n, m, years, day
INTEGER, ALLOCATABLE :: seed(:)
INTEGER, DIMENSION(1:8) :: dt_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed


!Ask the user to enter the number of years the program will be run
   WRITE(6,*) ''
   WRITE(6,*) 'Please enter the number of years to be tested for:'
   WRITE(6,*) '--> Question 2B: 1 year'
   WRITE(6,*) '--> Question 2C: 1&2&5 years'
   WRITE(6,*) '--> Question 2E: 20 years'
   WRITE(6,*) '--> Question 2F: 20 years'
   WRITE(6,*) ''
   READ(5,*) years


!Calculate number of trading days
   tr_days = years/dt

!Constants
   r = 0.04
   k = 0.2
   pi = 4.0*ATAN(1.0)
   S0 = 100.0
   realisations = 1000

   !Allocate the trading days and number of realisations to the portfolio array
      ALLOCATE(port(tr_days,realisations))

!Random Seed
!Sets up random numbers that reset each time the program is run
   CALL RANDOM_SEED(SIZE = i_seed)
   ALLOCATE(a_seed(1:i_seed))
   CALL RANDOM_SEED(GET = a_seed)
   CALL DATE_AND_TIME(values = dt_seed)
   a_seed(i_seed) = dt_seed(8); a_seed(1) = dt_seed(8)*dt_seed(7)*dt_seed(6)
   CALL RANDOM_SEED(PUT = a_seed)
   DEALLOCATE(a_seed)


!Open all the files to write the data to, before big do-loop
   OPEN(10,file='2B_rand_num.txt')
   OPEN(20,file='2C_1realisation.dat')
   OPEN(30,file='2C_1000realisations.dat')
   OPEN(40,file='2D_histogram.txt')
   OPEN(50,file='2E_mean_and_variance.dat')


!Set of big do-loop
   !Reset portfolio to zero
   port(:,:) = 0.0
   !Do-loop
   DO realisation=1, realisations
      S = 100.0
      time = 0.0
      DO day=1, tr_days
         time = time + dt


!-------------------------------------------------------------------------------
!Question 2B - (calculate one set of random numbers)
            !Call random numbers x1 and x2
            CALL RANDOM_NUMBER(x1)
            CALL RANDOM_NUMBER(x2)
            !Use the Box-Muller Transform to produce both random number producing
            !equations
            !Eqn(3): Random number using cosine
            z1 = SQRT((-2.0)*LOG(x1))*COS(2.0*pi*x2)
            !Eqn(4): Random number using sine
            z2 = SQRT((-2.0)*LOG(x1))*SIN(2.0*pi*x2)
            !Write one set of the random numbers to file '2B_rand_num.txt'
            WRITE(10,*) z1

!-------------------------------------------------------------------------------
!Question 2C - (calculate share price for 1 realisations, and 1000 realisations,
              !(for use in part 2E)
            !Calculate dW, random number based on time
            dW = SQRT(dt)*Z1
            !Eqn(1): Change in share price
            dS = S*((mu*dt) + (sigma*dW))
            !Update share (add 'change in share price')
            S = S + dS
            !Write share price and time for 1 realisations to file '2C_1realisation.dat'
            IF (realisation == 1) THEN
               WRITE(20,*) time, S
            END IF
            !Write share price and time for 1000 realisations to file '2C_1000realisations.dat'
            WRITE(30,*) time, S

!-------------------------------------------------------------------------------
!Question 2D (1, 2, and 5 years)
         !Portfolio t<=t*
         port(day,realisation) = (S*EXP(-(r*time))) - S0
         !Portfolio t>t*
         IF (S >= S0*(1.0+k)*EXP(r*time)) THEN
            port(day,realisation) = S0*k
         END IF
         !Write portfolio values, to text file '2D_histogram.txt'
         !to be plotted into a histogram
         !WRITE(40,*) port(day,realisation)
   END DO
END DO

WRITE(6,*) 'sup'

DO day=0, tr_days
    mean = SUM(port(day, :))/1000
    WRITE(40,*) mean
END DO

!-------------------------------------------------------------------------------
!Question 2E (20 years)


























!      DO n=1,tr_days
!               time = time + dt
      !Mean
!      mean = SUM(port(n,1000))/realisations
      !Variance
      !Set variance to a starting value of zero
!      variance = 0.0
!      DO m=1,realisations
!         variance = variance + ((port(n,m) - mean)**2.0)
!      END DO
!      variance = variance/(realisations - 1.0)
      !Write mean, variance and time to file '2E_mean_and_variance.dat'
!      WRITE(50,*) mean, variance, time
!      END DO
!-------------------------------------------------------------------------------
!Question 2F (20 years)



!-------------------------------------------------------------------------------
!Close all the files
   CLOSE(10)
   CLOSE(20)
   CLOSE(30)
   CLOSE(40)
   CLOSE(50)

END PROGRAM statistical_arbitrage
