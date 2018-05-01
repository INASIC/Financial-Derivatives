# What is in this Repository?

A statistical analysis of the recent price of Bitcoin compared to the British pound was performed, and a hedging strategy for a hypothetical stock presented in a report. A program that models a stochastic simulation for the hypothetical stock is provided as "statistical_arbitrage.f90", which uses the Black-Scholes model to simulate the predicted price movement of the share. Hedging strategies to maximise capital by trading this share are discussed and performed by this program, and their effectiveness discussed in the report.

# Access to the Report

In the case that the report PDF does not load on GitHub, access to it can be viewed on the following Google Doc:
https://docs.google.com/document/d/1dv7NzptODqQffby55s2Gbq3aBpOjWUZzqGSc7f93Odo/edit?usp=sharing.

# Compiling the Program

To compile the program, ensure that you have a FORTRAN compiler. If using the FORTRAN compiler gfortran then simply type in the terminal:
* gfortran -o sa statistical_arbitrage.f90 && ./sa
