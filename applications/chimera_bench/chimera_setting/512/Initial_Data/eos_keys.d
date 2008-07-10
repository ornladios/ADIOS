
cccccc   Equation of state parameters

!-----------------------------------------------------------------------
!         Cube grid
!
!      dgrid(i), tgrid(i), ygrid(i): log(rho), log(t), ye space is overlain with a uniform grid of 
!       'dgrid(i)' divisions per unit change in log(rho), 'tgrid(i)' divisions per unit change in 
!       log(t), and 'ygrid(i)' divisions per 0.5 change in ye. Equation of state, nuclear reaction
!       rate, and neutrino interaction variables at each radial zone are interpolated from values at
!       the nearest corners on this grid.
!
!      rhoes(k): The variables dgrid, tgrid, and ygrid are each 3 element arrays permitting different 
!       partitionings of log(rho), log(t), ye space in different density regimes delimited by rhoes. 
!       These different regimes are
!
!          regime 1:             rho < rhoes(1)
!          regime 2:        rhoes(1) < rho < rhoes(2)
!          regime 3:             rhoes(2) < rho
!
!      idty(j): the density regime (i.e., 1, 2, or 3) of radial zone j as given by the above 
!       inequalities.
!-----------------------------------------------------------------------

esgrid                       1      1.00000000E+01                      dgrid
esgrid                       2      3.00000000E+01                      dgrid
esgrid                       3      1.00000000E+01                      dgrid
esgrid                       1      2.00000000E+01                      tgrid
esgrid                       2      2.00000000E+01                      tgrid
esgrid                       3      2.00000000E+01                      tgrid
esgrid                       1      5.00000000E+01                      ygrid
esgrid                       2      5.00000000E+01                      ygrid
esgrid                       3      5.00000000E+01                      ygrid
esgrid                       1      3.00000000E+13                      rhoes
esgrid                       2      2.00000000E+14                      rhoes

!-----------------------------------------------------------------------
!         Equation of state identifier
!
!      eos_i - equation of state identifier
!
!          eos_i = 'L'  : Lattimer-Swesty equation of state is used.
!          eos_i = 'B'  : Cooperstein-BCK equation of state is used.
!-----------------------------------------------------------------------

eos_i                        L                                          eos_i

!-----------------------------------------------------------------------
!         Equation of state borders
!
!      eosrho: the border density between the LS EOS and the BCK EOS (/fm3).
!-----------------------------------------------------------------------

eosrho                              1.00000000E-08                      eosrho

!-----------------------------------------------------------------------
!         NSE flashing and deflashing controls
!
!      tnse: temperature at which material is flashed to nse (K)
!
!      tdnse: temperature at which material is deflashed from nse (K)
!-----------------------------------------------------------------------

tnse                                5.10560000E+09                      tnse
tdnse                               3.50000000E+09                      tdnse
