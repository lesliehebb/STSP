STSP
====

In and out of transit starspot modeling code

V4.4.2  2015-04-08

L & G fixed bug in combining solutions from the walkers, which should help speed up convergence.


### NOTE
The default mode should have these flags set in the top of the code, which improves speed during the MCMC solutions being run in bulk.

    #define QUIET 1
    #define QUIETMCMC 1
    #define ANYPRINTVIS 0
    
To get the optional visaulization outputs to work (action L), change them to this and recompile STSP:

    #define QUIET 0
    #define QUIETMCMC 0   // optional
    #define ANYPRINTVIS 1



V4.4.1  2014-09-09

Light curve file now has an extra column that identifies which starspots 
are interacting with the planet at each time step.  The column is
a binary number where:  
no spot-planet overlap = 0
spot0 overlaps with planet = +1
spot1 overlaps with planet = +2
spot k overlaps with planet = +2^k

To turn the number, n, to a string of 0/1 digits:
Loop over the number of spots (or the number of digits in the n)
For the ith digit:
  (1) k = n mod 2
  (2) spot_flag[i] = k   (this should be 1 or 0)
  (3) n = (n - k)/2   (subtract 1 and divide by 2)

V4.4  2014-05-23

Fixed bug related to combine 1 spot
Added functionality to injest a flattened light curve and to flat the models
Added functionality to use if fixed longitude mode

V4.33  2014-05-06

Change sigma value for brightness factor from 0.1 --> 0.001
that is used during seeded runs.

V4.32  2014-04-23

Puts a filecheck on the input file.  Bug fix for name of errstsp.txt filename.

V4.31  2014-03-21

Added bug fix that corrects problem with combining
across phi =0-->2pi boundary when in fixed latitude mode

V4.3

Added new fitted parameter, Brightness_ratio,
that corrects the normalization of the input light curve

V4.21 2014-01-14

Bug fix to program JRAD encountered in Baltimore when large spot & planet at very limb of star as spot rotating out
of view.

V4.2 of the starspot program

2013-12-18

Lastest starspot program.  Added two new functionalities
  -- allows for fixing the latitudes for the spots
  -- allows misaligned planet & star with lambda parameter

to compile

gcc -lm stsp.c -o stsp

to run

stsp  stsp_kepler17_starspot.in

New input file:

The input file for stsp has the following structure:

New input file:

The input file for stsp has the following structure:
##Input File Stucture
    #PLANET PROPERTIES
    1                      ; Number of planets
    132.792729             ; T0, epoch in days (time of the middle of first transit)  
    						  Better if this number is closer to zero.
    1.48571127             ; Planet Period      (days)
    0.0222153113           ; Depth of transit (Rp/Rs)^2         (Rplanet / Rstar )^ 2
    0.094349811            ; Duration (days) of transit   (physical duration of transit, not used)
    0.0147313109           ; Impact parameter  (0= planet cross over equator, not used)
    89.8535141             ; Inclination angle of orbit (90 deg = planet crosses over equator)
    0.0                    ; Lambda of orbit (0 deg = orbital axis along z-axis) - angle between spin axis of star and orbital axis of planet
    0.0                    ; ecosw  (currently place holders)
    0.0                    ; esinw  (currently place holders)
    #STAR PROPERTIES
    1.15976225             ; Mean Stellar density (Msun/Rsun^3)  (not used if nplanets = 0)
    12.1012                ; Stellar Rotation period (days)
    5801                   ; Stellar Temperature  (not used)
    0.0                    ; Stellar metallicity  (not used)
    32.0                   ; Tilt of the rotation axis of the star down from z-axis (degrees)
    0.83143 -0.68444 0.55014 -0.19879         ; Limb darkening (4 coefficients)
    100                    ; number of rings for limb darkening approximation
    #SPOT PROPERTIES
    6                      ; number of spots
    0.67                   ; fractional brightness (0.0= totally dark, 1.0=brightness of star)
    # FITTING properties
    kepler17.lc            ; lightcurve data file
    1201.000               ; start time to start fitting the light curve
    10.0                   ; duration of light curve to fit (days)
    # ACTION
    H                      ; M / S / T / L / H -- M: Affine invariant MCMC sampling starting with random values for all parameters (MCMC)
                           ; fM / fS / fT         S or s: Affine invariant MCMC seeded with initial values given in the next lines after this one (Seeded MCMC)
                                T: Affine invariant MCMC where all parameters for all chains are seeded (Totally seeded MCMC)
                                L or l: Just generate the light curve for a set of parameter values given in the next lines
                                H: Metropolis-Hastings MCMC sampling starting with random values for parameters
                           ; small s and l mean
    
    --- THE REMAINING LINES DEPEND ON THE VALUE OF THE ACTION PARAMETER --

    If the ACTION is set to M, then the following lines need to be random seed, scale factor, number of changes,
    #ACTION
    M              ; M= unseeded mcmc
    74384338       ; random seed
    1.25000        ; ascale
    40             ; number of chains
    5000           ; mcmc steps
    1              ; 1= combine one spot at a time, 0= combine all spots.
    
    If the ACTION is set to l or s, then the following lines need to give the radius (0.0-->1.0), latitude (in radians) and longitude (in radians) for each spot.
    For example, for 3 spots and to generate the light curve only, the following lines of the input file would read:
    # ACTION
    l
    0.220
    2.457
    3.829
    0.418
    1.955
    4.791
    0.157
    1.533
    0.021

    # ACTION
    L or S or T

    L
    kepler17_parambest.txt  ; file name for parameter file from a previous run that contains the parameters for all the objects

