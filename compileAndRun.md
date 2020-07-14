To compile

	gcc -O2 -lm stsp.c -o stsp

To run

	stsp  inputfile.in

The input file for stsp has the following structure:

Input File Stucture

#PLANET PROPERTIES
	1                      ; Number of planets
	132.792729             ; T0, time of the middle of first transit in days  (Better if this number is closer to zero)
	1.48571127             ; Planet Period      (days)
	0.0222153113           ; Depth of transit (Rp/Rs)^2         (Rplanet/Rstar)^2
	0.094349811            ; Duration (days) of transit   (physical duration of transit, not used)
	0.0147313109           ; Impact parameter  (0= planet cross over equator, use inclination angle instead)
	89.8535141             ; Inclination angle of orbit (90 deg = planet crosses over equator)
	0.0                    ; Lambda of orbit (0 deg = orbital axis along z-axis) - angle between spin axis of star and orbital axis of planet
	0.0                    ; ecosw  
	0.0                    ; esinw  
#STAR PROPERTIES
1.15976225             ; Mean Stellar density (Msun/Rsun^3)  (related to a/Rstar)
12.1012                ; Stellar Rotation period (days)
5801                   ; Stellar Temperature  (not used)
0.0                    ; Stellar metallicity  (not used)
32.0                   ; Tilt of the rotation axis of the star down from z-axis (degrees)
0.83143 -0.68444 0.55014 -0.19879         ; Limb darkening (4 coefficients)
100                    ; number of rings for limb darkening approximation
#SPOT PROPERTIES
6                      ; number of spots
0.67                   ; fractional brightness of spots (0.0= totally dark, 1.0=brightness of star)
# FITTING properties
kepler17.lc            ; lightcurve data file
1201.000               ; start time to start fitting the light curve
10.0                   ; duration of light curve to fit (days)
1.00002909234          ; real maximum of light curve data (corrected for noise), 0 -> STSP uses simple downfrommax
1                      ; is light curve flattened (to zero) outside of transits?
# ACTION
H                      ; M / S / T / L / H -- M: Affine invariant MCMC sampling starting with random values for all parameters (MCMC)
                       ; fM / fS / fT         S or s: Affine invariant MCMC seeded with initial values given in the next lines after this one (Seeded MCMC)
                            T: Affine invariant MCMC where all parameters for all chains are seeded (Totally seeded MCMC)
                            L or l: Just generate the light curve for a set of parameter values given in the next lines
                            H: Metropolis-Hastings MCMC sampling starting with random values for parameters
                       ; small s and l mean

--- THE REMAINING LINES DEPEND ON THE VALUE OF THE ACTION PARAMETER --

If the ACTION is set to m, then the following lines need to be random seed, scale factor, number of chains, number of steps/chain

# ACTION
   m              ; M= unseeded mcmc
   74384338       ; random seed
   1.25000        ; ascale
   40             ; number of chains
   5000           ; mcmc steps
   1              ; normalization factor, 0= use downfrommax normalization, 1= calculate brightness factor for every model

If the ACTION is set to s, then the following lines need to be the same 5 above + sigma for radius and angle variations and spot properties for Ns spots

# ACTION
   s               ; s= unseeded mcmc
   74384338        ; random seed
   1.25            ; ascale
   300             ; number of chains
   -14400          ; mcmc steps
   1               ; 0= use downfrommax normalization, 1= calculate brightness factor for every model
   0.002           ; sigma for radius variation
   0.01            ; sigma for angle variation
   0.0233321221            ; spot radius
   1.3235347387            ; theta
   4.9132315049            ; phi
   1                       ; brightness correction factor associated with this spot model

If the ACTION is set to l, then the following lines need to give the radius (0.0-->1.0), latitude (in radians) and longitude (in radians) for each spot.
For example, for 3 spots and to generate the light curve only, the following lines of the input file would read:

# ACTION
    l
    0.220		; spot 1 radius
    2.457		; spot 1 theta
    3.829		; spot 1 phi
    0.418		; spot 2 radius
    1.955		; spot 2 theta
    4.791		; spot 2 phi
    0.157		; spot 3 radius
    1.533		; spot 3 theta
    0.021		; spot 3 phi
    1			; brightness correction factor associated with this spot model

If the action is capitalized all the action parameters are in a text file. STSP expects the name of the text file with the Action parameters

# ACTION 

    L
    kepler17_parambest.txt  ; file name for parameter file from a previous run that contains the parameters for all the objects

