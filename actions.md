##Description of Action Calls to STSP

The STSP input file must specify an ``Action" to be implemented with the run of STSP and any parameters necessary to implement that Action.  We use a convention in which action calls with lower case letters mean that all the parameters needed to run the action are listed on subsequent lines in the input file.  Action calls with upper case letters imply that the parameters needed to run the action are listed in an external file with the name of the external file given in the input file.  


####Optimizer Action Calls to STSP

**[Action-L/l]** allows for generating a synthetic light curve based on a set of spot properties provided by the user.  The input properties are always listed in the order of $r_i$, $\theta_i$, and then $\phi_i$ for $i = 1 \ldots N_s$ spots.  The final necessary parameter is the normalization factor.

**[Action-M/m]** allows for running the affine invariant MCMC optimization algorithm with random initial conditions for all chains.  The input
parameters for this action are a random seed, the MCMC $a$ scale parameter, the number of chains, the number of steps, and finally the normalization factor.  All other actions that run the affine invariant MCMC optimization algorithm (Action-s/S and Action-T) require the five basic
parameters listed above plus additional parameters.

**[Action-S/s]** is a ``seeded" run of the affine invariant MCMC optimizer in which the user supplies a single set of spot properties which STSP uses to populate the starting conditions for all chains.  STSP draws random starting values for all spot parameters from a Gaussian distribution centered on the single set of properties provided by the user.  Like Action-m/M, Action-s/S takes the five basic parameters needed to run the affine-invariant MCMC, but it also takes two Gaussian sigma parameters that determine the width of the distributions for the spot radii and spot angles ($\theta_i$ and $\phi_i$).  Finally, it requires a set of spot properties for all $N_s$ spots, like Action-L/l.  

**[Action-T]** is a ``totally-seeded" run of the MCMC which populates all MCMC chains with a specific set of spot parameters.  This mode is used to continue runs that have previously been started, but need to be run for additional steps in order to achieve convergence.  This action requires the five basic affine-invariant MCMC parameters plus a set of spot parameters ($r_i$, $\theta_i$, $\phi_i$) for all $N_s$ spots and all chains.  This action requires an input file to implement, so there is no Action-t.

**[Action-U]** is a run of the MCMC in which some spots are fixed and others are allowed to vary.  This action is useful when the properties of spots in the path of the transit have already been defined, but additional spots are needed to also characterize the out-of-transit variability.

**[Action-H]** runs the single chain, Metropolis-Hastings version of the MCMC algorithm.  To implement this action, the user supplies a single set of spot parameters for a single chain and the number of steps to run.  STSP draws a set of random starting values for all spot parameters from a Gaussian distribution centered on the user supplied set of spot properties, thus this action also requires two sigma parameters to determine the width of the distributions from which the random starting values are drawn.

**[f]** modifier to an MCMC action to indicate that the out-of-transit variability has been removed from the light curve.  This action is applied to other actions in the following way: `fS', `fM', `fT', `fU'

####Diagnostic Action Calls to STSP

**[Action-I]**  takes a list of times corresponding to features in the observed light curve and reports the longitude and latitude on the star underneath the planet at those times.  This is useful for determining initial guesses for spot positions for seeded MCMC runs and requires
a file input with the list of input times in the same units as the T_0 epoch of the planet.

**[Action-A/a]** takes a defined spot configuration and calculates the fraction of the star's area probed by the planet during transit which is covered by those spots.  This action is usually implemented after running the optimizer and settling on a final starspot solution.



