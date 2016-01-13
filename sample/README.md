# Example Files

This directory contains example files which demonstrate the available "actions"
which can be used within `STSP`. 

#### Available Actions

| Action | Description | 
|--------|-------------|
|    L   | generate light curve from parameter file |
|    S   | seeded mcmc from parameter file | 
|    fS   | seeded mcmc from parameter file with fixed thetas |
|    fM   | unseeded mcmc with fixed thetas |
|    fs   | seeded mcmc from parameters with fixed thetas |
|    fT   | totally seeded mcmc from parameter file with fixed thetas |
|    H   | metropolis-hastings from seed parameters |
|    l (lowercase L)  | generate light curve from parameters |
|    M   | unseeded mcmc |
|    P   | plot chi squared as one spot varies (plot parameters are in plotdata function) |
|    s   | seeded mcmc from parameters |
|    T   | totally seeded mcmc from parameter file |
|    u   | partially seeded mcmc |

# Notes

### Limb-darkening parameters

STSP input files define the limb-darkening parameters of the host star using the "nonlinear" parameterization, such that: 

`I(mu)/I_0 = 1 - c_1*(1-mu^0.5) - c_2*(1-mu) - c_3*(1-mu^1.5) - c_4*(1-mu^2)`

rather than the commonly fit quadratic limb-darkening parameters: 

`I(mu)/I_0 = 1 - u_1*(1-mu) - u_2*(1-mu)^2`

If you've fit for (u1, u2) and need to transform into the non-linear parameters, use the following transformation: 

```
c_1 = 0
c_2 = u_1 + 2*u_2
c_3 = 0
c_4 = -u_2
```
