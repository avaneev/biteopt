This folder contains implementations of other optimization methods with the
same interface as BiteOpt. Using this methods is as simple as redefining
inherited class.

`cnmp.cpp` - CNelderMeadPlusOpt class that implements sequential Nelder-Mead
simplex method with "Z" parameter improvement (improves convergence in many
cases). Provides "stall count" information, but the plateau threshold value
unresearched.

`ccmaes.cpp` - Interface to CMA-ES optimization method, requires CMA-ES C
library which can be obtained here: [c-cmaes](https://github.com/CMA-ES/c-cmaes).
Does not provide "stall count" information.
