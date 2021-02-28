This folder contains implementations of other optimization methods with the
same interface as BiteOpt. Using this methods is as simple as redefining
inherited class.

`ccmaes.cpp` - Interface to CMA-ES optimization method, requires CMA-ES C
library which can be obtained here: [c-cmaes](https://github.com/CMA-ES/c-cmaes).
Does not provide "stall count" information.

`nmpopt.cpp` - CNelderMeadPlusOpt class that implements sequential
Nelder-Mead simplex method with "Z" parameter improvement, which the author
discovered, and which significantly increases the quality of this rather old
and ineffective optimization method.

The essence of improvement is that each optimization parameter is split into
several (Z) sub-parameters (e.g. 4). These sub-parameters are summed into a
single parameter value before objective function evaluation. While intutively
this should increase convergence time, because problem's dimensionality is
increased by a factor Z, in practice this does not happen: somehow Nelder-Mead
method manages to keep a fast convergence with such vastly increased parameter
space.
