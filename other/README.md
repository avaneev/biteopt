This folder contains implementations of other optimization methods with the
same interface as BiteOpt. Using this methods is as simple as redefining
inherited class.

`ccmaes.cpp` - Interface to CMA-ES optimization method, requires CMA-ES C
library which can be obtained here: [c-cmaes](https://github.com/CMA-ES/c-cmaes).
Does not provide "stall count" information.
