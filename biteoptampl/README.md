This is an implementation of BITEOPT solver using AMPL model parser.

Designed for derivative-free optimization of models with low to average
dimensionality. Global optimum is not guaranteed, but as test have shown
global optimum is often reached. Solver can be recommended for models with
complex objective function for which computation of derivative is impossible,
otherwise other solvers are recommended.

It supports both continuous and integer variables. Constraint support is
experimental, but is working on most problems. The solver will report if it
could not meet a constraint. While inequality constraints are strict,
equality constraints are evaluated within +/- 1e-6 tolerance.

Solver runs up to itmult\*2000\*n_var^1.75\*sqrt(depth) iterations (function
evaluations) per attempt, or lower if solver plateaus. So, for
higher-dimensional problems the budget can be quite high, it can be decreased
via the `itmult` parameter.

Available solver parameters:

    attcnt  - attempt count (default 10)
    depth   - solver's depth (default 9), expected value 1 to 32.
    itmult  - iteration number multiplier (default 1.0)
    nprob   - objective choice: 1 (default) = 1st

If you have any questions please write to aleksey.vaneev@gmail.com
