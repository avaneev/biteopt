This is an implementation of BiteOpt solver using AMPL model parser.

Designed for derivative-free optimization of models with low dimensionality.
Global optimum is not guaranteed, but as the tests have shown global optimum
is often reached. Solver can be recommended for models with a complex
non-linear objective function for which computation of derivative is
impossible, otherwise other solvers are recommended.

It supports both continuous and integer variables. Constraints satisfaction is
working for most problems (via penalties). The solver will report if it could
not meet constraints. Equality constraints are evaluated within the configured
tolerance (10<sup>-5</sup> by default), non-equality constraints are evaluated
strictly.

Solver runs for a dimension-dependent maximum number of iterations (function
evaluations) per attempt, which is a hard limit, or lower if solver plateaus.
For higher-dimensional problems this budget can be quite high, but it can be
decreased via the `itmult` parameter. Solver reports the number of attempts
where hard iteration limit was reached without reaching optimization plateau.
If all attempts reached this limit, increasing `itmult` is advisable.

The model `.mod` file should first be converted to the `.nl` format via
the `ampl -og` command.

Available solver parameters can be printed with the `-=` command line option.

If you have any questions please write to aleksey.vaneev@gmail.com
