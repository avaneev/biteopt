# BiteOpt - Stochastic Function Optimizer #
## Introduction ##

### CBiteOpt (biteopt.h) ###

BiteOpt stochastic function optimization class. This optimization strategy
is based on now outdated CBEOOptimizer and CBEOOptimizer2, CBEOOptimizerFan
stochastic derivative-less strategies. It uses an ordered list of several
current parameter vectors (called "fan elements") that are evolved towards
a lower cost. Highest-cost "fan element" in the list can be replaced with a
new solution if "fan element's" cost (plus some margin) is higher than that
of the new solution's. Having several "fan elements" allows the strategy
to space them apart from each other thus making them cover a larger
parameter search space collectively. Beside that, parameter randomization
and several "step in the right direction" operations are used that move the
parameter vector into position with a probabilistically lower objective
function value.

The benefit of this strategy is a relatively high robustness: it can
successfully optimize a wide range of test functions. Another benefit is a
low convergence time which depends on the complexity of the objective
function. Like many stochastic optimization strategies with fast
convergence, this strategy can't solve problems with narrow or rogue
optimums. Harder problems may require dozens of optimization attempts to
reach optimum.

### Notes ###

This strategy was tested on 100+ classic 2-dimensional and 30-dimensional
optimization problems and performed well. Due to its design this strategy may
be particularly good at improving an existing sub-optimal local solution.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10-20 times is a minimal requirement. This strategy is hugely
probabilistic and it depends on its initial state, which is selected randomly.
In most cases it is more efficient to rerun the optimization with a new random
seed than to wait for the optimization process to converge. Based on the
results of optimization of the test corpus, for 2-dimensional functions it is
reasonable to expect convergence in 2000 iterations (in a successful attempt),
for 10-dimensional functions it is reasonable to expect convergence in 10000
iterations (harder functions may require more iterations to converge). Most
classic 2-dimensional problems converge in 300 iterations.

The required number of optimization attempts is usually proportional to the
number of strongly competing minima in a function. Rogue optimums may not be
found by this strategy. A rogue optimum is an optimum that has a very small
area of descent and is placed apart from other competing minima. The
strategy favors minimum with a larger area of descent. The Damavandi test
function is a perfect example of the limitation of this strategy. In practice,
however, rogue optimums can be considered as undesired outliers that have an
unstable real-life performance due to existing parameter value tolerances.

To some degree this strategy is immune to noise in the objective function.
While this strategy was designed to be applied to continuous functions, it is
also immune to discontinuities to some degree, and it can solve problems that
utilize parameter value rounding. This strategy can't acceptably solve
high-dimensional problems that are implicitly or explicitly combinatorial.

Most hard constraints can be introduced by applying huge penalties to the
objective function. Even binary penalties like "if(x>1)cost+=x*10000" should
work acceptably in most cases.

The minimal and maximal allowed parameter values should be specified in a way
to cover a wider value range, in order to reduce boundary effects that may
greatly reduce convergence. If the target local or global minimum stands
very close to the parameter value boundaries this strategy may fail to
converge.

The FanSize parameter can be adjusted to increase quality of solutions,
especially in high-dimensional problems. FanSize controls robustness of the
strategy at the cost of convergence time. However, setting FanSize too high
may increase rejection rate (equivalent to problem's overspecification).

The optimizePlateau() function can be used to optimize functions that have a
complex landscape. This function works within the bounds of allocated
iteration limit, and performs re-initializations when the optimization process
reaches the plateau.

Strategy's "robustness" is a multi-factor non-formal estimation which includes
average convergence time, standard deviation of convergence time, the set of
functions the algorithm can solve successfully given randomized initial
conditions, in a given number of attempts. For global optimization robustness
also includes the lowest achieved cost.

Use the example.cpp program to see the basic usage example.

test2.cpp is a convergence test for all available functions. Performs many
optimization attempts on all functions. Prints various performance
information, including percentage of rejected attempts (rejection rate).

test3.cpp demonstrates the usage of the optimizePlateau() function.

test4.cpp is a convergence test for multi-dimensional functions.

## Development ##

While the basic algorithm of the strategy is finished, the built-in parameters
of the algorithm is an area of ongoing research. There are several things that
were discovered that may need to be addressed in the future:

1. The formula of dependence of FanSize on the number of dimensions may need
to be updated to better suit varying dimensionality of the problems.

2. The values of the CentProb and RandProb mutually affect optimization
rejection rate and convergence time. That is why it is reasonable to select
these values manually and optimize other parameters with these values fixed.

3. If high optimization attempt rejection rate is not problematic, the
parameters of the algorithm can be tuned to provide at least 20% lower
convergence time.

4. Parallelization of this algorithm is technically possible, but is
counter-productive (increases convergence time considerably). It is more
efficient to run several optimizers in parallel with different random seeds.

## Warning ##

When solving problems whose solutions are critical and may be
health-threatening, always use several optimization strategies (methods) to
find the optimal solution, do not rely on a single strategy when solving
multi-modal problems.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your project to the list of users.
