# BiteOpt - Stochastic Function Optimizer #
## Introduction ##

### CBiteOpt (biteopt.h) ###

BiteOpt stochastic optimization class. Implements a stochastic non-linear
bound-constrained derivative-free optimization strategy. It maintains a
cost-ordered population of previously evaluated solutions that are evolved
towards a lower cost. On every iteration, the highest-cost solution in the
list is unconditionally replaced with a new solution, and the list is
reordered. A population of solutions allows the strategy to space solution
vectors apart from each other thus making them cover a larger parameter
search space collectively. Beside that, parameter randomization and the
"step in the right direction" operation are used that move the solutions
into position with a probabilistically lower objective function value.

The benefit of this strategy is a relatively high robustness: it can
successfully optimize a wide range of 1-10 dimensional test functions.
Another benefit is a low convergence time which depends on the complexity
of the objective function. Like many stochastic optimization strategies
with fast convergence, this strategy can't solve problems with narrow or
rogue optimums. Harder problems may require dozens of optimization attempts
to reach optimum.

### Notes ###

This "black-box" strategy was tested on 300+ 1-10 dimensional optimization
problems and performed well. Due to its design this strategy may be
particularly good at improving an existing sub-optimal local solution. This
strategy offers a very fast convergence on 1-3 dimensional problems, moderate
convergence speed on 4-10 dimensional problems, and slow convergence speed
on >10 dimensional problems, usually much slower than the best competing
strategies. However, on 1-3 dimensional problems there is little competition
to this strategy available.

This strategy was compared with the results of this paper (on 241 published
non-convex problems): [Comparison of derivative-free optimization algorithms](http://archimedes.cheme.cmu.edu/?q=dfocomp)
This strategy was able to solve 64% of problems in 10 attempts, 2500
iterations each. For 2-dimensional problems, this strategy's success rate is
96%. For 3-9 dimensional problems the success rate is 66%, for 10-30
dimensional problems the success rate is only 28%.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10-20 times is a minimal requirement. This strategy is hugely
probabilistic and it depends on its initial state, which is selected randomly.
In most cases it is more efficient to rerun the optimization with a new random
seed than to wait for the optimization process to converge. Based on the
results of optimization of the test corpus, for 2-dimensional functions it is
reasonable to expect convergence in 1000 iterations (in a successful attempt),
for 10-dimensional functions it is reasonable to expect convergence in 5000
iterations (harder functions may require more iterations to converge). Most
classic 2-dimensional problems converge in 300 iterations or less, at 1e-6
precision.

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

Strategy's "robustness" is a multi-factor non-formal estimation which includes
average convergence time, standard deviation of convergence time, the set of
functions the algorithm can solve successfully given randomized initial
conditions, in a given number of attempts. For global optimization robustness
also includes the lowest achieved cost.

## Examples ##

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

1. Parallelization of this algorithm is technically possible, but may be
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
