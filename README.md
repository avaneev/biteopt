# BiteOpt - Stochastic Function Optimizer #
## Introduction ##

### CBiteOpt (biteopt.h) ###

BiteOpt stochastic optimization class. Implements a stochastic non-linear
bound-constrained derivative-free optimization strategy. It maintains a
cost-ordered population list of previously evaluated solutions that are
evolved towards a lower cost. On every iteration, the highest-cost solution in
the list can be replaced with a new solution, and the list reordered. A
population of solutions allows the strategy to space solution vectors apart
from each other thus making them cover a larger parameter search space
collectively. Beside that, parameter randomization and the "step in the
right direction" operation are used that move the solutions into position
with a probabilistically lower objective function value.

The benefit of this strategy is a relatively high robustness: it can
successfully optimize a wide range of multi-dimensional test functions.
Another benefit is a low convergence time which depends on the complexity
of the objective function. Like many stochastic optimization strategies
with fast convergence, this strategy can't solve problems with narrow or
rogue optimums. Hard (multi-modal) problems may require many optimization
attempts to reach optimum. The name "BiteOpt" is an acronym for "BITmask
Evolution OPTimization".

### CBiteOptDeep (biteoptd.h) ###

Deep stochastic optimization class. Based on an array of M CBiteOpt
objects. This "deep" strategy pushes the newly-obtained solution to the
next random CBiteOpt object which is then optimized. This strategy while
increasing the convergence time by a factor of M is able to solve even the
most noisy non-linear functions.

This strategy is most effective on stochastic functions or functions with
huge fluctuations near the global solution that are not very expensive to
calculate and that have a large iteration budget. Tests have shown that on
smooth functions that have many strongly competing minima this strategy
does not considerably increase the chance to find a global solution
relative to the CBiteOpt class, and still requires several runs at
different random seeds.

### Notes ###

The text below mainly addresses the CBiteOpt class.

This "black-box" strategy was tested on 400+ 1-10 dimensional optimization
problems and performed well, and it successfully solves even 600-dimensional
test problems found in some books. Due to its design this strategy may be
particularly good at improving an existing sub-optimal local solution.

This strategy was compared with the results of this paper (on 241 published C
non-convex problems, convex problems were not evaluated): [Comparison of derivative-free optimization algorithms](http://archimedes.cheme.cmu.edu/?q=dfocomp)
This strategy was able to solve 78% of non-convex problems in 10 attempts, 2500
iterations each. For 1-2 dimensional problems, this strategy's success rate is
99%. For 3-9 dimensional problems the success rate is 70%, for 10-30
dimensional problems the success rate is 86%, for >30 dimensional problems the
success rate is 38%. With a huge iteration budget (1 million) this strategy
solves 93% of problems. On a comparable test function suite and conditions
outlined at this page: [global_optimization](http://infinity77.net/global_optimization/multidimensional.html)
(excluding several ill-defined and overly simple functions, and including
several complex functions, use test2.cpp to run the test) this strategy's
success rate is >91% while the average number of objective function
evaluations is ~340.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10 times is a minimal general requirement. This strategy is hugely
probabilistic and it depends on its initial state, which is selected randomly.
In most cases it is more efficient to rerun the optimization with a new random
seed than to wait for the optimization process to converge. Based on the
results of optimization of the test corpus, for 2-dimensional functions it is
reasonable to expect convergence in 1000 iterations (in a successful attempt),
for 10-dimensional functions it is reasonable to expect convergence in 5000
iterations (harder functions may require more iterations to converge). Most
classic 2-dimensional problems converge in 400 iterations or less, at 1e-6
precision.

The required number of optimization attempts is usually proportional to the
number of strongly competing minima in a function. Rogue optimums may not be
found by this strategy. A rogue optimum is an optimum that has a very small
area of descent and is placed apart from other competing minima. The
strategy favors minimum with a larger area of descent. The Damavandi test
function is a perfect example of the limitation of this strategy. In practice,
however, rogue optimums can be considered as undesired outliers that rely on
unstable parameter values (if such parameters are used in real-world system
that has a certain parameter value precison, a system may leave the "rogue"
optimal regime easily).

To some degree this strategy is immune to noise in the objective function.
While this strategy was designed to be applied to continuous functions, it is
also immune to discontinuities to some degree, and it can solve problems that
utilize parameter value rounding (integer parameters). This strategy can't
acceptably solve high-dimensional problems that are implicitly or explicitly
combinatorial (e.g. Perm and Lennard-Jones atom clustering problems). Also
problems with many competing minima without a pronounced global gradient
(e.g. Bukin N.6) may not be solved acceptably as in most cases they require
exhaustive search.

Mixed integer programming can be achieved by using rounded parameter values in
the objective function while value constraints can be implemented as
penalties, in this way: constraint c1:x1+2.0\*x2-3.0\*x3<=0 can be used to
adjust objective function value: cost+=(c1<=0?0:100000+c1*9999), with 100000
penalty base and 9999 constraint scale chosen to assure no interaction with the
expected "normal" objective function values while providing a useful gradient.
Note that if the solution's value is equal or higher than the penalty base
it means either a feasible solution was not found or the chosen constraint
scale does not generate a useful gradient. See constr.cpp for an example of
constraint programming. constr2.cpp is an example of non-linear constraint
programming with both non-equalities and equalities.

The minimal and maximal allowed parameter values (bounds) should be specified
in a way to cover a wider value range, in order to reduce boundary effects
that may reduce convergence. It maybe beneficial to specify bounds in a way so
that the expected optimum is located at the center of the search space.

## Examples ##

Use the example.cpp program to see the basic usage example.

test2.cpp is a convergence test for all available functions. Performs many
optimization attempts on all functions. Prints various performance
information, including percentage of rejected attempts (rejection rate).

test4.cpp is a convergence test for multi-dimensional functions.

constr.cpp and constr2.cpp programs demonstrate use of constraint penalties.

constr3.cpp demonstrates use of the "deep" optimization strategy.

## Development ##

While the basic algorithm of the strategy is finished, the built-in parameters
of the algorithm is an area of ongoing research. There are several things that
were discovered that may need to be addressed in the future:

1. Parallelization of this algorithm is technically possible, but may be
counter-productive (increases convergence time considerably). It is more
efficient to run several optimizers in parallel with different random seeds.

2. The default population size formula 12+Dim*2 works well for most functions,
however some functions converge faster if a higher population size is used.
On the other hand, using an overly large population size may also increase
convergence time.

3. The method currently uses "short-cuts" which can be considered as "tricks"
which are non-universal. However, optimization of some functions benefits from
them greatly, they increase optimization success of test suites by several
percent.

## Users ##

This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your project to the list of users.
