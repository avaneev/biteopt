# BiteOpt - Stochastic Function Optimizer #

## Contents ##

* [Introduction](#introduction)
* [Comparison](#comparison)
* [Notes](#notes)
* [CBiteOpt (biteopt.h)](#cbiteopt-biteopth)
* [CBiteOptDeep (biteopt.h)](#cbiteoptdeep-biteopth)
* [Limitations](#limitations)
* [Constraint programming](#constraint-programming)
* [Convergence proof](#convergence-proof)
* [Tested uses](#tested-uses)
* [Examples](#examples)
* [Users](#users)
* [Method description](#method-description)

## Introduction ##

BITEOPT is a free open-source stochastic non-linear bound-constrained
derivative-free optimization strategy. The name "BiteOpt" is an acronym for
"BITmask Evolution OPTimization".

The benefit of this strategy is a relatively high robustness: it can
successfully optimize a wide range of multi-dimensional test functions.
Another benefit is a low convergence time which depends on the complexity
of the objective function. Hard (multi-modal) problems may require many
optimization attempts to reach optimum.

Instead of iterating through different "starting guesses" to find optimum
like in deterministic strategies, this strategy requires optimization attempts
with different random seeds. The stochastic nature of the method allows it to
automatically "fall" into different competing minima with each run. If there
are no competing minima in a function (or the true/global minimum is rogue and
cannot be detected), this strategy in absolute majority of runs will return
the same optimum.

## Comparison ##

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
solves 93% of problems (CBiteOptDeep solves 95% with M=25). On a comparable
test function suite and conditions outlined at this page:
[global_optimization](http://infinity77.net/global_optimization/multidimensional.html)
(excluding several ill-defined and overly simple functions, and including
several complex functions, use test2.cpp to run the test) this strategy's
success rate is >91% while the average number of objective function
evaluations is ~340.

## Notes ##

Strategy's hyper-parameters (probabilities) were pre-selected and should not
be changed.

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

## CBiteOpt (biteopt.h) ##

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

## CBiteOptDeep (biteopt.h) ##

Deep stochastic optimization class. Based on an array of M CBiteOpt
objects. This "deep" strategy pushes the newly-obtained solution to the
next CBiteOpt object which is then optimized. This strategy while
increasing the convergence time by a factor of about sqrt(M) is able to
solve even a very noisy non-linear functions.

This strategy is most effective on stochastic functions or functions with
huge fluctuations near the global solution that are not very expensive to
calculate and that have a large iteration budget. Tests have shown that on
smooth functions that have many strongly competing minima this strategy
increases the chance to find a global solution by a factor of sqrt(M)
relative to the CBiteOpt class, but still requires several runs at
different random seeds. However, the number of required runs in most cases
is lower by about sqrt(M). So, in many cases it's more efficient to
increase the iteration budget by a factor of sqrt(M) which will in turn
increase the chance to find a global optimum by a factor of sqrt(M) and
also reduce the number of optimization attempts by the same number.  In
practice, the chance to find a global optimum is increased even more
considerably with this "deep" strategy.

## Limitations ##

The required number of optimization attempts is usually proportional to the
number of strongly competing minima in a function. Rogue optimums may not be
found by this strategy. A rogue optimum is an optimum that has a very small,
almost undetectable area of descent and is placed apart from other competing
minima. The strategy favors minimum with a larger area of descent. The
Damavandi test function is a perfect example of the limitation of this
strategy (this test function is solved by this strategy, but requires a lot
of iterations). In practice, however, rogue optimums can be considered as
undesired outliers that rely on unstable parameter values (if such parameters
are used in real-world system that has a certain parameter value precision, a
system may leave the "rogue" optimal regime easily).

To some degree this strategy is immune to noise in the objective function.
While this strategy was designed to be applied to continuous functions, it is
also immune to discontinuities to some degree, and it can solve problems that
utilize parameter value rounding (integer parameters). This strategy can't
acceptably solve high-dimensional problems that are implicitly or explicitly
combinatorial (e.g. Perm and Lennard-Jones atom clustering problems). Also
problems with many competing minima without a pronounced global descent
towards global minimum (e.g. Bukin N.6 problem) may not be solved acceptably
as in most cases they require exhaustive search.

## Constraint programming ##

Mixed integer programming can be achieved by using rounded parameter values in
the objective function while value constraints can be implemented as
penalties, in this way: constraint c1:x1+2.0\*x2-3.0\*x3<=0 can be used to
adjust objective function value: cost+=(c1<=0?0:100000+c1*9999), with 100000
penalty base and 9999 constraint scale chosen to assure no interaction with the
expected "normal" objective function values while providing a useful gradient.
Note that if the solution's value is equal to or higher than the penalty base
it means either a feasible solution was not found or the chosen constraint
scale does not generate a useful gradient. See constr.cpp for an example of
constraint programming. constr2.cpp is an example of non-linear constraint
programming with both non-equalities and equalities.

The minimal and maximal allowed parameter values (bounds) should be specified
in a way to cover a wider value range, in order to reduce boundary effects
that may reduce convergence. It may be beneficial to specify bounds in a way
so that the expected optimum is located at the center of the search space.

## Convergence proof ##

Considering the structure of the strategy and the fact that on every
iteration only improving solutions are accepted into the population, with ever
decreasing upper bound on the objective function value, it is logically
impossible for the strategy to be divergent.

## Tested uses ##

This optimization strategy was tested for the following applications beside
synthetic benchmarking:

1. Hyperparameter optimization of complex non-linear black-box systems.

2. Non-linear least-squares problems, see calcHougen function in testfn.h
for an example problem.

## Examples ##

Use the example.cpp program to see the basic usage example of C++ interface.

example2.cpp program is an example of a simple C-like function
biteopt_minimize(). This is a minimization test for Hougen-Watson model for
reaction kinetics. Non-linear least squares problem.

    void biteopt_minimize( const int N, biteopt_func f, void* data,
        const double* lb, const double* ub, double* x, double* minf,
        const int iter, const int M = 1, const int attc = 10 )

    N     The number of parameters in an objective function.
    f     Objective function.
    data  Objective function's data.
    lb    Lower bounds of obj function parameters, should not be infinite.
    ub    Upper bounds of obj function parameters, should not be infinite.
    x     Minimizer.
    minf  Minimizer's value.
    iter  The number of iterations to perform in a single attempt.
          Corresponds to the number of obj function evaluations that are performed.
    M     Depth to use, 1 for plain CBiteOpt algorithm, >1 for CBiteOptDeep
          algorithm. Internally multiplies "iter" by sqrt(M).
    attc  The number of optimization attempts to perform.

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

2. The default population size formula 12+Dim\*2 works well for most
non-convex functions, however some functions converge better if a higher
population size is used. Tests have revealed that convex problems converge
better with a population as large as 12+Dim\*6, but possibly at the expense of
a longer convergence time.

3. The method currently uses "short-cuts" which can be considered as "tricks"
which are non-universal. However, optimization of some functions benefits from
them greatly, they increase optimization success of test suites by several
percent.

## Users ##

This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your project to the list of users.

## Method description ##

The algorithm consists of the following elements:

1. A cost-ordered population of previous solutions is maintained. A solution
is an independent parameter vector which is evolved towards a better solution.
On every iteration, the best solution is evolved.

![equation1](http://quicklatex.com/cache3/8b/ql_66f6b6e1379996e5dceb22a4238fb08b_l3.png)

Below, _i_ is either equal to rand(1, N) or in the range [1; N], depending on
the `AllpProb` probability. Probabilities are defined in the range [0; 1] and
in many instances in the code were replaced with simple resetting counters for
more efficiency. Parameter values are internally normalized to [0; 1] range
and, to stay in this range, are wrapped in a special manner before each
function evaluation. Algorithm's hyper-parameters (probabilities) were
pre-selected and should not be changed.

2. Depending on the `RandProb` probability, a single (or all) parameter value
randomization is performed using "bitmask inversion" operation.

![equation2](http://quicklatex.com/cache3/15/ql_0791b3be35d15c7e3a987e558dd2b615_l3.png)

![equation3](http://quicklatex.com/cache3/70/ql_fe4d99ed011962956503c80053908b70_l3.png)

![equation4](http://quicklatex.com/cache3/51/ql_c1f36edbe9027044ef6d120b4a2b8451_l3.png)

Plus, with `CentProb` probability the random "step in the right direction" operation is performed using the centroid vector, twice.

![equation5](http://quicklatex.com/cache3/3d/ql_faa3b8c6fd883fb228f68b258fb3b33d_l3.png)

![equation6](http://quicklatex.com/cache3/0f/ql_c9ea3d70b8eb29a1c0e112f9d15a5b0f_l3.png)

![equation7](http://quicklatex.com/cache3/24/ql_ccfec27c2b5a32348a59a84caa552924_l3.png)

![equation8](http://quicklatex.com/cache3/cf/ql_40565ca7f1625f7e0f28d06f3f6a8acf_l3.png)

With `RandProb2` probability an alternative randomization method is used.

![equation9](http://quicklatex.com/cache3/f7/ql_b2e88e9d9710ddf91a1561a021ee66f7_l3.png)

3. (Not together with N.2) the "step in the right direction" operation is
performed using the random previous solution, current best and worst
solutions. This is conceptually similar to Differential Evolution's "mutation"
operation.

![equation10](http://quicklatex.com/cache3/58/ql_a15224c3709fd3fda67760646e7b7e58_l3.png)

4. With `ScutProb` probability a "short-cut" parameter vector change operation
is performed.

![equation11](http://quicklatex.com/cache3/da/ql_3c8871a6c42a602b2bfac8c1c89affda_l3.png)

![equation12](http://quicklatex.com/cache3/67/ql_b434d4517882743a2528d2e10b7ff767_l3.png)

5. After each objective function evaluation, the highest-cost previous
solution is replaced using the upper bound cost constraint.
