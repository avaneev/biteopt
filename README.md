# BITEOPT - Derivative-Free Optimization Method #

## Contents ##

* [Introduction](#introduction)
* [Comparison](#comparison)
* [CBiteOpt (biteopt.h)](#cbiteopt-biteopth)
* [CBiteOptDeep (biteopt.h)](#cbiteoptdeep-biteopth)
* [Notes](#notes)
* [Limitations](#limitations)
* [Constraint programming](#constraint-programming)
* [Convergence proof](#convergence-proof)
* [Tested uses](#tested-uses)
* [Examples](#examples)
* [Users](#users)
* [Method description](#method-description)

## Introduction ##

BITEOPT is a free open-source stochastic non-linear bound-constrained
derivative-free optimization method (heuristic or strategy). The name
"BiteOpt" is an acronym for "BITmask Evolution OPTimization".

The benefit of this method is a relatively high robustness: it can
successfully optimize a wide range of multi-dimensional test functions.
Another benefit is a low convergence time which depends on the complexity
of the objective function. Hard (multi-modal) problems may require many
optimization attempts to reach optimum.

Instead of iterating through different "starting guesses" to find optimum
like in deterministic methods, this method requires optimization attempts
(runs) with different random seeds. The stochastic nature of the method allows
it to automatically "fall" into different competing minima with each run. If
there are no competing minima in a function (or the true/global minimum is
rogue and cannot be detected), this method in absolute majority of runs will
return the same optimum.

## Comparison ##

This "black-box" optimization algorithm was tested on 450+ 1-10 dimensional
optimization problems and performed well, and it successfully solves even
600-dimensional test problems found in some textbooks. But the main focus of
the method is to provide fast solutions for computationally expensive
"black-box" problems of medium dimensionality (up to 60).

This method was compared with the results of this paper (on 244 published C
non-convex smooth problems, convex and non-convex non-smooth problems were not
evaluated): [Comparison of derivative-free optimization algorithms](http://archimedes.cheme.cmu.edu/?q=dfocomp).
This method was able to solve 77% of non-convex smooth problems in 10
attempts, 2500 iterations each. It comes 2nd in the comparison on non-convex
smooth problems (see Fig.9 in the paper). With a huge iteration budget (up to 1
million) this method solves 95% of problems.

On a comparable test function suite and conditions outlined at this page:
[global_optimization](http://infinity77.net/global_optimization/multidimensional.html)
(excluding several ill-defined and overly simple functions, and including
several complex functions, use `test2.cpp` to run the test) this method's
success rate is >90% while the average number of objective function
evaluations is ~330.

At least in these comparisons, this method performs better than plain
CMA-ES which is also a well-performing stochastic optimization method. CMA-ES
outperforms CBiteOptDeep (e.g. with M=9) on synthetic function sets that
involve random coordinate axis rotations and offsets (e.g.
[BBOB suite](http://coco.gforge.inria.fr/)). Plain CBiteOpt lags far behind in
such benchmark arrangements. However, BiteOpt's development from its inception
was based on a wider selection of functions proposed by global optimization
researchers, without focus on synthetic transformations.

## CBiteOpt (biteopt.h) ##

BiteOpt optimization class. Implements a stochastic non-linear
bound-constrained derivative-free optimization method. It maintains a
cost-ordered population list of previously evaluated solutions that are
evolved towards a lower cost (objective function value). On every iteration,
the highest-cost solution in the list can be replaced with a new solution, and
the list reordered. A population of solutions allows the method to space
solution vectors apart from each other thus making them cover a larger
parameter search space collectively. Beside that, a range of parameter
randomization and the "step in the right direction" (Differential Evolution
"mutation") operations are used that move the solutions into positions with a
probabilistically lower objective function value.

## CBiteOptDeep (biteopt.h) ##

Deep optimization class. Based on an array of M CBiteOpt objects. This "deep"
method pushes the newly-obtained solution to the next CBiteOpt object which
is then optimized. This method while increasing the convergence time is able
to solve complex multi-modal functions.

This method is most effective on complex functions, possibly with noisy
fluctuations near the global solution that are not very expensive to calculate
and that have a large iteration budget. Tests have shown that on smooth
functions that have many strongly competing minima this "deep" method
considerably increases the chance to find a global solution relative to the
CBiteOpt class, but still requires several runs at different random seeds.
When using this method, the iteration budget increases but the number of
required optimization attempts usually decreases. In practice, it is not
always possible to predict the convergence time increase of the CBiteOptDeep
class, but increase does correlate to its M parameter. For sure, CBiteOptDeep
class often produces better solutions than the CBiteOpt class.

## Notes ##

Method's hyper-parameters (probabilities) were pre-selected and should not
be changed.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10 times is a minimal general requirement. This method is hugely
probabilistic and it depends on its initial state, which is selected randomly.
In most cases it is more efficient to rerun the optimization with a new random
seed than to wait for the optimization process to converge. Based on the
results of optimization of the test corpus, for 2-dimensional functions it is
reasonable to expect convergence in 1000 iterations (in a successful attempt),
for 10-dimensional functions it is reasonable to expect convergence in 5000
iterations (harder functions may require more iterations to converge). Most
classic 2-dimensional problems converge in 400 iterations or less, at 1e-6
precision.

## Limitations ##

The required number of optimization attempts is usually proportional to the
number of strongly competing minima in a function. Rogue optimums may not be
found by this method. A rogue optimum is an optimum that has a very small,
almost undetectable area of descent and is placed apart from other competing
minima. The method favors minimum with a larger area of descent. The
Damavandi test function is a perfect example of the limitation of this
method (this test function is solved by this method, but requires a lot
of iterations). In practice, however, rogue optimums can be considered as
undesired outliers that rely on unstable parameter values (if such parameters
are used in real-world system that has a certain parameter value precision, a
system may leave the "rogue" optimal regime easily).

To a small degree this method is immune to noise in the objective function.
While this method was designed to be applied to continuous functions, it is
immune to discontinuities, and it can solve problems that utilize parameter
value rounding (integer parameters). This method can't acceptably solve
high-dimensional problems that are implicitly or explicitly combinatorial
(e.g. Perm and Lennard-Jones atom clustering problems) as in such problems the
global descent vanishes at some point and the method is left with an
exponentially increasing number of local minima. Similarly, problems with many
competing minima without a pronounced global descent towards global minimum
(e.g. Bukin N.6 problem) may not be solved acceptably as in most cases they
require exhaustive search or a search involving knowledge of the structure of
the problem.

Difference between upper and lower parameter bound values affects precision of
the method. To increase precision, this difference should be kept as low as
possible: for example, [-100; 100] bounds should be used instead of
[-1000; 1000], if possible. Other than that, the bounds should be specified
in a way to cover a wider value range, in order to reduce boundary effects
that may reduce convergence.

Tests have shown, that in comparison to stochastic method like CMA-ES,
BiteOpt's convergence time varies more from attempt to attempt. For example,
on some problem CMA-ES's average convergence time may be 7000 iterations +/-
1400 while BiteOpt's may be 7000 +/- 3000. Such higher standard deviation
is mostly a negative property if only a single optimization attempt is
performed since it makes required budget unpredictable. But if several
attempts are performed, it is a positive property: it means that in some
optimization attempts BiteOpt converges faster and may find a better optimum
with the same iteration budget per attempt. Based on `test2.cpp`
(2-dimensional) and `test4.cpp` (10-dimensional) test corpuses, only about
1% of attempts require more than 3\*sigma iterations, 58% of attempts require
less than 0\*sigma. A typical probability distribution of percent of
attempts/sigma is as follows (discretized, not centered around 0 because it
deviates from standard distribution, average value corresponds to 0\*sigma):

![PDF plot](https://github.com/avaneev/biteopt/blob/master/attempt_pdf_plot.png)

## Constraint programming ##

Mixed integer programming can be achieved by using rounded parameter values in
the objective function. Note that using categorical variables may not be
effective, because they require an exhaustive search. Binary variables may be
used, in small quantities (otherwise the problem usually transforms into
a complex combinatorial problem).

Value constraints can be implemented as penalties, in this way: constraint
c1:x1+2.0\*x2-3.0\*x3<=0 can be used to adjust objective function value:
cost+=(c1<=0?0:100000+c1*9999), with 100000 penalty base (barrier) and 9999
constraint scale chosen to assure no interaction with the expected "normal"
objective function values while providing a useful gradient. Note that if the
solution's value is equal to or higher than the penalty base it means either a
feasible solution was not found or the chosen constraint scale does not
generate a useful gradient. See `constr.cpp` for an example of constraint
programming. `constr2.cpp` is an example of non-linear constraint programming
with both non-equalities and equalities. In practice, using quadratic
penalties may be more effective than adding the aforementioned barrier
constant.

It is not advisable to use constraints like (x1-round(x1)=0) commonly used
in model libraries to force integer or binary values, as such constraint
formulation does not provide a global descent. Instead, direct rounding should
be used on integer variables.

## Convergence proof ##

Considering the structure of the method and the fact that on every iteration
only improving solutions are accepted into the population, with ever
decreasing upper bound on the objective function value, it is logically
impossible for the method to be divergent. While it is strictly non-divergent,
the formal proof of ability of the method to converge is complicated, and
should be at least as good as partly random search and partly Differential
Evolution.

## Tested uses ##

This optimization method was tested for the following applications beside
synthetic benchmarking:

1. Hyperparameter optimization of complex non-linear black-box systems.
Namely, [AVIR](https://github.com/avaneev/avir) image resizing algorithm's
hyper-parameters, BiteOpt's own hyper-parameters, digital audio limiter
algorithm's parameters.

2. Non-linear least-squares problems, see calcHougen and calcOsborne functions
in `testfn.h` for example problems.

## Examples ##

Use the `example.cpp` program to see the basic usage example of C++ interface.

`example2.cpp` program is an example of a simple C-like function
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
          algorithm. Expected range is [1; 36]. Internally multiplies "iter"
          by sqrt(M).
    attc  The number of optimization attempts to perform.

`test2.cpp` is a convergence test for all available functions. Performs many
optimization attempts on all functions. Prints various performance
information, including percentage of rejected attempts (rejection rate).

`test3.cpp` is a convergence test for multi-dimensional functions with random
axis rotations and offsets.

`test4.cpp` is a convergence test for multi-dimensional functions without
randomization.

`constr.cpp` and `constr2.cpp` programs demonstrate use of constraint
penalties.

`constr3.cpp` demonstrates use of the "deep" optimization method.

## Development ##

While the basic algorithm of the method is finished, the built-in
hyper-parameters of the algorithm is an area of ongoing research. There are
several things that were discovered that may need to be addressed in the
future:

1. Parallelization of this algorithm is technically possible, but may be
counter-productive (increases convergence time considerably). It is more
efficient to run several optimizers in parallel with different random seeds.
Specifically saying, it is possible (tested to be working on some commits
before May 15, 2018) to generate a serie of candidate solutions, evaluate them
in parallel, and then update optimizer's state before generating a new batch
of candidate solutions. Later commits have changed the algorithm to a from
less suitable for such parallelization.

2. The default population size formula 10+Dim\*2 works well for most
non-convex functions, however some functions converge better if a higher
population size is used. Tests have revealed that convex problems converge
better with a population as large as 10+Dim\*6, but possibly at the expense of
a longer convergence time.

3. The method currently uses "short-cuts" which can be considered as "tricks"
which are non-universal. However, optimization of some functions benefits from
them greatly, they increase optimization success of test suites by several
percent. These "short-cuts" are not critically important to method's
convergence properties, but they increase convergence speed even for functions
that do not have minimum at a point where all arguments are equal.

4. The method uses resetting counters instead of direct probability
evaluation. This was done to reduce method's overhead (it was important to
keep overhead low so that optimization of method's own hyper-parameters does
not take too much time). In practice, resetting counters are equivalent to
`if( getRndValue() < Probability )` constructs, and actually provide slightly
better convergence properties, probably due to some "state automata" effect.

5. The method uses LCG pseudo-random number generator due to its efficiency.
The method was also tested with a more statistically-correct PRNG and the
difference turned out to be negligible.

## Users ##

This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your project to the list of users.

## Method description ##

The algorithm consists of the following elements:

1. A cost-ordered population of previous solutions is maintained. A solution
is an independent parameter vector which is evolved towards a better solution.
On every iteration, one of the 4 best solutions is evolved (best selection
allows method to be less sensitive to noise). At start, solution vectors
are initialized almost on hyper-box boundaries.

```math
x_\text{new}=x_\text{best}
```

Below, $`i`$ is either equal to $`rand(1, N)`$ or in the range $`[1; N]`$,
depending on the `AllpProb` probability. Probabilities are defined in the
range $`[0; 1]`$ and in many instances in the code were replaced with simple
resetting counters for more efficiency. Parameter values are internally
normalized to $`[0; 1]`$ range and, to stay in this range, are wrapped in a
special manner before each function evaluation. Algorithm's hyper-parameters
(probabilities) were pre-selected and should not be changed. Algorithm uses an
alike of state automata to switch between different probability values
depending on the candidate solution acceptance.

2. Depending on the `RandProb` probability, a single (or all) parameter value
randomization is performed using "bitmask inversion" operation.

![equation2](http://quicklatex.com/cache3/15/ql_0791b3be35d15c7e3a987e558dd2b615_l3.png)

![equation3](http://quicklatex.com/cache3/70/ql_fe4d99ed011962956503c80053908b70_l3.png)

![equation4](http://quicklatex.com/cache3/51/ql_c1f36edbe9027044ef6d120b4a2b8451_l3.png)

Plus, with `CentProb` probability the move around a random previous solution
is performed, utilizing a TPDF random value.

![equation5](http://quicklatex.com/cache3/3d/ql_faa3b8c6fd883fb228f68b258fb3b33d_l3.png)

![equation6](http://quicklatex.com/cache3/0f/ql_c9ea3d70b8eb29a1c0e112f9d15a5b0f_l3.png)

![equation7](http://quicklatex.com/cache3/24/ql_ccfec27c2b5a32348a59a84caa552924_l3.png)

![equation8](http://quicklatex.com/cache3/cf/ql_40565ca7f1625f7e0f28d06f3f6a8acf_l3.png)

With `RandProb2` probability an alternative randomization method is used
involving the current best solution and centroid vector.

![equation9](http://quicklatex.com/cache3/f7/ql_b2e88e9d9710ddf91a1561a021ee66f7_l3.png)

3. (Not together with N.2) the "step in the right direction" operation is
performed using the random previous solution, chosen best and worst
solutions, plus a difference of two random solutions. This is conceptually
similar to Differential Evolution's "mutation" operation. The used worst
solution is chosen from 3 worst solutions.

![equation10](http://quicklatex.com/cache3/58/ql_a15224c3709fd3fda67760646e7b7e58_l3.png)

4. With `ScutProb` probability a "short-cut" parameter vector change operation
is performed.

![equation11](http://quicklatex.com/cache3/da/ql_3c8871a6c42a602b2bfac8c1c89affda_l3.png)

![equation12](http://quicklatex.com/cache3/67/ql_b434d4517882743a2528d2e10b7ff767_l3.png)

5. After each objective function evaluation, the highest-cost previous
solution is replaced using the upper bound cost constraint.
