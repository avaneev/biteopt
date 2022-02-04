# BITEOPT - Derivative-Free Optimization Method #

## Introduction ##

BITEOPT is a free open-source stochastic non-linear bound-constrained
derivative-free optimization method (algorithm, heuristic, or strategy), for
global optimization. The name "BiteOpt" is an acronym for "BITmask Evolution
OPTimization".

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

BiteOpt uses self-optimization techniques making it objective
function-agnostic. In its inner workings, BiteOpt uses objective function
value's ranking, not the actual value. BiteOpt is a multi-faceted example of
a "probabilistic computing" system.

Python wrapper with `scipy.optimize`-alike API, currently in beta, is
available [here](https://scipybiteopt.readthedocs.io/en/latest/index.html).

## Comparison ##

This "black-box" optimization method was tested on 2000+ 1-60 dimensional
optimization problems and performed well, and it successfully solves even
600-dimensional test problems found in some textbooks. But the main focus of
the method is to provide fast solutions for computationally expensive
"black-box" problems of medium dimensionality (up to 60).

This method was compared with the results of this paper (on 244 published C
non-convex smooth problems, convex and non-convex non-smooth problems were not
evaluated): [Comparison of derivative-free optimization
algorithms](https://sahinidis.coe.gatech.edu/?q=dfocomp).
This method was able to solve 77% of non-convex smooth problems in 10
attempts, 2500 iterations each. It comes 2nd (very close to the 1st) in the
comparison on non-convex smooth problems (see Fig.9 in the paper). With a huge
iteration budget (up to 1 million) this method solves 97% of problems.

On a comparable test function suite and conditions outlined at this page:
[global_optimization](http://infinity77.net/global_optimization/multidimensional.html)
(excluding several ill-defined and overly simple functions, and including
several complex functions, use `test2.cpp` to run the test) this method's
attempt success rate is >94% (with 100% of functions solved) while the average
number of objective function evaluations is ~370.

At least in these comparisons, this method performs better than plain
CMA-ES which is also a well-performing stochastic optimization method. As of
version 2021.1, BiteOpt's "solvability" exceeds CMA-ES on synthetic function
sets that involve random coordinate axis rotations and offsets (e.g., [BBOB
suite](https://coco.gforge.inria.fr/)). BiteOptDeep (e.g., with M=6)
considerably outperforms CMA-ES in "solvability".

As a matter of sport curiosity, BiteOpt is able to solve, in reasonable time,
almost all functions proposed in classic academic literature on global
optimization. This is quite a feat for a derivative-free method (not to be
confused with large-scale analytic and gradient-based global optimization
methods). Of course, BiteOpt is capable of more than that. If you have a
reference to a function (with a known solution) published in literature that
BiteOpt can't solve, let the author know.

BiteOpt (state at commit 124) took 2nd place (1st by sum of ranks) in
[BBComp2018-1OBJ-expensive](https://www.ini.rub.de/PEOPLE/glasmtbl/projects/bbcomp/results/BBComp2018-1OBJ-expensive/summary.html)
competition track. Since the time of that commit the method improved in many
aspects, especially in low-dimensional convergence performance. Commit 124 can
be considered as "baseline effective" version of the method (it is also
maximally simple), with further commits implementing gradual improvements, but
also adding more complexity.

Also, BiteOpt (state at commit 256) took 2nd place (3rd by sum of ranks) in
[BBComp2019-1OBJ](https://www.ini.rub.de/PEOPLE/glasmtbl/projects/bbcomp/results/BBComp2019-1OBJ/summary.html)
competition track.

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

Since version 2021.1 BiteOpt uses a companion optimizer - SpherOpt - which
works independently and provides "reference points" to BiteOpt. Such companion
improves BiteOpt's convergence properties considerably, especially when the
parameter space is rotated. Since version 2021.15 BiteOpt uses an additional
companion optimizer - NMSeqOpt - which increases diversity of generated
solutions.

Since version 2021.3 BiteOpt became a self-optimizing method not requiring any
fune-tuning from the user nor the author.

## CBiteOptDeep (biteopt.h) ##

Deep optimization class. Based on an array of `M` CBiteOpt objects. This
"deep" method pushes the newly-obtained solution to the random CBiteOpt object
which is then optimized. This method, while increasing the convergence time,
is able to solve complex multi-modal functions.

This method is most effective on complex functions, possibly with noisy
fluctuations near the global solution, that are not very expensive to
calculate and that have a large iteration budget. Tests have shown that on
smooth functions that have many strongly competing minima this "deep" method
considerably increases the chance to find a global solution relative to the
CBiteOpt class, but still requires several runs with different random seeds.
When using this method, the required iteration budget usually increases by
a factor of M<sup>0.5</sup>, but the number of the required optimization
attempts usually decreases. In practice, it is not always possible to predict
the convergence time increase of the CBiteOptDeep class, but increase does
correlate to its `M` parameter. For some complex functions the use of
CBiteOptDeep even decreases convergence time. For sure, the CBiteOptDeep class
often produces better solutions than the CBiteOpt class.

## Notes ##

BiteOpt is a completely self-optimizing method. It does not feature
user-adjustable hyper-parameters. Even population size adjustments may not be
effective.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10 times is a minimal general requirement. The required number of
optimization attempts is usually proportional to the number of strongly
competing minima in a function. 

This method is hugely-probabilistic, and it depends on its initial state,
which is selected randomly. In most cases it is more efficient to rerun the
optimization with a new random seed than to wait for the optimization process
to converge. Based on the results of optimization of the test corpus, for
2-dimensional functions it is reasonable to expect convergence in 1000
iterations (in a successful attempt), for 10-dimensional functions it is
reasonable to expect convergence in 7000 iterations (harder functions may
require more iterations to converge). Most classic 2-dimensional problems
converge in 400 iterations or less, at 10<sup>-6</sup> precision.

Each run may generate an equally-usable candidate solution (not necessarily
having the least cost), so the researcher may select solution from any run
based on his/her own considerations. In this light, it may be incorrect to
assume that least-performing runs are "wasted". In practice, least-performing
runs may give more acceptable parameter values within the search space
compared to the best-performing runs.

Note that derivative-free optimization methods in general provide "asymptotic"
solutions for complex functions. Thus it is reasonable to assume that BiteOpt
gives an optimal solution with some implicit tolerance factor. Given a large
enough function evaluation budget, BiteOpt usually does find an optimal
solution which can be cross-checked with other solvers, but a solution of a
new unexplored function must be treated as "asymptotically optimal".

Also note that in some problem areas like [ESA GTOP](https://www.esa.int/gsp/ACT/projects/gtop/)
problem suite the attempt budget should as high as 2000 (beside using the
BiteOpt depth of at least 8). At the same time, iteration budget can be
kept reasonably moderate (100000). Despite a large attempt budget, on a 8-core
processor, this still allows one to get good solutions in a matter of minutes.

## Limitations ##

Rogue optimums may not be found by this method. A rogue optimum is an optimum
that has a very small, almost undetectable area of descent and is placed apart
from other competing minima. The method favors minimum with a larger area of
descent. The Damavandi test function is a perfect example of the limitation of
this method (this test function is solved by this method, but requires a lot
of iterations). In practice, however, rogue optimums can be considered as
undesired outliers that rely on unstable parameter values (if such parameters
are used in a real-world system that has a certain parameter value precision,
a system may leave the "rogue" optimal regime easily).

When the problem field requires one to locate such "rogue optimums", the best
approach is to use a magnitudes larger attempt budget (a so called "parallel
retry" approach).

To a small degree, this method is immune to noise in the objective function.
While this method was designed to be applied to continuous functions, it is
immune to discontinuities, and it can solve problems that utilize parameter
value rounding (integer parameters). This method can't acceptably solve
high-dimensional problems that are implicitly or explicitly combinatorial
(e.g., Perm, and Lennard-Jones atom clustering problems) as in such problems
the global descent vanishes at some point and the method is left with an
exponentially increasing number of local minima. Similarly, problems with many
competing minima without a pronounced global descent towards global minimum
(e.g., Bukin N.6 problem) may not be solved acceptably as in most cases they
require exhaustive search or a search involving knowledge of the structure of
the problem.

Difference between upper and lower parameter bound values should be specified
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
(2-dimensional) and `test3.cpp` (14-dimensional) test corpuses, less than
0.9% of attempts require more than 3\*sigma iterations, 54% of attempts
require less than the mean. A typical probability distribution of percent of
attempts/sigma is as follows (discretized, not centered around 0 because it
deviates from the standard distribution, the mean corresponds to 0\*sigma):

![PDF plot](https://github.com/avaneev/biteopt/blob/master/attempt_pdf_plot.png)

## Constraint Programming ##

Mixed integer programming can be achieved by using rounded parameter values in
the objective function. Note that using categorical variables may not be
effective, because they require combinatorial search. Binary variables may be
used, in small quantities (otherwise the problem usually transforms into
a combinatorial problem as well). While not very fast, BiteOpt is able to
solve binary combinatorial problems if the cost function is formulated as
a sum of differences between bit values and continuous variables in the range
[0; 1].

Equality and non-equality constraints can be implemented as penalties. The
author has found a general effective method to apply value constraints via
penalties. In the code below, `n_con` is the number of constraints,
`con_notmet` is the number of constraints not meeting tolerances, and the
`pn[]` is the array of linear penalty values for each constraint; a penalty
value should be set to 0 if it meets the tolerance. For derivative-free
methods, a suggested equality tolerance is 10<sup>-4</sup>, but a more common
10<sup>-6</sup> can be also used; lower values are not advised for use. Models
with up to 190 constraints, both equalities and non-equalities, were tested
with this method. In practice, on a large set of problems, this method finds a
feasible solution in up to 93% of cases.

	double pns = 0.0;
	int i;

	for( i = 0; i < n_con; i++ )
	{
		pns = pns * pow( 4.0, 1.0 / n_con ) + pn[ i ] +
			pn[ i ] * pn[ i ] + pn[ i ] * pn[ i ] * pn[ i ];
	}

    real_value = cost;
	cost += 1e10 * ( con_notmet + pns );

In essence, this method transforms each penalty value into a cubic penalty
value, places each penalty value into its own "stratum", and also applies a
"barrier value" that depends on the number of constraints not met. The barrier
value is suitably large for most practical constraint programming problems.

See `constr.cpp` for an example of constraint programming. `constr2.cpp` is an
example of non-linear constraint programming with both non-equalities and
equalities. To effectively solve constraint programming problems, the
CBiteOptDeep class should be used, with M=6 or higher.

It is not advisable to use constraints like (x1-round(x1)=0) commonly used
in model libraries to force integer or binary values, as such constraint
formulation does not provide a useful global gradient. Instead, direct
rounding should be used on integer variables.

## Convergence Proof ##

Considering the structure of the method and the fact that on every iteration
only improving solutions are accepted into the population, with ever
decreasing upper bound on the objective function value, it is logically
impossible for the method to be divergent. While it is strictly non-divergent,
the formal proof of ability of the method to converge is complicated, and
should be at least as good as partly random search and partly Differential
Evolution.

## Tested Uses ##

This optimization method was tested for the following applications beside
synthetic benchmarking:

* Hyper-parameter optimization of complex non-linear black-box systems.
Namely, [AVIR](https://github.com/avaneev/avir) image resizing algorithm's
hyper-parameters, digital audio limiter algorithm's parameters.

* Non-linear least-squares problems, see the calcHougen and calcOsborne
functions in the `testfn.h` file for example problems.

* BiteOptDeep was successfuly used for direct search of optimal short
symmetric FIR filters. Namely, in
[r8brain-free-src](https://github.com/avaneev/r8brain-free-src)
sample rate converter.

BiteOpt is also referenced in these research papers:

* [Password Strength Signaling: A Counter-Intuitive Defense Against Password
Cracking; Springer]
(https://link.springer.com/chapter/10.1007/978-3-030-90370-1_18)

* [The CIP2A�TOPBP1 axis safeguards chromosome stability and is a synthetic
lethal target for BRCA-mutated cancer; Nature Cancer]
(https://www.nature.com/articles/s43018-021-00266-w)

## Examples ##

Use the `example.cpp` program to see the basic usage example of C++ interface.

The `example2.cpp` program is a usage example of a simple C-like function
biteopt_minimize(). This is a minimization test for Hougen-Watson model for
reaction kinetics (non-linear least squares problem).

    int biteopt_minimize( const int N, biteopt_func f, void* data,
        const double* lb, const double* ub, double* x, double* minf,
        const int iter, const int M = 1, const int attc = 10,
        const int stopc = 0 )

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
    stopc Stopping criteria (convergence check). 0: off, 1: 64*N, 2: 128*N.

    This function returns the total number of function evaluations performed;
    useful if the "stopc>0" was used.

`test2.cpp` is a convergence test for all available functions. Performs many
optimization attempts on all functions. Prints various performance
information, including per-function and per-attempt success rates.

`test3.cpp` is a convergence test for multi-dimensional functions with random
axis rotations and offsets.

`test4.cpp` is a convergence test for multi-dimensional functions without
randomization.

`constr.cpp`, `constr2.cpp`, and `constr3.cpp` programs demonstrate the use of
constraint penalties.

## Development ##

There are several things that were discovered that may need to be addressed in
the future:

1. Parallelization of BiteOpt algorithm is technically possible, but may be
counter-productive (increases convergence time considerably). It is more
efficient to run several optimizers in parallel with different random seeds.
Specifically saying, it is possible (tested to be working on some code commits
before May 15, 2018) to generate series of candidate solutions, evaluate them
in parallel, and then update optimizer's state before generating a new batch
of candidate solutions. Later commits have changed the algorithm to a form
less suitable for such parallelization.

2. The method uses "short-cuts" which can be considered as "tricks"
(criticized in literature) which are non-universal, and reduce convergence
time out of proportion for many known test functions. These "short-cuts" are
not critically important to method's convergence properties, but they reduce
convergence time even for functions that do not have minimum at a point where
all arguments are equal. It just often happens that such "short-cuts" provide
useful "reference points" to the method. Removing these "short-cuts" will
increase average convergence time of the method, but in most cases won't
impact method's ability to find a global solution. "Short-cuts" are used only
in 4% of objective function evaluations on average.

3. The method uses LCG pseudo-random number generator due to its efficiency.
The method was also tested with a more statistically-correct PRNG and the
difference turned out to be negligible.

## Method's Philosophy ##

BiteOpt is an evolutionary optimization method. Unlike many established
optimization methods like CMA-ES where new populations are generated on each
iteration, with or without combining with the previous generation, BiteOpt
keeps and updates a single main population of solutions at any given time. A
new solution either replaces a worst solution or is discarded. In common terms
it means that population has some fixed "living space" which is only available
to the best fit (least cost) solutions. Structurally, this is similar to a
natural evolutionary environment which usually offers only a limited "living
space" to its members. Least fit members have little chance to stay in this
"living space".

BiteOpt uses several methods to generate new solutions, each method taking
various information from various internal populations. These methods are used
in a probabilistic manner without any predefined preference.

BiteOptDeep implements evolutionary method which can be seen in society and
nature: exchange of solutions between independent populations. Such exchange
allows to find better solutions in a teamwork of sufficiently diverse members,
it also reduces time (but not human-hours) to find a better solution. This
method is a model of Swiss presidency rotation (each independent population
represents an independent human).

The author did not originally employ results and reasoning available in papers
on Differential Evolution. Author's use of DE operation is based on
understanding that it provides an implicit gradient information. A candidate
solution is generated as a sum of best solution and a difference between a
random and the worst solution. Such difference between a random and the worst
solution generates a probabilistically correct step towards the minimum of a
function, relative to a better solution. Due to this understanding, it is
impossible to employ various DE variants in BiteOpt, only the difference
between high rank and low rank solutions generates a valuable information;
moreover, only a difference multiplied by a factor of 0.5 works in practice.
Also, since BiteOpt does not use random crossover in its DE-alike operations,
the used approach is closer to an intermix of Nelder-Mead ("reduction") and
DE (multi-vector "mutation").

BiteOpt is more like a stochastic meta-method, it is incorrect to assume it
leans towards some specific optimizer class: for example, it won't work
acceptably if only DE-alike solution generators are used by it. BiteOpt
encompasses Differential Evolution, Nelder-Mead, and author's original
SpherOpt, "bitmask inversion", and "bit mixing" solution generators. An
initial success with the "bitmask inversion" operation (coupled with a
stochastic "move" operation it looks quite a lot like a random search) was the
main driver for BiteOpt's further development.

## Method's Description ##

NOTE: as of version 2021.3 this topic is not yet up-to-date.

The algorithm consists of the following elements:

1. A cost-ordered population of previous solutions is maintained. A solution
is an independent parameter vector which is evolved towards a better solution.
On every iteration, the method utilizes a probabilistically-chosen candidate
solution generator. At start, solution vectors are initialized at the center
of the search space, using Gaussian sampling.

Parameter values are internally normalized to [0; 1] range and, to stay in
this range, are wrapped in a special manner before each function evaluation.
Algorithm uses an alike of probabilistic state automata (histograms) to switch
between algorithm flow paths, depending on the candidate solution acceptance.
In many instances candidate solution generators uses square of the random
variable: this has an effect of giving more weight to better solutions.

With some probability, an independent, algorithmically different, parallel
optimizer is engaged whose solution is evaluated for inclusion into the
population.

After each objective function evaluation, the highest-cost previous solution
is replaced using the upper bound cost constraint.

2. Depending on the `RandProb` probability, a single (or all) parameter value
randomization is performed using "bitmask inversion" operation (which is
approximately equivalent to `v=1-v` operation in normalized parameter space).
Below, _i_ is either equal to rand(1, N) or in the range [1; N], depending on
the `AllpProb` probability. `>>` is a bit shift-right operation, `MantSize` is
a constant equal to 54, `MantSizeSh` is a hyper-parameter that limits bit
shift operation range. Actual implementation is more complex as it uses
average of two such operations.

![equation](https://latex.codecogs.com/gif.latex?mask=(2^{MantSize}-1)\gg&space;\lfloor&space;rand(0\ldots1)^4\cdot&space;MantSizeSh\rfloor)

![equation](https://latex.codecogs.com/gif.latex?x_\text{new}[i]&space;=&space;\frac{\lfloor&space;x_\text{new}[i]\cdot&space;2^{MantSize}&space;\rfloor&space;\bigotimes&space;mask&space;}{2^{MantSize}})

Plus, with `CentProb` probability the move around a random previous solution
is performed, utilizing a TPDF random value. This operation is performed
twice.

![equation](https://latex.codecogs.com/gif.latex?x_\text{new}[i]=x_\text{new}[i]-rand_{TPDF}(-1\ldots1)\cdot&space;CentSpan\cdot&space;(x_\text{new}[i]-x_\text{rand}[i]))

With `RandProb2` probability an alternative randomization method is used
involving the best solution, centroid vector and a random solution.

![equation](https://latex.codecogs.com/gif.latex?x_\text{new}[i]=x_\text{new}[i]&plus;(-1)^{s}(x_\text{cent}[i]-x_\text{new}[i]),&space;\quad&space;i=1,\ldots,N,\\&space;\quad&space;s\in\{1,2\}=(\text{rand}(0\ldots1)<0.5&space;?&space;1:2))

3. (Not together with N.2) the "step in the right direction" operation is
performed using the random previous solution, chosen best and worst
solutions, plus a difference of two other random solutions. This is
conceptually similar to Differential Evolution's "mutation" operation. The
used worst solution is randomly chosen from 3 worst solutions.

![equation](https://latex.codecogs.com/gif.latex?x_\text{new}=x_\text{best}-\frac{(x_\text{worst}-x_\text{rand}-(x_\text{rand2}-x_\text{rand3}))}{2})

4. Alternatively, an "entropy bit mixing" method is used to create a candidate
solution. This method mixes (XORs) parameter values represented as raw bit
strings drawn from an odd number of parameter vectors. Probabilistically,
such composition creates a new random parameter vector with an overwhelming
number of bits common to better-performing solutions, and a fewer number of
bits without certainty.

5. With `ScutProb` probability a "short-cut" parameter vector change operation
is performed.

![equation](https://latex.codecogs.com/gif.latex?z=x_\text{new}[\text{rand}(1\ldots&space;N)])

![equation](https://latex.codecogs.com/gif.latex?x_\text{new}[i]=z,&space;\quad&space;i=1,\ldots,N)

## SMA-ES ##

This is a working optimization method called "SigMa Adaptation Evolution
Strategy". It has the same programmatic interface as the CBiteOpt class, so it
can be easily used in place of CBiteOpt.

SMA-ES is based on the same concept as CMA-ES, but performs vector sigma
adaptation. SMA-ES performs covariance matrix update like CMA-ES, but it
is a simple linear update using "leaky integrator" averaging filtering, not
adaptation. SMA-ES algorithm's operation is based on principles of control
signals.

The main difference to CMA-ES is that per-parameter sigmas are updated using
these elements:

1. Sigma auto-adapts due to weighted parameter covariance calculation. Better
fit solutions have more influence over expansion or contraction of the sigma.

2. SMA-ES approximates the "geometry" of the sample distribution. It ranges
from "spherical" to "needle" geometry (represented by a continuous `spc`
variable). When geometry is spherical, covariance matrix update filter is
tuned to an increased frequency (`CovUpdFast` instead of `CovUpdSlow`).

3. An asymmetry is introduced to the Gaussian sampling function, depending on
the centroid step size. Distribution is expanded in the direction of the step
and contracted in the opposite direction.

4. On every update, all per-parameter sigmas are contracted (multiplied) by
the `SigmaMulBase` coefficients, depending on sphericity. Additionally,
overly-contracted sigmas are expanded by the `SigmaMulExp` coefficient.

In overall, SMA-ES is a completely self-adaptive method. It has several
hyper-parameters that do not depend on problem's dimensionality.

Population size formula in SMA-ES is fixed to `13+Dims`: according to tests,
in average it suits all dimensionalities. Of course, particular problems may
converge better/faster with a lower or higher population size. The number of
objective function evaluations is twice the population size per sample
distribution update (best fit solutions enter the population): this aspect is
controlled via the `EvalFac` parameter, which adjusts method's overhead
with only a minor effect on convergence property. Method's typical
observed complexity is O(N<sup>1.6</sup>).

## SpherOpt ##

This is a "converging hyper-spheroid" optimization method. While it is not
as effective as, for example, CMA-ES, it also stands parameter space scaling,
offsetting and rotation well. Since version 2021.1 it is used as a companion
to BiteOpt, with excellent results.

This method is structurally similar to SMA-ES, but instead of Gaussian
sampling, SpherOpt selects random points on a hyper-spheroid (with a bit of
added jitter for lower dimensions). This makes the method very
computationally-efficient, but at the same time provides immunity to
coordinate axis rotations.

This method uses the same self-optimization technique as BiteOpt.

## NMSeqOpt ##

The CNMSeqOpt class implements sequential Nelder-Mead simplex method with
the "stall count" tracking. This optimizer is used as an additional parallel
optimizer in BiteOpt.

## DEOpt ##

The CDEOpt class implements a Differential Evolution-alike DFO solver, but in
the population-handling framework of BiteOpt. "best/3/bit". Mutation parameter
is fixed, equals to 0.25. Instead of a crossover, the method uses 1-bit
randomization with 50% probability. Population size is equal to 20\*Dims, by
default. Population is initialized with a Gaussian sampling.

## Thanks ##

The author would like to thank [Dietmar Wolz](https://github.com/dietmarwo)
for additional BiteOpt testing in the astrodynamics field. You may check out
the `fast-cma-es` for multi-threaded "parallel retry" wrapper for BiteOpt.

## Citing ##

```bibtex
@misc{biteopt2022,
    author = {Aleksey Vaneev},
    title = {{BITEOPT - Derivative-free optimization method}},
    note = {C++ source code, with description and examples},
    year = {2022},
    publisher = {GitHub},
    journal = {GitHub repository},
    howpublished = {Available at \url{https://github.com/avaneev/biteopt}},
}
```
