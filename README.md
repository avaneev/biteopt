# biteopt - "Bitmask Evolution" Function Optimizer #
## Introduction ##

### CBEOOptimizerFan (bitefan.h) ###

"Bitmask evolution" version "fan" optimization class. This strategy is
based on now outdated CBEOOptimizer and CBEOOptimizer2 stochastic
derivative-less strategies, and uses several current parameter vectors
("fan elements"). Highest cost "fan element" can be replaced with a new
solution if "fan element's" cost (plus some margin) is higher than that of
the new solution's. Having several "fan elements" allows parameter vectors
to be spaced apart from each other thus making them cover a larger
parameter search space collectively. The strategy was named as "bitmask
evolution", because at its core an operation of inversion of a random
segment of parameter value's lowest bits is used. Beside that, several
"step in the right direction" operations are used that move the solution
vector into position with a probably lower objective function value.

The benefit of this strategy is a relatively high robustness: it can
successfully optimize a wide range of test functions. Another benefit is a
low convergence time which depends on the complexity of the objective
function. Like many stochastic optimization strategies with fast
convergence, this strategy can't solve problems with narrow or rogue
optimums. Harder problems may require dozens of optimization attempts to
reach optimum.

### Notes ###

This strategy was tested on many classic 2-dimensional and 30-dimensional
optimization problems and performed well. Due to its design this strategy may
be particularly good at improving an existing sub-optimal local solution.

The FanSize parameter can be adjusted to increase quality of solutions,
especially in high dimensionality problems. By default, FanSize is equal to
the square of the number of dimensions divided by 3. FanSize controls
robustness of the strategy at the cost of convergence time. However, setting
FanSize too high may reduce convergence.

The minimal and maximal allowed parameter values should be specified in a way
to cover a wider value range, in order to reduce boundary effects that may
greatly reduce convergence. If the target local or global minimum stands
very close to the parameter value boundaries this strategy may fail to
converge.

It is usually necessary to run the optimization process several times with
different random seeds since the process may get stuck in a local minimum.
Running 10-20 times is a minimal requirement. The optimizePlateau() function
can be used to optimize functions that have a complex landscape. This function
works within the bounds of allocated iteration limit, and performs
re-initializations when the optimization process reaches the plateau.

Most hard constraints can be introduced by applying huge penalties to the
objective function. Even binary penalties like "if(x>1)cost+=10000" should
work acceptably in most cases.

Strategy's "robustness" is a multi-factor non-formal estimation which includes
average convergence time, standard deviation of convergence time, the set of
functions the method can solve successfully given randomized initial
conditions, in a given number of attempts. For global optimization robustness
also includes the lowest achieved cost.

Use the test.cpp program to see the basic usage example.

test2.cpp is a convergence test for all available functions. Performs many
optimization attempts on all functions. Prints various performance
information, including percentage of rejected attempts.

test3.cpp demonstrates the usage of the optimizePlateau() function.

test4.cpp is a convergence test for multi-dimensional functions.

## Development ##

While the basic algorithm of the strategy is finished, the built-in parameters
of the algorithm is an area of ongoing research. There are several things that
were discovered that may need to be addressed in the future:

1. The CostMult parameter must probably depend on the FanSize parameter. For
example, at FanSize=8 the best CostMult which gives a faster convergence is
around 2.5, at FanSize=40 the best CostMult is around 1.2. This is
understandable, because when FanSize is higher, the difference between the
lowest cost and highest cost "fan element" is naturally higher due to a larger
pool of "fan elements".

2. The formula of dependence of FanSize on the number of dimensions may need
to be updated to better suit varying dimensionality of the problems.

3. When the algorithm reaches convergence point, parameter vector used in
"non-bitmask inversion" round of function evaluation changes very little.
Which probably means that the "bitmask inversion" operation's frequency should
be increased when the convergence point approaches. Also the range of bits
in this operation should probably be shifted towards lower bits. This should
reduce convergence time when the optimum's area was correctly found, but may
also increase optimization rejection rate.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your project to the list of users.
