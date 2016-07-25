# biteopt - "Bitmask Evolution" Optimizer #
## Introduction ##

"Bitmask evolution" optimization class. Implements a very simple
evolutionary optimization method (strategy) which involves inversion of a
random segment of parameter value's lowest bits at each step. Additionally
includes a crossing-over operation which in some cases improves convergence
considerably. In some cases crossing-over reduces convergence, but only
slightly. For more robustness it is possible to assign several internal
values to each optimization parameter.

This strategy was tested on several classic 2-parameter optimization
problems and it performed fairly well. Global problems (with multiple local
minima) may not be handled well by this strategy, but in practice this
strategy strives to provide "minimum among minima" nevertheless.

The CBEOOptimizer2 class is a further evolution of this strategy. Additionally
includes the "step in the right direction" operation and a history of best
solutions. This version provides a quite fast convergence time, a very small
code size and minimal memory requirement. The only drawback is that this
strategy requires 3 instead of 2 random number generator calls per
parameter on each step.

The CBEOHive class implements "hive" optimization strategy which utilizes
several CBEOOptimizer objects in parallel exchanging solutions between
them. This can offer a benefit for some "hard" functions by reducing the
number of required function evaluations by a factor of 2. For other
functions there may be a negative benefit. CBEOHive is best used with
low (0.1-0.2) crossing-over probabilities and ValuesPerParam=3 or 4.

Use the test.cpp program to see the basic usage example.

test2.cpp is a more complex test which performs optimization of several
functions and calculates the average convergence time. Can be used with
CBEOOptimizer, CBEOOptimizer2 and CBEOHive classes.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your product to the list of users.
