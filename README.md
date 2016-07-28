# biteopt - "Bitmask Evolution" Optimizer #
## Introduction ##

### CBEOOptimizer ###

"Bitmask evolution" optimization class. Implements a very simple
evolutionary optimization method (strategy) which involves inversion of a
random segment of parameter value's lowest bits at each step. Additionally
includes a crossing-over operation which in some cases improves convergence
considerably. For more robustness it is possible to assign several internal
values to each optimization parameter. This strategy is associated with
a very small code size and minimal memory requirement.

### CBEOOptimizer2 ###

The CBEOOptimizer2 class is a further evolution of this strategy. Additionally
includes the "step in the right direction" operation. This version provides a
faster convergence time.

### CBEOOptimizerFan ###

This strategy is based on the CBEOOptimizer2 strategy, but uses several
current parameter vectors ("fan elements"). Any parameter vector can be
replaced with a new solution if parameter vector's cost is higher than that
of the new solution's and the "distance" of the new solution is not
considerably low. The "distance" constraint allows parameter vectors to be
spaced apart from each other thus making them cover a larger parameter search
space collectively. The "fan elements" are used unevenly: some are used more
frequently than the others.

The benefit of this strategy is increased robustness: it can optimize
successfully a wider range of functions. Another benefit is a considerably
decreased convergence time in deeper optimizations.

This strategy is associated with a high overhead per function evaluation.
Due to this fact, for simple functions and not deep optimization it may be
more beneficial to use the CBEOOptimizer2 class.

### * * * ###

All these strategies were tested on several classic 2-parameter optimization
problems and performed fairly well. Global problems (with multiple local
minima) may not be handled well by these strategies, but in practice these
strategies strive to provide the "minimum among minima" nevertheless.

Use the test.cpp program to see the basic usage example.

test2.cpp is a more complex test which performs optimization of several
functions and calculates the average convergence time. Can be used with
CBEOOptimizer, CBEOOptimizer2 and CBEOOptimizerFan classes.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your product to the list of users.
