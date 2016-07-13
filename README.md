# biteopt - Bitmask Evolution Optimizer #
## Introduction ##

Bitmask evolution optimization class. Implements a very simple evolution
strategy which involves inversion of a segment of parameter value's bits
at each step. Additionally includes crossing over operation which in some
cases improves convergence considerably. In some cases crossing over
reduces convergence, but only slightly. For more robustness it is possible
to assign several internal values to each optimization parameter.

This strategy was tested on several classic 2-parameter optimization
problems and it performed fairly well. Global (multiple local minima)
problems may not be handled well by this strategy, but in practice this
strategy strives to provide "minimum among minimums" nevertheless.

Use the test.cpp program to see the usage example.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your product to the list of users.
