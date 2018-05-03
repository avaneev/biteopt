// High-dimensional convergence test of the BiteOpt strategy.
//
// On each line:
// AI - average number of iterations taken in successful attempts.
// RI - std.dev of the number of iterations in successful attempts.
// At - average required number of attempts to converge.
// C - minimal objective function value detected in all successful attempts.
// RjC - minimal objective function value detected in all rejected attempts.
// PowerSum_10 - the name of the test function and dimensionality.

#include <stdio.h>
#include "tester.h"

int main()
{
	CTester Tester;

	rnd.init( 1 );
	Tester.init( 10, OptCorpusNDRotOfs, 0.01, 15, 70000, true, true );

	Tester.run();
}
