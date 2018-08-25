// High-dimensional convergence test of the BiteOpt method.
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

	Tester.init( 14, OptCorpusNDRotOfs, 0.01, 150, 50000, true, false, true );

	Tester.run();
}
