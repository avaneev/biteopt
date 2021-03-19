// High-dimensional convergence test of the BiteOpt method.
//
// On each line:
// I - average number of iterations taken in successful attempts.
// R - std.dev of the number of iterations in successful attempts.
// A - average required number of attempts to converge.
// C - minimal objective function value detected in all attempts.
// RC - average objective function value detected in all rejected attempts.
// PowerSum_10 - the name of the test function and dimensionality.

#include "tester.h"
#include "testfn.h"

int main()
{
	CTester Tester;

	Tester.init( 0.01, 20, 70000, true );
	Tester.addCorpus( 14, OptCorpusND, false, false );

	Tester.run();
}
