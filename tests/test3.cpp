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

const CTestFn* TestFn[] = {
	&TestFnSchwefel222, &TestFnRosenbrock, &TestFnStretchedV, NULL };

int main()
{
	CTester Tester;

//	Tester.init( 0.01, 26, 2500000, true );
//	Tester.addCorpus( 56, TestFn, true, true );

	Tester.init( 0.01, 26, 170000, true );
	Tester.addCorpus( 14, OptCorpusNDRotOfs, true, true );
	Tester.addCorpus( 14, OptCorpusNDRotOfs, true, false );
//	Tester.addCorpus( 14, OptCorpusND, true, false );
//	Tester.addCorpus( 14, OptCorpusNDRotOfsSol, true, true );

	Tester.run();
}
