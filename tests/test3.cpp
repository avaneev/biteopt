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

#if defined( _WIN32 ) || defined( _WIN64 )
	#include <windows.h>
	#define USEPERF 1
#endif // defined( _WIN32 ) || defined( _WIN64 )

int main()
{
	CTester Tester;

	Tester.init( 0.01, 35, 14000, true );
	Tester.addCorpus( 14, OptCorpusNDRotOfs, true, false );
//	Tester.addCorpus( 14, OptCorpusNDRotOfsSol, true, true );

	#if USEPERF
	LARGE_INTEGER Freq;
	QueryPerformanceFrequency( &Freq );
	LARGE_INTEGER t1;
	QueryPerformanceCounter( &t1 );
	#endif // USEPERF

	Tester.run();

	#if USEPERF
	LARGE_INTEGER t2;
	QueryPerformanceCounter( &t2 );
	printf("time: %f s\n", ( t2.QuadPart - t1.QuadPart ) /
		(double) Freq.QuadPart );
	#endif // USEPERF
}
