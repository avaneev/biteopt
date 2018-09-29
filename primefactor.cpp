// Prime-number factoring example. Uses BiteOptDeep with a huge (5000) depth
// factor. This is purely a proof-of-concept example, may not work for
// factoring of huge numbers, so far works for 17-digit numbers. Since this
// method is probabilistic, various variables should be changed to obtain
// solutions - iteration number, depth, attempt count, MinMaxFactor values,
// Split. The resulting "minf" is equal to 0 when a precise solution was
// found.
//
// Currently method uses just an even-odd check as a primality test. Using a
// more elaborative and precise test may produce better convergence property.

#include <stdio.h>
#include <stdint.h>
#include "biteopt.h"

const int64_t val = 160000003LL * 148030163LL; // n=pq
const int Split = 4;
const int N = Split * 2;
const double lb[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
const double ub[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
const double MinFactor1 = 0.0;
const double MaxFactor1 = 400000000.0;
const double Factor1Diff = ( MaxFactor1 - MinFactor1 ) / Split;
const double MinFactor2 = 0.0;
const double MaxFactor2 = 400000000.0;
const double Factor2Diff = ( MaxFactor2 - MinFactor2 ) / Split;
bool DoPrint = false;

double fn( int ND, const double* x, void* data )
{
	int64_t a1 = 0;
	int64_t a2 = 0;
	double cost = 0.0;
	int i;

	for( i = 0; i < Split; i++ )
	{
		a1 += (int64_t) ( MinFactor1 + x[ i ] * Factor1Diff + 0.5 );
		a2 += (int64_t) ( MinFactor2 + x[ Split + i ] * Factor2Diff + 0.5 );
	}

	if(( a1 & 1 ) == 0 )
	{
		cost += 1e10;
	}

	if(( a2 & 1 ) == 0 )
	{
		cost += 1e10;
	}

	int64_t d = val - a1 * a2;

	if( d < 0 )
	{
		d = -d;
	}

	if( DoPrint )
	{
		printf( "%lli\n", a1 );
		printf( "%lli\n", a2 );
	}

	return( d + cost );
}

int main()
{
	double x[ N ];
	double minf;
	biteopt_minimize( N, fn, NULL, lb, ub, x, &minf, 2000000, 5000, 10 );

	printf( "minf = %.10g\n", minf );

	DoPrint = true;
	fn( N, x, NULL );
}
