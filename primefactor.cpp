// Prime-number factoring example. Uses BiteOptDeep with a huge (5000) depth
// factor.

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "biteopt.h"

const int N = 16;
double lb[ N ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
double ub[ N ] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
int64_t val = 20000003LL * 10000079LL;
bool DoPrint = false;

double fn( int ND, const double* x, void* data )
{
	int64_t a1 = 0;
	int64_t a2 = 0;
	double cost = 0.0;
	int i;

	for( i = 0; i < 8; i++ )
	{
		a1 += (int) ( x[ i ] * 25000000.0 + 0.5 );
		a2 += (int) ( x[ 8 + i ] * 25000000.0 + 0.5 );
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
	biteopt_minimize( N, fn, NULL, lb, ub, x, &minf, 2000000, 5000, 3 );

	printf( "minf = %.10g\n", minf );

	DoPrint = true;
	fn( N, x, NULL );
}
