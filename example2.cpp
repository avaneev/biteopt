// Minimization test for Hougen-Watson model for reaction kinetics.
// Non-linear least squares problem.

#include <stdio.h>
#include "biteopt.h"

const int N = 5;
double lb[ N ] = { 0.01, 0.01, 0.01, 0.01, 0.01 };
double ub[ N ] = { 2.0, 2.0, 2.0, 2.0, 2.0 };

double fn( int ND, const double* x, void* data )
{
	const double x1[ 13 ] =
		{ 470, 285, 470, 470, 470, 100, 100, 470, 100, 100, 100, 285, 285 };

	const double x2[ 13 ] =
		{ 300, 80, 300, 80, 80, 190, 80, 190, 300, 300, 80, 300, 190 };

	const double x3[ 13 ] =
		{ 10, 10, 120, 120, 10, 10, 65, 65, 54, 120, 120, 10, 120 };

	const double rate[ 13 ] =
		{ 8.55,3.79,4.82,0.02,2.75,14.39,2.54,4.35,13,8.5,0.05,11.32,3.13 };

	double s = 0.0;
	int i;

	for( i = 0; i < 13; i++ )
	{
		double v = (x[0]*x2[i]-x3[i]/x[4])/
			(1.0+x[1]*x1[i]+x[2]*x2[i]+x[3]*x3[i]);

		double d = rate[i]-v;
		s += d * d;
	}

	return( s );
}

int main()
{
	double x[ N ];
	double minf;
	int c = biteopt_minimize( N, fn, NULL, lb, ub, x, &minf, 20000, 1, 1, 1 );

	printf( "evals_done = %i\n", c );
	printf( "minf = %.10g\n", minf );

	int i;

	for( i = 0; i < N; i++ )
	{
		printf( "x[%i] = %.8g\n", i, x[ i ]);
	}
}
