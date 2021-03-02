// Non-linearly constrained problem:
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2031.htm

#include <stdio.h>
#include "biteopt.h"

const int N = 4;

class CTestOpt : public CBiteOptDeep
{
public:
	CTestOpt()
	{
		updateDims( ::N, 6 );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.0;
		p[ 1 ] = 0.0;
		p[ 2 ] = -0.55;
		p[ 3 ] = -0.55;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 1200.0;
		p[ 1 ] = 1200.0;
		p[ 2 ] = 0.55;
		p[ 3 ] = 0.55;
	}

	static double penalty( const double v )
	{
		return( v <= 0.0 ? 0.0 : 100000 + v * 9999 );
	}

	static double penalty0( const double v )
	{
		return( fabs( v ) <= 0.0001 ? 0.0 : 100000 + fabs( v ) * 9999 );
	}

	virtual double optcost( const double* const p )
	{
		double cost = 3.0*p[1]+1e-6*pow(p[0],3.0)+2.0*p[1]+
			2e-6/3.0*pow(p[1],3.0);

		cost += penalty( p[2]-p[3]-0.55 );
		cost += penalty( p[3]-p[2]-0.55 );
		cost += penalty0( 1000*(sin(-p[2]-0.25)+sin(-p[3]-0.25))+894.8-p[0] );
		cost += penalty0( 1000*(sin(p[2]-0.25)+sin(p[2]-p[3]-0.25))+894.8-p[1] );
		cost += penalty0( 1000*(sin(p[3]-0.25)+sin(p[3]-p[2]-0.25))+1294.8 );

		return( cost );
	}
};

int main()
{
	CBiteRnd rnd;
	rnd.init( 1 ); // Needs to be seeded with different values on each run.

	CTestOpt opt;
	opt.init( rnd );

	int i;

	for( i = 0; i < 200000; i++ )
	{
		opt.optimize( rnd );
	}

	printf( "BestCost: %f\n", opt.getBestCost() );

	for( i = 0; i < N; i++ )
	{
		printf( "x[%i] = %.10g\n", i, opt.getBestParams()[ i ]);
	}

	return( 0 );
}
