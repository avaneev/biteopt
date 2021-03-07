// Non-linearly constrained problem:
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2031.htm

#include <stdio.h>
#include "biteopt.h"

const int N = 4;

class CTestOpt : public CBiteOptDeep
{
public:
	int con_notmet;

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

	void penalty( double& pn, const double v )
	{
		if( v > 1e-8 )
		{
			pn = pn * 1.30 + ( v + v * v + v * v * v ) * 1e4;
			con_notmet++;
		}
	}

	void penalty0( double& pn, double v )
	{
		v = fabs( v );

		if( v > 1e-8 )
		{
			pn = pn * 1.25 + ( v + v * v ) * 1e4;
			con_notmet++;
		}
	}

	virtual double optcost( const double* const p )
	{
		double cost = 3.0*p[1]+1e-6*pow(p[0],3.0)+2.0*p[1]+
			2e-6/3.0*pow(p[1],3.0);

		con_notmet = 0;
		double pn = 0.0;

		penalty( pn, p[2]-p[3]-0.55 );
		penalty( pn, p[3]-p[2]-0.55 );
		penalty0( pn, 1000*(sin(-p[2]-0.25)+sin(-p[3]-0.25))+894.8-p[0] );
		penalty0( pn, 1000*(sin(p[2]-0.25)+sin(p[2]-p[3]-0.25))+894.8-p[1] );
		penalty0( pn, 1000*(sin(p[3]-0.25)+sin(p[3]-p[2]-0.25))+1294.8 );

		if( con_notmet > 0 )
		{
			cost += pn;
		}

		return( cost );
	}
};

int main()
{
	CBiteRnd rnd;
	rnd.init( 1 ); // Needs to be seeded with different values on each run.

	CTestOpt opt;
	opt.updateDims( N, 6 );
	opt.init( rnd );

	int i;

	for( i = 0; i < 200000; i++ )
	{
		opt.optimize( rnd );
	}

	const double minf = opt.optcost( opt.getBestParams() );

	printf( "BestCost: %f\n", minf );
	printf( "Constraints not met: %i\n", opt.con_notmet );

	for( i = 0; i < N; i++ )
	{
		printf( "x[%i] = %.10g\n", i, opt.getBestParams()[ i ]);
	}

	return( 0 );
}
