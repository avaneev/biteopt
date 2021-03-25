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

	double penalty( const double v )
	{
		if( v > 1e-8 )
		{
			con_notmet++;
			return( v );
		}

		return( 0.0 );
	}

	double penalty0( double v )
	{
		v = fabs( v );

		if( v > 1e-4 )
		{
			con_notmet++;
			return( v );
		}

		return( 0.0 );
	}

	virtual double optcost( const double* const p )
	{
		double cost = 3.0*p[1]+1e-6*pow(p[0],3.0)+2.0*p[1]+
			2e-6/3.0*pow(p[1],3.0);

		const int n_con = 5;
		double pn[ n_con ];
		con_notmet = 0;

		pn[ 0 ] = penalty( p[2]-p[3]-0.55 );
		pn[ 1 ] = penalty( p[3]-p[2]-0.55 );
		pn[ 2 ] = penalty0( 1000*(sin(-p[2]-0.25)+sin(-p[3]-0.25))+894.8-p[0] );
		pn[ 3 ] = penalty0( 1000*(sin(p[2]-0.25)+sin(p[2]-p[3]-0.25))+894.8-p[1] );
		pn[ 4 ] = penalty0( 1000*(sin(p[3]-0.25)+sin(p[3]-p[2]-0.25))+1294.8 );

		double pns = 0.0;
		int i;

		for( i = 0; i < n_con; i++ )
		{
			pns = pns * pow( 4.0, 1.0 / n_con ) + pn[ i ] * pn[ i ];
		}

		cost += 1e8 * ( con_notmet * con_notmet + pns );

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

	for( i = 0; i < 500000; i++ )
	{
		if( opt.optimize( rnd ) > N * 64 )
		{
			printf( "Finished at iter %i\n", i + 1 );
			break;
		}
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
