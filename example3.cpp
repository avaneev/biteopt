// ZDT4 2-objective minimization example.

#include <stdio.h>
#include "biteopt.h"
#include "bitehv2.h"

const int D = 10; // ZDT4 parameters.

class CTestOpt : public CBiteOpt
{
public:
	RankSumFrontHV hv_solver;

	CTestOpt()
		: hv_solver( 1.0, 300.0 ) // Reference point
	{
		updateDims( D );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.0;

		int i;

		for( i = 1; i < D; i++ )
		{
			p[ i ] = -5.0;
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 1.0;

		int i;

		for( i = 1; i < D; i++ )
		{
			p[ i ] = 5.0;
		}
	}

	virtual double optcost( const double* const p )
	{
		// ZDT4:
		// g = 1 + 10*(n-1) + sum_{i=2..n}(x_i^2 - 10*cos(4*pi*x_i))
		// min over [-5,5]^29 is at x_i=0 giving g=1, hence f2=1-sqrt(f1).
		double g=1.0 + 10.0*(D-1);
		for(int j=1;j<D;j++)
			g += p[j]*p[j] - 10.0*cos(4.0*3.14159265358979324*p[j]);

		double f1 = p[0];
		double f2 = g*(1.0-sqrt(f1/g));

		if( f1 < 1.0 && f2 < 300.0 )
		{
			hv_solver.addPoint( f1, f2, p, D );
		}

		return( hv_solver.hypervolume() * ( -1.0 + f1 * f2 ));
	}
};

int main()
{
	CBiteRnd rnd;
	rnd.init( 1 ); // Needs to be seeded with different values on each run.

	CTestOpt opt;
	opt.init( rnd );

	const double ohv = 899.0 / 3.0; // Optimal ZDT4 hypervolume at ref.point.
	double hvn; // Normalized hypervolume.
	int i;

	for( i = 0; i < 1000000; i++ )
	{
		opt.optimize( rnd );

		hvn = opt.hv_solver.hypervolume() / ohv;

		if( hvn > 0.99995 )
		{
			break;
		}
	}

	printf( "IterCount: %i\n", i );
	printf( "hypervolume: %f\n", opt.hv_solver.hypervolume() );

	printf( "normalized hv: %.4f\n", hvn );
	printf( "Front archive size: %i\n", opt.hv_solver.getFront().size() );

	return( 0 );
}
