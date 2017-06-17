// Example of the optimizePlateau() function use.

#include <stdio.h>
#include "bitefan.h"

class CTestOpt : public CBEOOptimizerFan
{
public:
	CTestOpt()
	{
		updateDims( 2 );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = -15;
		p[ 1 ] = -3;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = -5;
		p[ 1 ] = 3;
	}

	virtual double optcost( const double* const p ) const
	{
		const double x = p[ 0 ];
		const double y = p[ 1 ];

		return( 100 * sqrt( fabs( y - 0.01 * x * x )) + 0.01 * fabs( x + 10 ));
	}
};

int main()
{
	CBEORnd rnd;
	rnd.init( 1 ); // Needs to be seeded with different values on each run.

	CTestOpt opt;
	opt.optimizePlateau( rnd, 0.001, 100, 15000 );

	printf( "BestCost: %f\n", opt.getBestCost() );
	printf( "%f ", opt.getBestParams()[ 0 ]);
	printf( "%f\n", opt.getBestParams()[ 1 ]);

	return( 0 );
}
