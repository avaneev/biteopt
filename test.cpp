#include <stdio.h>
#include "biteopt.h"

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

/**
 * Optimization test class.
 */

class CTestOpt : public CBEOOptimizer< 2, 4 >
{
public:
	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = -10;
		p[ 1 ] = -10;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 10;
		p[ 1 ] = 10;
	}

	virtual double optcost( const double* const p ) const
	{
		const double x = p[ 0 ];
		const double y = p[ 1 ];

		return( sqr( x + 2 * y - 7 ) + sqr( 2 * x + y - 5 ));
//		return( 0.26 * ( x * x + y * y ) - 0.48 * x * y );
//		return( x * x + y * y );
//		return( 2 * x * x - 1.05 * sqr( sqr( x )) + sqr( sqr( sqr( x ))) / 6 + x * y + y * y );
	}
};

int main()
{
	CBEORnd rnd;
	rnd.init( 0 ); // Needs to be seeded with different values on each run.

	CTestOpt opt;
	opt.init( rnd );

	int i;

	for( i = 0; i < 10000; i++ )
	{
		opt.optimize( rnd );

		if( opt.getBestCost() < 0.000001 )
		{
			break;
		}
	}

	printf( "IterCount: %i\n", i );
	printf( "BestCost: %f\n", opt.getBestCost() );
	printf( "%f ", opt.getBestParams()[ 0 ]);
	printf( "%f\n", opt.getBestParams()[ 1 ]);

	return( 0 );
}
