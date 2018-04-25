// Constrained problem:
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page506.htm

#include <stdio.h>
#include "biteopt.h"

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

class CTestOpt : public CBiteOpt
{
public:
	CTestOpt()
	{
		updateDims( 13 );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.0;
		p[ 1 ] = 0.0;
		p[ 2 ] = 0.0;
		p[ 3 ] = 0.0;
		p[ 4 ] = 0.0;
		p[ 5 ] = 0.0;
		p[ 6 ] = 0.0;
		p[ 7 ] = 0.0;
		p[ 8 ] = 0.0;
		p[ 9 ] = 0.0;
		p[ 10 ] = 0.0;
		p[ 11 ] = 0.0;
		p[ 12 ] = 0.0;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 1.0;
		p[ 1 ] = 1.0;
		p[ 2 ] = 1.0;
		p[ 3 ] = 1.0;
		p[ 4 ] = 1.0;
		p[ 5 ] = 1.0;
		p[ 6 ] = 1.0;
		p[ 7 ] = 1.0;
		p[ 8 ] = 1.0;
		p[ 9 ] = 100.0;
		p[ 10 ] = 100.0;
		p[ 11 ] = 100.0;
		p[ 12 ] = 1.0;
	}

	static double penalty( const double v )
	{
		return( v <= 0.0 ? 0.0 : 100000 + v * 9999 );
	}

	virtual double optcost( const double* const p )
	{
		double s1 = 0.0;
		double s2 = 0.0;
		double s3 = 0.0;
		int i;

		for( i = 0; i < 4; i++ )
		{
			s1 += p[ i ];
			s2 += sqr( p[ i ]);
		}

		for( i = 4; i < 13; i++ )
		{
			s3 += p[ i ];
		}

		double cost = 5.0*s1-5.0*s2-s3;

		cost += penalty( 2.0*p[0]+2.0*p[1]+p[9]+p[10]-10.0 );
		cost += penalty( 2.0*p[0]+2.0*p[2]+p[9]+p[11]-10.0 );
		cost += penalty( 2.0*p[1]+2.0*p[2]+p[10]+p[11]-10.0 );
		cost += penalty( -8.0*p[0]+p[9] );
		cost += penalty( -8.0*p[1]+p[10] );
		cost += penalty( -8.0*p[2]+p[11] );
		cost += penalty( -2.0*p[3]-p[4]+p[9] );
		cost += penalty( -2.0*p[5]-p[6]+p[10] );
		cost += penalty( -2.0*p[7]-p[8]+p[11] );

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

	for( i = 0; i < 15000; i++ )
	{
		opt.optimize( rnd );

		if( opt.getBestCost() < -14.999999 )
		{
			break;
		}
	}

	printf( "IterCount: %i\n", i );
	printf( "BestCost: %f\n", opt.getBestCost() );

	for( i = 0; i < 13; i++ )
	{
		printf( "x[%i] = %f\n", i, opt.getBestParams()[ i ]);
	}

	return( 0 );
}
