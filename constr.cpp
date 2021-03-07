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
	int con_notmet;

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

	void penalty( double& pn, const double v )
	{
		if( v > 1e-8 )
		{
			pn = pn * 1.30 + ( v + v * v + v * v * v ) * 1e4;
			con_notmet++;
		}
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

		con_notmet = 0;
		double pn = 0.0;

		penalty( pn, 2.0*p[0]+2.0*p[1]+p[9]+p[10]-10.0 );
		penalty( pn, 2.0*p[0]+2.0*p[2]+p[9]+p[11]-10.0 );
		penalty( pn, 2.0*p[1]+2.0*p[2]+p[10]+p[11]-10.0 );
		penalty( pn, -8.0*p[0]+p[9] );
		penalty( pn, -8.0*p[1]+p[10] );
		penalty( pn, -8.0*p[2]+p[11] );
		penalty( pn, -2.0*p[3]-p[4]+p[9] );
		penalty( pn, -2.0*p[5]-p[6]+p[10] );
		penalty( pn, -2.0*p[7]-p[8]+p[11] );

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
	opt.init( rnd );

	int i;

	for( i = 0; i < 30000; i++ )
	{
		opt.optimize( rnd );

		if( opt.getBestCost() < -14.999999 )
		{
			break;
		}
	}

	const double minf = opt.optcost( opt.getBestParams() );

	printf( "IterCount: %i\n", i );
	printf( "BestCost: %f\n", minf );
	printf( "Constraints not met: %i\n", opt.con_notmet );

	for( i = 0; i < 13; i++ )
	{
		printf( "x[%i] = %f\n", i, opt.getBestParams()[ i ]);
	}

	return( 0 );
}
