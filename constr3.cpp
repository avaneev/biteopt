// Constrained problem:
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2424.htm

#include <stdio.h>
#include <string.h>
#include "biteopt.h"

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

const int N = 10;

class CTestOpt : public CBiteOptDeep
{
public:
	virtual void getMinValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = -10.0;
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = 10.0;
		}
	}

	static double penalty( const double v )
	{
		return( v <= 0.0001 ? 0.0 : 100000 + v * 9999 );
	}

	virtual double optcost( const double* const p )
	{
		double cost = sqr(p[0])+sqr(p[1])+p[0]*p[1]-14*p[0]-16*p[1]+sqr(p[2]-10)+
			4*sqr(p[3]-5)+sqr(p[4]-3)+2*sqr(p[5]-1)+5*sqr(p[6])+
			7*sqr(p[7]-11)+2*sqr(p[8]-10)+sqr(p[9]-7)+45;

		cost += penalty( 4*p[0]+5*p[1]-3*p[6]+9*p[7]-105 );
		cost += penalty( 10*p[0]-8*p[1]-17*p[6]+2*p[7] );
		cost += penalty( -8*p[0]+2*p[1]+5*p[8]-2*p[9]-12 );
		cost += penalty( 3*sqr(p[0]-2)+4*sqr(p[1]-3)+2*sqr(p[2])-7*p[3]-120 );
		cost += penalty( 5*sqr(p[0])+8*p[1]+sqr(p[2]-6)-2*p[3]-40 );
		cost += penalty( 0.5*sqr(p[0]-8)+2*sqr(p[1]-4)+3*sqr(p[4])-p[5]-30 );
		cost += penalty( sqr(p[0])+2*sqr(p[1]-2)-2*p[0]*p[1]+14*p[4]-6*p[5] );
		cost += penalty( -3*p[0]+6*p[1]+12*sqr(p[8]-8)-7*p[9] );

		return( cost );
	}
};

int main()
{
	CBiteRnd rnd;
	rnd.init( 1 ); // Needs to be seeded with different values on each run.

	int i;

	CTestOpt opt;
	opt.updateDims( N, 15 );
	opt.init( rnd );

	for( i = 0; i < 250000; i++ )
	{
		opt.optimize( rnd );
	}

	printf( "BestCost: %f\n", opt.getBestCost() );

	for( i = 0; i < N; i++ )
	{
		printf( "x[%i] = %f\n", i, opt.getBestParams()[ i ]);
	}

	// Optimum provided by function's source.

	double p[ 10 ];
	p[ 0 ] = 2.171996;
	p[ 1 ] = 2.363683;
	p[ 2 ] = 8.773926;
	p[ 3 ] = 5.095984;
	p[ 4 ] = 0.9906548;
	p[ 5 ] = 1.430574;
	p[ 6 ] = 1.321644;
	p[ 7 ] = 9.828726;
	p[ 8 ] = 8.280092;
	p[ 9 ] = 8.375927;
	printf( "Source opt: %f\n",opt.optcost(p));

	return( 0 );
}
