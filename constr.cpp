// Constrained problem:
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page506.htm

#include <stdio.h>
#include "biteopt.h"

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

static const double tol = 1e-15;
const int N = 13;

class CTestOpt : public CBiteOpt
{
public:
	int con_notmet;
	double real_value;

	CTestOpt()
	{
		updateDims( N );
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

	double penalty( double v )
	{
		if( v > tol )
		{
			con_notmet++;
			return( v - tol );
		}

		return( 0.0 );
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

		const int n_con = 9;
		double pn[ n_con ];
		con_notmet = 0;

		pn[ 0 ] = penalty( 2.0*p[0]+2.0*p[1]+p[9]+p[10]-10.0 );
		pn[ 1 ] = penalty( 2.0*p[0]+2.0*p[2]+p[9]+p[11]-10.0 );
		pn[ 2 ] = penalty( 2.0*p[1]+2.0*p[2]+p[10]+p[11]-10.0 );
		pn[ 3 ] = penalty( -8.0*p[0]+p[9] );
		pn[ 4 ] = penalty( -8.0*p[1]+p[10] );
		pn[ 5 ] = penalty( -8.0*p[2]+p[11] );
		pn[ 6 ] = penalty( -2.0*p[3]-p[4]+p[9] );
		pn[ 7 ] = penalty( -2.0*p[5]-p[6]+p[10] );
		pn[ 8 ] = penalty( -2.0*p[7]-p[8]+p[11] );

		real_value = cost;

		if( con_notmet > 0 )
		{
			const double ps = 1.0 + 1.0 / n_con;
			double pns = 0.0;

			for( i = 0; i < n_con; i++ )
			{
				const double v = pn[ i ];
				pns = ( 1.0 + pns ) * ps + ( v + v * v + v * v * v ) * 0.33333;
			}

			cost += 1e10 * pns;
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

	for( i = 0; i < 500000; i++ )
	{
		if( opt.optimize( rnd ) > N * 128 )
		{
			break;
		}
	}

	opt.optcost( opt.getBestParams() );

	printf( "Finished at iteration %i\n", i + 1 );
	printf( "Objective = %.8g\n", opt.real_value );
	printf( "Constraints not met: %i\n", opt.con_notmet );

	for( i = 0; i < N; i++ )
	{
		printf( "x[%i] = %f\n", i, opt.getBestParams()[ i ]);
	}

	return( 0 );
}
