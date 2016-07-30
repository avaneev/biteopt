#include <stdio.h>
#include "biteopt.h"
#include "biteopt2.h"
#include "bitefan.h"
//#include "CNMSimplexSolver.h"

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

#if !defined( M_PI )
	#define M_PI 3.14159265358979324
#endif // !defined( M_PI )

const int ParamCount = 2;
const int FnCount = 8;
CBEORnd rnd;

/**
 * Optimization test class.
 */

class CTestOpt : public CBEOOptimizerFan< ParamCount >
{
public:
	int fn;

	virtual void getMinValues( double* const p ) const
	{
		if( fn == 9 )
		{
			p[ 0 ] = -15;
			p[ 1 ] = -3;
		}
		else
		{
			p[ 0 ] = -10 + rnd.getRndValue() * 5.5;
			p[ 1 ] = -10 + rnd.getRndValue() * 5.5;
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		if( fn == 9 )
		{
			p[ 0 ] = -5;
			p[ 1 ] = 3;
		}
		else
		{
			p[ 0 ] = 10 - rnd.getRndValue() * 5.5;
			p[ 1 ] = 10 - rnd.getRndValue() * 5.5;
		}
	}

	virtual double optcost( const double* const p ) const
	{
		const double x = p[ 0 ];
		const double y = p[ 1 ];

		if( fn == 1 ) return( sqr( x + 2 * y - 7 ) + sqr( 2 * x + y - 5 ));
		if( fn == 2 ) return( 0.26 * ( x * x + y * y ) - 0.48 * x * y );
		if( fn == 3 ) return( x * x + y * y );
		if( fn == 4 ) return( sqr( sin( 3 * M_PI * x )) + sqr( x - 1 ) * ( 1 + sqr( sin( 3 * M_PI * y ))) + sqr( y - 1 ) * ( 1 + sqr( sin( 2 * M_PI * y ))));
		if( fn == 5 ) return( 0.5 + ( sqr( sin( x * x - y * y )) - 0.5 ) / sqr( 1 + 0.001 * ( x * x + y * y )));
		if( fn == 6 ) return( -20 * exp( -0.2 * sqrt( 0.5 * ( x * x + y * y ))) - exp( 0.5 * ( cos( 2 * M_PI * x ) + cos( 2 * M_PI * y ))) + 2.71828182845904524 + 20 );
		if( fn == 7 ) return( 100 * sqr( y - x * x ) + sqr( x - 1 ));
		if( fn == 8 ) return( sqr( 1.5 - x + x * y ) + sqr( 2.25 - x + x * y * y ) + sqr( 2.625 - x + x * y * y * y ));
		if( fn == 9 ) return( 100 * sqrt( fabs( y - 0.01 * x * x )) + 0.01 * fabs( x + 10 ));
		return( 2 * x * x - 1.05 * sqr( sqr( x )) + sqr( sqr( sqr( x ))) / 6 + x * y + y * y );
	}
};

	inline __declspec( naked ) uint64_t rdtsc()
	{
		__asm
		{
			rdtsc
			ret
		}
	}

int main()
{
	CTestOpt opt;
	rnd.init( 0 );

	double ItAvg = 0.0;
	double ItRtAvg = 0.0;
	double RjAvg = 0.0;
	int k;

	uint64_t tc = 0;

	for( k = 0; k < FnCount; k++ )
	{
		opt.fn = k;
		const int IterCount = 10000;
		int Iters[ IterCount ];
		int Rejects[ IterCount ];
		int AvgIter = 0;
		int AvgRej = 0;
		double AvgCost = 0.0;
		double AvgP1 = 0.0;
		double AvgP2 = 0.0;
		int Rej = 0;
		int j;

		for( j = 0; j < IterCount; j++ )
		{
			double Params[ ParamCount ];
			double minv[ ParamCount ];
			double maxv[ ParamCount ];
			opt.getMinValues( minv );
			opt.getMaxValues( maxv );
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = minv[ i ] + ( maxv[ i ] - minv[ i ]) *
					rnd.getRndValue();
			}

/*			i = 10000;
			AvgCost += vox :: solveNMSimplex( opt, ParamCount, Params, true,
				0.000001, &i );

			AvgP1 += Params[ 0 ];
			AvgP2 += Params[ 1 ];
*/
			opt.init( rnd, Params );

			opt.optimize( rnd );

			double PrevBestCost = opt.getBestCost();
			int PrevBestCostCount = 0;
			const int InnerIterCount = 10000;

			for( i = 0; i < InnerIterCount; i++ )
			{
				if( opt.getBestCost() < 0.001 )
				{
					break;
				}

				if( PrevBestCost > opt.getBestCost() )
				{
					PrevBestCost = opt.getBestCost();
					PrevBestCostCount = 0;
				}
				else
				{
					PrevBestCostCount++;

					if( PrevBestCostCount == 10000 )
					{
						i = InnerIterCount;
						break;
					}
				}

				const uint64_t t1 = rdtsc();
				opt.optimize( rnd );
				tc += rdtsc() - t1;
			}

			if( i == InnerIterCount )
			{
				Rej++;
			}

			AvgCost += opt.getBestCost();
			AvgP1 += opt.getBestParams()[ 0 ];
			AvgP2 += opt.getBestParams()[ 1 ];

			Iters[ j ] = i;
			AvgIter += i;
			Rejects[ j ] = Rej;
			AvgRej += Rej;
		}

		AvgCost /= IterCount;
		AvgP1 /= IterCount;
		AvgP2 /= IterCount;

		double Avg = (double) AvgIter / IterCount;
		double RMS = 0.0;
		double Avg2 = (double) AvgRej / IterCount;
		double RMS2 = 0.0;

		for( j = 0; j < IterCount; j++ )
		{
			double v = Iters[ j ] - Avg;
			RMS += v * v;
			v = Rejects[ j ] - Avg2;
			RMS2 += v * v;
		}

		RMS = sqrt( RMS / IterCount );
		RMS2 = sqrt( RMS2 / IterCount );

		printf( "AvgIt:%6.1f RMSIt:%6.1f AvgRj:%5.1f "
			"RMSRj:%5.1f Cost:%10.8f%6.3f%6.3f\n", Avg,
			RMS, Avg2, RMS2, AvgCost, AvgP1, AvgP2 );

		ItAvg += Avg;
		RjAvg += Avg2;
		ItRtAvg += RMS / Avg;
	}

	ItAvg /= FnCount;
	ItRtAvg /= FnCount;
	RjAvg /= FnCount;
	printf( "ItAvg: %.1f\n", ItAvg );
	printf( "ItRtAvg: %.3f\n", ItRtAvg );
	printf( "RjAvg: %.1f\n", RjAvg );
	printf( "%llu\n", tc );

	return( 0 );
}
