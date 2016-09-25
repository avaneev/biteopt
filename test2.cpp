#include <stdio.h>
#include "testopt.h"
//#include "CNMSimplexSolver.h"

const int FnCount = 13;
const double CostThreshold = 0.001;
const int IterCount = 10000;
const int InnerIterCount = 10000;

int main()
{
	rnd.init( 0 );

	CTestOpt opt;

	double ItAvg = 0.0;
	double ItRtAvg = 0.0;
	double RjAvg = 0.0;
	int k;

	uint64_t tc = 0;

	for( k = 0; k < FnCount; k++ )
	{
		opt.fn = k;
		int Iters[ IterCount ];
		double AvgIter = 0;
		double AvgCost = 0.0;
		double AvgP1 = 0.0;
		double AvgP2 = 0.0;
		int Rej = 0;
		int j;

		for( j = 0; j < IterCount; j++ )
		{
			int i;
/*			double Params[ ParamCount ];
			double minv[ ParamCount ];
			double maxv[ ParamCount ];
			opt.getMinValues( minv );
			opt.getMaxValues( maxv );

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = minv[ i ] + ( maxv[ i ] - minv[ i ]) *
					rnd.getRndValue();
			}

			i = 10000;
			AvgCost += vox :: solveNMSimplex( opt, ParamCount, Params, true,
				0.000001, &i );

			AvgP1 += Params[ 0 ];
			AvgP2 += Params[ 1 ];
*/
			opt.sign1 = ( rnd.getRndValue() < 0.5 ? 1.0 : -1.0 );
			opt.sign2 = ( rnd.getRndValue() < 0.5 ? 1.0 : -1.0 );
			opt.init( rnd );
			opt.optimize( rnd );

			double PrevBestCost = opt.getBestCost();
			int PrevBestCostCount = 0;

			for( i = 0; i < InnerIterCount; i++ )
			{
				if( opt.getBestCost() < CostThreshold )
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

				const uint64_t t1 = __rdtsc();
				opt.optimize( rnd );
				tc += __rdtsc() - t1;
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
		}

		AvgCost /= IterCount;
		AvgP1 /= IterCount;
		AvgP2 /= IterCount;

		const double Avg = AvgIter / IterCount;
		double RMS = 0.0;

		for( j = 0; j < IterCount; j++ )
		{
			const double v = Iters[ j ] - Avg;
			RMS += v * v;
		}

		RMS = sqrt( RMS / IterCount );

		printf( "AvgIt:%6.1f RMSIt:%6.1f Rj:%5i Cost:%10.8f%6.3f%6.3f\n",
			Avg, RMS, Rej, AvgCost, AvgP1, AvgP2 );

		ItAvg += Avg;
		RjAvg += Rej;
		ItRtAvg += RMS / Avg;
	}

	ItAvg /= FnCount;
	ItRtAvg /= FnCount;
	RjAvg /= FnCount;
	printf( "ItAvg: %.1f (avg convergence time)\n", ItAvg );
	printf( "ItRtAvg: %.3f (avg ratio of std.deviation and average)\n", ItRtAvg );
	printf( "RjAvg: %.1f (avg number of rejects, out of %i)\n", RjAvg,
		InnerIterCount );

	printf( "%llu\n", tc );

	return( 0 );
}
