#include <stdio.h>
#include "testopt.h"

const int FnCount = 9;
const double CostThreshold = 0.001;
const int IterCount = 10000;
const int InnerIterCount = 10000;
const int FanParamCount = 7;
const int FanIterCount = 2000;

static double roundp( const double x )
{
	return( floor( x * 1000000.0 + 0.5 ) / 1000000.0 );
}

class CFanOpt : public CBEOOptimizerFan< FanParamCount, 3 >
{
public:
	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.7;
		p[ 1 ] = 0.1;
		p[ 2 ] = 0.3;
		p[ 3 ] = 0.3;
		p[ 4 ] = 0.3;
		p[ 5 ] = 0.3;
		p[ 6 ] = 0.0;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 3.5;
		p[ 1 ] = 1.0;
		p[ 2 ] = 1.0;
		p[ 3 ] = 2.5;
		p[ 4 ] = 3.5;
		p[ 5 ] = 3.5;
		p[ 6 ] = 1.0;
	}

	virtual double optcost( const double* const p ) const
	{
		rnd.init( 0 );

		CTestOpt opt;
		opt.CostMult = roundp( p[ 0 ]);
		opt.BestMult = roundp( p[ 1 ]);
		opt.HistMult = roundp( p[ 2 ]);
		opt.HistMult2 = roundp( p[ 3 ]);
		opt.PrevMult = roundp( p[ 4 ]);
		opt.CentMult = roundp( p[ 5 ]);
		opt.CentOffs = roundp( p[ 6 ]);

		double ItAvg = 0.0;
		double RMSAvg = 0.0;
		double ItRtAvg = 0.0;
		double RjAvg = 0.0;
		int k;

		for( k = 0; k < FnCount; k++ )
		{
			opt.fn = k;
			int Iters[ IterCount ];
			double AvgIter = 0.0;
			double IterAvgCost = 0.0;
			double AvgP1 = 0.0;
			double AvgP2 = 0.0;
			int Rej = 0;
			int j;

			for( j = 0; j < IterCount; j++ )
			{
				int i;

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

					opt.optimize( rnd );
				}

				if( i == InnerIterCount )
				{
					Rej++;
				}

				IterAvgCost += opt.getBestCost();
				AvgP1 += opt.getBestParams()[ 0 ];
				AvgP2 += opt.getBestParams()[ 1 ];

				Iters[ j ] = i;
				AvgIter += i;
			}

			IterAvgCost /= IterCount;
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

			ItAvg += Avg;
			RMSAvg += RMS;
			RjAvg += Rej;
			ItRtAvg += RMS / Avg;
		}

		ItAvg /= FnCount;
		RMSAvg /= FnCount;
		ItRtAvg /= FnCount;
		RjAvg /= FnCount;

		return( ItAvg * pow( ItRtAvg, 0.5 ));
	}
};

int main()
{
	CBEORnd rnd2;
	rnd2.init( 1 );

	double Params[ FanParamCount ];
	Params[ 0 ] = 2.040037;
	Params[ 1 ] = 0.712701;
	Params[ 2 ] = 0.710949;
	Params[ 3 ] = 1.544876;
	Params[ 4 ] = 1.374690;
	Params[ 5 ] = 1.278461;
	Params[ 6 ] = 0.755409;

	CFanOpt opt;
	opt.init( rnd2, Params );
	int i;

	for( i = 0; i < FanIterCount; i++ )
	{
		opt.optimize( rnd2 );

		printf( "%f\n", opt.getBestCost() );
		int j;

		for( j = 0; j < FanParamCount; j++ )
		{
			Params[ j ] = roundp( opt.getBestParams()[ j ]);
		}

		printf( "CostMult = %.6f;\n", Params[ 0 ]);
		printf( "BestMult = %.6f;\n", Params[ 1 ]);
		printf( "HistMult = %.6f;\n", Params[ 2 ]);
		printf( "HistMult2 = %.6f;\n", Params[ 3 ]);
		printf( "PrevMult = %.6f;\n", Params[ 4 ]);
		printf( "CentMult = %.6f;\n", Params[ 5 ]);
		printf( "CentOffs = %.6f;\n", Params[ 6 ]);
	}

	return( 0 );
}
