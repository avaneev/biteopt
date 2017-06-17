#include <stdio.h>
#include "tester.h"

const int FanParamCount = 6;
const int FanIterCount = 4000;

static double roundp( const double x )
{
	return( floor( x * 100000000.0 + 0.5 ) / 100000000.0 );
}

class CFanOpt : public CBEOOptimizerFan< FanParamCount, 3, 16 >
{
public:
	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.7;
		p[ 1 ] = 0.1;
		p[ 2 ] = 0.0;
		p[ 3 ] = 0.3;
		p[ 4 ] = 0.3;
		p[ 5 ] = 0.0;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 5.0;
		p[ 1 ] = 1.0;
		p[ 2 ] = 1.0;
		p[ 3 ] = 1.0;
		p[ 4 ] = 2.0;
		p[ 5 ] = 1.5;
	}

	virtual double optcost( const double* const p ) const
	{
		rnd.init( 0 );

//		CTester< 2 > Tester;
		CTester< 10 > Tester;
		Tester.opt -> CostMult = roundp( p[ 0 ]);
		Tester.opt -> BestMult = roundp( p[ 1 ]);
		Tester.opt -> HistMult = roundp( p[ 2 ]);
		Tester.opt -> PrevMult = roundp( p[ 3 ]);
		Tester.opt -> CentMult = roundp( p[ 4 ]);
		Tester.opt -> CentOffs = roundp( p[ 5 ]);
//		Tester.init( OptCorpus2D, 0.001, 5000, 10000, true, false );
		Tester.init( OptCorpusND, 0.01, 120, 150000, false, false );

		Tester.run();

		return( Tester.ItAvg * pow( Tester.ItRtAvg, 0.5 ) *
			( 1.0 + Tester.RjAvg * 50.0 ));
	}
};

int main()
{
	CBEORnd rnd2;
	rnd2.init( 1 );

	double Params[ FanParamCount ];
	Params[ 0 ] = 1.36232993;
	Params[ 1 ] = 0.63217411;
	Params[ 2 ] = 0.55474765;
	Params[ 3 ] = 0.43400424;
	Params[ 4 ] = 1.31470436;
	Params[ 5 ] = 0.61169491;

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

		printf( "CostMult = %.8f;\n", Params[ 0 ]);
		printf( "BestMult = %.8f;\n", Params[ 1 ]);
		printf( "HistMult = %.8f;\n", Params[ 2 ]);
		printf( "PrevMult = %.8f;\n", Params[ 3 ]);
		printf( "CentMult = %.8f;\n", Params[ 4 ]);
		printf( "CentOffs = %.8f;\n", Params[ 5 ]);
	}

	return( 0 );
}
