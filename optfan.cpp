#include <stdio.h>
#include "tester.h"

const int FanParamCount = 6;
const int FanIterCount = 4000;

static double roundp( const double x )
{
	return( floor( x * 100000000.0 + 0.5 ) / 100000000.0 );
}

class CFanOpt : public CBEOOptimizerFan
{
public:
	CFanOpt()
	{
		updateDims( FanParamCount, 16 );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.1;
		p[ 1 ] = 0.1;
		p[ 2 ] = 0.0;
		p[ 3 ] = 0.0;
		p[ 4 ] = 0.3;
		p[ 5 ] = 0.0;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 4.0;
		p[ 1 ] = 1.0;
		p[ 2 ] = 1.5;
		p[ 3 ] = 1.5;
		p[ 4 ] = 3.0;
		p[ 5 ] = 1.5;
	}

	virtual double optcost( const double* const p ) const
	{
		rnd.init( 0 );

		CTester Tester;
		Tester.opt -> CostMult = roundp( p[ 0 ]);
		Tester.opt -> BestMult = roundp( p[ 1 ]);
		Tester.opt -> HistMult = roundp( p[ 2 ]);
		Tester.opt -> PrevMult = roundp( p[ 3 ]);
		Tester.opt -> CentMult = roundp( p[ 4 ]);
		Tester.opt -> CentOffs = roundp( p[ 5 ]);
//		Tester.init( 2, OptCorpus2D, 0.001, 5000, 10000, true, false );
		Tester.init( 10, OptCorpusND, 0.01, 120, 20000, false, false );

		Tester.run();

		return( Tester.Score );
	}
};

int main()
{
	CBEORnd rnd2;
	rnd2.init( 1 );

	double Params[ FanParamCount ];
	Params[ 0 ] = 1.19002383;
	Params[ 1 ] = 0.63504273;
	Params[ 2 ] = 0.57471040;
	Params[ 3 ] = 0.42249913;
	Params[ 4 ] = 1.26880243;
	Params[ 5 ] = 0.60540170;

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
