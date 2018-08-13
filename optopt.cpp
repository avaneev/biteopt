// BiteOpt's algorithm's hyper-parameter optimization code.

#include <stdio.h>
#include "tests/tester.h"
//#include "tester2.h"

const int FanParamCount = 6;
const int FanIterCount = 4000;

static double roundp( const double x )
{
	return( floor( x * 100000000.0 + 0.5 ) / 100000000.0 );
}

class CFanOpt : public CBiteOpt
{
public:
	CFanOpt()
	{
		updateDims( FanParamCount );
	}

	virtual void getMinValues( double* const p ) const
	{
		p[ 0 ] = 0.0;
		p[ 1 ] = 0.0;
		p[ 2 ] = 0.0;
		p[ 3 ] = 0.0;
		p[ 4 ] = 0.0;
		p[ 5 ] = 0.0;
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 1.0;
		p[ 1 ] = 1.0;
		p[ 2 ] = 1.0;
		p[ 3 ] = 3.0;
		p[ 4 ] = 1.0;
		p[ 5 ] = 1.0;
	}

	virtual double optcost( const double* const p )
	{
		CTester Tester;
//		CTester2 Tester;
		Tester.opt -> RandProb[ 0 ] = roundp( p[ 0 ]);
		Tester.opt -> RandProb[ 1 ] = roundp( p[ 1 ]);
		Tester.opt -> CentProb[ 0 ] = 1.0;
		Tester.opt -> CentProb[ 1 ] = roundp( p[ 2 ]);
		Tester.opt -> CentSpan = roundp( p[ 3 ]);
		Tester.opt -> AllpProb[ 0 ] = roundp( p[ 4 ]);
		Tester.opt -> AllpProb[ 1 ] = 1.0;
		Tester.opt -> ScutProb = 0.11;
		Tester.opt -> RandProb2[ 0 ] = roundp( p[ 5 ]);
		Tester.opt -> RandProb2[ 1 ] = 1.0;

		rnd.init( 0 );
		Tester.init( 2, TestCorpusAll, 0.000001, 60, 2000, false, false );

		Tester.run();

		return( Tester.Score );
	}
};

int main()
{
	CBiteRnd rnd2;
	rnd2.init( 1 );

	double Params[ FanParamCount ];
	Params[ 0 ] = 4.54665746;
	Params[ 1 ] = 0.64152578;
	Params[ 2 ] = 0.55839563;

	CFanOpt opt;
	opt.init( rnd2/*, Params*/ );
	int i;

	for( i = 0; i < FanIterCount; i++ )
	{
		opt.optimize( rnd2 );

		printf( "%i\n// Cost=%f\n", i, opt.getBestCost() );

		int j;

		for( j = 0; j < FanParamCount; j++ )
		{
			Params[ j ] = roundp( opt.getBestParams()[ j ]);
		}

		printf( "RandProb[ 0 ] = %.8f;\n", Params[ 0 ]);
		printf( "RandProb[ 1 ] = %.8f;\n", Params[ 1 ]);
		printf( "CentProb[ 0 ] = %.8f;\n", 1.0 );
		printf( "CentProb[ 1 ] = %.8f;\n", Params[ 2 ]);
		printf( "CentSpan = %.8f;\n", Params[ 3 ]);
		printf( "AllpProb[ 0 ] = %.8f;\n", Params[ 4 ]);
		printf( "AllpProb[ 1 ] = %.8f;\n", 1.0 );
		printf( "ScutProb = %.8f;\n", 0.11 );
		printf( "RandProb2[ 0 ] = %.8f;\n", Params[ 5 ]);
		printf( "RandProb2[ 1 ] = %.8f;\n", 1.0 );
	}

	return( 0 );
}
