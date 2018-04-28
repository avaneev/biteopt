// BiteOpt's algorithm's hyper-parameter optimization code.

#include <stdio.h>
#include "tester.h"
//#include "tester2.h"

const int FanParamCount = 3;
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
	}

	virtual void getMaxValues( double* const p ) const
	{
		p[ 0 ] = 1.0;
		p[ 1 ] = 1.0;
		p[ 2 ] = 3.0;
	}

	virtual double optcost( const double* const p )
	{
		CTester Tester;
//		CTester2 Tester;
		Tester.opt -> MinxMult = 0.5;
		Tester.opt -> RandProb = roundp( p[ 0 ]);
		Tester.opt -> CentProb = roundp( p[ 1 ]);
		Tester.opt -> CentSpan = roundp( p[ 2 ]);
		Tester.opt -> AllpProb = 0.2;
		Tester.opt -> ScutProb = 0.11;
		Tester.opt -> RandProb2 = 0.25;

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

		printf( "// Cost=%f\n", opt.getBestCost() );

		int j;

		for( j = 0; j < FanParamCount; j++ )
		{
			Params[ j ] = roundp( opt.getBestParams()[ j ]);
		}

		printf( "MinxMult = %.8f;\n", 0.5 );
		printf( "RandProb = %.8f;\n", Params[ 0 ]);
		printf( "CentProb = %.8f;\n", Params[ 1 ]);
		printf( "CentSpan = %.8f;\n", Params[ 2 ]);
		printf( "AllpProb = %.8f;\n", 0.2 );
		printf( "ScutProb = %.8f;\n", 0.11 );
		printf( "RandProb2 = %.8f;\n", 0.25 );
	}

	return( 0 );
}
