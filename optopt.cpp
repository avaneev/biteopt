// BiteOpt's algorithm's hyper-parameter optimization code.

#include <stdio.h>
#include "tests/tester.h"

const int FanParamCount = 14;
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
		p[ 6 ] = 0.0;
		p[ 7 ] = 0.0;
		p[ 8 ] = 0.0;
		p[ 9 ] = 0.0;
		p[ 10 ] = 10.0;
		p[ 11 ] = 10.0;
		p[ 12 ] = 11.0;
		p[ 13 ] = 1.0;
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
		p[ 8 ] = 3.0;
		p[ 9 ] = 3.0;
		p[ 10 ] = 96.999;
		p[ 11 ] = 96.999;
		p[ 12 ] = 16.0;
		p[ 13 ] = 3.0;
	}

	virtual double optcost( const double* const p )
	{
		CTester Tester;
		Tester.opt -> RandProb[ 0 ] = roundp( p[ 0 ]);
		Tester.opt -> RandProb[ 1 ] = roundp( p[ 1 ]);
		Tester.opt -> RandProb2[ 0 ] = roundp( p[ 2 ]);
		Tester.opt -> RandProb2[ 1 ] = roundp( p[ 3 ]);
		Tester.opt -> AllpProb[ 0 ] = roundp( p[ 4 ]);
		Tester.opt -> AllpProb[ 1 ] = roundp( p[ 5 ]);
		Tester.opt -> CentProb[ 0 ] = roundp( p[ 6 ]);
		Tester.opt -> CentProb[ 1 ] = roundp( p[ 7 ]);
		Tester.opt -> CentSpan[ 0 ] = roundp( p[ 8 ]);
		Tester.opt -> CentSpan[ 1 ] = roundp( p[ 9 ]);
		Tester.opt -> ScutProb = 0.06;
		Tester.opt -> MantSizeSh = roundp( p[ 10 ]);
		Tester.opt -> MantSizeSh2 = roundp( p[ 11 ]);
		Tester.opt -> PopSizeBase = roundp( p[ 12 ]);
		Tester.opt -> PopSizeMult = roundp( p[ 13 ]);

		// Run low-dimensional and 14-dimensional test corpuses.

		Tester.init( 0.000001, 70, 2000, false );
		Tester.addCorpus( 2, TestCorpusAll, false, false );
		Tester.addCorpus( 2, OptCorpusNDRotOfs, true, false );
		Tester.run();

		double a1 = Tester.ItAvg2l10n;
		double b1 = Tester.AtAvg;

		Tester.init( 0.01, 35, 14000, false );
		Tester.addCorpus( 14, OptCorpusNDRotOfsSol, true, false );
		Tester.run();

		double a2 = Tester.ItAvg2l10n;
		double b2 = Tester.AtAvg;

		// Apply weighting to obtained statistics and calculate score.

		double Score = ( 0.65 * a1 + 0.35 * a2 ) * ( 0.75 * b1 + 0.25 * b2 );

		return( Score );
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
		printf( "RandProb2[ 0 ] = %.8f;\n", Params[ 2 ]);
		printf( "RandProb2[ 1 ] = %.8f;\n", Params[ 3 ]);
		printf( "AllpProb[ 0 ] = %.8f;\n", Params[ 4 ]);
		printf( "AllpProb[ 1 ] = %.8f;\n", Params[ 5 ]);
		printf( "CentProb[ 0 ] = %.8f;\n", Params[ 6 ]);
		printf( "CentProb[ 1 ] = %.8f;\n", Params[ 7 ]);
		printf( "CentSpan[ 0 ] = %.8f;\n", Params[ 8 ]);
		printf( "CentSpan[ 1 ] = %.8f;\n", Params[ 9 ]);
		printf( "ScutProb = %.8f;\n", 0.06 );
		printf( "MantSizeSh = %.8f;\n", Params[ 10 ]);
		printf( "MantSizeSh2 = %.8f;\n", Params[ 11 ]);
		printf( "PopSizeBase = %.8f;\n", Params[ 12 ]);
		printf( "PopSizeMult = %.8f;\n", Params[ 13 ]);
	}

	return( 0 );
}
