// BiteOpt's algorithm's self hyper-parameter optimization code.

#include <stdio.h>
#include "tests/tester.h"

const int FanParamCount = 21;
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
		p[ 10 ] = 0.06;
		p[ 11 ] = 10.0;
		p[ 12 ] = 10.0;
		p[ 13 ] = 11.0;
		p[ 14 ] = 1.95;
		p[ 15 ] = 0.05;
		p[ 16 ] = 0.05;
		p[ 17 ] = 0.05;
		p[ 18 ] = 0.05;
		p[ 19 ] = 2.0;
		p[ 20 ] = 2.0;
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
		p[ 10 ] = 0.08;
		p[ 11 ] = 96.999;
		p[ 12 ] = 96.999;
		p[ 13 ] = 16.0;
		p[ 14 ] = 2.05;
		p[ 15 ] = 0.333;
		p[ 16 ] = 0.333;
		p[ 17 ] = 1.0;
		p[ 18 ] = 1.0;
		p[ 19 ] = 12.0;
		p[ 20 ] = 30.0;
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
		Tester.opt -> ScutProb = roundp( p[ 10 ]);
		Tester.opt -> MantSizeSh = roundp( p[ 11 ]);
		Tester.opt -> MantSizeSh2 = roundp( p[ 12 ]);
		Tester.opt -> PopSizeBase = roundp( p[ 13 ]);
		Tester.opt -> PopSizeMult = roundp( p[ 14 ]);
		Tester.opt -> ParOptProb[ 0 ] = roundp( p[ 15 ]);
		Tester.opt -> ParOptProb[ 1 ] = roundp( p[ 16 ]);
		Tester.opt -> EntmProb[ 0 ] = roundp( p[ 17 ]);
		Tester.opt -> EntmProb[ 1 ] = roundp( p[ 18 ]);
		Tester.opt -> getParOpt() -> CentPow = roundp( p[ 19 ]);
		Tester.opt -> getParOpt() -> RadPow = roundp( p[ 20 ]);

		// Run low-dimensional and 14-dimensional test corpuses.

		Tester.init( 0.000001, 60, 2000, false );
		Tester.addCorpus( 2, TestCorpusAll, false, false );
		Tester.addCorpus( 2, OptCorpusNDRotOfs, true, false );
		Tester.run();

		double a1 = Tester.ItAvgl10n;
		double b1 = Tester.AtAvg;

		Tester.init( 0.01, 10, 12000, false );
		Tester.addCorpus( 14, OptCorpusNDRotOfsSol, true, true );
		Tester.run();

		double a2 = Tester.ItAvgl10n;
		double b2 = Tester.AtAvg / 1.15;

		Tester.init( 0.01, 10, 10000, false );
		Tester.addCorpus( 11, OptCorpusNDRotOfsSol, true, false );
		Tester.run();

		a2 += Tester.ItAvgl10n;
		b2 += Tester.AtAvg;

		Tester.init( 0.01, 10, 8000, false );
		Tester.addCorpus( 10, OptCorpusNDRotOfsSol, true, true );
		Tester.run();

		a2 += Tester.ItAvgl10n;
		b2 += Tester.AtAvg / 1.15;

		Tester.init( 0.01, 10, 7000, false );
		Tester.addCorpus( 8, OptCorpusNDRotOfsSol, true, false );
		Tester.run();

		a2 += Tester.ItAvgl10n;
		b2 += Tester.AtAvg;

		Tester.init( 0.01, 10, 5000, false );
		Tester.addCorpus( 6, OptCorpusNDRotOfsSol, true, false );
		Tester.run();

		a2 += Tester.ItAvgl10n;
		b2 += Tester.AtAvg;

		Tester.init( 0.01, 10, 3000, false );
		Tester.addCorpus( 4, OptCorpusNDRotOfsSol, true, true );
		Tester.run();

		a2 += Tester.ItAvgl10n;
		b2 += Tester.AtAvg / 1.15;

		a2 /= 6.0;
		b2 /= 6.0;

		// Apply weighting to obtained statistics and calculate score.

		double Score = ( 0.55 * a1 + 0.45 * a2 ) * ( 0.55 * b1 + 0.45 * b2 );

		return( Score );
	}
};

int main()
{
	CBiteRnd rnd2;
	rnd2.init( 1 );

	CFanOpt opt;
	opt.init( rnd2 );
	int i;

	for( i = 0; i < FanIterCount; i++ )
	{
		opt.optimize( rnd2 );

		printf( "%i\n// Cost=%f\n", i, opt.getBestCost() );

		static const char* const HyperNames[ FanParamCount ] = {
			"RandProb[ 0 ]",
			"RandProb[ 1 ]",
			"RandProb2[ 0 ]",
			"RandProb2[ 1 ]",
			"AllpProb[ 0 ]",
			"AllpProb[ 1 ]",
			"CentProb[ 0 ]",
			"CentProb[ 1 ]",
			"CentSpan[ 0 ]",
			"CentSpan[ 1 ]",
			"ScutProb",
			"MantSizeSh",
			"MantSizeSh2",
			"PopSizeBase",
			"PopSizeMult",
			"ParOptProb[ 0 ]",
			"ParOptProb[ 1 ]",
			"EntmProb[ 0 ]",
			"EntmProb[ 1 ]",
			"ParOpt.CentPow",
			"ParOpt.RadPow",
		};

		int j;

		for( j = 0; j < FanParamCount; j++ )
		{
			const double v = roundp( opt.getBestParams()[ j ]);

			printf( "%s = %.8f;\n", HyperNames[ j ], v );
		}
	}
}
