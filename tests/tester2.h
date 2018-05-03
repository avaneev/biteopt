//$ nocpp

#include <stdio.h>
#include <string.h>
#include "biteopt.h"
#include "aset.h"

/**
 * Function test corpus class.
 */

class CTester2
{
public:
	/**
	 * Function optimizer class.
	 */

	class CTestOpt : public CBiteOptDeep
	{
	public:
		const CFunc* Func;

		virtual void getMinValues( double* const p ) const
		{
			int i;

			for( i = 0; i < Func -> Dims; i++ )
			{
				p[ i ] = Func -> MinValues[ i ];
			}
		}

		virtual void getMaxValues( double* const p ) const
		{
			int i;

			for( i = 0; i < Func -> Dims; i++ )
			{
				p[ i ] = Func -> MaxValues[ i ];
			}
		}

		virtual double optcost( const double* const p )
		{
			return( (*Func -> CalcFunc)( (int*) &Func -> Dims, (double*) p ));
		}
	};

	CTestOpt* opt; ///< Optimizer.
		///<
	double Score; ///< Optimization score.
		///<

	CTester2()
		: opt( new CTestOpt() )
	{
	}

	~CTester2()
	{
		delete opt;
	}

	/**
	 * Function runs the test. On return, Score variable will be updated.
	 */

	void run( const bool DoPrint = false )
	{
		const int AttemptCount = 10;
		const int IterCount = 2500;
		const int RoundCount = 3;
		CBiteRnd rnd;
		rnd.init( 1 );
		Score = 0.0;
		double Score2 = 0.0;
		int l;

		for( l = 0; l < RoundCount; l++ )
		{
			int SuccessCount = 0;
			int FuncCount = 0;
			double sums = 0.0;
			int k;

			for( k = 0; k < corpus_size; k++ )
			{
//				if( corpus[ k ] -> Dims > 2 )
//					continue;
//				if( corpus[ k ] -> Dims < 3 || corpus[ k ] -> Dims > 9 )
//					continue;
//				if( corpus[ k ] -> Dims < 10 || corpus[ k ] -> Dims > 30 )
//					continue;
//				if( corpus[ k ] -> Dims < 31 )
//					continue;

				FuncCount++;
				opt -> Func = corpus[ k ];
				opt -> updateDims( opt -> Func -> Dims, 1 );

				double mv1 = opt -> Func -> OptValue *
					( opt -> Func -> OptValue < 0 ? 0.99 : 1.01 );

				double mv2 = opt -> Func -> OptValue + 0.01;
				double ExpMinValue = ( mv1 > mv2 ? mv1 : mv2 );
				double BestCost = 1e100;
				double WorstCost = -1e100;
				int j;

				for( j = 0; j < AttemptCount; j++ )
				{
					opt -> init( rnd );
					int i;

					for( i = 0; i < IterCount; i++ )
					{
						opt -> optimize( rnd );

						if( opt -> getBestCost() < BestCost )
						{
							BestCost = opt -> getBestCost();

							if( BestCost <= ExpMinValue )
							{
								WorstCost = BestCost;
								j = AttemptCount;
								break;
							}
						}
					}

					if( opt -> getBestCost() > WorstCost )
					{
						WorstCost = opt -> getBestCost();
					}
				}

				const double s = 1.0 -
					fabs( BestCost - opt -> Func -> OptValue ) /
					( fabs( WorstCost - opt -> Func -> OptValue ) + 1e-10 );

				const double d = BestCost - opt -> Func -> OptValue;
				char sm;

				if( BestCost <= ExpMinValue )
				{
					sm = ' ';
					SuccessCount++;
				}
				else
				{
					sm = '*';
				}

				if( DoPrint )
				{
					printf( "%-15s %3i %c %15.9f %15.9f\n",
						opt -> Func -> Name,
						opt -> Func -> Dims, sm, BestCost,
						opt -> Func -> OptValue );
				}

				sums += s;
			}

			Score += sums / FuncCount;
			Score2 += -100.0 * SuccessCount / FuncCount;
		}

		Score /= RoundCount;
		Score2 /= RoundCount;

		if( DoPrint )
		{
			printf( "Success: %.4f%%\n", -Score2 );
			printf( "Score: %.4f\n", Score );
		}

		Score = Score2;
	}
};
