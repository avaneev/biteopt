//$ nocpp

#include <stdio.h>
#include <emmintrin.h>
#include "biteopt.h"
#include "testfn.h"

CBiteRnd rnd;

/**
 * Function test corpus class.
 */

class CTester
{
public:
	/**
	 * Function optimizer class.
	 */

	class CTestOpt : public CBiteOpt
	{
	public:
		const CTestFn* fn; ///< Test function.
			///<
		int Dims; ///< Dimensions in the function.
			///<
		double* minv; ///< Minimal parameter values.
			///<
		double* maxv; ///< Maximal parameter values.
			///<
		double optv; ///< Optimal value.
			///<
		double* signs; ///< Signs to apply to function parameters.
			///<
		double* tp; ///< Temporary parameter storage.
			///<

		CTestOpt()
			: minv( NULL )
			, maxv( NULL )
			, signs( NULL )
			, tp( NULL )
		{
		}

		~CTestOpt()
		{
			delete[] minv;
			delete[] maxv;
			delete[] signs;
			delete[] tp;
		}

		void updateDims( const int aDims, const int aFanSize = 0 )
		{
			Dims = aDims;
			delete[] minv;
			delete[] maxv;
			delete[] signs;
			delete[] tp;
			minv = new double[ Dims ];
			maxv = new double[ Dims ];
			signs = new double[ Dims ];
			tp = new double[ Dims ];
			CBiteOpt :: updateDims( Dims, aFanSize );
		}

		virtual void getMinValues( double* const p ) const
		{
			int i;

			for( i = 0; i < Dims; i++ )
			{
				p[ i ] = minv[ i ];
			}
		}

		virtual void getMaxValues( double* const p ) const
		{
			int i;

			for( i = 0; i < Dims; i++ )
			{
				p[ i ] = maxv[ i ];
			}
		}

		virtual double optcost( const double* const p ) const
		{
			int i;

			for( i = 0; i < Dims; i++ )
			{
				tp[ i ] = p[ i ] * signs[ i ];
			}

			return( (*fn -> CalcFunc)( tp, Dims ));
		}
	};

	int DefDims; ///< The default number of dimensions to use.
		///<
	CTestOpt* opt; ///< Optimizer.
		///<
	double ItAvg; ///< Average convergence time after run().
		///<
	double ItRtAvg; ///< Average ratio of std.dev and average after run().
		///<
	double RjAvg; ///< Average number of rejects.
		///<
	uint64_t tc; ///< Processor clock ticks used in evaluation.
		///<
	double Score; ///< Optimization score.
		///<

	CTester()
		: opt( new CTestOpt() )
	{
	}

	~CTester()
	{
		delete opt;
	}

	/**
	 * Function initalizes the tester.
	 *
 	 * @param aDefDims The number of dimensions in each function with variable
	 * number of dimensions.
	 * @param Corpus NULL-limited list of test functions.
	 * @param Threshold Objective function value threshold - stop condition.
	 * @param IterCount The number of attempts to solve a function to perform.
	 * @param InnerIterCount The maximal number of solver iterations to
	 * perform.
	 * @param DoRandomize "True" if randomization of value range and value
	 * sign should be performed.
	 * @param DoPrint "True" if results should be printed to "stdout".
	 */

	void init( const int aDefDims, const CTestFn** Corpus,
		const double aThreshold, const int aIterCount,
		const int aInnerIterCount, const bool aDoRandomize,
		const bool aDoPrint )
	{
		DefDims = aDefDims;
		CostThreshold = aThreshold;
		IterCount = aIterCount;
		InnerIterCount = aInnerIterCount;
		DoRandomize = aDoRandomize;
		DoPrint = aDoPrint;
		FnCount = 0;
		Funcs = Corpus;

		while( true )
		{
			if( *Corpus == NULL )
			{
				break;
			}

			FnCount++;
			Corpus++;
		}
	}

	/**
	 * Function runs the test. On return, the ItAvg, ItRtAvg and tc variables
	 * will be updated.
	 */

	void run()
	{
		ItAvg = 0.0;
		ItRtAvg = 0.0;
		RjAvg = 0.0;
		tc = 0;
		int* Iters = new int[ IterCount ];
		int k;

		for( k = 0; k < FnCount; k++ )
		{
			_mm_empty();
			double AvgIter = 0;
			double AvgCost = 0.0;
			double AvgRjCost = 0.0;
			int Rej = 0;
			int j;
			int i;

			opt -> fn = Funcs[ k ];
			const int Dims = ( opt -> fn -> Dims == 0 ?
				DefDims : opt -> fn -> Dims );

			opt -> updateDims( Dims );

			if( !DoRandomize )
			{
				for( i = 0; i < Dims; i++ )
				{
					opt -> signs[ i ] = 1.0;
				}
			}

			for( j = 0; j < IterCount; j++ )
			{
				if( opt -> fn -> ParamFunc != NULL )
				{
					opt -> optv = (*opt -> fn -> ParamFunc)(
						opt -> minv, opt -> maxv, Dims );
				}
				else
				{
					for( i = 0; i < Dims; i++ )
					{
						opt -> minv[ i ] = opt -> fn -> RangeMin;
						opt -> maxv[ i ] = opt -> fn -> RangeMax;
					}

					opt -> optv = opt -> fn -> OptValue;
				}

				if( DoRandomize )
				{
					for( i = 0; i < Dims; i++ )
					{
						opt -> signs[ i ] =
							( rnd.getRndValue() < 0.5 ? 1.0 : -1.0 );
					}
				}

				_mm_empty();
				opt -> init( rnd );
				i = 0;

				while( true )
				{
					const uint64_t t1 = __rdtsc();
					opt -> optimize( rnd );
					tc += __rdtsc() - t1;

					i++;

					if( opt -> getBestCost() - opt -> optv < CostThreshold )
					{
						AvgCost += opt -> getBestCost();
						Iters[ j ] = i;
						AvgIter += i;
						break;
					}

					if( i == InnerIterCount )
					{
						AvgRjCost += opt -> getBestCost();
						Rej++;
						Iters[ j ] = -1;
						break;
					}
				}
			}

			AvgCost /= ( IterCount - Rej );
			AvgRjCost /= Rej;
			double Avg;
			double RMS;

			if( Rej >= IterCount )
			{
				Avg = 1000000.0;
				RMS = 1000000.0;
			}
			else
			{
				Avg = AvgIter / ( IterCount - Rej );
				RMS = 0.0;

				for( j = 0; j < IterCount; j++ )
				{
					if( Iters[ j ] >= 0 )
					{
						const double v = Iters[ j ] - Avg;
						RMS += v * v;
					}
				}

				RMS = sqrt( RMS / ( IterCount - Rej ));
			}

			ItAvg += Avg;
			ItRtAvg += RMS / Avg;
			RjAvg += (double) Rej / IterCount;

			if( DoPrint )
			{
				printf( "AIt:%6.0f RIt:%6.0f Rj:%5.2f%% C:%11.8f RjC:%7.4f "
					"%s_%i\n", Avg, RMS, 100.0 * Rej / IterCount, AvgCost,
					AvgRjCost, opt -> fn -> Name, Dims );
			}
		}

		_mm_empty();
		ItAvg /= FnCount;
		ItRtAvg /= FnCount;
		RjAvg /= FnCount;
//		Score = fabs( ItRtAvg - 0.160 ) * 50000.0 + ItAvg *
//			( 1.0 + RjAvg * 100.0 );
		Score = RjAvg * 100.0 + fabs( ItAvg - 360.0/*360.0*/ ) * 0.1 +
			fabs( ItRtAvg - 0.250 ) * 50;

		if( DoPrint )
		{
			printf( "ItAvg: %.1f (avg convergence time)\n", ItAvg );
			printf( "ItRtAvg: %.6f (avg ratio of std.dev and average)\n",
				ItRtAvg );

			printf( "RjAvg: %.2f%% (avg percentage of rejects)\n",
				RjAvg * 100.0 );

			printf( "Score: %f\n", Score );
			printf( "Ticks: %llu\n", tc );
		}

		delete[] Iters;
	}

protected:
	const CTestFn** Funcs; ///< Test functions corpus.
		///<
	int FnCount; ///< Test function count.
		///<
	double CostThreshold; ///< Cost threshold (finish criteria).
		///<
	int IterCount; ///< Iteration count.
		///<
	int InnerIterCount; ///< Inner iteration count (the number of
		///< optimization calls).
	bool DoRandomize; ///< Randomize argument signs.
		///<
	bool DoPrint; ///< Print results to stdout.
		///<
};
