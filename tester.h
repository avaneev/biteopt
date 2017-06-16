//$ nocpp

#include <stdio.h>
#include "bitefan.h"
#include "testfn.h"

CBEORnd rnd;

/**
 * Function test corpus class.
 *
 * @tparam Dims The number of dimensions in each function.
 */

template< int Dims >
class CTester
{
public:
	/**
	 * Function optimizer class.
	 */

	class CTestOpt : public CBEOOptimizerFan< Dims >
	{
	public:
		const CTestFn* fn; ///< Test function.
			///<
		bool DoRandomize; ///< Randomize signs and function range.
			///<
		double signs[ Dims ]; ///< Signs to apply to function parameters.
			///<

		virtual void getMinValues( double* const p ) const
		{
			int i;

			if( DoRandomize )
			{
				for( i = 0; i < Dims; i++ )
				{
					p[ i ] = fn -> RangeMin *
						( 0.5 + rnd.getRndValue() * 0.5 );
				}
			}
			else
			{
				for( i = 0; i < Dims; i++ )
				{
					p[ i ] = fn -> RangeMin;
				}
			}
		}

		virtual void getMaxValues( double* const p ) const
		{
			int i;

			if( DoRandomize )
			{
				for( i = 0; i < Dims; i++ )
				{
					p[ i ] = fn -> RangeMax *
						( 0.5 + rnd.getRndValue() * 0.5 );
				}
			}
			else
			{
				for( i = 0; i < Dims; i++ )
				{
					p[ i ] = fn -> RangeMax;
				}
			}
		}

		virtual double optcost( const double* const p ) const
		{
			double pp[ Dims ];
			int i;

			for( i = 0; i < Dims; i++ )
			{
				pp[ i ] = p[ i ] * signs[ i ];
			}

			return( (*fn -> Calc)( pp, Dims ));
		}
	};

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
	 * @param Corpus NULL-limited list of test functions.
	 * @param Threshold Objective function value threshold - stop condition.
	 * @param IterCount The number of attempts to solve a function to perform.
	 * @param InnerIterCount The maximal number of solver iterations to
	 * perform.
	 * @param DoRandomize "True" if randomization of value range and value
	 * sign should be performed.
	 * @param DoPrint "True" if results should be printed to "stdout".
	 */

	void init( const CTestFn** Corpus, const double aThreshold,
		const int aIterCount, const int aInnerIterCount,
		const bool aDoRandomize, const bool aDoPrint )
	{
		CostThreshold = aThreshold;
		IterCount = aIterCount;
		InnerIterCount = aInnerIterCount;
		DoRandomize = aDoRandomize;
		DoPrint = aDoPrint;
		FnCount = 0;

		while( true )
		{
			const CTestFn* const fn = *Corpus;

			if( fn == NULL )
			{
				break;
			}

			if( fn -> Dims == Dims || fn -> Dims == 0 )
			{
				Funcs[ FnCount ] = fn;
				FnCount++;
			}

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
		int k;

		for( k = 0; k < FnCount; k++ )
		{
			int Iters[ IterCount ];
			double AvgIter = 0;
			double AvgCost = 0.0;
			double AvgRjCost = 0.0;
			int Rej = 0;
			int j;
			int i;

			opt -> fn = Funcs[ k ];
			opt -> DoRandomize = DoRandomize;

			if( !DoRandomize )
			{
				for( i = 0; i < Dims; i++ )
				{
					opt -> signs[ i ] = 1.0;
				}
			}

			for( j = 0; j < IterCount; j++ )
			{
				if( DoRandomize )
				{
					for( i = 0; i < Dims; i++ )
					{
						opt -> signs[ i ] =
							( rnd.getRndValue() < 0.5 ? 1.0 : -1.0 );
					}
				}

				opt -> init( rnd );
				i = 0;

				while( true )
				{
					const uint64_t t1 = __rdtsc();
					opt -> optimize( rnd );
					tc += __rdtsc() - t1;

					i++;

					if( opt -> getBestCost() -
						opt -> fn -> OptValue < CostThreshold )
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

			if( DoPrint )
			{
				printf( "AIt:%6.0f RIt:%6.0f Rj:%5.2f%% C:%11.8f RjC:%7.4f %s\n",
					Avg, RMS, 100.0 * Rej / IterCount, AvgCost, AvgRjCost,
					opt -> fn -> Name );
			}

			ItAvg += Avg;
			ItRtAvg += RMS / Avg;
			RjAvg += (double) Rej / IterCount;
		}

		ItAvg /= FnCount;
		ItRtAvg /= FnCount;
		RjAvg /= FnCount;

		if( DoPrint )
		{
			printf( "ItAvg: %.1f (avg convergence time)\n", ItAvg );
			printf( "ItRtAvg: %.3f (avg ratio of std.dev and average)\n",
				ItRtAvg );

			printf( "RjAvg: %.2f%% (avg percentage of rejects)\n",
				RjAvg * 100.0 );

			printf( "Ticks: %llu\n", tc );
		}
	}

protected:
	static const int FnCountMax = 50; ///< Funcs array capacity.
		///<
	const CTestFn* Funcs[ FnCountMax ]; ///< Test functions corpus.
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
