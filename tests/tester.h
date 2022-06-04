//$ nocpp

#include <conio.h>
#include <stdio.h>
#include "tester_types.h"
#include "../biteopt.h"
#include "../spheropt.h"
#include "../smaesopt.h"
#include "../nmsopt.h"
#include "../deopt.h"
//#include "../other/ccmaes.h"

#define OPT_CLASS CBiteOpt//CDEOpt//CSMAESOpt//CSpherOpt//CBiteOptDeep//CNMSeqOpt//CCMAESOpt//
#define OPT_DIMS_PARAMS Dims // updateDims() parameters.
//#define OPT_PLATEAU_MUL 256 // Uncomment to enable plateau check.
//#define EVALBINS 1
#define OPT_STATS 0 // Set to 1 to enable histogram statistics output.
#define OPT_TIME 0 // Set to 1 to evaluate timings.

#if 0
	#define OPT_THREADS 1
	#include "../../../libvox/Sources/Core/CWorkerThreadPool.h"
	using namespace vox;
	CSyncObject StatsSync;
#else // OPT_THREADS
	#define OPT_THREADS 0
#endif // OPT_THREADS

#if defined( _WIN32 )
	#include <windows.h>
	#define OPT_PERF
#endif // defined( _WIN32 )

/**
 * Summary optimization statistics class.
 */

class CSumStats
{
public:
	int SumIt; ///< The overall number of performed function evaluations in
		///< successful attempts.
		///<
	int SumItAll; ///< The overall number of performed function evaluations in
		///< all attempts.
		///<
	double SumIt_l10n; ///< = sum( log( It ) / log( 10 )) in completed
		///< attempts.
		///<
	int SumItImpr; ///< Sum of improving iterations across successful
		///< attempts.
		///<
	int SumItImprAll; ///< Sum of improving iterations across all attempts.
		///<
	double SumRjCost; ///< Summary unbiased costs detected in all rejected
		///< attempts.
		///<
	double SumRMS; ///< = sum( RMS ).
		///<
	double SumRMS_l10n; ///< = sum( log( RMS / N ) / log( 10 )).
		///<
	int SumRMSCount; ///< The number of elements summed in the SumRMS.
		///<
	int TotalAttempts; ///< The overall number of function evaluation attempts.
		///< Must be set externally.
		///<
	int ComplAttempts; ///< The overall number of successful function
		///< attempts.
		///<
	int ComplFuncs; ///< The total number of optimized functions.
		///<

	#if defined( EVALBINS )
		static const int binspan = 50;
		static const int binc = 1 + binspan * 2;
		double bins[ binc ];
		double RMSSpan;
	#endif // defined( EVALBINS )

	#if OPT_STATS
		int SumSelsMap[ CBiteOpt :: MaxHistCount * 8 ]; ///< Sums of
			///< individual choices associated with histograms over all
			///< successful attempts.
			///<
		double SumAvgSels[ CBiteOpt :: MaxHistCount ]; ///< Sum of average
			///< histogram choices over all successful attempts.
			///<
		double SumDevSels[ CBiteOpt :: MaxHistCount ]; ///< Sum of average
			///< histogram choice deviations (squared) over all successful
			///< attempts.
			///<
	#endif // OPT_STATS

	#if OPT_TIME
		double TimeOpt; ///< Overall time spent in optimize() function,
			///< seconds.
		double TimeFunc; ///< Overall time spent in objective function,
			///< seconds.
	#endif // OPT_TIME

	void clear()
	{
		SumIt = 0;
		SumItAll = 0;
		SumIt_l10n = 0.0;
		SumItImpr = 0;
		SumItImprAll = 0;
		SumRjCost = 0.0;
		SumRMS = 0.0;
		SumRMS_l10n = 0.0;
		SumRMSCount = 0;
		ComplAttempts = 0;
		ComplFuncs = 0;

		#if defined( EVALBINS )

		RMSSpan = 3.0;

		int k;

		for( k = 0; k < binc; k++ )
		{
			bins[ k ] = 0.0;
		}

		#endif // defined( EVALBINS )

		#if OPT_STATS
		memset( SumSelsMap, 0, sizeof( SumSelsMap ));
		memset( SumAvgSels, 0, sizeof( SumAvgSels ));
		memset( SumDevSels, 0, sizeof( SumDevSels ));
		#endif // OPT_STATS

		#if OPT_TIME
		TimeOpt = 0.0;
		TimeFunc = 0.0;
		#endif // OPT_TIME
	}
};

/**
 * Function optimization statistics class, across all attempts.
 */

class CFuncStats
{
public:
	double MinCost; ///< Minimal cost detected in all attempts.
		///<
	double SumRjCost; ///< Summary rejects cost.
		///<
	int ComplAttempts; ///< The number of completed attempts.
		///<
	int SumItCompl; ///< Sum of iterations in completed attempts.
		///<

	CFuncStats()
	{
		clear();
	}

	void clear()
	{
		MinCost = 1e300;
		SumRjCost = 0.0;
		ComplAttempts = 0;
		SumItCompl = 0;
	}
};

/**
 * Function optimizer class.
 */

class CTestOpt : public OPT_CLASS
#if OPT_THREADS
	, public CWorkerThread
#endif // OPT_THREADS
{
public:
	const CTestFn* fn; ///< Test function.
		///<
	double CostThreshold; ///< Successful optimization attempt cost threshold.
		///<
	CSumStats* SumStats; ///< Pointer to summary statistics object.
		///<
	CFuncStats* FuncStats; ///< Pointer to function's statistics
		///< object.
		///<
	CBiteRnd rnd; ///< Random number generator.
		///<
	int Dims; ///< Dimensions in the function.
		///<
	int MaxIters; ///< Maximal number of iterations per attempt.
		///<
	int Index; ///< Attempt index.
		///<
	int* Iters; ///< Pointer to array containing finishing iteration counts
		///< for each attempt.
		///<
	double* minv; ///< Minimal parameter values.
		///<
	double* maxv; ///< Maximal parameter values.
		///<
	double optv; ///< Optimal value.
		///<
	double* shifts; ///< Shifts to apply to function parameters.
		///<
	double* signs; ///< Signs or scales to apply to function parameters.
		///<
	double** rots; ///< Rotation matrix to apply to function parameters,
		///< should be provided externally if DoRandomizeAll == true.
		///<
	double* tp; ///< Temporary parameter storage.
		///<
	double* tp2; ///< Temporary parameter storage 2.
		///<
	bool DoRandomize; ///< Apply value randomization.
		///<
	bool DoRandomizeAll; ///< Apply all randomizations, "false" only
		///< shift and scale randomizations will be applied.
		///<

	#if OPT_STATS
	int SumSels[ CBiteOpt :: MaxHistCount ]; ///< Sum of histogram choices.
		///<
	int* Sels[ CBiteOpt :: MaxHistCount ]; ///< Histogram choices at each
		///< optimization step.
		///<
	int SelAlloc; ///< The number of items allocated in each Sels element.
		///<
	#endif // OPT_STATS

	#if OPT_THREADS && OPT_TIME
		double funct; ///< Function evaluation timing.
			///<
	#endif // OPT_THREADS && OPT_TIME

	CTestOpt()
		: Dims( 0 )
		, minv( NULL )
		, maxv( NULL )
		, shifts( NULL )
		, signs( NULL )
		, tp( NULL )
		, tp2( NULL )
	#if OPT_STATS
		, SelAlloc( 0 )
	#endif // OPT_STATS
	{
		#if OPT_STATS
		memset( Sels, 0, sizeof( Sels ));
		#endif // OPT_STATS
	}

	virtual ~CTestOpt()
	{
		delete[] minv;
		delete[] maxv;
		delete[] shifts;
		delete[] signs;
		delete[] tp;
		delete[] tp2;

		#if OPT_STATS
		int i;

		for( i = 0; i < getHistCount(); i++ )
		{
			delete[] Sels[ i ];
		}
		#endif // OPT_STATS
	}

	void updateDims( const int aDims )
	{
		if( Dims == aDims )
		{
			return;
		}

		Dims = aDims;
		delete[] minv;
		delete[] maxv;
		delete[] shifts;
		delete[] signs;
		delete[] tp;
		delete[] tp2;
		minv = new double[ Dims ];
		maxv = new double[ Dims ];
		shifts = new double[ Dims ];
		signs = new double[ Dims ];
		tp = new double[ Dims ];
		tp2 = new double[ Dims ];

		OPT_CLASS :: updateDims( OPT_DIMS_PARAMS );
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

	virtual double optcost( const double* const p )
	{
		#if OPT_THREADS && OPT_TIME
			TClock t1( CSystem :: getClock() );
		#endif // OPT_THREADS && OPT_TIME

		if( !DoRandomize )
		{
			const double fv = (*fn -> CalcFunc)( p, Dims );

			#if OPT_THREADS && OPT_TIME
				funct += CSystem :: getClockDiffSec( t1 );
			#endif // OPT_THREADS && OPT_TIME

			return( fv );
		}

		int i;

		for( i = 0; i < Dims; i++ )
		{
			tp2[ i ] = p[ i ] * signs[ i ] + shifts[ i ];
		}

		if( !DoRandomizeAll )
		{
			const double fv = (*fn -> CalcFunc)( tp2, Dims );

			#if OPT_THREADS && OPT_TIME
				funct += CSystem :: getClockDiffSec( t1 );
			#endif // OPT_THREADS && OPT_TIME

			return( fv );
		}

		for( i = 0; i < Dims; i++ )
		{
			tp[ i ] = 0.0;
			int j;

			for( j = 0; j < Dims; j++ )
			{
				tp[ i ] += rots[ i ][ j ] * tp2[ j ];
			}
		}

		const double fv = (*fn -> CalcFunc)( tp, Dims );

		#if OPT_THREADS && OPT_TIME
			funct += CSystem :: getClockDiffSec( t1 );
		#endif // OPT_THREADS && OPT_TIME

		return( fv );
	}

	#if OPT_STATS
	double calcDevSel( const int h, const int c, const double avg,
		int* const SumSelsMap )
	{
		int* const s = Sels[ h ];
		double sd = 0.0;
		int i;

		for( i = 0; i < c; i++ )
		{
			const double d = s[ i ] - avg;
			SumSelsMap[ s[ i ]]++;
			sd += d * d;
		}

		return( sd / c );
	}
	#endif // OPT_STATS

	void performOpt()
	{
		init( rnd );

		#if OPT_STATS
		int k;

		if( SelAlloc < MaxIters )
		{
			SelAlloc = MaxIters;

			for( k = 0; k < getHistCount(); k++ )
			{
				delete[] Sels[ k ];
				Sels[ k ] = new int[ MaxIters ];
			}
		}

		memset( SumSels, 0, sizeof( SumSels ));
		#endif // OPT_STATS

		int i = 0;
		int ImprIters = 0;

		#if OPT_THREADS && OPT_TIME
			double optt = 0.0;
			funct = 0.0;
		#endif // OPT_THREADS && OPT_TIME

		while( true )
		{
			#if OPT_THREADS && OPT_TIME
				TClock t1( CSystem :: getClock() );
			#endif // OPT_THREADS && OPT_TIME

			const int sc = optimize( rnd );

			#if OPT_THREADS && OPT_TIME
				optt += CSystem :: getClockDiffSec( t1 );
			#endif // OPT_THREADS && OPT_TIME

			#if OPT_STATS
				CBiteOptHistBase** const h = getHists();

				for( k = 0; k < getHistCount(); k++ )
				{
					const int s = h[ k ] -> getSel();
					Sels[ k ][ i ] = s;
					SumSels[ k ] += s;
				}
			#endif // OPT_STATS

			if( sc == 0 )
			{
				ImprIters++;
			}

			i++;

			if( getBestCost() <= CostThreshold )
			{
				#if OPT_STATS
				double DevSels[ CBiteOpt :: MaxHistCount ];
				int SumSelsMap[ CBiteOpt :: MaxHistCount * 8 ];
				memset( SumSelsMap, 0, sizeof( SumSelsMap ));

				for( k = 0; k < getHistCount(); k++ )
				{
					DevSels[ k ] = calcDevSel( k, i,
						(double) SumSels[ k ] / i, SumSelsMap + k * 8 );
				}
				#endif // OPT_STATS

				#if OPT_THREADS
					VOXSYNC( StatsSync );

					#if OPT_TIME
					SumStats -> TimeOpt += optt;
					SumStats -> TimeFunc += funct;
					#endif // OPT_TIME
				#endif // OPT_THREADS

				if( getBestCost() < FuncStats -> MinCost )
				{
					FuncStats -> MinCost = getBestCost();
				}

				Iters[ Index ] = i;
				FuncStats -> ComplAttempts++;
				FuncStats -> SumItCompl += i;
				SumStats -> ComplAttempts++;
				SumStats -> SumIt += i;
				SumStats -> SumItAll += i;
				SumStats -> SumItImpr += ImprIters;
				SumStats -> SumItImprAll += ImprIters;
				SumStats -> SumIt_l10n += log( (double) i / Dims ) /
					log( 10.0 );

				#if OPT_STATS
				for( k = 0; k < getHistCount(); k++ )
				{
					int j;

					for( j = 0; j < 8; j++ )
					{
						SumStats -> SumSelsMap[ k * 8 + j ] +=
							SumSelsMap[ k * 8 + j ];
					}

					SumStats -> SumAvgSels[ k ] +=
						(double) SumSels[ k ] / i;

					SumStats -> SumDevSels[ k ] += DevSels[ k ];
				}
				#endif // OPT_STATS

				break;
			}

			#if defined( OPT_PLATEAU_MUL )
			if( sc > Dims * OPT_PLATEAU_MUL || i == MaxIters )
			#else // defined( OPT_PLATEAU_MUL )
			if( i == MaxIters )
			#endif // defined( OPT_PLATEAU_MUL )
			{
				#if OPT_THREADS
					VOXSYNC( StatsSync );

					#if OPT_TIME
					SumStats -> TimeOpt += optt;
					SumStats -> TimeFunc += funct;
					#endif // OPT_TIME
				#endif // OPT_THREADS

				if( getBestCost() < FuncStats -> MinCost )
				{
					FuncStats -> MinCost = getBestCost();
				}

				Iters[ Index ] = -1;
				FuncStats -> SumRjCost += getBestCost();
				SumStats -> SumItAll += i;
				SumStats -> SumItImprAll += ImprIters;
				SumStats -> SumRjCost += getBestCost() - optv;

				break;
			}
		}
	}

#if OPT_THREADS
protected:
	virtual ecode performWork()
	{
		performOpt();

		VOXRET;
	}
#endif // OPT_THREADS
};

/**
 * Function test corpus class.
 */

class CTester
{
public:
	int DefDims; ///< The default number of dimensions to use.
		///<
	CSumStats SumStats; ///< Summary statistics.
		///<
	double SuccessAt; ///< Average Success attempts, available after run().
		///<
	double AvgIt; ///< Average Iters, available after run().
		///<
	double AvgIt_l10n; ///< Average Iters/ln(10), available after run().
		///<
	double AvgRMS; ///< Average RMS, available after run().
		///<
	double AvgRMS_l10n; ///< Average RMS/ln(10), available after run().
		///<
	double AvgRjCost; ///< Average reject cost, available after run().
		///<

	CTester()
		: RMCacheCount( 0 )
	{
	}

	/**
	 * Function initalizes the tester. Then the addCorpus() function should be
	 * called at least once.
	 *
	 * @param Threshold Objective function value threshold - stop condition.
	 * @param IterCount The number of attempts to solve a function to perform.
	 * @param InnerIterCount The maximal number of solver iterations to
	 * perform.
	 * @param DoPrint "True" if results should be printed to "stdout".
	 */

	void init( const double aThreshold, const int aIterCount,
		const int aInnerIterCount, const bool aDoPrint )
	{
		CostThreshold = aThreshold;
		IterCount = aIterCount;
		InnerIterCount = aInnerIterCount;
		DoPrint = aDoPrint;
		FnCount = 0;
	}

	/**
	 * Function adds a function corpus to the tester.
	 *
 	 * @param aDefDims The number of dimensions in each function with variable
	 * number of dimensions.
	 * @param Corpus NULL-limited list of test functions.
	 * @param DoRandomize "True" if randomization of parameter shifts and
	 * rotations should be performed.
	 * @param DoRandomizeAll "True" if all randomizations should be applied,
	 * otherwise only shifts and scales will be applied.
	 */

	void addCorpus( const int aDefDims, const CTestFn** Corpus,
		const bool aDoRandomize, const bool aDoRandomizeAll )
	{
		while( FnCount < MaxFuncs )
		{
			if( *Corpus == NULL )
			{
				break;
			}

			Funcs[ FnCount ] = *Corpus;
			FuncData[ FnCount ].DefDims = aDefDims;
			FuncData[ FnCount ].DoRandomize = aDoRandomize;
			FuncData[ FnCount ].DoRandomizeAll = aDoRandomizeAll;
			FnCount++;
			Corpus++;
		}
	}

	/**
	 * Function runs the test. On return, SumStats object will be updated.
	 */

	void run()
	{
		#if defined( OPT_PERF )
		LARGE_INTEGER Freq;
		QueryPerformanceFrequency( &Freq );
		LARGE_INTEGER t1;
		QueryPerformanceCounter( &t1 );
		#endif // defined( OPT_PERF )

		SumStats.clear();
		SumStats.TotalAttempts = FnCount * IterCount;

		CTestOpt* opt;
		int* Iters = new int[ IterCount ];
		int k;

		#if OPT_THREADS
		CWorkerThreadPool Threads;

		for( k = 0; k < CSystem :: getProcessorCount(); k++ )
		{
			Threads.add( new CTestOpt() );
		}
		#else // OPT_THREADS
			opt = new CTestOpt();
		#endif // OPT_THREADS

		for( k = 0; k < FnCount; k++ )
		{
			const CTestFnData* const fndata = &FuncData[ k ];
			const int Dims = ( Funcs[ k ] -> Dims == 0 ?
				fndata -> DefDims : Funcs[ k ] -> Dims );

			CFuncStats FuncStats;
			int j;

			for( j = 0; j < IterCount; j++ )
			{
				#if OPT_THREADS
				VOXERRSKIP( Threads.getIdleThread( opt ));
				#endif // OPT_THREADS

				opt -> rnd.init( k + j * 10000 );
				opt -> updateDims( Dims );
				opt -> fn = Funcs[ k ];
				opt -> DoRandomize = fndata -> DoRandomize;
				opt -> DoRandomizeAll = fndata -> DoRandomizeAll;
				opt -> MaxIters = InnerIterCount;
				opt -> SumStats = &SumStats;
				opt -> FuncStats = &FuncStats;
				opt -> Index = j;
				opt -> Iters = Iters;
				opt -> optv = opt -> fn -> OptValue;

				int i;

				if( opt -> fn -> ParamFunc != NULL )
				{
					(*opt -> fn -> ParamFunc)( opt -> minv, opt -> maxv, Dims,
						&opt -> optv );
				}
				else
				{
					for( i = 0; i < Dims; i++ )
					{
						opt -> minv[ i ] = opt -> fn -> RangeMin;
						opt -> maxv[ i ] = opt -> fn -> RangeMax;
					}
				}

				if( CostThreshold > 0.0 )
				{
					opt -> CostThreshold = opt -> optv + CostThreshold;
				}
				else
				{
					const double d = fabs( opt -> optv ) * -CostThreshold;

					opt -> CostThreshold =
						opt -> optv + ( d < 0.01 ? 0.01 : d );
				}

				if( fndata -> DoRandomize )
				{
					if( fndata -> DoRandomizeAll )
					{
						getRotationMatrix( Dims, opt -> rots, opt -> rnd );
					}

					for( i = 0; i < Dims; i++ )
					{
						const double d =
							( opt -> maxv[ i ] - opt -> minv[ i ]) * 0.5;

						opt -> shifts[ i ] = d *
							( opt -> rnd.get() - 0.5 ) * 2.0;

						opt -> signs[ i ] = 1.0 +
							opt -> rnd.get() * 0.5;

						opt -> minv[ i ] -= d * 2.5;
						opt -> maxv[ i ] += d * 2.5;
					}
				}

				#if OPT_THREADS
				opt -> start();
				#else // OPT_THREADS
				opt -> performOpt();
				#endif // OPT_THREADS
			}

			#if OPT_THREADS
			VOXERRSKIP( Threads.waitAllForFinish() );
			#endif // OPT_THREADS

			const double MinCost = FuncStats.MinCost;
			const double AvgRjCost0 = FuncStats.SumRjCost /
				( IterCount - FuncStats.ComplAttempts );

			double Avg = 0.0;
			double RMS = 0.0;

			if( FuncStats.ComplAttempts > 0 )
			{
				SumStats.ComplFuncs++;
				Avg = (double) FuncStats.SumItCompl / FuncStats.ComplAttempts;

				RMS = 0.0;
				int RMSCount = 0;

				for( j = 0; j < IterCount; j++ )
				{
					if( Iters[ j ] >= 0 )
					{
						const double v = Iters[ j ] - Avg;
						RMS += v * v;
						RMSCount++;
					}
				}

				if( RMSCount > 0 && RMS > 0.0 )
				{
					RMS = sqrt( RMS / RMSCount );
					SumStats.SumRMS += RMS;
					SumStats.SumRMS_l10n += log( RMS / Dims ) / log( 10.0 );
					SumStats.SumRMSCount++;
				}

				#if defined( EVALBINS )
				for( j = 0; j < IterCount; j++ )
				{
					if( Iters[ j ] >= 0 )
					{
						const double v = ( Iters[ j ] - Avg ) / RMS /
							SumStats.RMSSpan;

						int z = (int) (( v + 1.0 ) * SumStats.binspan + 0.5 );

						if( z < 0 )
						{
							z = 0;
						}
						else
						if( z >= SumStats.binc )
						{
							z = SumStats.binc - 1;
						}

						SumStats.bins[ z ]++;
					}
				}
				#endif // defined( EVALBINS )
			}

			const double At = 1.0 / ( (double) FuncStats.ComplAttempts /
				IterCount );

			if( DoPrint )
			{
				printf( "I:%6.0f R:%5.0f A:%5.2f C:%11.5f RC:%11.5f %s_%i%c\n",
					Avg, RMS, At, MinCost, AvgRjCost0, opt -> fn -> Name,
					Dims, ( fndata -> DoRandomize ? 'r' : ' ' ));
			}

			if( _kbhit() && _getch() == 27 )
			{
				break;
			}
		}

		#if defined( EVALBINS )
		for( k = 0; k < SumStats.binc; k++ )
		{
			printf( "%2.2f\t%2.2f\n", SumStats.RMSSpan *
				( k - SumStats.binspan ) / SumStats.binspan,
				100.0 * SumStats.bins[ k ] / SumStats.ComplAttempts );
		}
		#endif // defined( EVALBINS )

		const double SuccessFn = 100.0 * SumStats.ComplFuncs / FnCount;
		SuccessAt = 100.0 * SumStats.ComplAttempts / SumStats.TotalAttempts;
		AvgIt = (double) SumStats.SumIt / SumStats.ComplAttempts;
		AvgIt_l10n = SumStats.SumIt_l10n / SumStats.ComplAttempts;
		AvgRMS = SumStats.SumRMS / SumStats.SumRMSCount;
		AvgRMS_l10n = SumStats.SumRMS_l10n / SumStats.SumRMSCount;
		AvgRjCost = SumStats.SumRjCost /
			( SumStats.TotalAttempts - SumStats.ComplAttempts );

		if( DoPrint )
		{
			#if defined( EVALBINS )
			printf( ">=%.0f-sigma: %.2f%%\n", SumStats.RMSSpan,
				100.0 * SumStats.bins[ SumStats.binc - 1 ] /
				SumStats.ComplAttempts );
			#endif // defined( EVALBINS )

			printf( "Func count: %i, Attempts: %i, MaxIters/Attempt: %i\n",
				FnCount, IterCount, InnerIterCount );

			printf( "GoodItersAvg: %.2f%% (percent of improving iterations "
				"in all attempts)\n",
				100.0 * SumStats.SumItImprAll / SumStats.SumItAll );

			printf( "Func success: %.2f%%\n", SuccessFn );
			printf( "Attempts success: %.2f%%\n", SuccessAt );
			printf( "AvgIt: %.3f (avg iterations across all successful "
				"attempts)\n", AvgIt );

			printf( "AvgIt_l10n: %.3f (avg log10(it/N) across all successful "
				"attempts)\n", AvgIt_l10n );

			printf( "AvgRMS: %.1f (avg std.dev of convergence time)\n",
				AvgRMS );

			printf( "AvgRMS_l10n: %.1f (avg log10(std.dev/N) of convergence time)\n",
				AvgRMS_l10n );

			printf( "AvgRjCost: %.6f (avg unbiased reject cost)\n",
				AvgRjCost );

			#if OPT_TIME
			printf( "TimeOpt: %.3f (time in optimize)\n", SumStats.TimeOpt );
			printf( "TimeFunc: %.3f (time in objfunc)\n", SumStats.TimeFunc );
			printf( "Overhead: %.1f%% (optimize overhead)\n", 100.0 - 100.0 *
				SumStats.TimeFunc / SumStats.TimeOpt );
			#endif // OPT_TIME

			#if defined( OPT_PERF )
			LARGE_INTEGER t2;
			QueryPerformanceCounter( &t2 );
			printf("time: %.3f s\n", ( t2.QuadPart - t1.QuadPart ) /
				(double) Freq.QuadPart );
			#endif // defined( OPT_PERF )

			#if OPT_STATS
			CBiteOptHistBase** const h = opt -> getHists();
			const char** const hnames = opt -> getHistNames();

			printf( "\nmin\tmax\tbias\tcount\t"
				"sel0\tsel1\tsel2\tsel3\thist\n" );

			for( k = 0; k < opt -> getHistCount(); k++ )
			{
				const int cc = h[ k ] -> getChoiceCount();
				const double avg = SumStats.SumAvgSels[ k ] /
					SumStats.ComplAttempts;

				const double std = sqrt( SumStats.SumDevSels[ k ] /
					SumStats.ComplAttempts );

				const double v0 = ( avg - std ) / ( cc - 1 );
				const double v1 = ( avg + std ) / ( cc - 1 );

				printf( "%6.3f\t%6.3f\t%6.3f\t%i",
					v0, v1, avg / ( cc - 1 ), cc );

				int j;

				for( j = 0; j < 4; j++ )
				{
					if( SumStats.SumSelsMap[ k * 8 + j ] > 0 )
					{
						printf( "\t%-5.1f", 100.0 *
							SumStats.SumSelsMap[ k * 8 + j ] /
							SumStats.SumIt );
					}
					else
					{
						printf( "\t     " );
					}
				}

				printf( "\t%s\n", hnames[ k ]);
			}
			#endif // OPT_THREADS
		}

		delete[] Iters;

		#if !OPT_THREADS
		delete opt;
		#endif // !OPT_THREADS
	}

protected:
	static const int MaxFuncs = 1000; ///< The maximal number of functions
		///< possible to add to the tester.
		///<

	/**
	 * Structure holds auxilliary function data.
	 */

	struct CTestFnData
	{
		int DefDims; ///< The number of dimensions to use if function
			///< supports any number of dimensions.
			///<
		bool DoRandomize; ///< Perform parameter shift randomization.
			///<
		bool DoRandomizeAll; ///< Use all randomization techniques, "false"
			///< if only shift and scale.
			///<
	};

	const CTestFn* Funcs[ MaxFuncs ]; ///< Test functions corpus.
		///<
	CTestFnData FuncData[ MaxFuncs ]; ///< Test function aux data.
		///<
	int FnCount; ///< Test function count.
		///<
	double CostThreshold; ///< Cost threshold (finish criteria).
		///<
	int IterCount; ///< Iteration count.
		///<
	int InnerIterCount; ///< Inner iteration count (the number of
		///< optimization calls).
	bool DoPrint; ///< Print results to stdout.
		///<

	/**
	 * Structure that holds cached random rotation matrices for a specified
	 * dimensionality.
	 */

	struct CRotMatCacheItem
	{
		static const int EntryCount = 2048; ///< The number of cache entries.
			///<
		int Dims; ///< Dimension count.
			///<
		double* rmbuf[ EntryCount ]; ///< Rotation matrix buffer pointers.
			///<
		double** rm[ EntryCount ]; ///< Row-wise rotation matrices pointers.
			///<

		CRotMatCacheItem()
			: Dims( 0 )
		{
			memset( rmbuf, 0, sizeof( rmbuf ));
			memset( rm, 0, sizeof( rm ));
		}

		~CRotMatCacheItem()
		{
			int i;

			for( i = 0; i < EntryCount; i++ )
			{
				delete[] rmbuf[ i ];
				delete[] rm[ i ];
			}
		}
	};

	static const int RMCacheSize = 8; ///< The maximal number of different
		///< dimensionalities supported by cache in the same run.
		///<
	int RMCacheCount; ///< The number of cache entires actually used.
		///<
	CRotMatCacheItem RMCache[ RMCacheSize ]; ///< Rotation matrix cache for
		///< different dimensionalities.
		///<

	/**
	 * Function creates random rotation matrix. Code from BBOB competition.
	 */

	void makeRotationMatrix( double** B, const int _DIM, CBiteRnd& rrnd )
	{
		double prod;
		int i, j, k;

		for( i = 0; i < _DIM; i++ )
		{
			for( j = 0; j < _DIM; j++ )
			{
				double unif = rrnd.get();

				if( unif == 0.0 )
				{
					unif = 1e-99;
				}

				double unif2 = rrnd.get();

				if( unif2 == 0.0 )
				{
					unif2 = 1e-99;
				}

				B[i][j] = sqrt(-2.0*log(unif)) * cos(2.0*M_PI*unif2);

				if( B[i][j] == 0 )
				{
					B[i][j] = 1e-99;
				}
			}
		}

		for( i = 0; i < _DIM; i++ )
		{
			for( j = 0; j < i; j++ )
			{
				prod = 0.0;

				for( k = 0; k < _DIM; k++ )
				{
					prod += B[k][i] * B[k][j];
				}

				for( k = 0; k < _DIM; k++ )
				{
					B[k][i] -= prod * B[k][j];
				}
			}

			prod = 0.0;

			for( k = 0; k < _DIM; k++ )
			{
				prod += B[k][i] * B[k][i];
			}

			prod = sqrt( prod );

			for( k = 0; k < _DIM; k++ )
			{
				B[k][i] /= prod;
			}
		}
	}

	/**
	 * Function creates or returns a cached random rotation matrix.
	 *
	 * @param _DIM Dimensionality.
	 * @param[out] B Row-wise matrix pointers.
	 * @param rrnd Random number generator.
	 */

	void getRotationMatrix( const int _DIM, double**& B, CBiteRnd& rrnd )
	{
		CRotMatCacheItem* Cache = NULL;
		int i;

		for( i = 0; i < RMCacheCount; i++ )
		{
			if( RMCache[ i ].Dims == _DIM )
			{
				Cache = &RMCache[ i ];
				break;
			}
		}

		if( Cache == NULL )
		{
			if( RMCacheCount == RMCacheSize )
			{
				printf( "rmcache overflow\n" );
				Cache = &RMCache[ 0 ];
			}
			else
			{
				Cache = &RMCache[ RMCacheCount ];
				Cache -> Dims = _DIM;
				RMCacheCount++;
			}
		}

		const int ri = rrnd.getInt( CRotMatCacheItem :: EntryCount );

		if( Cache -> rm[ ri ] == NULL )
		{
			Cache -> rm[ ri ] = new double*[ _DIM ];
			Cache -> rmbuf[ ri ] = new double[ _DIM * _DIM ];

			for( i = 0; i < _DIM; i++ )
			{
				Cache -> rm[ ri ][ i ] = Cache -> rmbuf[ ri ] + i * _DIM;
			}

			makeRotationMatrix( Cache -> rm[ ri ], _DIM, rrnd );
		}

		B = Cache -> rm[ ri ];
	}
};
