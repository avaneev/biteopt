//$ nocpp

#include <conio.h>
#include <stdio.h>
#include "../biteopt.h"
#include "../spheropt.h"
#include "../smaesopt.h"
//#include "../other/nmpopt.h"
//#include "../other/ccmaes.h"

#include "testfn.h"

#define OPT_CLASS CBiteOpt//CBiteOptDeep//CSpherOpt//CSMAESOpt//CNelderMeadPlusOpt//CCMAESOpt//
#define OPT_DIMS_PARAMS Dims
//#define EVALBINS 1

#if 0
	#define OPT_THREADS 1
	#include "/libvox/Sources/Core/CWorkerThreadPool.h"
	using namespace vox;
	CSyncObject StatsSync;
#else // OPT_THREADS
	#define OPT_THREADS 0
#endif // OPT_THREADS

/**
 * Summary optimization statistics class.
 */

class CSumStats
{
public:
	double SumIt_l10n; ///< = sum( log( It ) / log( 10 )) in completed
		///< attempts.
		///<
	double SumRjCost; ///< Summary costs detected in all attempts (successful
		///< attempts include cost prior to success).
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
	int SumIters; ///< The overall number of performed function evaluations
		///< in successful attempts.
		///<
	int SumImprIters; ///< Sum of improving iterations across sucessful
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

	void clear()
	{
		SumIt_l10n = 0.0;
		SumRjCost = 0.0;
		SumRMS_l10n = 0.0;
		SumRMSCount = 0;
		ComplAttempts = 0;
		SumIters = 0;
		SumImprIters = 0;
		ComplFuncs = 0;

		#if defined( EVALBINS )

		RMSSpan = 3.0;

		int k;

		for( k = 0; k < binc; k++ )
		{
			bins[ k ] = 0.0;
		}

		#endif // defined( EVALBINS )
	}
};

/**
 * Function optimization statistics class, across all attempts.
 */

class CFuncStats
{
public:
	double MinCost; ///< Minimal cost detected in successes.
		///<
	double MinRjCost; ///< Minimal cost detected in rejects.
		///<
	int ComplAttempts; ///< The number of completed attempts.
		///<
	int SumComplIters; ///< Sum of iterations in completed attempts.
		///<

	CFuncStats()
	{
		clear();
	}

	void clear()
	{
		MinCost = 1e300;
		MinRjCost = 1e300;
		ComplAttempts = 0;
		SumComplIters = 0;
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

	CTestOpt()
		: minv( NULL )
		, maxv( NULL )
		, shifts( NULL )
		, signs( NULL )
		, tp( NULL )
		, tp2( NULL )
	{
	}

	virtual ~CTestOpt()
	{
		delete[] minv;
		delete[] maxv;
		delete[] shifts;
		delete[] signs;
		delete[] tp;
		delete[] tp2;
	}

	void updateDims( const int aDims )
	{
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
		if( !DoRandomize )
		{
			return( (*fn -> CalcFunc)( p, Dims ));
		}

		int i;

		for( i = 0; i < Dims; i++ )
		{
			tp2[ i ] = p[ i ] * signs[ i ] + shifts[ i ];
		}

		if( !DoRandomizeAll )
		{
			return( (*fn -> CalcFunc)( tp2, Dims ));
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

		return( (*fn -> CalcFunc)( tp, Dims ));
	}

	void performOpt()
	{
		init( rnd );

		int i = 0;
		int ImprIters = 0;
		double PrevCost = 0.0;

		while( true )
		{
			if( optimize( rnd ) == 0 )
			{
				ImprIters++;
			}

			i++;

			if( getBestCost() - optv < CostThreshold )
			{
				#if OPT_THREADS
					VOXSYNC( StatsSync );
				#endif // OPT_THREADS

				if( getBestCost() < FuncStats -> MinCost )
				{
					FuncStats -> MinCost = getBestCost();
				}

				Iters[ Index ] = i;
				FuncStats -> ComplAttempts++;
				FuncStats -> SumComplIters += i;
				SumStats -> ComplAttempts++;
				SumStats -> SumIters += i;
				SumStats -> SumImprIters += ImprIters;
				SumStats -> SumIt_l10n += log( (double) i / Dims ) /
					log( 10.0 );

				SumStats -> SumRjCost += PrevCost;

				break;
			}

			PrevCost = getBestCost() - optv;

			if( i == MaxIters )
			{
				#if OPT_THREADS
					VOXSYNC( StatsSync );
				#endif // OPT_THREADS

				if( getBestCost() < FuncStats -> MinRjCost )
				{
					FuncStats -> MinRjCost = getBestCost();
				}

				Iters[ Index ] = -1;
				SumStats -> SumRjCost += PrevCost;

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
	double AvgRMS; ///< Average RMS, available after run().
		///<
	double AvgIt; ///< Average Iters, available after run().
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
			int i;
			int j;

			for( j = 0; j < IterCount; j++ )
			{
				#if OPT_THREADS
				VOXERRSKIP( Threads.getIdleThread( opt ));
				#endif // OPT_THREADS

				opt -> updateDims( Dims );
				opt -> fn = Funcs[ k ];
				opt -> CostThreshold = CostThreshold;
				opt -> DoRandomize = fndata -> DoRandomize;
				opt -> DoRandomizeAll = fndata -> DoRandomizeAll;
				opt -> MaxIters = InnerIterCount;
				opt -> SumStats = &SumStats;
				opt -> FuncStats = &FuncStats;
				opt -> Index = j;
				opt -> Iters = Iters;
				opt -> rnd.init( k + j * 10000 );

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
							( opt -> rnd.getRndValue() - 0.5 ) * 2.0;

						opt -> minv[ i ] -= d * 2.5;
						opt -> maxv[ i ] += d * 2.5;

						opt -> signs[ i ] = 1.0 +
						opt -> rnd.getRndValue() * 0.5;
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

			const double MinCost = ( FuncStats.ComplAttempts == 0 ?
				1.0 / FuncStats.ComplAttempts : FuncStats.MinCost );

			const double MinRjCost = ( FuncStats.ComplAttempts == IterCount ?
				1.0 / ( IterCount - FuncStats.ComplAttempts ) :
				FuncStats.MinRjCost );

			double Avg = 0.0;
			double RMS = 0.0;

			if( FuncStats.ComplAttempts > 0 )
			{
				SumStats.ComplFuncs++;
				Avg = (double) FuncStats.SumComplIters /
					FuncStats.ComplAttempts;

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
				printf( "AI:%6.0f RI:%5.0f At:%5.2f C:%13.8f RjC:%7.4f "
					"%s_%i%c\n", Avg, RMS, At, MinCost, MinRjCost,
					opt -> fn -> Name, Dims,
					( fndata -> DoRandomize ? 'r' : ' ' ));
//				printf( "C:%20.13f %20.13f %s_%i\n",
//					MinRjCost, opt -> optv, opt -> fn -> Name, Dims );
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

		SuccessAt = 100.0 * SumStats.ComplAttempts / SumStats.TotalAttempts;
		AvgRMS = SumStats.SumRMS_l10n / SumStats.SumRMSCount;
		AvgIt = SumStats.SumIt_l10n / SumStats.ComplAttempts;
		AvgRjCost = SumStats.SumRjCost * FnCount / SumStats.TotalAttempts;

		if( DoPrint )
		{
			#if defined( EVALBINS )
			printf( ">=%.0f-sigma: %.2f%%\n", SumStats.RMSSpan,
				100.0 * SumStats.bins[ SumStats.binc - 1 ] /
				SumStats.ComplAttempts );
			#endif // defined( EVALBINS )

			printf( "Func count: %i, Attempts: %i, MaxIters/Attempt: %i\n",
				FnCount, IterCount, InnerIterCount );

			printf( "GoodItersAvg: %.2f%% (percent of improving "
				"iterations in successful attempts)\n", 100.0 *
				SumStats.SumImprIters / SumStats.SumIters );

			printf( "Attempts success: %.2f%%\n", SuccessAt );
			printf( "AvgRMS_l10n: %.1f (avg log10(std.dev/N) of convergence time)\n",
				AvgRMS );

			printf( "AvgIt_l10n: %.3f (avg log10(it/N) across all successful "
				"attempts)\n", AvgIt );

			printf( "AvgRjCost: %.6f (avg reject cost)\n", AvgRjCost );
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
				double unif = rrnd.getRndValue();

				if( unif == 0.0 )
				{
					unif = 1e-99;
				}

				double unif2 = rrnd.getRndValue();

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

		const int ri = (int) ( rrnd.getRndValue() *
			CRotMatCacheItem :: EntryCount );

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
