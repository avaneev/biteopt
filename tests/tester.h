//$ nocpp

#include <stdio.h>
#include "../biteopt.h"
#include "../spheropt.h"
//#include "../smaesopt.h"
//#include "../other/nmpopt.h"
//#include "../other/ccmaes.h"

CBiteRnd rnd;

#include "testfn.h"

#define OPT_CLASS CBiteOpt//CBiteOptDeep//CSMAESOpt//CSpherOpt//CNelderMeadPlusOpt//CCMAESOpt//
#define OPT_DIMS_PARAMS Dims
//#define EVALBINS 1

/**
 * Function test corpus class.
 */

class CTester
{
public:
	/**
	 * Function optimizer class.
	 */

	class CTestOpt : public OPT_CLASS
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
	};

	int DefDims; ///< The default number of dimensions to use.
		///<
	CTestOpt* opt; ///< Optimizer.
		///<
	double ItAvg; ///< Average of convergence time after run() across
		///< functions.
		///<
	double ItAvg_l10n; ///< = log( AtAvg / N ) / log( 10 ).
		///<
	double ItAvg2; ///< Average of convergence time after run() across all
		///< attempts.
		///<
	double ItAvg2_l10n; ///< = avg( log( At / N ) / log( 10 )).
		///<
	double RMSAvg; ///< Std.dev of convergence time after run().
		///<
	double RMSAvg_l10n; ///< = log( RMSAvg / N ) / log( 10 ).
		///<
	double ItRtAvg; ///< Average ratio of std.dev and average after run().
		///<
	double RjAvg; ///< Average number of rejects.
		///<
	double AtAvg; ///< Average number of attempts.
		///<
	double CostAvg; ///< Average achieved cost among all functions, including
		///< in successful and rejected attempts, zero-based.
		///<
	double Success; ///< Success rate.

	CTester()
		: opt( new CTestOpt() )
		, RMCacheCount( 0 )
	{
	}

	~CTester()
	{
		delete opt;
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
	 * Function runs the test. On return, the ItAvg, RMSAvg, ItRtAvg, RjAvg,
	 * AtAvg class variables will be updated.
	 */

	void run()
	{
		ItAvg = 0.0;
		ItAvg_l10n = 0.0;
		ItAvg2 = 0.0;
		int ItAvg2Count = 0;
		ItAvg2_l10n = 0.0;
		RMSAvg = 0.0;
		RMSAvg_l10n = 0.0;
		ItRtAvg = 0.0;
		int AvgCount = 0;
		RjAvg = 0.0;
		AtAvg = 0.0;
		CostAvg = 0.0;
		double RejTotal = 0.0;
		int ComplTotal = 0;
		int* Iters = new int[ IterCount ];
		int k;
		double GoodIters = 0.0;
		int GoodItersCount = 0;

		#if defined( EVALBINS )
		const int binspan = 50;
		const int binc = 1 + binspan * 2;
		double bins[ binc ];
		const double RMSSpan = 3.0;

		for( k = 0; k < binc; k++ )
		{
			bins[ k ] = 0.0;
		}
		#endif // defined( EVALBINS )

		for( k = 0; k < FnCount; k++ )
		{
			const CTestFnData* const fndata = &FuncData[ k ];
			double AvgIter = 0.0;
			double MinCost = 1e300; // Minimal cost detected in successes.
			double MinRjCost = 1e300; // Minimal cost detected in rejects.
			double MinCost2 = 1e300; // Minimal cost detected in successes.
			double AvgRjCost = 0.0; // Average cost detected in rejects.
			int Rej = 0; // The number of rejected attempts.
			int j;
			int i;

			opt -> fn = Funcs[ k ];
			const int Dims = ( opt -> fn -> Dims == 0 ?
				fndata -> DefDims : opt -> fn -> Dims );

			opt -> updateDims( Dims );
			opt -> DoRandomize = fndata -> DoRandomize;
			opt -> DoRandomizeAll = fndata -> DoRandomizeAll;

			for( j = 0; j < IterCount; j++ )
			{
				rnd.init( k + j * 10000 );

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
						getRotationMatrix( Dims, opt -> rots, rnd );
					}

					for( i = 0; i < Dims; i++ )
					{
						const double d =
							( opt -> maxv[ i ] - opt -> minv[ i ]) * 0.5;

						opt -> shifts[ i ] = d *
							( rnd.getRndValue() - 0.5 ) * 2.0;

						opt -> minv[ i ] -= d * 2.5;
						opt -> maxv[ i ] += d * 2.5;

						opt -> signs[ i ] = 1.0 + rnd.getRndValue() * 0.5;
					}
				}

				opt -> init( rnd );
				i = 0;
				int impriters = 0;

				while( true )
				{
					if( opt -> optimize( rnd ) == 0 )
					{
						impriters++;
					}

					i++;

					if( opt -> getBestCost() - opt -> optv < CostThreshold )
					{
						if( opt -> getBestCost() < MinCost )
						{
							MinCost = opt -> getBestCost();
						}

						if( opt -> getBestCost() - opt -> optv < MinCost2 )
						{
							MinCost2 = opt -> getBestCost() - opt -> optv;
						}

						GoodIters += (double) impriters / i;
						GoodItersCount++;
						ComplTotal++;
						const int itc = i + opt -> getInitEvals();
						Iters[ j ] = itc;
						AvgIter += itc;
						ItAvg2 += itc;
						ItAvg2_l10n += log( (double) itc / Dims ) /
							log( 10.0 );

						ItAvg2Count++;
						break;
					}

					if( i == InnerIterCount )
					{
						if( opt -> getBestCost() < MinRjCost )
						{
							MinRjCost = opt -> getBestCost();
						}

						AvgRjCost += opt -> getBestCost() - opt -> optv;
						GoodIters += (double) impriters / i;
						GoodItersCount++;
						Rej++;
						Iters[ j ] = -1;
						break;
					}
				}
			}

			MinCost = ( Rej >= IterCount ?
				1.0 / ( Rej - IterCount ) : MinCost );

			MinRjCost = ( Rej == 0 ? 1.0 / Rej : MinRjCost );
			double Avg;
			double RMS;

			if( Rej >= IterCount )
			{
				Avg = 0.0;
				RMS = 0.0;
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

				ItAvg += Avg;
				ItAvg_l10n += log( Avg / Dims ) / log( 10.0 );
				RMSAvg += RMS;
				RMSAvg_l10n += ( RMS == 0.0 ? 10.0 :
					log( RMS / Dims ) / log( 10.0 ));

				ItRtAvg += RMS / Avg;
				AvgCount++;

				#if defined( EVALBINS )
				for( j = 0; j < IterCount; j++ )
				{
					if( Iters[ j ] >= 0 )
					{
						const double v = ( Iters[ j ] - Avg ) / RMS / RMSSpan;
						int z = (int) (( v + 1.0 ) * binspan + 0.5 );

						if( z < 0 )
						{
							z = 0;
						}
						else
						if( z >= binc )
						{
							z = binc - 1;
						}

						bins[ z ]++;
					}
				}
				#endif // defined( EVALBINS )
			}

			const double Rj = (double) Rej / IterCount;
			RjAvg += Rj;
			const double At = 1.0 / ( 1.0 - (double) Rej / IterCount );
			AvgRjCost = ( Rej > 0 ? AvgRjCost / Rej : 0.0 );
			CostAvg += ( Rej >= IterCount ? AvgRjCost :
				( Rej == 0 ? MinCost2 : ( AvgRjCost * Rej +
				MinCost2 * ( IterCount - Rej )) / IterCount ));

			RejTotal += Rej;

			if( DoPrint )
			{
				printf( "AI:%6.0f RI:%5.0f At:%5.2f C:%13.8f RjC:%7.4f "
					"%s_%i%c\n", Avg, RMS, At, MinCost,
					MinRjCost, opt -> fn -> Name, Dims,
					( fndata -> DoRandomize ? 'r' : ' ' ));
//				printf( "C:%20.13f %20.13f %s_%i\n",
//					MinRjCost, opt -> optv, opt -> fn -> Name, Dims );
			}
		}

		#if defined( EVALBINS )
		for( k = 0; k < binc; k++ )
		{
			printf( "%2.2f\t%2.2f\n", RMSSpan * ( k - binspan ) / binspan,
				bins[ k ] / ComplTotal * 100.0 );
		}
		#endif // defined( EVALBINS )

		ItAvg /= AvgCount;
		ItAvg_l10n /= AvgCount;
		RMSAvg /= AvgCount;
		RMSAvg_l10n /= AvgCount;
		ItRtAvg /= AvgCount;
		ItAvg2 /= ItAvg2Count;
		ItAvg2_l10n /= ItAvg2Count;
		RjAvg /= FnCount;
		AtAvg = 1.0 / ( 1.0 - (double) RejTotal / IterCount / FnCount );
		CostAvg /= FnCount;
		Success = 100.0 * ComplTotal / FnCount / IterCount;

		if( DoPrint )
		{
			#if defined( EVALBINS )
			printf( ">=%.0f-sigma: %.2f%%\n", RMSSpan,
				bins[ binc - 1 ] / ComplTotal * 100.0 );
			#endif // defined( EVALBINS )

			printf( "GoodItersAvg: %.2f%% (avg percentage of improving "
				"iterations in all attempts)\n",
				GoodIters / GoodItersCount * 100.0 );

			printf( "Success: %.2f%%\n", Success );
			printf( "ItAvg: %.1f (avg convergence time)\n", ItAvg );
			printf( "ItAvg_l10n: %.3f (avg log10(it/N))\n", ItAvg_l10n );
			printf( "ItAvg2: %.1f (avg convergence time across all "
				"successful attempts)\n", ItAvg2 );

			printf( "ItAvg2_l10n: %.3f (avg log10(it/N) across all successful "
				"attempts)\n", ItAvg2_l10n );

			printf( "RMSAvg: %.1f (avg std.dev of convergence time)\n",
				RMSAvg );

			printf( "ItRtAvg: %.6f (avg ratio of std.dev and average)\n",
				ItRtAvg );

			printf( "RjAvg: %.2f%% (avg percentage of rejects)\n",
				RjAvg * 100.0 );

			printf( "AtAvg: %.3f (avg number of attempts)\n", AtAvg );
			printf( "CostAvg: %.6f (avg cost)\n", CostAvg );
		}

		delete[] Iters;
	}

protected:
	static const int MaxFuncs = 500; ///< The maximal number of functions
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
		///< dimensionalities supported by cache.
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

				double unif2 = rnd.getRndValue();

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
