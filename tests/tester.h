//$ nocpp

#include <stdio.h>
#include <emmintrin.h>
#include "../biteopt.h"

CBiteRnd rnd;

#include "testfn.h"

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
		double* rotsbuf; ///< rots buffer.
			///<
		double** rots; ///< Rotation matrix to apply to function parameters.
			///<
		double* shifts; ///< Shifts to apply to function parameters.
			///<
		double* signs; ///< Signs or scales to apply to function parameters.
			///<
		double* tp; ///< Temporary parameter storage.
			///<
		bool DoRandomize; ///< Apply value randomization.
			///<

		CTestOpt()
			: minv( NULL )
			, maxv( NULL )
			, rotsbuf( NULL )
			, rots( NULL )
			, shifts( NULL )
			, signs( NULL )
			, tp( NULL )
		{
		}

		~CTestOpt()
		{
			delete[] minv;
			delete[] maxv;
			delete[] rotsbuf;
			delete[] rots;
			delete[] shifts;
			delete[] signs;
			delete[] tp;
		}

		void updateDims( const int aDims )
		{
			Dims = aDims;
			delete[] minv;
			delete[] maxv;
			delete[] rotsbuf;
			delete[] rots;
			delete[] shifts;
			delete[] signs;
			delete[] tp;
			minv = new double[ Dims ];
			maxv = new double[ Dims ];
			rotsbuf = new double[ Dims * Dims ];
			rots = new double*[ Dims ];
			shifts = new double[ Dims ];
			signs = new double[ Dims ];
			tp = new double[ Dims ];

			int i;

			for( i = 0; i < Dims; i++ )
			{
				rots[ i ] = rotsbuf + i * Dims;
			}

			CBiteOpt :: updateDims( Dims );
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

/*			for( i = 0; i < Dims; i++ )
			{
				tp[ i ] = p[ i ] * signs[ i ] + shifts[ i ];
			}
*/
			for( i = 0; i < Dims; i++ )
			{
				tp[ i ] = 0.0;
				int j;

				for( j = 0; j < Dims; j++ )
				{
					tp[ i ] += rots[ i ][ j ] *
						( p[ i ] * signs[ i ] + shifts[ i ]);
				}
			}

			return( (*fn -> CalcFunc)( tp, Dims ));
		}
	};

	int DefDims; ///< The default number of dimensions to use.
		///<
	CTestOpt* opt; ///< Optimizer.
		///<
	double ItAvg; ///< Average of convergence time after run().
		///<
	double RMSAvg; ///< Std.dev of convergence time after run().
		///<
	double ItRtAvg; ///< Average ratio of std.dev and average after run().
		///<
	double RjAvg; ///< Average number of rejects.
		///<
	double AtAvg; ///< Average number of attempts.
		///<
	double Score; ///< Optimization score.
		///<
	double Success; ///< Success rate.

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
	 * @param DoRandomize "True" if randomization of value shifts should be
	 * performed.
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
	 * Function runs the test. On return, the ItAvg, RMSAvg, ItRtAvg, RjAvg,
	 * AtAvg class variables will be updated.
	 */

	void run()
	{
		ItAvg = 0.0;
		RMSAvg = 0.0;
		ItRtAvg = 0.0;
		RjAvg = 0.0;
		AtAvg = 0.0;
		int RejTotal = 0;
		int ComplTotal = 0;
		int* Iters = new int[ IterCount ];
		int k;
		double GoodIters = 0.0;
		int GoodItersCount = 0;

		const int binspan = 50;
		const int binc = 1 + binspan * 2;
		double bins[ binc ];
		const double RMSSpan = 3.0;

		for( k = 0; k < binc; k++ )
		{
			bins[ k ] = 0.0;
		}

		for( k = 0; k < FnCount; k++ )
		{
			_mm_empty();
			double AvgIter = 0;
			double MinCost = 1e300; // Minimal cost detected in successes.
			double MinRjCost = 1e300; // Minimal cost detected in rejects.
			int Rej = 0; // The number of rejected attempts.
			int j;
			int i;

			opt -> fn = Funcs[ k ];
			const int Dims = ( opt -> fn -> Dims == 0 ?
				DefDims : opt -> fn -> Dims );

			opt -> updateDims( Dims );
			opt -> DoRandomize = DoRandomize;

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
					makeRotationMatrix( opt -> rots, Dims, rnd );

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

				_mm_empty();
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

						GoodIters += (double) impriters / i;
						GoodItersCount++;
						ComplTotal++;
						Iters[ j ] = i + opt -> getInitEvals();
						AvgIter += i + opt -> getInitEvals();
						break;
					}

					if( i == InnerIterCount )
					{
						if( opt -> getBestCost() < MinRjCost )
						{
							MinRjCost = opt -> getBestCost();
						}

						GoodIters += (double) impriters / i;
						GoodItersCount++;
						Rej++;
						Iters[ j ] = -1;
						break;
					}
				}
			}

			MinCost = ( Rej >= IterCount ?
				1.0 / ( Rej - IterCount) : MinCost );

			MinRjCost = ( Rej == 0 ? 1.0 / Rej : MinRjCost );
			double Avg;
			double RMS;

			if( Rej >= IterCount )
			{
				Avg = 1.0;
				RMS = 1.0;
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
			}

			ItAvg += Avg;
			RMSAvg += RMS;
			ItRtAvg += RMS / Avg;
			const double Rj = (double) Rej / IterCount;
			RjAvg += Rj;
			const double At = 1.0 / ( 1.0 - (double) Rej / IterCount );
			AtAvg += At;
			RejTotal += Rej;

			if( DoPrint )
			{
				printf( "AI:%6.0f RI:%5.0f At:%5.2f C:%13.10f RjC:%7.4f "
					"%s_%i\n", Avg, RMS, At, MinCost,
					MinRjCost, opt -> fn -> Name, Dims );
			}
		}

/*		for( k = 0; k < binc; k++ )
		{
			printf( "%2.2f\t", bins[ k ] / ComplTotal * 100.0 );
		}

		printf( "\n" );

		for( k = 0; k < binc; k++ )
		{
			printf( "%2.2f\t", RMSSpan * ( k - binspan ) / binspan );
		}

		printf( "\n" );
*/
		_mm_empty();
		ItAvg /= FnCount;
		RMSAvg /= FnCount;
		ItRtAvg /= FnCount;
		RjAvg /= FnCount;
		AtAvg = 1.0 / ( 1.0 - (double) RejTotal / IterCount / FnCount );
		Score = ( AtAvg - 1.0 ) * 100.0 +
			fabs( ItAvg - 334.0 ) * 0.1;
		Success = 100.0 * ComplTotal / FnCount / IterCount;

//		Score = -GoodIters / GoodItersCount * 100.0 /
//			(( AtAvg - 1.0 ) * 100.0 ) / ItAvg;

		if( DoPrint )
		{
			printf( ">=%.0f-sigma: %.2f%%\n", RMSSpan,
				bins[ binc - 1 ] / ComplTotal * 100.0 );

			printf( "GoodItersAvg: %.2f%% (avg percentage of improving "
				"iterations in all attempts)\n",
				GoodIters / GoodItersCount * 100.0 );

			printf( "Success: %.2f%%\n", Success );
			printf( "ItAvg: %.1f (avg convergence time)\n", ItAvg );
			printf( "RMSAvg: %.1f (avg std.dev of convergence time)\n",
				RMSAvg );

			printf( "ItRtAvg: %.6f (avg ratio of std.dev and average)\n",
				ItRtAvg );

			printf( "RjAvg: %.2f%% (avg percentage of rejects)\n",
				RjAvg * 100.0 );

			printf( "AtAvg: %.3f (avg number of attempts)\n", AtAvg );
			printf( "Score: %f\n", Score );
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
};
