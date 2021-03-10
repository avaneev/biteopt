//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBiteOpt and CBiteOptDeep classes.
 *
 * @section license License
 *
 * Copyright (c) 2016-2021 Aleksey Vaneev
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * @version 2021.14
 */

#ifndef BITEOPT_INCLUDED
#define BITEOPT_INCLUDED

#include "spheropt.h"

/**
 * BiteOpt optimization class. Implements a stochastic non-linear
 * bound-constrained derivative-free optimization method.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CBiteOpt : public CBiteOptBase< int64_t >
{
public:
	typedef int64_t ptype; ///< Parameter value storage type.
		///<

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: HistCount( 0 )
	{
		setParPopCount( 4 );

		addHist( ParPopPHist );
		addHist( ParPopHist[ 0 ]);
		addHist( ParPopHist[ 1 ]);
		addHist( ParPopHist[ 2 ]);
		addHist( PopChangeHist );
		addHist( ParPopUpdHist );
		addHist( ScutHist );
		addHist( MethodHist );
		addHist( DrawHist );
		addHist( M1Hist );
		addHist( M1SHist[ 0 ]);
		addHist( M1SHist[ 1 ]);
		addHist( MinSolPwrHist[ 0 ]);
		addHist( MinSolPwrHist[ 1 ]);
		addHist( MinSolPwrHist[ 2 ]);
		addHist( MinSolMulHist[ 0 ]);
		addHist( MinSolMulHist[ 1 ]);
		addHist( MinSolMulHist[ 2 ]);
		addHist( Gen1AllpHist );
		addHist( Gen1MoveHist );
		addHist( Gen1MoveAsyncHist );
		addHist( Gen1MoveDEHist );
		addHist( Gen1MoveSpanHist );
		addHist( Gen4RedFacHist );
		addHist( Gen4MixFacHist );
		addHist( Gen5BinvHist );
		addHist( *ParOpt.getHists()[ 0 ]);
		addHist( *ParOpt.getHists()[ 1 ]);
		addHist( *ParOpt.getHists()[ 2 ]);
		addHist( *ParOpt.getHists()[ 3 ]);
	}

	static const int MaxHistCount = 64; ///< The maximal number of histograms
		///< that can be added (for static arrays).
		///<

	/**
	 * Function returns a pointer to an array of histograms in use.
	 */

	CBiteOptHistBase** getHists()
	{
		return( Hists );
	}

	/**
	 * Function returns a pointer to an array of histogram names.
	 */

	static const char** getHistNames()
	{
		static const char* HistNames[] = {
			"ParPopPHist",
			"ParPopHist[ 0 ]",
			"ParPopHist[ 1 ]",
			"ParPopHist[ 2 ]",
			"PopChangeHist",
			"ParPopUpdHist",
			"ScutHist",
			"MethodHist",
			"DrawHist",
			"M1Hist",
			"M1SHist[ 0 ]",
			"M1SHist[ 1 ]",
			"MinSolPwrHist[ 0 ]",
			"MinSolPwrHist[ 1 ]",
			"MinSolPwrHist[ 2 ]",
			"MinSolMulHist[ 0 ]",
			"MinSolMulHist[ 1 ]",
			"MinSolMulHist[ 2 ]",
			"Gen1AllpHist",
			"Gen1MoveHist",
			"Gen1MoveAsyncHist",
			"Gen1MoveDEHist",
			"Gen1MoveSpanHist",
			"Gen4RedFacHist",
			"Gen4MixFacHist",
			"Gen5BinvHist",
			"ParOpt.CentPowHist",
			"ParOpt.RadPowHist",
			"ParOpt.EvalFacHist",
			"ParOpt.PopChangeHist"
		};

		return( HistNames );
	}

	/**
	 * Function returns the number of histograms in use.
	 */

	int getHistCount() const
	{
		return( HistCount );
	}

	/**
	 * Function updates dimensionality of *this object. Function does nothing
	 * if dimensionality has not changed since the last call. This function
	 * should be called at least once before calling the init() function.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0 or negative, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 :
			13 + aParamCount * 4 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		ParOpt.Owner = this;
		ParOpt.updateDims( aParamCount );
	}

	/**
	 * Function initializes *this optimizer. Does not perform objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid for initial population vectors.
	 * @param InitRadius Initial radius, multiplier relative to the default
	 * sigma value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();
		updateDiffValues( true );

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		const double sd = 0.25 * InitRadius;
		int i;
		int j;

		if( InitParams == NULL )
		{
			for( j = 0; j < PopSize; j++ )
			{
				ptype* const p = PopParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					p[ i ] = wrapParamInt( rnd, getGaussianInt(
						rnd, sd, IntMantMult >> 1 ));
				}
			}
		}
		else
		{
			ptype* const p0 = PopParams[ 0 ];

			for( i = 0; i < ParamCount; i++ )
			{
				p0[ i ] = wrapParamInt( rnd,
					(ptype) (( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ]));
			}

			for( j = 1; j < PopSize; j++ )
			{
				ptype* const p = PopParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					p[ i ] = wrapParamInt( rnd, getGaussianInt(
						rnd, sd, p0[ i ]));
				}
			}
		}

		updateCentroid();

		ParamCountRnd = ParamCount * rnd.getRawScaleInv();
		ParamCntr = (int) ( rnd.getUniformRaw() * ParamCountRnd );
		AllpProbDamp = (int) ( CBiteRnd :: getRawScale() * 1.8 / ParamCount );
		PrevSelMethod = MethodHist.selectRandom( rnd );
		CentUpdateCtr = 0;

		for( i = 0; i < HistCount; i++ )
		{
			Hists[ i ] -> reset( rnd );
		}

		ParOpt.init( rnd, InitParams, InitRadius );

		DoInitEvals = true;
		InitEvalIndex = 0;
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @param PushOpt Optimizer where the recently obtained solution should be
	 * "pushed", used for deep optimization algorithm.
	 * @return The number of non-improving iterations so far. A high value
	 * means optimizer has reached an optimization plateau. The suggested
	 * threshold value is ParamCount * 64. When this value was reached, the
	 * probability of plateau is high. This value, however, should not be
	 * solely relied upon when considering a stopping criteria: a hard
	 * iteration limit should be always used as in some cases convergence time
	 * may be very high with small, but frequent improving steps. This value
	 * is best used to allocate iteration budget between optimization attempts
	 * more efficiently.
	 */

	int optimize( CBiteRnd& rnd, CBiteOpt* const PushOpt = NULL )
	{
		int i;

		if( DoInitEvals )
		{
			const ptype* const p = PopParams[ InitEvalIndex ];

			for( i = 0; i < ParamCount; i++ )
			{
				NewValues[ i ] = getRealValue( p, i );
			}

			const double NewCost = optcost( NewValues );
			sortPop( NewCost, InitEvalIndex );
			updateBestCost( NewCost, NewValues );

			InitEvalIndex++;

			if( InitEvalIndex == PopSize )
			{
				for( i = 0; i < ParPopCount; i++ )
				{
					ParPops[ i ] -> copy( *this );
				}

				DoInitEvals = false;
			}

			return( 0 );
		}

		bool DoEval = true;
		double NewCost;
		int SelMethod;

		static const int ScutProbLim =
			(int) ( 0.1 * CBiteRnd :: getRawScale() ); // Short-cut
			// probability limit, in raw scale.

		bool DoScut = false;

		if( rnd.getUniformRaw() < ScutProbLim )
		{
			if( select( ScutHist, rnd ))
			{
				DoScut = true;
			}
		}

		if( DoScut )
		{
			SelMethod = -1;

			// Parameter value short-cuts, they considerably reduce
			// convergence time for some functions while not severely
			// impacting performance for other functions.

			i = (int) ( rnd.getUniformRaw() * ParamCountRnd );

			const double r = rnd.getRndValue();
			const double r2 = r * r;

			// "Same-value parameter vector" short-cut.

			const int si = (int) ( r2 * r2 * CurPopSize );
			const ptype* const rp = getParamsOrdered( si );

			const double v = getRealValue( rp, i );

			for( i = 0; i < ParamCount; i++ )
			{
				TmpParams[ i ] = (ptype) (( v - MinValues[ i ]) /
					DiffValues[ i ]);
			}
		}
		else
		{
			const int SelDraw = select( DrawHist, rnd );

			if( SelDraw == 0 )
			{
				SelMethod = selectForce( MethodHist, PrevSelMethod );
			}
			else
			if( SelDraw == 1 )
			{
				SelMethod = select( MethodHist, rnd );
			}
			else
			{
				SelMethod = selectRandom( MethodHist, rnd );
			}

			if( SelMethod == 0 )
			{
				generateSol2( rnd );
			}
			else
			if( SelMethod == 1 )
			{
				if( select( M1Hist, rnd ))
				{
					if( select( M1SHist[ 0 ], rnd ))
					{
						generateSol3( rnd );
					}
					else
					{
						generateSol2b( rnd );
					}
				}
				else
				{
					if( select( M1SHist[ 1 ], rnd ))
					{
						generateSol4( rnd );
					}
					else
					{
						generateSol5( rnd );
					}
				}
			}
			else
			if( SelMethod == 2 )
			{
				generateSol1( rnd );
			}
			else
			{
				ParOpt.optimize( rnd, &NewCost, NewValues );
				DoEval = false;

				for( i = 0; i < ParamCount; i++ )
				{
					TmpParams[ i ] = (ptype) (( NewValues[ i ] -
						MinValues[ i ]) / DiffValues[ i ]);
				}

				updateBestCost( NewCost, NewValues );
			}
		}

		// Evaluate objective function with new parameters.

		if( DoEval )
		{
			// Wrap parameter values so that they stay in the [0; 1] range.

			for( i = 0; i < ParamCount; i++ )
			{
				TmpParams[ i ] = wrapParamInt( rnd, TmpParams[ i ]);
				NewValues[ i ] = getRealValue( TmpParams, i );
			}

			NewCost = optcost( NewValues );

			updateBestCost( NewCost, NewValues );
		}

		if( !isAcceptedCost( NewCost ))
		{
			// Upper bound cost constraint check failed, reject this solution.

			applyHistsDecr();

			if( SelMethod >= 0 )
			{
				PrevSelMethod = MethodHist.selectRandom( rnd );
			}

			StallCount++;

			if( CurPopSize < PopSize )
			{
				if( select( PopChangeHist, rnd ))
				{
					// Increase population size on fail.

					const double r = rnd.getRndValue();

					incrCurPopSize( CurPopSize1 -
						(int) ( r * r * CurPopSize ));
				}
			}
		}
		else
		{
			applyHistsIncr();

			if( SelMethod >= 0 )
			{
				PrevSelMethod = SelMethod;
			}

			if( NewCost == PopCosts[ CurPopSize1 ])
			{
				StallCount++;
			}
			else
			{
				StallCount = 0;
			}

			updatePop( NewCost, TmpParams, false, false );

			if( PushOpt != NULL && PushOpt != this &&
				!PushOpt -> DoInitEvals )
			{
				PushOpt -> updatePop( NewCost, TmpParams, false, true );
				PushOpt -> updateParPop( NewCost, TmpParams );
			}

			if( CurPopSize > PopSize / 3 )
			{
				if( select( PopChangeHist, rnd ))
				{
					// Decrease population size on success.

					decrCurPopSize();
				}
			}
		}

		// Diverging populations technique.

		if( select( ParPopUpdHist, rnd ))
		{
			updateParPop( NewCost, TmpParams );
		}

		CentUpdateCtr++;

		if( CentUpdateCtr >= CurPopSize * 32 )
		{
			// Update centroids of parallel populations that use running
			// average to reduce error accumulation.

			CentUpdateCtr = 0;

			for( i = 0; i < ParPopCount; i++ )
			{
				ParPops[ i ] -> updateCentroid();
			}
		}

		return( StallCount );
	}

protected:
	CBiteOptHistBase* Hists[ MaxHistCount ]; ///< Pointers to histogram
		///< objects, for indexed access in some cases.
		///<
	int HistCount; ///< The number of histograms in use.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	int AllpProbDamp; ///< Damped Allp probability, in raw PRNG scale. Applied
		///< for higher dimensions as the "all parameter" randomization is
		///< ineffective in the higher dimensions.
		///<
	CBiteOptHist< 2 > ParPopPHist; ///< Parallel population use probability
		///< histogram.
		///<
	CBiteOptHist< 4 > ParPopHist[ 3 ]; ///< Parallel population
		///< histograms for solution generators (template's Count parameter
		///< should match ParPopCount).
		///<
	CBiteOptHist< 2 > PopChangeHist; ///< Population size change
		///< histogram.
		///<
	CBiteOptHist< 2 > ParPopUpdHist; ///< Parallel population update
		///< histogram.
		///<
	CBiteOptHistBinary ScutHist; ///< Short-cut method's histogram.
		///<
	CBiteOptHist< 4 > MethodHist; ///< Population generator method
		///< histogram.
		///<
	CBiteOptHist< 3 > DrawHist; ///< Method draw histogram.
		///<
	CBiteOptHist< 2 > M1Hist; ///< Method 1's sub-method histogram.
		///<
	CBiteOptHist< 2 > M1SHist[ 2 ]; ///< Method 1's sub-sub-method histograms.
		///<
	CBiteOptHist< 4 > MinSolPwrHist[ 3 ]; ///< Index of least-cost
		///< population, power factor.
		///<
	CBiteOptHist< 4 > MinSolMulHist[ 3 ]; ///< Index of least-cost
		///< population, multiplier.
		///<
	CBiteOptHistBinary Gen1AllpHist; ///< Generator method 1's Allp
		///< histogram.
		///<
	CBiteOptHistBinary Gen1MoveHist; ///< Generator method 1's Move
		///< histogram.
		////<
	CBiteOptHist< 2 > Gen1MoveAsyncHist; ///< Generator method 1's Move
		///< async histogram.
		///<
	CBiteOptHist< 2 > Gen1MoveDEHist; ///< Generator method 1's Move DE
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen1MoveSpanHist; ///< Generator method 1's Move span
		///< histogram.
		///<
	CBiteOptHist< 3 > Gen4RedFacHist; ///< Generator method 4's RedFac
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen4MixFacHist; ///< Generator method 4's mixing
		///< count histogram.
		///<
	CBiteOptHistBinary Gen5BinvHist; ///< Generator method 5's random
		///< inversion technique histogram.
		///<
	int PrevSelMethod; ///< Previously successfully used method; contains
		///< random method index if optimization was not successful.
		///<
	int CentUpdateCtr; ///< Centroid update counter.
		///<
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
		///<
	int InitEvalIndex; ///< Current initial population index.
		///<

	/**
	 * Parallel optimizer class.
	 */

	class CParOpt : public CSpherOpt
	{
	public:
		CBiteOpt* Owner; ///< Owner object.

		virtual void getMinValues( double* const p ) const
		{
			Owner -> getMinValues( p );
		}

		virtual void getMaxValues( double* const p ) const
		{
			Owner -> getMaxValues( p );
		}

		virtual double optcost( const double* const p )
		{
			return( Owner -> optcost( p ));
		}
	};

	CParOpt ParOpt; ///< Parallel optimizer.
		///<

	/**
	 * Function adds a histogram to the Hists list.
	 *
	 * @param h Histogram object to add.
	 */

	void addHist( CBiteOptHistBase& h )
	{
		Hists[ HistCount ] = &h;
		HistCount++;
	}

	/**
	 * Function updates a selected parallel population.
	 *
	 * @param NewCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 */

	void updateParPop( const double NewCost, const ptype* const UpdParams )
	{
		const int p = getMaxDistParPop( NewCost, UpdParams );

		if( p >= 0 )
		{
			ParPops[ p ] -> updatePop( NewCost, UpdParams, true, true );
		}
	}

	/**
	 * Function selects a parallel population to use for solution generation.
	 * With certain probability, *this object's own population will be
	 * returned instead of parallel population.
	 *
	 * @param gi Solution generator index (0-2).
	 * @param rnd PRNG object.
	 */

	CBiteOptPop& selectParPop( const int gi, CBiteRnd& rnd )
	{
		if( select( ParPopPHist, rnd ))
		{
			return( *ParPops[ select( ParPopHist[ gi ], rnd )]);
		}

		return( *this );
	}

	/**
	 * Function returns a dynamically-selected minimal population index, used
	 * in some solution generation methods.
	 *
	 * @param gi Solution generator index (0-2).
	 * @param rnd PRNG object.
	 * @param ps Population size.
	 */

	int getMinSolIndex( const int gi, CBiteRnd& rnd, const int ps )
	{
		static const double pp[ 4 ] = { 0.05, 0.125, 0.25, 0.5 };
		const double r = ps * pow( rnd.getRndValue(),
			ps * pp[ select( MinSolPwrHist[ gi ], rnd )]);

		static const double rm[ 4 ] = { 0.0, 0.125, 0.25, 0.5 };

		return( (int) ( r * rm[ select( MinSolMulHist[ gi ], rnd )]));
	}

	/**
	 * The original "bitmask inversion" solution generator. Most of the time
	 * adjusts only a single parameter of the very best solution, yet manages
	 * to produce excellent "reference points".
	 */

	void generateSol1( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const CBiteOptPop& ParPop = selectParPop( 0, rnd );

		memcpy( Params, ParPop.getParamsOrdered(
			getMinSolIndex( 0, rnd, ParPop.getCurPopSize() )),
			ParamCount * sizeof( Params[ 0 ]));

		// Select a single random parameter or all parameters for further
		// operations.

		int i;
		int a;
		int b;
		bool DoAllp = false;

		if( rnd.getUniformRaw() < AllpProbDamp )
		{
			if( select( Gen1AllpHist, rnd ))
			{
				DoAllp = true;
			}
		}

		if( DoAllp )
		{
			a = 0;
			b = ParamCount - 1;
		}
		else
		{
			a = ParamCntr;
			b = ParamCntr;
			ParamCntr = ( ParamCntr == 0 ? ParamCount : ParamCntr ) - 1;
		}

		// Bitmask inversion operation, works as a "driver" of optimization
		// process.

		const double r1 = rnd.getRndValue();
		const double r12 = r1 * r1;
		const int ims = (int) ( r12 * r12 * 48.0 );
		const int64_t imask = ( ims > 63 ? 0 : IntMantMask >> ims );

		const double r2 = rnd.getRndValue();
		const int im2s = (int) ( r2 * r2 * 96.0 );
		const int64_t imask2 = ( im2s > 63 ? 0 : IntMantMask >> im2s );

		const int si1 = (int) ( r1 * r12 * CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		for( i = a; i <= b; i++ )
		{
			Params[ i ] = (( Params[ i ] ^ imask ) +
				( rp1[ i ] ^ imask2 )) >> 1;
		}

		if( select( Gen1MoveHist, rnd ))
		{
			const double r3 = rnd.getRndValue();
			const int si2 = (int) ( r3 * r3 * CurPopSize );
			const ptype* const rp2 = getParamsOrdered( si2 );

			if( select( Gen1MoveDEHist, rnd ))
			{
				// Apply a DE-based move.

				const ptype* const rp3 = getParamsOrdered(
					CurPopSize1 - si2 );

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] -= (( rp3[ i ] - rp2[ i ]) >> 1 );
				}
			}
			else
			{
				if( select( Gen1MoveAsyncHist, rnd ))
				{
					a = 0;
					b = ParamCount - 1;
				}

				// Random move around random previous solution vector.

				static const double SpanMults[ 4 ] = {
					0.5 * CBiteRnd :: getRawScaleInv(),
					1.5 * CBiteRnd :: getRawScaleInv(),
					2.0 * CBiteRnd :: getRawScaleInv(),
					2.5 * CBiteRnd :: getRawScaleInv()
				};

				const double m = SpanMults[ select( Gen1MoveSpanHist, rnd )];
				const double m1 = rnd.getTPDFRaw() * m;
				const double m2 = rnd.getTPDFRaw() * m;

				for( i = a; i <= b; i++ )
				{
					Params[ i ] -= (ptype) (( Params[ i ] - rp2[ i ]) * m1 );
					Params[ i ] -= (ptype) (( Params[ i ] - rp2[ i ]) * m2 );
				}
			}
		}
	}

	/**
	 * The "Digital Evolution"-based solution generator.
	 */

	void generateSol2( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int si1 = getMinSolIndex( 1, rnd, CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );
		const ptype* const rp3 = getParamsOrdered( CurPopSize1 - si1 );

		const double r2 = rnd.getRndValue();
		const int si2 = 1 + (int) ( r2 * CurPopSize1 );
		const ptype* const rp2 = getParamsOrdered( si2 );

		const double r4 = rnd.getRndValue();
		const int si4 = (int) ( r4 * r4 * CurPopSize );
		const ptype* const rp4 = getParamsOrdered( si4 );
		const ptype* const rp5 = getParamsOrdered( CurPopSize1 - si4 );

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// The "step in the right direction" (Differential Evolution
			// "mutation") operation towards the best (minimal) and away from
			// the worst (maximal) parameter vector, plus a difference of two
			// random vectors.

			Params[ i ] = rp1[ i ] - ((( rp3[ i ] - rp2[ i ]) +
				( rp5[ i ] - rp4[ i ])) >> 1 );
		}
	}

	/**
	 * An alternative "Digital Evolution"-based solution generator.
	 */

	void generateSol2b( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		// Select worst and a random previous solution from the ordered list,
		// apply offsets to reduce sensitivity to noise.

		const int si1 = getMinSolIndex( 1, rnd, CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		const double r2 = rnd.getRndValue();
		const int si2 = (int) ( r2 * r2 * CurPopSize );
		const ptype* const rp2 = getParamsOrdered( si2 );

		const ptype* const rp3 = getParamsOrdered( CurPopSize1 - si2 );

		// Select two more previous solutions to be used in the mix.

		const double r4 = rnd.getRndValue();
		const int si4 = (int) ( r4 * r4 * CurPopSize );
		const ptype* const rp4 = getParamsOrdered( si4 );
		const ptype* const rp5 = getParamsOrdered( CurPopSize1 - si4 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rp1[ i ] - ((( rp3[ i ] - rp2[ i ]) +
				( rp5[ i ] - rp4[ i ])) >> 1 );
		}
	}

	/**
	 * Alternative randomized solution generator, works well for convex
	 * functions. Uses the very best solution and a random previous solution.
	 * "mp * mp" is equivalent of giving more weight to better solutions.
	 */

	void generateSol3( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		if( NeedCentUpdate )
		{
			updateCentroid();
		}

		const ptype* const MinParams = getParamsOrdered(
			getMinSolIndex( 2, rnd, CurPopSize ));

		const ptype* const cp = getCentroid();

		const double r1 = rnd.getRndValue();
		const int si1 = (int) ( r1 * r1 * CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const int64_t m1 = rnd.getBit();
			const int64_t m2 = 1 - m1;

			Params[ i ] = cp[ i ] * m1 +
				( MinParams[ i ] + ( MinParams[ i ] - rp1[ i ])) * m2;
		}
	}

	/**
	 * "Entropy bit mixing"-based solution generator. Performs crossing-over
	 * of an odd number (this is important) of random solutions via XOR
	 * operation. Slightly less effective than the DE-based mixing, but makes
	 * the optimization method more diverse overall.
	 */

	void generateSol4( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		CBiteOptPop& ParPop = selectParPop( 1, rnd );

		static const double RedFacs[ 3 ] = { 0.75, 1.0, 1.5 };

		int UseSize;
		const ptype** const UseParams = ParPop.getSparsePopParams(
			UseSize, RedFacs[ select( Gen4RedFacHist, rnd )]);

		const int km = 3 + ( select( Gen4MixFacHist, rnd ) << 1 );
		int k;
		int i;

		for( k = 0; k < km; k++ )
		{
			const double r1 = rnd.getRndValue();
			const int si1 = (int) ( r1 * r1 * UseSize );
			const ptype* const rp1 = UseParams[ si1 ];

			if( k == 0 )
			{
				memcpy( Params, rp1, ParamCount * sizeof( Params[ 0 ]));
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] ^= rp1[ i ];
				}
			}
		}
	}

	/**
	 * A novel "Randomized bit crossing-over" candidate solution generation
	 * method. Effective, but on its own cannot stand coordinate system
	 * offsets, converges slowly. Completely mixes bits of two
	 * randomly-selected solutions, plus changes 1 random bit.
	 *
	 * This method is fundamentally similar to biological DNA crossing-over.
	 */

	void generateSol5( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const CBiteOptPop& ParPop = selectParPop( 2, rnd );

		const double r1 = rnd.getRndValue();
		const ptype* const CrossParams1 = ParPop.getParamsOrdered(
			(int) ( r1 * r1 * ParPop.getCurPopSize() ));

		const double r2 = rnd.getRndValue();
		const ptype* const CrossParams2 = getParamsOrdered(
			(int) ( r2 * r2 * CurPopSize ));

		const bool UseInv = select( Gen5BinvHist, rnd );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit mixing mask.

			const int64_t crpl = rnd.getUniformRaw2();

			int64_t v1 = CrossParams1[ i ];
			int64_t v2 = ( UseInv ? ~CrossParams2[ i ] : CrossParams2[ i ]);

			if( rnd.getBit() )
			{
				const int b = (int) ( rnd.getRndValue() * IntMantBits );
				const int64_t m = ~( 1LL << b );
				const int64_t bv = (int64_t) rnd.getBit() << b;

				v1 &= m;
				v2 &= m;
				v1 |= bv;
				v2 |= bv;
			}

			Params[ i ] = ( v1 & crpl ) | ( v2 & ~crpl );
		}
	}
};

/**
 * Deep optimization class. Based on an array of M CBiteOpt objects. This
 * "deep" method pushes the newly-obtained solution to the next CBiteOpt
 * object which is then optimized.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CBiteOptDeep : public CBiteOptInterface
{
public:
	CBiteOptDeep()
		: ParamCount( 0 )
		, OptCount( 0 )
		, Opts( NULL )
	{
	}

	virtual ~CBiteOptDeep()
	{
		deleteBuffers();
	}

	virtual const double* getBestParams() const
	{
		return( BestOpt -> getBestParams() );
	}

	virtual double getBestCost() const
	{
		return( BestOpt -> getBestCost() );
	}

	/**
	 * Function returns a pointer to an array of histograms in use by the
	 * current CBiteOpt object.
	 */

	CBiteOptHistBase** getHists()
	{
		return( CurOpt -> getHists() );
	}

	/**
	 * Function returns a pointer to an array of histogram names.
	 */

	static const char** getHistNames()
	{
		return( CBiteOpt :: getHistNames() );
	}

	/**
	 * Function returns the number of histograms in use by the current
	 * CBiteOpt object.
	 */

	int getHistCount() const
	{
		return( CurOpt -> getHistCount() );
	}

	/**
	 * Function updates dimensionality of *this object. Function does nothing
	 * if dimensionality has not changed since the last call. This function
	 * should be called at least once before calling the init() function.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param M The number of CBiteOpt objects. This number depends on the
	 * complexity of the problem, if the default value does not produce a good
	 * solution, it should be increased together with the iteration count.
	 * Minimal value is 1, in this case a plain CBiteOpt optimization will be
	 * performed.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int M = 6,
		const int PopSize0 = 0 )
	{
		if( aParamCount == ParamCount && M == OptCount )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		OptCount = M;
		Opts = new CBiteOptWrap*[ OptCount ];

		int i;

		for( i = 0; i < OptCount; i++ )
		{
			Opts[ i ] = new CBiteOptWrap( this );
			Opts[ i ] -> updateDims( aParamCount, PopSize0 );
		}
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 * @param InitRadius Initial radius, relative to the default value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		int i;

		for( i = 0; i < OptCount; i++ )
		{
			Opts[ i ] -> init( rnd, InitParams, InitRadius );
		}

		BestOpt = Opts[ 0 ];
		CurOpt = Opts[ 0 ];
		StallCount = 0;
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return The number of non-improving iterations so far. The plateau
	 * threshold value is ParamCount * 64.
	 */

	int optimize( CBiteRnd& rnd )
	{
		if( OptCount == 1 )
		{
			StallCount = Opts[ 0 ] -> optimize( rnd );

			return( StallCount );
		}

		CBiteOptWrap* PushOpt;

		if( OptCount == 2 )
		{
			PushOpt = Opts[ CurOpt == Opts[ 0 ]];
		}
		else
		{
			while( true )
			{
				PushOpt = Opts[ (int) ( rnd.getRndValue() * OptCount )];

				if( PushOpt != CurOpt )
				{
					break;
				}
			}
		}

		const int sc = CurOpt -> optimize( rnd, PushOpt );

		if( CurOpt -> getBestCost() <= BestOpt -> getBestCost() )
		{
			BestOpt = CurOpt;
		}

		if( sc == 0 )
		{
			StallCount = 0;
		}
		else
		{
			StallCount++;
			CurOpt = PushOpt;
		}

		return( StallCount );
	}

protected:
	/**
	 * Wrapper class for CBiteOpt class.
	 */

	class CBiteOptWrap : public CBiteOpt
	{
	public:
		CBiteOptDeep* Owner; ///< Owner object.
			///<

		CBiteOptWrap( CBiteOptDeep* const aOwner )
			: Owner( aOwner )
		{
		}

		virtual void getMinValues( double* const p ) const
		{
			Owner -> getMinValues( p );
		}

		virtual void getMaxValues( double* const p ) const
		{
			Owner -> getMaxValues( p );
		}

		virtual double optcost( const double* const p )
		{
			return( Owner -> optcost( p ));
		}
	};

	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int OptCount; ///< The total number of optimization objects in use.
		///<
	CBiteOptWrap** Opts; ///< Optimization objects.
		///<
	CBiteOptWrap* BestOpt; ///< Optimizer that contains the best solution.
		///<
	CBiteOptWrap* CurOpt; ///< Current optimization object index.
		///<
	int StallCount; ///< The number of iterations without improvement.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	void deleteBuffers()
	{
		if( Opts != NULL )
		{
			int i;

			for( i = 0; i < OptCount; i++ )
			{
				delete Opts[ i ];
			}

			delete[] Opts;
			Opts = NULL;
		}
	}
};

/**
 * Objective function.
 */

typedef double (*biteopt_func)( int N, const double* x,
	void* func_data );

/**
 * Wrapper class for the biteopt_minimize() function.
 */

class CBiteOptMinimize : public CBiteOptDeep
{
public:
	int N; ///< The number of dimensions in objective function.
	biteopt_func f; ///< Objective function.
	void* data; ///< Objective function's data.
	const double* lb; ///< Parameter's lower bounds.
	const double* ub; ///< Parameter's lower bounds.

	virtual void getMinValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = lb[ i ];
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = ub[ i ];
		}
	}

	virtual double optcost( const double* const p )
	{
		return( (*f)( N, p, data ));
	}
};

/**
 * Function performs minimization using the CBiteOpt or CBiteOptDeep
 * algorithm.
 *
 * @param N The number of parameters in an objective function.
 * @param f Objective function.
 * @param data Objective function's data.
 * @param lb Lower bounds of obj function parameters, should not be infinite.
 * @param ub Upper bounds of obj function parameters, should not be infinite.
 * @param[out] x Minimizer.
 * @param[out] minf Minimizer's value.
 * @param iter The number of iterations to perform in a single attempt.
 * Corresponds to the number of obj function evaluations that are performed.
 * @param M Depth to use, 1 for plain CBiteOpt algorithm, >1 for CBiteOptDeep
 * algorithm. Expected range is [1; 36]. Internally multiplies "iter" by
 * sqrt(M). 
 * @param attc The number of optimization attempts to perform.
 */

inline void biteopt_minimize( const int N, biteopt_func f, void* data,
	const double* lb, const double* ub, double* x, double* minf,
	const int iter, const int M = 1, const int attc = 10 )
{
	CBiteOptMinimize opt;
	opt.N = N;
	opt.f = f;
	opt.data = data;
	opt.lb = lb;
	opt.ub = ub;
	opt.updateDims( N, M );

	CBiteRnd rnd;
	rnd.init( 1 );

	const int useiter = (int) ( iter * sqrt( (double) M ));
	int k;

	for( k = 0; k < attc; k++ )
	{
		opt.init( rnd );

		int i;

		for( i = 0; i < useiter; i++ )
		{
			opt.optimize( rnd );
		}

		if( k == 0 || opt.getBestCost() <= *minf )
		{
			memcpy( x, opt.getBestParams(), N * sizeof( x[ 0 ]));
			*minf = opt.getBestCost();
		}
	}
}

#endif // BITEOPT_INCLUDED
