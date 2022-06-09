//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBiteOpt and CBiteOptDeep classes.
 *
 * @section license License
 *
 * Copyright (c) 2016-2022 Aleksey Vaneev
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
 */

#ifndef BITEOPT_INCLUDED
#define BITEOPT_INCLUDED

#define BITEOPT_VERSION "2022.17"

#include "spheropt.h"
#include "nmsopt.h"

/**
 * BiteOpt optimization class. Implements a stochastic non-linear
 * bound-constrained derivative-free optimization method.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CBiteOpt : public CBiteOptBase< int64_t >
{
public:
	typedef int64_t ptype; ///< Parameter value storage type (should be a
		///< signed integer type, same as CBiteOptBase template parameter).
		///<

	CBiteOpt()
	{
		setParPopCount( 4 );

		addHist( MethodHist, "MethodHist" );
		addHist( M1Hist, "M1Hist" );
		addHist( M1AHist, "M1AHist" );
		addHist( M1BHist, "M1BHist" );
		addHist( M1BAHist, "M1BAHist" );
		addHist( M1BBHist, "M1BBHist" );
		addHist( M2Hist, "M2Hist" );
		addHist( M2BHist, "M2BHist" );
		addHist( PopChangeIncrHist, "PopChangeIncrHist" );
		addHist( PopChangeDecrHist, "PopChangeDecrHist" );
		addHist( ParOpt2Hist, "ParOpt2Hist" );
		addHist( ParPopPHist[ 0 ], "ParPopPHist[ 0 ]" );
		addHist( ParPopPHist[ 1 ], "ParPopPHist[ 1 ]" );
		addHist( ParPopPHist[ 2 ], "ParPopPHist[ 2 ]" );
		addHist( ParPopPHist[ 3 ], "ParPopPHist[ 3 ]" );
		addHist( AltPopPHist, "AltPopPHist" );
		addHist( AltPopHist[ 0 ], "AltPopHist[ 0 ]" );
		addHist( AltPopHist[ 1 ], "AltPopHist[ 1 ]" );
		addHist( AltPopHist[ 2 ], "AltPopHist[ 2 ]" );
		addHist( AltPopHist[ 3 ], "AltPopHist[ 3 ]" );
		addHist( MinSolPwrHist[ 0 ], "MinSolPwrHist[ 0 ]" );
		addHist( MinSolPwrHist[ 1 ], "MinSolPwrHist[ 1 ]" );
		addHist( MinSolPwrHist[ 2 ], "MinSolPwrHist[ 2 ]" );
		addHist( MinSolPwrHist[ 3 ], "MinSolPwrHist[ 3 ]" );
		addHist( MinSolMulHist[ 0 ], "MinSolMulHist[ 0 ]" );
		addHist( MinSolMulHist[ 1 ], "MinSolMulHist[ 1 ]" );
		addHist( MinSolMulHist[ 2 ], "MinSolMulHist[ 2 ]" );
		addHist( MinSolMulHist[ 3 ], "MinSolMulHist[ 3 ]" );
		addHist( Gen1AllpHist, "Gen1AllpHist" );
		addHist( Gen1MoveHist, "Gen1MoveHist" );
		addHist( Gen1MoveAsyncHist, "Gen1MoveAsyncHist" );
		addHist( Gen1MoveSpanHist, "Gen1MoveSpanHist" );
		addHist( Gen4MixFacHist, "Gen4MixFacHist" );
		addHist( Gen7PowFacHist, "Gen7PowFacHist" );
		addHist( Gen8NumHist, "Gen8NumHist" );
		addHist( Gen8SpanHist, "Gen8SpanHist" );
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
			10 + aParamCount * 3 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		ParOpt.Owner = this;
		ParOpt.updateDims( aParamCount, 11 + aPopSize / 3 );
		ParOptPop.initBuffers( aParamCount, aPopSize );

		ParOpt2.Owner = this;
		ParOpt2.updateDims( aParamCount, aPopSize * 4 / 3 );
		ParOpt2Pop.initBuffers( aParamCount, aPopSize );

		OldPop.initBuffers( aParamCount, aPopSize );
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

		resetCommonVars( rnd );

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
					p[ i ] = wrapParam( rnd, getGaussianInt(
						rnd, sd, IntMantMult >> 1 ));
				}
			}
		}
		else
		{
			ptype* const p0 = PopParams[ 0 ];

			for( i = 0; i < ParamCount; i++ )
			{
				p0[ i ] = wrapParam( rnd,
					(ptype) (( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ]));
			}

			for( j = 1; j < PopSize; j++ )
			{
				ptype* const p = PopParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					p[ i ] = wrapParam( rnd,
						getGaussianInt( rnd, sd, p0[ i ]));
				}
			}
		}

		updateCentroid();

		AllpProbDamp = 1.8 * ParamCountI;
		CentUpdateCtr = 0;

		ParOpt.init( rnd, InitParams, InitRadius );
		ParOpt2.init( rnd, InitParams, InitRadius );
		UseParOpt = 0;

		ParOptPop.resetCurPopPos();
		ParOpt2Pop.resetCurPopPos();
		OldPop.resetCurPopPos();

		DoInitEvals = true;
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
	 * threshold value is ParamCount * 128. When this value was reached, the
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
			const ptype* const p = PopParams[ CurPopPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				NewValues[ i ] = getRealValue( p, i );
			}

			const double NewCost = optcost( NewValues );
			sortPop( NewCost, CurPopPos );
			updateBestCost( NewCost, NewValues );

			CurPopPos++;

			if( CurPopPos == PopSize )
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
		const int SelMethod = select( MethodHist, rnd );

		if( SelMethod == 0 )
		{
			generateSol2( rnd );
		}
		else
		if( SelMethod == 1 )
		{
			if( select( M1Hist, rnd ))
			{
				if( select( M1AHist, rnd ))
				{
					generateSol2c( rnd );
				}
				else
				{
					generateSol2b( rnd );
				}
			}
			else
			{
				const int SelM1B = select( M1BHist, rnd );

				if( SelM1B == 0 )
				{
					if( select( M1BAHist, rnd ))
					{
						generateSol4( rnd );
					}
					else
					{
						generateSol5b( rnd );
					}
				}
				else
				if( SelM1B == 1 )
				{
					if( select( M1BBHist, rnd ))
					{
						generateSol5( rnd );
					}
					else
					{
						generateSol7( rnd );
					}
				}
				else
				{
					generateSol6( rnd );
				}
			}
		}
		else
		if( SelMethod == 2 )
		{
			if( select( M2Hist, rnd ))
			{
				generateSol1( rnd );
			}
			else
			{
				if( select( M2BHist, rnd ))
				{
					generateSol3( rnd );
				}
				else
				{
					generateSol8( rnd );
				}
			}
		}
		else
		{
			DoEval = false;
			CBiteOptPop* UpdPop;

			if( UseParOpt == 1 )
			{
				// Re-assign optimizer 2 based on comparison of its
				// efficiency with optimizer 1.

				UseParOpt = select( ParOpt2Hist, rnd );
			}

			if( UseParOpt == 0 )
			{
				const int sc = ParOpt.optimize( rnd, &NewCost, NewValues );

				if( sc > 0 )
				{
					UseParOpt = 1; // On stall, select optimizer 2.
				}

				if( sc > ParamCount * 64 )
				{
					ParOpt.init( rnd, getBestParams(), 0.5 );
					ParOptPop.resetCurPopPos();
				}

				UpdPop = &ParOptPop;
			}
			else
			{
				const int sc = ParOpt2.optimize( rnd, &NewCost, NewValues );

				if( sc > 0 )
				{
					UseParOpt = 0; // On stall, select optimizer 1.
				}

				if( sc > ParamCount * 16 )
				{
					ParOpt2.init( rnd, getBestParams(), 1.0 );
					ParOpt2Pop.resetCurPopPos();
				}

				UpdPop = &ParOpt2Pop;
			}

			for( i = 0; i < ParamCount; i++ )
			{
				TmpParams[ i ] = (ptype) (( NewValues[ i ] -
					MinValues[ i ]) * DiffValuesI[ i ]);
			}

			UpdPop -> updatePop( NewCost, TmpParams, false, true );
		}

		if( DoEval )
		{
			// Evaluate objective function with new parameters, if the
			// solution was not provided by the parallel optimizer.
			// Wrap parameter values so that they stay in the [0; 1] range.

			for( i = 0; i < ParamCount; i++ )
			{
				TmpParams[ i ] = wrapParam( rnd, TmpParams[ i ]);
				NewValues[ i ] = getRealValue( TmpParams, i );
			}

			NewCost = optcost( NewValues );
		}

		updateBestCost( NewCost, NewValues );

		if( !isAcceptedCost( NewCost ))
		{
			// Upper bound cost constraint check failed, reject this solution.

			applyHistsDecr( rnd );

			StallCount++;

			if( CurPopSize < PopSize )
			{
				if( select( PopChangeIncrHist, rnd ))
				{
					// Increase population size on fail.

					incrCurPopSize( CurPopSize1 -
						rnd.getSqrInt( CurPopSize ));
				}
			}
		}
		else
		{
			if( NewCost == PopCosts[ CurPopSize1 ])
			{
				StallCount++;
			}
			else
			{
				StallCount = 0;
			}

			if( rnd.get() < ParamCountI )
			{
				OldPop.updatePop( PopCosts[ CurPopSize1 ],
					PopParams[ CurPopSize1 ], false, true );
			}

			const int p = updatePop( NewCost, TmpParams, true, false );
			applyHistsIncr( rnd, 1.0 - p * CurPopSizeI );

			if( PushOpt != NULL && PushOpt != this &&
				!PushOpt -> DoInitEvals && NewCost > PopCosts[ 0 ])
			{
				PushOpt -> updatePop( NewCost, TmpParams, true, true );
				PushOpt -> updateParPop( NewCost, TmpParams );
			}

			if( CurPopSize > PopSize / 2 )
			{
				if( select( PopChangeDecrHist, rnd ))
				{
					// Decrease population size on success.

					decrCurPopSize();
				}
			}
		}

		// "Diverging populations" technique.

		updateParPop( NewCost, TmpParams );

		CentUpdateCtr++;

		if( CentUpdateCtr >= CurPopSize * 8 )
		{
			// Update centroids of populations that use running average, to
			// reduce error accumulation.

			CentUpdateCtr = 0;

			updateCentroid();

			for( i = 0; i < ParPopCount; i++ )
			{
				ParPops[ i ] -> updateCentroid();
			}
		}

		return( StallCount );
	}

protected:
	double AllpProbDamp; ///< Damped Allp probability. Applied for higher
		///< dimensions as the "all parameter" randomization is ineffective in
		///< the higher dimensions.
		///<
	CBiteOptHist< 4 > MethodHist; ///< Population generator 4-method
		///< histogram.
		///<
	CBiteOptHist< 2 > M1Hist; ///< Method 1's sub-method histogram.
		///<
	CBiteOptHist< 2 > M1AHist; ///< Method 1's sub-sub-method A histogram.
		///<
	CBiteOptHist< 3 > M1BHist; ///< Method 1's sub-sub-method B histogram.
		///<
	CBiteOptHist< 2 > M1BAHist; ///< Method 1's sub-sub-method A2 histogram.
		///<
	CBiteOptHist< 2 > M1BBHist; ///< Method 1's sub-sub-method B2 histogram.
		///<
	CBiteOptHist< 2 > M2Hist; ///< Method 2's sub-method histogram.
		///<
	CBiteOptHist< 2 > M2BHist; ///< Method 2's sub-sub-method B histogram.
		///<
	CBiteOptHist< 2 > PopChangeIncrHist; ///< Population size change increase
		///< histogram.
		///<
	CBiteOptHist< 2 > PopChangeDecrHist; ///< Population size change decrease
		///< histogram.
		///<
	CBiteOptHist< 2 > ParOpt2Hist; ///< Parallel optimizer 2 use
		///< histogram.
		///<
	CBiteOptHist< 2 > ParPopPHist[ 4 ]; ///< Parallel population use
		///< probability histogram.
		///<
	CBiteOptHist< 2 > AltPopPHist; ///< Alternative population use
		///< histogram.
		///<
	CBiteOptHist< 2 > AltPopHist[ 4 ]; ///< Alternative population type use
		///< histogram.
		///<
	CBiteOptHist< 4 > MinSolPwrHist[ 4 ]; ///< Index of least-cost
		///< population, power factor.
		///<
	CBiteOptHist< 4 > MinSolMulHist[ 4 ]; ///< Index of least-cost
		///< population, multiplier.
		///<
	CBiteOptHist< 2 > Gen1AllpHist; ///< Generator method 1's Allp
		///< histogram.
		///<
	CBiteOptHist< 2 > Gen1MoveHist; ///< Generator method 1's Move
		///< histogram.
		///<
	CBiteOptHist< 2 > Gen1MoveAsyncHist; ///< Generator method 1's Move
		///< async histogram.
		///<
	CBiteOptHist< 4 > Gen1MoveSpanHist; ///< Generator method 1's Move span
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen4MixFacHist; ///< Generator method 4's mixing
		///< count histogram.
		///<
	CBiteOptHist< 4 > Gen7PowFacHist; ///< Generator method 7's Power
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen8NumHist; ///< Generator method 8's NumSols
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen8SpanHist; ///< Generator method 8's random span
		///< histogram.
		///<
	int CentUpdateCtr; ///< Centroid update counter.
		///<
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
		///<
	CBiteOptPop OldPop; ///< Population of older solutions, updated
		///< probabilistically.
		///<

	/**
	 * Parallel optimizer class.
	 */

	template< class T >
	class CParOpt : public T
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

	CParOpt< CSpherOpt > ParOpt; ///< Parallel optimizer.
		///<
	CBiteOptPop ParOptPop; ///< Population of parallel optimizer's solutions.
		///< Includes only its solutions.
		///<
	CParOpt< CNMSeqOpt > ParOpt2; ///< Parallel optimizer2.
		///<
	CBiteOptPop ParOpt2Pop; ///< Population of parallel optimizer 2's
		///< solutions. Includes only its solutions.
		///<
	int UseParOpt; ///< Parallel optimizer currently being in use.
		///<

	/**
	 * Function updates an appropriate parallel population.
	 *
	 * @param NewCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 */

	void updateParPop( const double NewCost, const ptype* const UpdParams )
	{
		const int p = getMinDistParPop( NewCost, UpdParams );

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
	 * @param gi Solution generator index (0-3).
	 * @param rnd PRNG object.
	 */

	CBiteOptPop& selectParPop( const int gi, CBiteRnd& rnd )
	{
		if( select( ParPopPHist[ gi ], rnd ))
		{
			return( *ParPops[ rnd.getInt( ParPopCount )]);
		}

		return( *this );
	}

	/**
	 * Function selects an alternative, parallel optimizer's, population, to
	 * use in some solution generators.
	 *
	 * @param gi Solution generator index (0-3).
	 * @param rnd PRNG object.
	 */

	CBiteOptPop& selectAltPop( const int gi, CBiteRnd& rnd )
	{
		if( select( AltPopPHist, rnd ))
		{
			if( select( AltPopHist[ gi ], rnd ))
			{
				if( ParOptPop.getCurPopPos() >= CurPopSize )
				{
					return( ParOptPop );
				}
			}
			else
			{
				if( ParOpt2Pop.getCurPopPos() >= CurPopSize )
				{
					return( ParOpt2Pop );
				}
			}
		}

		return( *this );
	}

	/**
	 * Function returns a dynamically-selected minimal population index, used
	 * in some solution generation methods.
	 *
	 * @param gi Solution generator index (0-3).
	 * @param rnd PRNG object.
	 * @param ps Population size.
	 */

	int getMinSolIndex( const int gi, CBiteRnd& rnd, const int ps )
	{
		static const double pp[ 4 ] = { 0.05, 0.125, 0.25, 0.5 };
		const double r = ps * pow( rnd.get(),
			ps * pp[ select( MinSolPwrHist[ gi ], rnd )]);

		static const double rm[ 4 ] = { 0.0, 0.125, 0.25, 0.5 };

		return( (int) ( r * rm[ select( MinSolMulHist[ gi ], rnd )]));
	}

	/**
	 * The original "bitmask inversion with random move" solution generator.
	 * Most of the time adjusts only a single parameter of a better solution,
	 * yet manages to produce excellent "reference points".
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

		if( rnd.get() < AllpProbDamp )
		{
			if( select( Gen1AllpHist, rnd ))
			{
				DoAllp = true;
			}
		}

		if( DoAllp )
		{
			a = 0;
			b = ParamCount;
		}
		else
		{
			a = rnd.getInt( ParamCount );
			b = a + 1;
		}

		// Bitmask inversion operation, works as the main "driver" of
		// optimization process.

		const double r1 = rnd.get();
		const double r12 = r1 * r1;
		const int ims = (int) ( r12 * r12 * 48.0 );
		const ptype imask = (ptype) ( IntMantMask >> ims );

		const int im2s = rnd.getSqrInt( 96 );
		const ptype imask2 = (ptype) ( im2s > 63 ? 0 : IntMantMask >> im2s );

		const int si1 = (int) ( r1 * r12 * CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		for( i = a; i < b; i++ )
		{
			Params[ i ] = (( Params[ i ] ^ imask ) +
				( rp1[ i ] ^ imask2 )) >> 1;
		}

		if( select( Gen1MoveHist, rnd ))
		{
			const int si2 = rnd.getSqrInt( CurPopSize );
			const ptype* const rp2 = getParamsOrdered( si2 );

			if( select( Gen1MoveAsyncHist, rnd ))
			{
				a = 0;
				b = ParamCount;
			}

			// Random move around a random previous solution vector.

			static const double SpanMults[ 4 ] = { 0.5, 1.5, 2.0, 2.5 };

			const double m = SpanMults[ select( Gen1MoveSpanHist, rnd )];
			const double m1 = rnd.getTPDF() * m;
			const double m2 = rnd.getTPDF() * m;

			for( i = a; i < b; i++ )
			{
				Params[ i ] += (ptype) (( rp2[ i ] - Params[ i ]) * m1 );
				Params[ i ] += (ptype) (( rp2[ i ] - Params[ i ]) * m2 );
			}
		}
	}

	/**
	 * The "Differential Evolution"-based solution generator.
	 */

	void generateSol2( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int si1 = getMinSolIndex( 1, rnd, CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );
		const ptype* const rp3 = getParamsOrdered( CurPopSize1 - si1 );

		const int si2 = 1 + rnd.getInt( CurPopSize1 );
		const ptype* const rp2 = getParamsOrdered( si2 );

		const int si4 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp4 = getParamsOrdered( si4 );
		const ptype* const rp5 = getParamsOrdered( CurPopSize1 - si4 );

		// The "step in the right direction" (Differential Evolution
		// "mutation") operation towards the best (minimal) and away from
		// the worst (maximal) parameter vector, plus a difference of two
		// random vectors.

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rp1[ i ] + ((( rp2[ i ] - rp3[ i ]) +
				( rp4[ i ] - rp5[ i ])) >> 1 );
		}
	}

	/**
	 * An alternative "Differential Evolution"-based solution generator.
	 */

	void generateSol2b( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		// rand/2/none DE-alike mutation.

		const int si1 = getMinSolIndex( 2, rnd, CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		const int si2 = rnd.getInt( CurPopSize );
		const ptype* const rp2 = getParamsOrdered( si2 );
		const ptype* const rp3 = getParamsOrdered( CurPopSize1 - si2 );

		const CBiteOptPop& AltPop = selectAltPop( 0, rnd );

		const int si4 = rnd.getInt( CurPopSize );
		const ptype* const rp4 = AltPop.getParamsOrdered( si4 );
		const ptype* const rp5 = AltPop.getParamsOrdered( CurPopSize1 - si4 );

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rp1[ i ] + (( rp2[ i ] - rp3[ i ]) +
				( rp4[ i ] - rp5[ i ]));
		}
	}

	/**
	 * "Differential Evolution"-based solution generator, first implemented in
	 * the CDEOpt class.
	 */

	void generateSol2c( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		memset( Params, 0, ParamCount * sizeof( Params[ 0 ]));

		const int si1 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		const int PairCount = 3;
		const int pc = 1 + 2 * PairCount;
		int PopIdx[ pc ];
		PopIdx[ 0 ] = si1;

		int pp = 1;
		int i;
		int j;

		if( CurPopSize1 <= pc )
		{
			while( pp < pc )
			{
				PopIdx[ pp ] = rnd.getInt( CurPopSize );
				pp++;
			}
		}
		else
		{
			while( pp < pc )
			{
				const int sii = rnd.getInt( CurPopSize );

				for( j = 0; j < pp; j++ )
				{
					if( PopIdx[ j ] == sii )
					{
						break;
					}
				}

				if( j == pp )
				{
					PopIdx[ pp ] = sii;
					pp++;
				}
			}
		}

		for( j = 0; j < PairCount; j++ )
		{
			const ptype* const rp2 = getParamsOrdered( PopIdx[ 1 + j * 2 ]);
			const ptype* const rp3 = getParamsOrdered( PopIdx[ 2 + j * 2 ]);

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] += rp2[ i ] - rp3[ i ];
			}

			const int k = rnd.getInt( ParamCount );
			const int b = rnd.getInt( IntMantBits );

			Params[ k ] += ( (ptype) rnd.getBit() << b ) -
				( (ptype) rnd.getBit() << b );
		}

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rp1[ i ] + ( Params[ i ] >> 1 );
		}
	}

	/**
	 * "Centroid mix with DE" solution generator, works well for convex
	 * functions. For DE operation, uses a better solution and a random
	 * previous solution.
	 */

	void generateSol3( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const ptype* const MinParams = getParamsOrdered(
			getMinSolIndex( 3, rnd, CurPopSize ));

		const ptype* const cp = getCentroid();

		const int si1 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = ( rnd.getBit() ? cp[ i ] :
				MinParams[ i ] + ( MinParams[ i ] - rp1[ i ]));
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

		CBiteOptPop& AltPop = selectAltPop( 1, rnd );
		CBiteOptPop& ParPop = selectParPop( 1, rnd );

		int UseSize[ 2 ];
		UseSize[ 0 ] = CurPopSize;
		UseSize[ 1 ] = ParPop.getCurPopSize();

		const ptype** UseParams[ 2 ];
		UseParams[ 0 ] = AltPop.getPopParams();
		UseParams[ 1 ] = ParPop.getPopParams();

		const int km = 5 + ( select( Gen4MixFacHist, rnd ) << 1 );

		int p = rnd.getBit();
		int si1 = rnd.getSqrInt( UseSize[ p ]);
		const ptype* rp1 = UseParams[ p ][ si1 ];

		memcpy( Params, rp1, ParamCount * sizeof( Params[ 0 ]));

		int k;

		for( k = 1; k < km; k++ )
		{
			p = rnd.getBit();
			si1 = rnd.getSqrInt( UseSize[ p ]);
			rp1 = UseParams[ p ][ si1 ];

			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] ^= rp1[ i ];
			}
		}
	}

	/**
	 * A novel "Randomized bit crossing-over" candidate solution generation
	 * method. Effective, but on its own cannot stand coordinate system
	 * offsets, converges slowly. Completely mixes bits of two
	 * randomly-selected solutions, plus changes 1 random bit.
	 *
	 * This method is similar to a biological DNA crossing-over, but on a
	 * single-bit scale.
	 */

	void generateSol5( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const CBiteOptPop& ParPop = selectParPop( 2, rnd );

		const int si1 = rnd.getSqrInt( ParPop.getCurPopSize() );
		const ptype* const CrossParams1 = ( rnd.getBit() ?
			ParPop.getParamsOrdered( si1 ) :
			ParPop.getParamsOrdered( ParPop.getCurPopSize() - 1 - si1 ));

		const CBiteOptPop& AltPop = selectAltPop( 2, rnd );

		const ptype* const CrossParams2 = AltPop.getParamsOrdered(
			rnd.getSqrInt( CurPopSize ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit mixing mask.

			const ptype crpl = (ptype) ( rnd.getRaw() & IntMantMask );

			const ptype v1 = CrossParams1[ i ];
			const ptype v2 = CrossParams2[ i ];

			Params[ i ] = ( v1 & crpl ) | ( v2 & ~crpl );

			if( rnd.getBit() )
			{
				// Randomize a single bit, with 50% probability.

				const int b = rnd.getInt( IntMantBits );

				Params[ i ] += ( (ptype) rnd.getBit() << b ) -
					( (ptype) rnd.getBit() << b );
			}
		}
	}

	/**
	 * "Randomized parameter cross-over" solution generator. Similar to the
	 * "randomized bit cross-over", but works with the whole parameter values.
	 */

	void generateSol5b( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;
		const ptype* CrossParams[ 2 ];

		const CBiteOptPop& ParPop = selectParPop( 3, rnd );

		CrossParams[ 0 ] = ParPop.getParamsOrdered(
			rnd.getSqrInt( ParPop.getCurPopSize() ));

		const CBiteOptPop& AltPop = selectAltPop( 3, rnd );

		if( rnd.getBit() )
		{
			CrossParams[ 1 ] = AltPop.getParamsOrdered(
				CurPopSize1 - rnd.getSqrInt( CurPopSize ));
		}
		else
		{
			CrossParams[ 1 ] = AltPop.getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));
		}

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = CrossParams[ rnd.getBit()][ i ];
		}
	}

	/**
	 * A short-cut solution generator. Parameter value short-cuts: they
	 * considerably reduce convergence time for some functions while not
	 * severely impacting performance for other functions.
	 */

	void generateSol6( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const double r = rnd.getSqr();
		const int si = (int) ( r * r * CurPopSize );
		const double v = getRealValue( getParamsOrdered( si ),
			rnd.getInt( ParamCount ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = (ptype) (( v - MinValues[ i ]) * DiffValuesI[ i ]);
		}
	}

	/**
	 * A solution generator that randomly combines solutions from the main
	 * and "old" populations. Conceptually, it can be called a
	 * weighted-random crossover that combines solutions from diverse
	 * sources.
	 */

	void generateSol7( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const bool UseOldPop = ( OldPop.getCurPopPos() > 2 );
		static const double p[ 4 ] = { 1.5, 1.75, 2.0, 2.25 };
		const double pwr = p[ select( Gen7PowFacHist, rnd )];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const double rv = pow( rnd.get(), pwr );

			if( UseOldPop && rnd.getBit() && rnd.getBit() )
			{
				Params[ i ] = OldPop.getParamsOrdered(
					(int) ( rv * OldPop.getCurPopPos() ))[ i ];
			}
			else
			{
				Params[ i ] = getParamsOrdered(
					(int) ( rv * CurPopSize ))[ i ];
			}
		}
	}

	/**
	 * Solution generator that is DE-alike in its base. It calculates a
	 * centroid of a number of best solutions, and then applies "mutation"
	 * operation between the centroid and the solutions, using a random
	 * multiplier. This generator is similar to the "move" operation of
	 * generator 1.
	 */

	void generateSol8( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int NumSols = 5 + select( Gen8NumHist, rnd );
		const ptype* rp[ 8 ];

		// Calculate centroid of a number of selected solutions.

		int si0 = rnd.getSqrInt( CurPopSize );
		const ptype* rp0 = getParamsOrdered( si0 );
		rp[ 0 ] = rp0;
		memcpy( Params, rp0, ParamCount * sizeof( Params[ 0 ]));

		int j;
		int i;

		for( j = 1; j < NumSols; j++ )
		{
			si0 = rnd.getSqrInt( CurPopSize );
			rp0 = getParamsOrdered( si0 );
			rp[ j ] = rp0;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] += rp0[ i ];
			}
		}

		const double m = 1.0 / NumSols;

		for( i = 0; i < ParamCount; i++ )
		{
			NewValues[ i ] = Params[ i ] * m; // Centroid.
			Params[ i ] = (ptype) NewValues[ i ];
		}

		static const double Spans[ 4 ] = { 1.5, 2.5, 3.5, 4.5 };
		const double gm = Spans[ select( Gen8SpanHist, rnd )] * sqrt( m );

		for( j = 0; j < NumSols; j++ )
		{
			const double r = rnd.getGaussian() * gm;
			rp0 = rp[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] += (ptype) (( NewValues[ i ] - rp0[ i ]) * r );
			}
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

	const char** getHistNames() const
	{
		return( CurOpt -> getHistNames() );
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
				PushOpt = Opts[ rnd.getInt( OptCount )];

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
			CurOpt = PushOpt;
			StallCount++;
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
	CBiteOptWrap* CurOpt; ///< Current optimizer object.
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

typedef double( *biteopt_func )( int N, const double* x, void* func_data );

/**
 * Wrapper class for the biteopt_minimize() function.
 */

class CBiteOptMinimize : public CBiteOptDeep
{
public:
	int N; ///< The number of dimensions in objective function.
	biteopt_func f; ///< Objective function.
	void* data; ///< Objective function's data.
	const double* lb; ///< Parameters' lower bounds.
	const double* ub; ///< Parameters' upper bounds.

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
		return(( *f )( N, p, data ));
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
 * @param stopc Stopping criteria (convergence check). 0: off, 1: 128*N,
 * 2: 256*N.
 * @param rf Random number generator function; 0: use the default BiteOpt
 * PRNG. Note that the external RNG should be seeded externally.
 * @param rdata Data pointer to pass to the "rf" function.
 * @return The total number of function evaluations performed; useful if the
 * "stopc" was used.
 */

inline int biteopt_minimize( const int N, biteopt_func f, void* data,
	const double* lb, const double* ub, double* x, double* minf,
	const int iter, const int M = 1, const int attc = 10,
	const int stopc = 0, biteopt_rng rf = 0, void* rdata = 0 )
{
	CBiteOptMinimize opt;
	opt.N = N;
	opt.f = f;
	opt.data = data;
	opt.lb = lb;
	opt.ub = ub;
	opt.updateDims( N, M );

	CBiteRnd rnd;
	rnd.init( 1, rf, rdata );

	const int sct = ( stopc <= 0 ? 0 : 128 * N * stopc );
	const int useiter = (int) ( iter * sqrt( (double) M ));
	int evals = 0;
	int k;

	for( k = 0; k < attc; k++ )
	{
		opt.init( rnd );

		int i;

		for( i = 0; i < useiter; i++ )
		{
			const int sc = opt.optimize( rnd );

			if( sct > 0 && sc >= sct )
			{
				evals++;
				break;
			}
		}

		evals += i;

		if( k == 0 || opt.getBestCost() <= *minf )
		{
			memcpy( x, opt.getBestParams(), N * sizeof( x[ 0 ]));
			*minf = opt.getBestCost();
		}
	}

	return( evals );
}

#endif // BITEOPT_INCLUDED
