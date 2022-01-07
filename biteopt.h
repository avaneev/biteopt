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
 *
 * @version 2022.1
 */

#ifndef BITEOPT_INCLUDED
#define BITEOPT_INCLUDED

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
		addHist( PopChangeHist, "PopChangeHist" );
		addHist( ParOpt2Hist, "ParOpt2Hist" );
		addHist( ParPopPHist[ 0 ], "ParPopPHist[ 0 ]" );
		addHist( ParPopPHist[ 1 ], "ParPopPHist[ 1 ]" );
		addHist( ParPopPHist[ 2 ], "ParPopPHist[ 2 ]" );
		addHist( ParPopHist[ 0 ], "ParPopHist[ 0 ]" );
		addHist( ParPopHist[ 1 ], "ParPopHist[ 1 ]" );
		addHist( ParPopHist[ 2 ], "ParPopHist[ 2 ]" );
		addHist( AltPopPHist, "AltPopPHist" );
		addHist( AltPopHist[ 0 ], "AltPopHist[ 0 ]" );
		addHist( AltPopHist[ 1 ], "AltPopHist[ 1 ]" );
		addHist( AltPopHist[ 2 ], "AltPopHist[ 2 ]" );
		addHist( MinSolPwrHist[ 0 ], "MinSolPwrHist[ 0 ]" );
		addHist( MinSolPwrHist[ 1 ], "MinSolPwrHist[ 1 ]" );
		addHist( MinSolPwrHist[ 2 ], "MinSolPwrHist[ 2 ]" );
		addHist( MinSolMulHist[ 0 ], "MinSolMulHist[ 0 ]" );
		addHist( MinSolMulHist[ 1 ], "MinSolMulHist[ 1 ]" );
		addHist( MinSolMulHist[ 2 ], "MinSolMulHist[ 2 ]" );
		addHist( Gen1AllpHist, "Gen1AllpHist" );
		addHist( Gen1MoveHist, "Gen1MoveHist" );
		addHist( Gen1MoveAsyncHist, "Gen1MoveAsyncHist" );
		addHist( Gen1MoveSpanHist, "Gen1MoveSpanHist" );
		addHist( Gen4MixFacHist, "Gen4MixFacHist" );
		addHist( Gen5BinvHist, "Gen5BinvHist" );
		addHist( *ParOpt.getHists()[ 0 ], "ParOpt.CentPowHist" );
		addHist( *ParOpt.getHists()[ 1 ], "ParOpt.RadPowHist" );
		addHist( *ParOpt.getHists()[ 2 ], "ParOpt.EvalFacHist" );
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
			13 + aParamCount * 3 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		ParOpt.Owner = this;
		ParOpt.updateDims( aParamCount );
		ParOptPop.initBuffers( aParamCount, aPopSize );

		ParOpt2.Owner = this;
		ParOpt2.updateDims( aParamCount );
		ParOpt2Pop.initBuffers( aParamCount, aPopSize );
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

		AllpProbDamp = 1.8 / ParamCount;
		CentUpdateCtr = 0;

		ParOpt.init( rnd, InitParams, InitRadius );
		ParOpt2.init( rnd, InitParams, InitRadius );
		UseParOpt = 0;

		ParOptPop.resetCurPopPos();
		ParOpt2Pop.resetCurPopPos();

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
					generateSol2b( rnd );
				}
				else
				{
					generateSol3( rnd );
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
						generateSol4b( rnd );
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
						generateSol5b( rnd );
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
			generateSol1( rnd );
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
					ParOpt.init( rnd, getBestParams() );
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
					ParOpt2.init( rnd, getBestParams(), -1.0 );
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
				if( select( PopChangeHist, rnd ) == 0 )
				{
					// Increase population size on fail.

					incrCurPopSize( CurPopSize1 -
						(int) ( rnd.getRndValueSqr() * CurPopSize ));
				}
			}
		}
		else
		{
			applyHistsIncr( rnd );

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

			if( CurPopSize > PopSize / 2 )
			{
				if( select( PopChangeHist, rnd ) == 1 )
				{
					// Decrease population size on success.

					decrCurPopSize();
				}
			}
		}

		// "Diverging populations" technique.

		updateParPop( NewCost, TmpParams );

		CentUpdateCtr++;

		if( CentUpdateCtr >= CurPopSize * 32 )
		{
			// Update centroids of parallel populations that use running
			// average, to reduce error accumulation.

			CentUpdateCtr = 0;

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
	CBiteOptHist< 2 > PopChangeHist; ///< Population size change
		///< histogram.
		///<
	CBiteOptHist< 2 > ParOpt2Hist; ///< Parallel optimizer 2 use
		///< histogram.
		///<
	CBiteOptHist< 2 > ParPopPHist[ 3 ]; ///< Parallel population use
		///< probability histogram.
		///<
	CBiteOptHist< 4 > ParPopHist[ 3 ]; ///< Parallel population
		///< histograms for solution generators (template's Count parameter
		///< should match ParPopCount).
		///<
	CBiteOptHist< 2 > AltPopPHist; ///< Alternative population use
		///< histogram.
		///<
	CBiteOptHist< 2 > AltPopHist[ 3 ]; ///< Alternative population type use
		///< histogram.
		///<
	CBiteOptHist< 4 > MinSolPwrHist[ 3 ]; ///< Index of least-cost
		///< population, power factor.
		///<
	CBiteOptHist< 4 > MinSolMulHist[ 3 ]; ///< Index of least-cost
		///< population, multiplier.
		///<
	CBiteOptHist< 2 > Gen1AllpHist; ///< Generator method 1's Allp
		///< histogram.
		///<
	CBiteOptHist< 2 > Gen1MoveHist; ///< Generator method 1's Move
		///< histogram.
		////<
	CBiteOptHist< 2 > Gen1MoveAsyncHist; ///< Generator method 1's Move
		///< async histogram.
		///<
	CBiteOptHist< 4 > Gen1MoveSpanHist; ///< Generator method 1's Move span
		///< histogram.
		///<
	CBiteOptHist< 4 > Gen4MixFacHist; ///< Generator method 4's mixing
		///< count histogram.
		///<
	CBiteOptHist< 2 > Gen5BinvHist; ///< Generator method 5's random
		///< inversion technique histogram.
		///<
	int CentUpdateCtr; ///< Centroid update counter.
		///<
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
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
	 * @param gi Solution generator index (0-2).
	 * @param rnd PRNG object.
	 */

	CBiteOptPop& selectParPop( const int gi, CBiteRnd& rnd )
	{
		if( select( ParPopPHist[ gi ], rnd ))
		{
			return( *ParPops[ select( ParPopHist[ gi ], rnd )]);
		}

		return( *this );
	}

	/**
	 * Function selects an alternative, parallel optimizer's, population, to
	 * use in some solution generators.
	 *
	 * @param gi Solution generator index (0-2).
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

		if( rnd.getRndValue() < AllpProbDamp )
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
			a = (int) ( rnd.getRndValue() * ParamCount );
			b = a + 1;
		}

		// Bitmask inversion operation, works as the main "driver" of
		// optimization process.

		const double r1 = rnd.getRndValue();
		const double r12 = r1 * r1;
		const int ims = (int) ( r12 * r12 * 48.0 );
		const ptype imask = (ptype) ( ims > 63 ? 0 : IntMantMask >> ims );

		const int im2s = (int) ( rnd.getRndValueSqr() * 96.0 );
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
			const int si2 = (int) ( rnd.getRndValueSqr() * CurPopSize );
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
				Params[ i ] -= (ptype) (( Params[ i ] - rp2[ i ]) * m1 );
				Params[ i ] -= (ptype) (( Params[ i ] - rp2[ i ]) * m2 );
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

		const int si2 = 1 + (int) ( rnd.getRndValue() * CurPopSize1 );
		const ptype* const rp2 = getParamsOrdered( si2 );

		const int si4 = (int) ( rnd.getRndValueSqr() * CurPopSize );
		const ptype* const rp4 = getParamsOrdered( si4 );
		const ptype* const rp5 = getParamsOrdered( CurPopSize1 - si4 );

		// The "step in the right direction" (Differential Evolution
		// "mutation") operation towards the best (minimal) and away from
		// the worst (maximal) parameter vector, plus a difference of two
		// random vectors.

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
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

		const int si2 = (int) ( rnd.getRndValueSqr() * CurPopSize );
		const ptype* const rp2 = getParamsOrdered( si2 );

		const ptype* const rp3 = getParamsOrdered( CurPopSize1 - si2 );

		// Select two more previous solutions to be used in the mix.

		const CBiteOptPop& AltPop = selectAltPop( 0, rnd );

		const int si4 = (int) ( rnd.getRndValueSqr() * CurPopSize );
		const ptype* const rp4 = AltPop.getParamsOrdered( si4 );
		const ptype* const rp5 = AltPop.getParamsOrdered( CurPopSize1 - si4 );

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rp1[ i ] - ((( rp3[ i ] - rp2[ i ]) +
				( rp5[ i ] - rp4[ i ])) >> 1 );
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
			getMinSolIndex( 2, rnd, CurPopSize ));

		if( NeedCentUpdate )
		{
			updateCentroid();
		}

		const ptype* const cp = getCentroid();

		const int si1 = (int) ( rnd.getRndValueSqr() * CurPopSize );
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

		CBiteOptPop& AltPop = selectAltPop( 2, rnd );
		CBiteOptPop& ParPop = selectParPop( 1, rnd );

		int UseSize[ 2 ];
		UseSize[ 0 ] = CurPopSize;
		UseSize[ 1 ] = ParPop.getCurPopSize();

		const ptype** UseParams[ 2 ];
		UseParams[ 0 ] = AltPop.getPopParams();
		UseParams[ 1 ] = ParPop.getPopParams();

		const int km = 3 + ( select( Gen4MixFacHist, rnd ) << 1 );

		int p = rnd.getBit();
		int si1 = (int) ( rnd.getRndValueSqr() * UseSize[ p ]);
		const ptype* rp1 = UseParams[ p ][ si1 ];

		memcpy( Params, rp1, ParamCount * sizeof( Params[ 0 ]));

		int k;

		for( k = 1; k < km; k++ )
		{
			p = rnd.getBit();
			si1 = (int) ( rnd.getRndValueSqr() * UseSize[ p ]);
			rp1 = UseParams[ p ][ si1 ];

			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] ^= rp1[ i ];
			}
		}
	}

	/**
	 * Solution generator similar to generateSol4, but uses solutions from the
	 * main population only, and includes "crossover" approach first
	 * implemented in the generateSol5b() function.
	 */

	void generateSol4b( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int km = 3 + ( select( Gen4MixFacHist, rnd ) << 1 );

		int si1 = (int) ( rnd.getRndValueSqr() * CurPopSize );
		const ptype* rp1 = getParamsOrdered( si1 );

		memcpy( Params, rp1, ParamCount * sizeof( Params[ 0 ]));

		int k;

		for( k = 1; k < km; k++ )
		{
			si1 = (int) ( rnd.getRndValueSqr() * CurPopSize );
			int si2 = (int) ( rnd.getRndValueSqr() * CurPopSize );

			const ptype* CrossParams[ 2 ];
			CrossParams[ 0 ] = getParamsOrdered( si1 );
			CrossParams[ 1 ] = getParamsOrdered( si2 );

			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] ^= CrossParams[ rnd.getBit()][ i ];
			}
		}
	}

	/**
	 * A novel "Randomized bit crossing-over" candidate solution generation
	 * method. Effective, but on its own cannot stand coordinate system
	 * offsets, converges slowly. Completely mixes bits of two
	 * randomly-selected solutions, plus changes 1 random bit.
	 *
	 * This method is fundamentally similar to a biological DNA crossing-over.
	 */

	void generateSol5( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const CBiteOptPop& ParPop = selectParPop( 2, rnd );

		const ptype* const CrossParams1 = ParPop.getParamsOrdered(
			(int) ( rnd.getRndValueSqr() * ParPop.getCurPopSize() ));

		const CBiteOptPop& AltPop = selectAltPop( 1, rnd );

		const ptype* const CrossParams2 = AltPop.getParamsOrdered(
			(int) ( rnd.getRndValueSqr() * CurPopSize ));

		const bool UseInv = select( Gen5BinvHist, rnd );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit mixing mask.

			const ptype crpl = (ptype) ( rnd.getUniformRaw2() & IntMantMask );

			ptype v1 = CrossParams1[ i ];
			ptype v2 = ( UseInv && rnd.getBit() ?
				~CrossParams2[ i ] : CrossParams2[ i ]);

			if( rnd.getBit() )
			{
				const int b = (int) ( rnd.getRndValue() * IntMantBits );

				const ptype m = ~( (ptype) 1 << b );
				const ptype bv = (ptype) rnd.getBit() << b;

				v1 &= m;
				v2 &= m;
				v1 |= bv;
				v2 |= bv;
			}

			Params[ i ] = ( v1 & crpl ) | ( v2 & ~crpl );
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

		const CBiteOptPop& ParPop = selectParPop( 2, rnd );

		CrossParams[ 0 ] = ParPop.getParamsOrdered(
			(int) ( rnd.getRndValueSqr() * ParPop.getCurPopSize() ));

		const CBiteOptPop& AltPop = selectAltPop( 1, rnd );

		CrossParams[ 1 ] = AltPop.getParamsOrdered(
			(int) ( rnd.getRndValueSqr() * CurPopSize ));

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

		const double r = rnd.getRndValueSqr();
		const int si = (int) ( r * r * CurPopSize );
		const double v = getRealValue( getParamsOrdered( si ),
			(int) ( rnd.getRndValue() * ParamCount ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = (ptype) (( v - MinValues[ i ]) * DiffValuesI[ i ]);
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
 * @param stopc Stopping criteria (convergence check). 0: off, 1: 64*N,
 * 2: 128*N.
 * @return The total number of function evaluations performed; useful if the
 * "stopc" was used.
 */

inline int biteopt_minimize( const int N, biteopt_func f, void* data,
	const double* lb, const double* ub, double* x, double* minf,
	const int iter, const int M = 1, const int attc = 10,
	const int stopc = 0 )
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

	const int sct = ( stopc <= 0 ? 0 : 64 * N * stopc );
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
