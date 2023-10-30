//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBiteOpt and CBiteOptDeep classes.
 *
 * @section license License
 *
 * Copyright (c) 2016-2023 Aleksey Vaneev
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

#define BITEOPT_VERSION "2023.6"

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

	CBiteOpt()
	{
		setParPopCount( 4 );

		addSel( MethodSel, "MethodSel" );
		addSel( M1Sel, "M1Sel" );
		addSel( M1ASel, "M1ASel" );
		addSel( M1BSel, "M1BSel" );
		addSel( M1CSel, "M1CSel" );
		addSel( M2Sel, "M2Sel" );
		addSel( M2BSel, "M2BSel" );
		addSel( PopChangeIncrSel, "PopChangeIncrSel" );
		addSel( PopChangeDecrSel, "PopChangeDecrSel" );
		addSel( ParOpt2Sel, "ParOpt2Sel" );
		addSel( ParPopPSel[ 0 ], "ParPopPSel[ 0 ]" );
		addSel( ParPopPSel[ 1 ], "ParPopPSel[ 1 ]" );
		addSel( ParPopPSel[ 2 ], "ParPopPSel[ 2 ]" );
		addSel( ParPopPSel[ 3 ], "ParPopPSel[ 3 ]" );
		addSel( AltPopPSel, "AltPopPSel" );
		addSel( AltPopSel[ 0 ], "AltPopSel[ 0 ]" );
		addSel( AltPopSel[ 1 ], "AltPopSel[ 1 ]" );
		addSel( AltPopSel[ 2 ], "AltPopSel[ 2 ]" );
		addSel( AltPopSel[ 3 ], "AltPopSel[ 3 ]" );
		addSel( MinSolPwrSel[ 0 ], "MinSolPwrSel[ 0 ]" );
		addSel( MinSolPwrSel[ 1 ], "MinSolPwrSel[ 1 ]" );
		addSel( MinSolPwrSel[ 2 ], "MinSolPwrSel[ 2 ]" );
		addSel( MinSolPwrSel[ 3 ], "MinSolPwrSel[ 3 ]" );
		addSel( MinSolMulSel[ 0 ], "MinSolMulSel[ 0 ]" );
		addSel( MinSolMulSel[ 1 ], "MinSolMulSel[ 1 ]" );
		addSel( MinSolMulSel[ 2 ], "MinSolMulSel[ 2 ]" );
		addSel( MinSolMulSel[ 3 ], "MinSolMulSel[ 3 ]" );
		addSel( Gen1AllpSel, "Gen1AllpSel" );
		addSel( Gen1MoveAsyncSel, "Gen1MoveAsyncSel" );
		addSel( Gen1MoveSpanSel, "Gen1MoveSpanSel" );
		addSel( Gen2ModeSel, "Gen2ModeSel" );
		addSel( Gen2bModeSel, "Gen2bModeSel" );
		addSel( Gen2cModeSel, "Gen2cModeSel" );
		addSel( Gen2dModeSel, "Gen2dModeSel" );
		addSel( Gen3ModeSel, "Gen3ModeSel" );
		addSel( Gen4MixFacSel, "Gen4MixFacSel" );
		addSel( Gen5bModeSel, "Gen5bModeSel" );
		addSel( Gen7PowFacSel, "Gen7PowFacSel" );
		addSel( Gen8ModeSel, "Gen8ModeSel" );
		addSel( Gen8NumSel, "Gen8NumSel" );
		addSel( Gen8SpanSel[ 0 ], "Gen8SpanSel[ 0 ]" );
		addSel( Gen8SpanSel[ 1 ], "Gen8SpanSel[ 1 ]" );
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
			9 + aParamCount * 3 );

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
		initCommonVars( rnd );

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
			const ptype* const Params = getCurParams();

			for( i = 0; i < ParamCount; i++ )
			{
				NewValues[ i ] = getRealValue( Params, i );
			}

			NewCost = optcost( NewValues );
			updateBestCost( NewCost, NewValues,
				updatePop( NewCost, Params, false ));

			if( CurPopPos == PopSize )
			{
				updateCentroid();

				for( i = 0; i < ParPopCount; i++ )
				{
					ParPops[ i ] -> copy( *this );
				}

				DoInitEvals = false;
			}

			return( 0 );
		}

		DoEval = true;

		const int SelMethod = select( MethodSel, rnd );

		if( SelMethod == 0 )
		{
			generateSol2( rnd );
		}
		else
		if( SelMethod == 1 )
		{
			const int SelM1 = select( M1Sel, rnd );

			if( SelM1 == 0 )
			{
				const int SelM1A = select( M1ASel, rnd );

				if( SelM1A == 0 )
				{
					generateSol2b( rnd );
				}
				else
				if( SelM1A == 1 )
				{
					generateSol2c( rnd );
				}
				else
				{
					generateSol2d( rnd );
				}
			}
			else
			if( SelM1 == 1 )
			{
				if( select( M1BSel, rnd ))
				{
					generateSol4( rnd );
				}
				else
				{
					generateSol5b( rnd );
				}
			}
			else
			if( SelM1 == 2 )
			{
				if( select( M1CSel, rnd ))
				{
					generateSol5( rnd );
				}
				else
				{
					generateSol10( rnd );
				}
			}
			else
			{
				generateSol6( rnd );
			}
		}
		else
		if( SelMethod == 2 )
		{
			if( select( M2Sel, rnd ))
			{
				generateSol1( rnd );
			}
			else
			{
				const int SelM2B = select( M2BSel, rnd );

				if( SelM2B == 0 )
				{
					generateSol3( rnd );
				}
				else
				if( SelM2B == 1 )
				{
					generateSol7( rnd );
				}
				else
				if( SelM2B == 2 )
				{
					generateSol8( rnd );
				}
				else
				{
					generateSol9( rnd );
				}
			}
		}
		else
		{
			generateSolPar( rnd );
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

		const int p = updatePop( NewCost, TmpParams, true );

		if( p > CurPopSize1 )
		{
			// Upper bound cost constraint check failed, reject this solution.

			applySelsDecr( rnd );

			StallCount++;

			if( DoEval && CurPopSize < PopSize )
			{
				if( select( PopChangeIncrSel, rnd ))
				{
					// Increase population size on fail.

					incrCurPopSize();
				}
			}
		}
		else
		{
			updateBestCost( NewCost, NewValues, p );
			applySelsIncr( rnd, 1.0 - p * CurPopSizeI );

			StallCount = 0;

			if( rnd.get() < ParamCountI )
			{
				OldPop.updatePop( *getObjPtr( PopParams[ CurPopSize1 ]),
					PopParams[ CurPopSize1 ], false );
			}

			if( PushOpt != NULL && PushOpt != this &&
				!PushOpt -> DoInitEvals && p > 0 )
			{
				PushOpt -> updatePop( NewCost, TmpParams, true );
				PushOpt -> updateParPop( NewCost, TmpParams );
			}

			if( DoEval && CurPopSize > PopSize / 2 )
			{
				if( select( PopChangeDecrSel, rnd ))
				{
					// Decrease population size on success.

					decrCurPopSize();
				}
			}
		}

		// "Diverging populations" technique.

		updateParPop( NewCost, TmpParams );

		return( StallCount );
	}

protected:
	CBiteSel< 4 > MethodSel; ///< Population generator 4-method selector.
	CBiteSel< 4 > M1Sel; ///< Method 1's sub-method selector.
	CBiteSel< 3 > M1ASel; ///< Method 1's sub-sub-method A selector.
	CBiteSel< 2 > M1BSel; ///< Method 1's sub-sub-method B selector.
	CBiteSel< 2 > M1CSel; ///< Method 1's sub-sub-method C selector.
	CBiteSel< 2 > M2Sel; ///< Method 2's sub-method selector.
	CBiteSel< 4 > M2BSel; ///< Method 2's sub-sub-method B selector.
	CBiteSel< 2 > PopChangeIncrSel; ///< Population size change increase
		///< selector.
	CBiteSel< 2 > PopChangeDecrSel; ///< Population size change decrease
		///< selector.
	CBiteSel< 2 > ParOpt2Sel; ///< Parallel optimizer 2 use selector.
	CBiteSel< 2 > ParPopPSel[ 4 ]; ///< Parallel population use
		///< probability selectors.
	CBiteSel< 2 > AltPopPSel; ///< Alternative population use selector.
	CBiteSel< 2 > AltPopSel[ 4 ]; ///< Alternative population type use
		///< selectors.
	CBiteSel< 4 > MinSolPwrSel[ 4 ]; ///< Power factor selectors, for
		///< least-cost population index selection.
	CBiteSel< 4 > MinSolMulSel[ 4 ]; ///< Multiplier selectors, for
		///< least-cost population index selection.
	CBiteSel< 2 > Gen1AllpSel; ///< Generator method 1's Allp selector.
	CBiteSel< 2 > Gen1MoveAsyncSel; ///< Generator method 1's Move async
		///< selector.
	CBiteSel< 4 > Gen1MoveSpanSel; ///< Generator method 1's Move span
		///< selector.
	CBiteSel< 2 > Gen2ModeSel; ///< Generator method 2's Mode selector.
	CBiteSel< 2 > Gen2bModeSel; ///< Generator method 2b's Mode selector.
	CBiteSel< 2 > Gen2cModeSel; ///< Generator method 2c's Mode selector.
	CBiteSel< 2 > Gen2dModeSel; ///< Generator method 2d's Mode selector.
	CBiteSel< 4 > Gen3ModeSel; ///< Generator method 3's Mode selector.
	CBiteSel< 4 > Gen4MixFacSel; ///< Generator method 4's mixing count
		///< selector.
	CBiteSel< 2 > Gen5bModeSel; ///< Generator method 5b's Mode selector.
	CBiteSel< 4 > Gen7PowFacSel; ///< Generator method 7's Power selector.
	CBiteSel< 2 > Gen8ModeSel; ///< Generator method 8's mode selector.
	CBiteSel< 4 > Gen8NumSel; ///< Generator method 8's NumSols selector.
	CBiteSel< 4 > Gen8SpanSel[ 2 ]; ///< Generator method 8's random span
		///< selectors.
	CBitePop OldPop; ///< Population of older solutions, updated
		///< probabilistically.
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
	bool DoEval; ///< Temporary variable which equals to "true" if the
		///< newly-generated solution should be evaluated via the optcost()
		///< function.
	double NewCost; ///< Temporary variable that receives objective function's
		///< value (cost).

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
	CBitePop ParOptPop; ///< Population of parallel optimizer's solutions.
		///< Includes only its solutions.
	CParOpt< CNMSeqOpt > ParOpt2; ///< Parallel optimizer2.
	CBitePop ParOpt2Pop; ///< Population of parallel optimizer 2's solutions.
		///< Includes only its solutions.
	int UseParOpt; ///< Parallel optimizer currently being in use.

	/**
	 * Function updates an appropriate parallel population.
	 *
	 * @param UpdCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 */

	void updateParPop( const double UpdCost, const ptype* const UpdParams )
	{
		const int p = getMinDistParPop( UpdCost, UpdParams );

		if( p >= 0 )
		{
			ParPops[ p ] -> updatePop( UpdCost, UpdParams, true );
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

	CBitePop& selectParPop( const int gi, CBiteRnd& rnd )
	{
		if( select( ParPopPSel[ gi ], rnd ))
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

	CBitePop& selectAltPop( const int gi, CBiteRnd& rnd )
	{
		if( select( AltPopPSel, rnd ))
		{
			if( select( AltPopSel[ gi ], rnd ))
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
		const double r = ps * rnd.getPow( ps *
			pp[ select( MinSolPwrSel[ gi ], rnd )]);

		static const double rm[ 4 ] = { 0.0, 0.125, 0.25, 0.5 };

		return( (int) ( r * rm[ select( MinSolMulSel[ gi ], rnd )]));
	}

	/**
	 * The original "bitmask inversion with random move" solution generator.
	 * Most of the time adjusts only a single parameter of a better solution,
	 * yet manages to produce excellent "reference points".
	 */

	void generateSol1( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const CBitePop& ParPop = selectParPop( 0, rnd );

		copyParams( Params, ParPop.getParamsOrdered(
			getMinSolIndex( 0, rnd, ParPop.getCurPopSize() )));

		// Select a single random parameter or all parameters for further
		// operations.

		int a;
		int b;
		bool DoAllp = false;

		if( rnd.get() < 1.8 * ParamCountI )
		{
			if( select( Gen1AllpSel, rnd ))
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

		const int si1 = (int) ( r1 * r12 * ParPop.getCurPopSize() );
		const ptype* const rp1 = ParPop.getParamsOrdered( si1 );
		int i;

		for( i = a; i < b; i++ )
		{
			Params[ i ] = (( Params[ i ] ^ imask ) +
				( rp1[ i ] ^ imask2 )) >> 1;
		}

		if( rnd.get() < 1.0 - ParamCountI )
		{
			const ptype* const rp2 = getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			if( rnd.get() < sqrt( ParamCountI ))
			{
				if( select( Gen1MoveAsyncSel, rnd ))
				{
					a = 0;
					b = ParamCount;
				}
			}

			// Random move around a random previous solution vector.

			static const double SpanMults[ 4 ] = { 0.5, 1.5, 2.0, 2.5 };

			const double m = SpanMults[ select( Gen1MoveSpanSel, rnd )];
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
	 * The "Differential Evolution"-based solution generator. Note that
	 * compared to a usual DE, this generator does not use crossover, and
	 * it uses one, or an average of two best solutions as the base.
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

		const int Mode = select( Gen2ModeSel, rnd );
		int i;

		if( Mode == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + ((( rp2[ i ] - rp3[ i ]) +
					( rp4[ i ] - rp5[ i ])) >> 1 );
			}
		}
		else
		{
			const ptype* const rp1b = getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = (( rp1[ i ] + rp1b[ i ]) +
					( rp2[ i ] - rp3[ i ]) + ( rp4[ i ] - rp5[ i ])) >> 1;
			}
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

		const CBitePop& AltPop = selectAltPop( 0, rnd );

		const int si4 = rnd.getInt( CurPopSize );
		const ptype* const rp4 = AltPop.getParamsOrdered( si4 );
		const ptype* const rp5 = AltPop.getParamsOrdered( CurPopSize1 - si4 );

		const int Mode = select( Gen2bModeSel, rnd );
		int i;

		if( Mode == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + (( rp2[ i ] - rp3[ i ]) +
					( rp4[ i ] - rp5[ i ]));
			}
		}
		else
		{
			const ptype* const rp1b = getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = (( rp1[ i ] + rp1b[ i ]) >> 1 ) +
					( rp2[ i ] - rp3[ i ]) + ( rp4[ i ] - rp5[ i ]);
			}
		}
	}

	/**
	 * "Differential Evolution"-based solution generator, almost an exact
	 * replica of the CDEOpt optimizer.
	 */

	void generateSol2c( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;
		zeroParams( Params );

		const int si1 = rnd.getPowInt( 4.0, CurPopSize / 2 );
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

		const ptype* const rp2 = getParamsOrdered( PopIdx[ 1 ]);
		const ptype* const rp3 = getParamsOrdered( PopIdx[ 2 ]);
		const ptype* const rp4 = getParamsOrdered( PopIdx[ 3 ]);
		const ptype* const rp5 = getParamsOrdered( PopIdx[ 4 ]);
		const ptype* const rp6 = getParamsOrdered( PopIdx[ 5 ]);
		const ptype* const rp7 = getParamsOrdered( PopIdx[ 6 ]);

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = ( rp2[ i ] - rp3[ i ]) + ( rp4[ i ] - rp5[ i ]) +
				( rp6[ i ] - rp7[ i ]);
		}

		if( rnd.getBit() && rnd.getBit() )
		{
			const int k = rnd.getInt( ParamCount );

			// Produce sparsely-random bit-strings.

			const ptype v1 = (ptype) ( rnd.getRaw() & rnd.getRaw() &
				rnd.getRaw() & rnd.getRaw() & rnd.getRaw() & IntMantMask );

			const ptype v2 = (ptype) ( rnd.getRaw() & rnd.getRaw() &
				rnd.getRaw() & rnd.getRaw() & rnd.getRaw() & IntMantMask );

			Params[ k ] += v1 - v2; // Apply in TPDF manner.
		}

		const int Mode = select( Gen2cModeSel, rnd );

		if( Mode == 0 )
		{
			int si2 = si1 + rnd.getBit() * 2 - 1;

			if( si2 < 0 )
			{
				si2 = 1;
			}

			const ptype* const rp1b = getParamsOrdered( si2 );

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = ( rp1[ i ] + rp1b[ i ] + Params[ i ]) >> 1;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + ( Params[ i ] >> 1 );
			}
		}
	}

	/**
	 * An alternative "Differential Evolution"-based solution generator that
	 * uses "OldPop" population.
	 */

	void generateSol2d( CBiteRnd& rnd )
	{
		if( OldPop.getCurPopPos() < 3 )
		{
			generateSol2c( rnd );
			return;
		}

		ptype* const Params = TmpParams;

		const ptype* const rp1 = getParamsOrdered(
			rnd.getSqrInt( CurPopSize ));

		const ptype* const rp2 = getParamsOrdered(
			rnd.getInt( CurPopSize ));

		const ptype* const rp3 = OldPop.getParamsOrdered(
			rnd.getInt( OldPop.getCurPopPos() ));

		const int Mode = select( Gen2dModeSel, rnd );
		int i;

		if( Mode == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + (( rp2[ i ] - rp3[ i ]) >> 1 );
			}
		}
		else
		{
			const ptype* const rp1b = getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = (( rp1[ i ] + rp1b[ i ]) +
					( rp2[ i ] - rp3[ i ])) >> 1;
			}
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

		const ptype* const rp1 = getParamsOrdered(
			getMinSolIndex( 3, rnd, CurPopSize ));

		const ptype* const rp2 = getParamsOrdered(
			rnd.getSqrInt( CurPopSize ));

		const int Mode = select( Gen3ModeSel, rnd );
		int i;

		if( Mode == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + ( rp1[ i ] - rp2[ i ]);
			}
		}
		else
		{
			static const double CentProb[ 4 ] = { 0.0, 0.25, 0.5, 0.75 };
			const double p = CentProb[ Mode ];

			const ptype* const cp = getCentroid();

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = ( rnd.get() < p ? cp[ i ] :
					rp1[ i ] + ( rp1[ i ] - rp2[ i ]));
			}
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

		const CBitePop* UsePops[ 2 ];
		UsePops[ 0 ] = &selectAltPop( 1, rnd );
		UsePops[ 1 ] = &selectParPop( 1, rnd );

		int UseSize[ 2 ];
		UseSize[ 0 ] = CurPopSize;
		UseSize[ 1 ] = UsePops[ 1 ] -> getCurPopSize();

		const int km = 5 + ( select( Gen4MixFacSel, rnd ) << 1 );

		int p = rnd.getBit();
		const ptype* rp1 = UsePops[ p ] -> getParamsOrdered(
			rnd.getSqrInt( UseSize[ p ]));

		copyParams( Params, rp1 );
		int k;

		for( k = 1; k < km; k++ )
		{
			p = rnd.getBit();
			rp1 = UsePops[ p ] -> getParamsOrdered(
				rnd.getSqrInt( UseSize[ p ]));

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

		const CBitePop& ParPop = selectParPop( 2, rnd );

		const int si1 = rnd.getSqrInt( ParPop.getCurPopSize() );
		const ptype* const CrossParams1 = ParPop.getParamsOrdered( si1 );

		const CBitePop& AltPop = selectAltPop( 2, rnd );

		const ptype* const CrossParams2 = AltPop.getParamsOrdered(
			rnd.getSqrInt( CurPopSize ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit-mixing mask.

			const ptype crpl = (ptype) ( rnd.getRaw() & IntMantMask );

			Params[ i ] = ( CrossParams1[ i ] & crpl ) |
				( CrossParams2[ i ] & ~crpl );

			// Randomize a single bit, with 50% probability.

			const int b = rnd.getInt( IntMantBits );

			Params[ i ] += ( (ptype) rnd.getBit() << b ) -
				( (ptype) rnd.getBit() << b );
		}
	}

	/**
	 * "Randomized parameter cross-over" solution generator. Similar to the
	 * "randomized bit cross-over", but works with the whole parameter values.
	 */

	void generateSol5b( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;
		const ptype* CrossParams[ 4 ];

		const CBitePop& ParPop = selectParPop( 3, rnd );

		CrossParams[ 0 ] = ParPop.getParamsOrdered(
			rnd.getSqrInt( ParPop.getCurPopSize() ));

		const CBitePop& AltPop = selectAltPop( 3, rnd );

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

		const int Mode = select( Gen5bModeSel, rnd );
		int i;

		if( Mode == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = CrossParams[ rnd.getBit() ][ i ];
			}
		}
		else
		{
			CrossParams[ 2 ] = ParPop.getParamsOrdered(
				rnd.getSqrInt( ParPop.getCurPopSize() ));

			CrossParams[ 3 ] = AltPop.getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = CrossParams[ rnd.getBit() << 1 |
					rnd.getBit() ][ i ];
			}
		}
	}

	/**
	 * A short-cut solution generator. Parameter value short-cuts: they
	 * considerably reduce convergence time for some functions while not
	 * severely impacting performance for other functions.
	 *
	 * Can use variation with randomization between two values, and a slight
	 * move towards real 0.
	 */

	void generateSol6( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const double r = rnd.getPow( 4.0 );
		const int si = (int) ( r * CurPopSize );

		double v[ 2 ];
		v[ 0 ] = getRealValue( getParamsOrdered( si ),
			rnd.getInt( ParamCount ));

		if( rnd.getBit() )
		{
			v[ 1 ] = getRealValue( getParamsOrdered( si ),
				rnd.getInt( ParamCount ));
		}
		else
		{
			v[ 1 ] = v[ 0 ];
		}

		const double m = 1.0 - r * r;
		v[ 0 ] *= m; // Move towards real 0, useful for some functions.
		v[ 1 ] *= m;

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = (ptype) (( v[ rnd.getBit() ] - MinValues[ i ]) *
				DiffValuesI[ i ]);
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
		const double pwr = p[ select( Gen7PowFacSel, rnd )];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const double rv = rnd.getPow( pwr );

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

		const int Mode = select( Gen8ModeSel, rnd );
		const int NumSols = 5 + select( Gen8NumSel, rnd );
		const ptype* rp[ 8 ];

		// Calculate centroid of a number of selected solutions.

		const ptype* rp0 = getParamsOrdered( rnd.getSqrInt( CurPopSize ));
		rp[ 0 ] = rp0;
		copyParams( Params, rp0 );

		int j;
		int i;

		for( j = 1; j < NumSols; j++ )
		{
			rp0 = getParamsOrdered( rnd.getSqrInt( CurPopSize ));
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

		// Apply "move" operations in one of two modes.

		if( Mode == 0 )
		{
			static const double Spans[ 4 ] = { 1.5, 2.5, 3.5, 4.5 };
			const double gm = Spans[ select( Gen8SpanSel[ Mode ], rnd )] *
				sqrt( m );

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
		else
		{
			static const double Spans[ 4 ] = { 0.5, 1.5, 2.5, 3.5 };
			const double gm = Spans[ select( Gen8SpanSel[ Mode ], rnd )];

			for( j = 0; j < NumSols; j++ )
			{
				const double r = rnd.getGaussian() * gm;
				rp0 = rp[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] += (ptype) (( Params[ i ] - rp0[ i ]) * r );
				}
			}
		}
	}

	/**
	 * A "water drain" solution generator: makes a fixed-multiplier step from
	 * a better random solution 1 towards or away from worse random solution
	 * 2. Moderately efficient on its own.
	 */

	void generateSol9( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int si1 = rnd.getInt( CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		const int si2 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp2 = getParamsOrdered( CurPopSize1 - si2 );
		int i;

		// Such overall sign inversion seems unuseful, but has benefits in
		// practice.

		if( rnd.getBit() )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] - (( rp2[ i ] - rp1[ i ]) >> 1 ) *
					( 1 - 2 * rnd.getBit() );
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rp1[ i ] + (( rp2[ i ] - rp1[ i ]) >> 1 ) *
					( 1 - 2 * rnd.getBit() );
			}
		}
	}

	/**
	 * Solution generator based on SpherOpt's converging hyper-spheroid.
	 */

	void generateSol10( CBiteRnd& rnd )
	{
		ptype* const Params = TmpParams;

		const int si1 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp1 = getParamsOrdered( si1 );

		const int si2 = rnd.getSqrInt( CurPopSize );
		const ptype* const rp2 = getParamsOrdered( CurPopSize1 - si2 );
		int i;

		// Calculate centroid.

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = ( rp1[ i ] + rp2[ i ]) >> 1;
		}

		// Calculate radius.

		double Radius = 0.0;

		for( i = 0; i < ParamCount; i++ )
		{
			const ptype v1 = rp1[ i ] - Params[ i ];
			const ptype v2 = rp2[ i ] - Params[ i ];
			Radius += (double) v1 * v1 + 0.45 * v2 * v2;
		}

		// Select a point on a hyper-spheroid.

		double s2 = 1e-300;

		for( i = 0; i < ParamCount; i++ )
		{
			NewValues[ i ] = rnd.get() - 0.5;
			s2 += NewValues[ i ] * NewValues[ i ];
		}

		// Add hyper-spheroid-based offset to the centroid.

		const double d = sqrt( Radius / s2 );

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] += (ptype) ( NewValues[ i ] * d );
		}
	}

	/**
	 * Solution generator that obtains solution from an independently-running
	 * parallel optimizer.
	 */

	void generateSolPar( CBiteRnd& rnd )
	{
		DoEval = false;
		CBitePop* UpdPop;

		if( UseParOpt == 1 )
		{
			// Re-assign optimizer 2 based on comparison of its
			// efficiency with optimizer 1.

			UseParOpt = select( ParOpt2Sel, rnd );
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

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = (ptype) (( NewValues[ i ] -
				MinValues[ i ]) * DiffValuesI[ i ]);
		}

		UpdPop -> updatePop( NewCost, TmpParams, false );
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
	 * Function returns a pointer to an array of selectors in use by the
	 * current CBiteOpt object.
	 */

	CBiteSelBase** getSels()
	{
		return( CurOpt -> getSels() );
	}

	/**
	 * Function returns a pointer to an array of selector names.
	 */

	const char** getSelNames() const
	{
		return( CurOpt -> getSelNames() );
	}

	/**
	 * Function returns the number of selectors in use by the current CBiteOpt
	 * object.
	 */

	int getSelCount() const
	{
		return( CurOpt -> getSelCount() );
	}

	/**
	 * Function updates dimensionality of *this object. Function does nothing
	 * if dimensionality has not changed since the last call. This function
	 * should be called at least once before calling the init() function.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param M The number of CBiteOpt objects. This number depends on the
	 * complexity of the objective function, if the default value does not
	 * produce a good solution, it should be increased together with the
	 * iteration count. Minimal value is 1, in this case a plain CBiteOpt
	 * optimization will be performed.
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
	int OptCount; ///< The total number of optimization objects in use.
	CBiteOptWrap** Opts; ///< Optimization objects.
	CBiteOptWrap* BestOpt; ///< Optimizer that contains the best solution.
	CBiteOptWrap* CurOpt; ///< Current optimizer object.
	int StallCount; ///< The number of iterations without improvement.

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
		memcpy( p, lb, N * sizeof( p[ 0 ]));
	}

	virtual void getMaxValues( double* const p ) const
	{
		memcpy( p, ub, N * sizeof( p[ 0 ]));
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
 * @param f_minp If non-zero, a pointer to the stopping value: optimization
 * will stop when this objective value is reached.
 * @return The total number of function evaluations performed; useful if the
 * "stopc" and/or "*f_minp" were used.
 */

inline int biteopt_minimize( const int N, biteopt_func f, void* data,
	const double* lb, const double* ub, double* x, double* minf,
	const int iter, const int M = 1, const int attc = 10,
	const int stopc = 0, biteopt_rng rf = 0, void* rdata = 0,
	double* f_minp = 0 )
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

		bool IsFinished = false;
		int i;

		for( i = 0; i < useiter; i++ )
		{
			const int sc = opt.optimize( rnd );

			if( f_minp != 0 && opt.getBestCost() <= *f_minp )
			{
				evals++;
				IsFinished = true;
				break;
			}

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

		if( IsFinished )
		{
			break;
		}
	}

	return( evals );
}

#endif // BITEOPT_INCLUDED
