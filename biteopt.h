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
 * @version 2021.7
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

class CBiteOpt : public CBiteOptBase
{
public:
	/**
	 * Constructor.
	 */

	CBiteOpt()
		: MantMult( 1ULL << MantSize )
		, MantMultI( 1.0 / ( 1ULL << MantSize ))
		, IntParams( NULL )
	{
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
		const int aPopSize = ( PopSize0 > 0 ?
			PopSize0 : 7 + aParamCount * 3 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();
		initBaseBuffers( aParamCount, aPopSize );

		IntParams = new uint64_t[ aParamCount ];
		Params = CurParams[ aPopSize ];

		ParOpt.Owner = this;
		ParOpt.updateDims( aParamCount );
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
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
		updateDiffValues();

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		const double sd = 0.25 * InitRadius;
		int i;
		int j;

		if( InitParams == NULL )
		{
			resetCentroid();

			for( j = 0; j < PopSize; j++ )
			{
				double* const p = CurParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					const double v = wrapParam( rnd,
						0.5 + getGaussian( rnd ) * sd );

					p[ i ] = v;
					CentParams[ i ] += v;
				}
			}
		}
		else
		{
			double* const p0 = CurParams[ 0 ];

			for( i = 0; i < ParamCount; i++ )
			{
				const double v = wrapParam( rnd,
					( InitParams[ i ] - MinValues[ i ]) / DiffValues[ i ]);

				p0[ i ] = v;
				CentParams[ i ] = v;
			}

			for( j = 1; j < PopSize; j++ )
			{
				double* const p = CurParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					const double v = wrapParam( rnd,
						p0[ i ] + getGaussian( rnd ) * sd );

					p[ i ] = v;
					CentParams[ i ] += v;
				}
			}
		}

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] *= PopSizeI;
		}

		ParamCountRnd = (double) ParamCount / rnd.getRawScale();
		ParamCntr = (int) ( rnd.getUniformRaw() * ParamCountRnd );

		ParOpt.init( rnd, InitParams, InitRadius );

		ParPopHist.reset();
		ParPopMHist.reset();
		SelParPopM = 0;

		ScutHist.reset();
		MethodHist.reset();
		DrawHist.reset();
		D3Hist.reset();
		Gen1AllpHist.reset();
		Gen1CentHist.reset();
		Gen1SpanHist.reset();

		const double AllpProbDamp = ( ParamCount < 3 ? 1.0 :
			2.0 / ParamCount ); // Allp probability damping. Applied for
			// higher dimensions as the "all parameter" randomization is
			// ineffective in higher dimensions.

		AllpProbs[ 0 ] = (int) ( 0.6 * AllpProbDamp *
			CBiteRnd :: getRawScale() );

		AllpProbs[ 1 ] = (int) ( 0.9 * AllpProbDamp *
			CBiteRnd :: getRawScale() );

		PrevSelMethod = MethodHist.selectRandom( rnd );
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
	 * threshold value is ParamCount * 16. When this value was reached, the
	 * probability of plateau is high. This value, however, should not be
	 * solely relied upon when considering a stopping criteria: a hard
	 * iteration limit should be always used as in some cases convergence time
	 * may be very high with small, but frequent improving steps. This value
	 * is best used to allocate iteration budget between optimization attempts
	 * more efficiently.
	 */

	int optimize( CBiteRnd& rnd, CBiteOpt* const PushOpt = NULL )
	{
		double NewCost;
		int i;

		if( DoInitEvals )
		{
			const double* const p = CurParams[ InitEvalIndex ];

			for( i = 0; i < ParamCount; i++ )
			{
				NewParams[ i ] = getRealValue( p, i );
			}

			NewCost = optcost( NewParams );
			sortPop( NewCost, InitEvalIndex );
			updateBestCost( NewCost, NewParams );

			InitEvalIndex++;

			if( InitEvalIndex == PopSize )
			{
				for( i = 0; i < ParPopCount; i++ )
				{
					ParPops[ i ].copy( *this );
				}

				DoInitEvals = false;
			}

			return( 0 );
		}

		const int SelParPop = ParPopHist.select( rnd );
		UseParPop = &ParPops[ SelParPop ];
		WasParPopUsed = false;

		bool DoEval = true;
		int SelMethod = -1;
		int SelDraw;
		int SelD3;

		static const int ScutProbs[ 2 ] = {
			(int) ( 0.03 * CBiteRnd :: getRawScale() ),
			(int) ( 0.09 * CBiteRnd :: getRawScale() )
		}; // Short-cut probability range, in raw scale.

		int SelScut = ScutHist.select( rnd );

		if( rnd.getUniformRaw() < ScutProbs[ SelScut ])
		{
			// Parameter value short-cuts, they considerably reduce
			// convergence time for some functions while not severely
			// impacting performance for other functions.

			i = (int) ( rnd.getUniformRaw() * ParamCountRnd );

			const double r = rnd.getRndValue();
			const double r2 = r * r;

			if( getBit( rnd ))
			{
				// "Centroid offset" short-cut.

				const int si = (int) ( r2 * r2 * PopSize );
				const double* const rp = getParamsOrdered( si );

				const double v = getRealValue( rp, i ) -
					getRealValue( CentParams, i );

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = ( getRealValue( rp, i ) - v -
						MinValues[ i ]) / DiffValues[ i ];
				}
			}
			else
			{
				// "Same-value parameter vector" short-cut.

				const int si = (int) ( r * r2 * PopSize );
				const double* const rp = getParamsOrdered( si );

				const double v = getRealValue( rp, i );

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = ( v - MinValues[ i ]) / DiffValues[ i ];
				}
			}
		}
		else
		{
			SelScut = getBit( rnd );
			SelDraw = DrawHist.select( rnd );
			SelD3 = -1;

			if( SelDraw == 0 )
			{
				SelMethod = PrevSelMethod;
			}
			else
			if( SelDraw == 1 )
			{
				SelMethod = MethodHist.selectRandom( rnd );
			}
			else
			{
				SelMethod = MethodHist.select( rnd );
			}

			if( SelMethod == 0 )
			{
				ParOpt.optimize( rnd, &NewCost, Params );
				DoEval = false;
			}
			else
			{
				mp = rnd.getRndValue();
				mp2 = mp * mp;
				mpi = (int) ( mp * mp2 * 4 );

				if( SelMethod == 1 )
				{
					generateSol1( rnd );
				}
				else
				if( SelMethod == 2 )
				{
					generateSol2( rnd );
				}
				else
				{
					SelD3 = D3Hist.select( rnd );

					if( SelD3 == 0 )
					{
						generateSol3( rnd );
					}
					else
					if( SelD3 == 1 )
					{
						generateSol4( rnd );
					}
					else
					{
						generateSol5( rnd );
					}
				}
			}
		}

		// Evaluate objective function with new parameters.

		if( DoEval )
		{
			// Wrap parameter values so that they stay in the [0; 1] range.

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = wrapParam( rnd, Params[ i ]);
				NewParams[ i ] = getRealValue( Params, i );
			}

			NewCost = optcost( NewParams );

			updateBestCost( NewCost, NewParams );
		}
		else
		{
			updateBestCost( NewCost, Params, true );
		}

		if( PushOpt != NULL && PushOpt != this &&
			!PushOpt -> DoInitEvals )
		{
			PushOpt -> updatePop( NewCost, Params, true );
		}

		if( NewCost > CurCosts[ PopSize1 ])
		{
			// Upper bound cost constraint check failed, reject this solution.

			if( WasParPopUsed )
			{
				ParPopHist.decr( SelParPop );
				ParPopMHist.decr( SelParPopM );
			}

			ScutHist.decr( SelScut );

			if( SelMethod >= 0 )
			{
				PrevSelMethod = MethodHist.selectRandom( rnd );

				MethodHist.decr( SelMethod );
				DrawHist.decr( SelDraw );

				if( SelD3 < 0 )
				{
					if( SelMethod == 1 )
					{
						Gen1AllpHist.decr( SelGen1Allp );
						Gen1CentHist.decr( SelGen1Cent );

						if( SelGen1Span >= 0 )
						{
							Gen1SpanHist.decr( SelGen1Span );
						}
					}
				}
				else
				{
					D3Hist.decr( SelD3 );
				}
			}

			StallCount++;
		}
		else
		{
			if( WasParPopUsed )
			{
				ParPopHist.incr( SelParPop );
				ParPopMHist.incr( SelParPopM );
			}

			ScutHist.incr( SelScut );

			if( SelMethod >= 0 )
			{
				PrevSelMethod = SelMethod;

				MethodHist.incr( SelMethod );
				DrawHist.incr( SelDraw );

				if( SelD3 < 0 )
				{
					if( SelMethod == 1 )
					{
						Gen1AllpHist.incr( SelGen1Allp );
						Gen1CentHist.incr( SelGen1Cent );

						if( SelGen1Span >= 0 )
						{
							Gen1SpanHist.incr( SelGen1Span );
						}
					}
				}
				else
				{
					D3Hist.incr( SelD3 );
				}
			}

			if( NewCost == CurCosts[ PopSize1 ])
			{
				StallCount++;
			}
			else
			{
				StallCount = 0;
			}

			updatePop( NewCost, Params );
		}

		// Diverging populations technique.

		SelParPopM = ParPopMHist.select( rnd );
		double d = ParPops[ 0 ].getDistanceSqr( Params );
		int p = 0;

		for( i = 1; i < ParPopCount; i++ )
		{
			const double nd = ParPops[ i ].getDistanceSqr( Params );

			if( SelParPopM == 0 )
			{
				if( nd > d )
				{
					d = nd;
					p = i;
				}
			}
			else
			{
				if( nd < d )
				{
					d = nd;
					p = i;
				}
			}
		}

		ParPops[ p ].updatePop( NewCost, Params, true );

		return( StallCount );
	}

protected:
	static const int MantSize = 54; ///< Mantissa size of the bitmask
		///< operations.
		///<
	static const uint64_t MantSizeMask = ( 1ULL << MantSize ) - 1; ///< Mask
		///< that corresponds to mantissa.
		///<
	static const int ParPopCount = 3; ///< The number of parallel populations.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double MantMultI; ///< =1/MantMult.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	int AllpProbs[ 2 ]; ///< Generator method 1's Allp probability range,
		///< in raw scale.
		///<
	CBiteOptPop ParPops[ ParPopCount ]; ///< Parallel populations.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<
	uint64_t* IntParams; ///< Temporary integer value parameter buffer.
		///<
	CBiteOptHist< 3, 3, 1 > ParPopHist; ///< Parallel population histogram
		///< (count must match ParPopCoun).
		///<
	CBiteOptHistBinary ParPopMHist; ///< Parallel population update
		///< method histogram.
		///<
	int SelParPopM; ///< Parallel population update method.
		///<
	CBiteOptPop* UseParPop; ///< Parallel population to use in some solution
		///< generation methods.
		///<
	bool WasParPopUsed; ///< "True" if UseParPop was used in some solution
		///< generation.
		///<
	CBiteOptHistBinary ScutHist; ///< Short-cut method's histogram.
		///<
	CBiteOptHist< 4, 4, 2 > MethodHist; ///< Population generator method
		///< histogram.
		///<
	CBiteOptHist< 3, 6, 1 > DrawHist; ///< Method draw histogram.
		///<
	CBiteOptHist< 3, 6, 1 > D3Hist; ///< Draw method 3's histogram.
		///<
	CBiteOptHistBinary Gen1AllpHist; ///< Generator method 1's Allp
		///< histogram.
		///<
	CBiteOptHistBinary Gen1CentHist; ///< Generator method 1's Cent
		///< histogram.
		///<
	CBiteOptHistBinary Gen1SpanHist; ///< Generator method 1's Cent
		///< histogram.
		///<
	int PrevSelMethod; ///< Previously successfully used method; contains
		///< random method index if optimization was not successful.
		///<
	int SelGen1Allp; ///< Generator method 1's Allp selector.
		///<
	int SelGen1Cent; ///< Generator method 1's Cent selector.
		///<
	int SelGen1Span; ///< Generator method 1's Span selector.
		///<
	double mp, mp2; ///< Temporary variables.
		///<
	int mpi; ///< Temporary variable.
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

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] IntParams;
	}

	/**
	 * The original "bitmask inversion" solution generator. Most of the time
	 * adjusts only a single parameter of the very best solution, yet manages
	 * to produce excellent "reference points".
	 */

	void generateSol1( CBiteRnd& rnd )
	{
		WasParPopUsed = true;
		memcpy( Params, UseParPop -> getParamsOrdered( mpi ),
			ParamCount * sizeof( Params[ 0 ]));

		// Select a single random parameter or all parameters for further
		// operations.

		int i;
		int a;
		int b;

		SelGen1Allp = Gen1AllpHist.select( rnd );

		if( rnd.getUniformRaw() < AllpProbs[ SelGen1Allp ])
		{
			a = 0;
			b = ParamCount - 1;
		}
		else
		{
			SelGen1Allp = getBit( rnd );

			a = ParamCntr;
			b = ParamCntr;
			ParamCntr = ( ParamCntr == 0 ? ParamCount : ParamCntr ) - 1;
		}

		// Bitmask inversion operation, works as a "driver" of optimization
		// process.

		const int imasks = (int) ( mp2 * mp2 * MantSize );
		const uint64_t imask = ( imasks > 63 ? 0 : MantSizeMask >> imasks );

		const double rr = rnd.getRndValue();
		const int imask2s = (int) ( rr * rr * MantSize * 2.0 );
		const uint64_t imask2 = ( imask2s > 63 ? 0 : MantSizeMask >> imask2s );

		const int si0 = (int) ( mp * mp2 * PopSize );
		const double* const rp0 = getParamsOrdered( si0 );

		for( i = a; i <= b; i++ )
		{
			const uint64_t v1 = (uint64_t) ( Params[ i ] * MantMult );
			const uint64_t v2 = (uint64_t) ( rp0[ i ] * MantMult );
			uint64_t v0 = (( v1 ^ imask ) + ( v2 ^ imask2 )) >> 1;
			Params[ i ] = v0 * MantMultI;
		}

		static const int CentProbs[ 2 ] = {
			(int) ( 0.6 * CBiteRnd :: getRawScale() ),
			(int) ( 0.95 * CBiteRnd :: getRawScale() )
		}; // "Centroid move" probability range.

		SelGen1Cent = Gen1CentHist.select( rnd );
		SelGen1Span = -1;

		if( rnd.getUniformRaw() < CentProbs[ SelGen1Cent ])
		{
			// Random move around random previous solution vector.

			static const double SpanMults[ 2 ] = {
				1.5 * CBiteRnd :: getRawScaleInv(),
				3.0 * CBiteRnd :: getRawScaleInv()
			};

			SelGen1Span = Gen1SpanHist.select( rnd );
			const double m = SpanMults[ SelGen1Span ];
			const double m1 = rnd.getTPDFRaw() * m;
			const double m2 = rnd.getTPDFRaw() * m;

			const int si = (int) ( mp2 * PopSize );
			const double* const rp1 = getParamsOrdered( si );

			for( i = a; i <= b; i++ )
			{
				Params[ i ] -= ( Params[ i ] - rp1[ i ]) * m1;
				Params[ i ] -= ( Params[ i ] - rp1[ i ]) * m2;
			}
		}
		else
		{
			SelGen1Cent = getBit( rnd );
		}
	}

	/**
	 * The original "Digital Evolution"-based solution generator.
	 */

	void generateSol2( CBiteRnd& rnd )
	{
		// Select worst and a random previous solution from the ordered list,
		// apply offsets to reduce sensitivity to noise.

		const double* const MinParams = getParamsOrdered( mpi );
		const int si0 = mpi + (int) ( mp * ( PopSize1 - mpi ));
		const double* const OrigParams = getParamsOrdered( si0 );
		const double* const MaxParams = getParamsOrdered( PopSize1 - mpi );

		// Select two more previous solutions to be used in the mix.

		const double r = rnd.getRndValue();
		const int si1 = (int) ( r * r * PopSize );
		const double* const rp1 = getParamsOrdered( si1 );
		const double* const rp2 = getParamsOrdered( PopSize1 - si1 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// The "step in the right direction" (Differential Evolution
			// "mutation") operation towards the best (minimal) and away from
			// the worst (maximal) parameter vector, plus a difference of two
			// random vectors.

			Params[ i ] = MinParams[ i ] -
				(( MaxParams[ i ] - OrigParams[ i ]) -
				( rp1[ i ] - rp2[ i ])) * 0.5;
		}
	}

	/**
	 * Alternative randomized solution generator, works well for convex
	 * functions. Uses the very best solution and a random previous solution.
	 * "mp * mp" is equivalent of giving more weight to better solutions.
	 */

	void generateSol3( CBiteRnd& rnd )
	{
		WasParPopUsed = true;
		const double* const MinParams = UseParPop -> getParamsOrdered( mpi );
		const double* const cp = UseParPop -> getCentroid();
		const int si1 = (int) ( mp2 * PopSize );
		const double* const rp1 = getParamsOrdered( si1 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const double m1 = getBit( rnd );
			const double m2 = 1.0 - m1;

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
		int i;
		int k;

		for( k = 0; k < 7; k++ )
		{
			const double r = rnd.getRndValue();
			const int si = (int) ( r * r * PopSize );
			const double* const rp = getParamsOrdered( si );

			if( k == 0 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					IntParams[ i ] = (uint64_t) ( rp[ i ] * MantMult );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					IntParams[ i ] ^= (uint64_t) ( rp[ i ] * MantMult );
				}
			}
		}

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = IntParams[ i ] * MantMultI;
		}
	}

	/**
	 * A novel "Randomized bit crossing-over" candidate solution generation
	 * method. Effective, but on its own cannot stand coordinate system
	 * offsets, converges slowly. Completely mixes bits of two existing
	 * solutions, plus changes 1 random bit.
	 *
	 * This method is fundamentally similar to biological DNA crossing-over.
	 */

	void generateSol5( CBiteRnd& rnd )
	{
		WasParPopUsed = true;
		const double* const CrossParams1 =
			UseParPop -> getParamsOrdered( (int) ( mp2 * PopSize ));

		const double p2 = rnd.getRndValue();
		const double* const CrossParams2 =
			getParamsOrdered( (int) ( p2 * p2 * PopSize ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit mixing mask.

			const uint64_t crpl = ( rnd.getUniformRaw() |
				( (uint64_t) rnd.getUniformRaw() << rnd.getRawBitCount() )) &
				( rnd.getUniformRaw() |
				( (uint64_t) rnd.getUniformRaw() << rnd.getRawBitCount() ));
				// Produce sparsely-random mask.

			uint64_t v1 = (uint64_t) ( CrossParams1[ i ] * MantMult );
			uint64_t v2 = (uint64_t) ( CrossParams2[ i ] * MantMult );

			if( getBit( rnd ))
			{
				const int b = (int) ( rnd.getRndValue() * MantSize );
				const uint64_t m = ~( 1ULL << b );
				const uint64_t bv = (uint64_t) getBit( rnd ) << b;

				v1 &= m;
				v2 &= m;
				v1 |= bv;
				v2 |= bv;
			}

			Params[ i ] = (( v1 & crpl ) | ( v2 & ~crpl )) * MantMultI;
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

	void updateDims( const int aParamCount, const int M = 8,
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
	 * threshold value is ParamCount * 16.
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
				const double r = rnd.getRndValue();
				PushOpt = Opts[ (int) ( r * OptCount )];

				if( PushOpt != CurOpt )
				{
					break;
				}
			}
		}

		const int sc = CurOpt -> optimize( rnd, PushOpt );

		if( BestOpt -> getBestCost() >= CurOpt -> getBestCost() )
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
