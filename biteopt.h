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
 * @version 2021.4
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
		ParOpt.CentPow = 7.0;
		ParOpt.RadPow = 18.0;
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
			PopSize0 : 13 + aParamCount * 2 );

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
	 * @param InitRadius Initial radius, relative to the default value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();
		resetCentroid();
		updateDiffValues();

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		const double sd = 0.25 * InitRadius;
		int i;
		int j;

		if( InitParams == NULL )
		{
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
				CentParams[ i ] += v;
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

		ScutCntr = rnd.getRndValue();
		AllpCntr = rnd.getRndValue();
		CentCntr = rnd.getRndValue();
		ParamCountRnd = (double) ParamCount / rnd.getRawScale();
		ParamCntr = (int) ( rnd.getUniformRaw() * ParamCountRnd );
		AllpProbDamp = 2.0 / ParamCount;

		ParOpt.init( rnd, InitParams, InitRadius );

		ScutHist.reset();
		MethodHist.reset();
		DrawHist.reset();
		D3Hist.reset();
		Gen1Hist.reset();

		PrevSelMethod = MethodHist.selectRandom( rnd );
		mp = 0.0;
		mp2 = 0.0;
		mpi = 0;
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
			insertPopOrder( NewCost, InitEvalIndex, InitEvalIndex );
			updateBestCost( NewCost, NewParams );

			InitEvalIndex++;

			if( InitEvalIndex == PopSize )
			{
				DoInitEvals = false;
			}

			return( 0 );
		}

		bool DoEval = true;
		int SelMethod = -1;
		int SelDraw = -1;
		int SelM3 = -1;
		SelGen1 = -1;

		static const double ScutProbs[ 2 ] = { 0.012, 0.091 };
		int SelScut = ScutHist.select( rnd );
		ScutCntr += ScutProbs[ SelScut & 1 ];

		if( ScutCntr >= 1.0 )
		{
			ScutCntr -= 1.0;

			// Parameter value short-cuts, they considerably reduce
			// convergence time for some functions while not severely
			// impacting performance for other functions.
			//
			// Reuses any previously generated "mp" value.

			i = (int) ( rnd.getUniformRaw() * ParamCountRnd );

			if(( SelScut >> 1 ) == 0 )
			{
				// "Centroid offset" short-cut.

				const int si = (int) ( mp2 * mp2 * PopSize );
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

				const int si = (int) ( mp * mp2 * PopSize );
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
			SelScut = -1;
			SelDraw = DrawHist.select( rnd );

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

			mp = rnd.getRndValue();
			mp2 = mp * mp;
			mpi = (int) ( mp * mp2 * 4 );

			if( SelMethod == 0 )
			{
				ParOpt.optimize( rnd, &NewCost, Params );
				DoEval = false;
			}
			else
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
				SelM3 = D3Hist.select( rnd );

				if( SelM3 == 0 )
				{
					generateSol3( rnd );
				}
				else
				if( SelM3 == 1 )
				{
					generateSol4( rnd );
				}
				else
				{
					generateSol5( rnd );
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

		}

		if( PushOpt != NULL && PushOpt != this &&
			!PushOpt -> DoInitEvals )
		{
			PushOpt -> updatePop( NewCost, Params );
		}

		const int sH = PopOrder[ PopSize1 ];

		if( NewCost >= CurCosts[ sH ])
		{
			// Upper bound cost constraint check failed, reject this solution.

			if( SelMethod >= 0 )
			{
				PrevSelMethod = MethodHist.selectRandom( rnd );
			}

			ScutHist.decr( SelScut );
			MethodHist.decr( SelMethod );
			DrawHist.decr( SelDraw );
			D3Hist.decr( SelM3 );
			Gen1Hist.decr( SelGen1 );

			StallCount++;
		}
		else
		{
			if( SelMethod >= 0 )
			{
				PrevSelMethod = SelMethod;
			}

			if( DoEval )
			{
				updateBestCost( NewCost, NewParams );
			}
			else
			{
				updateBestCost( NewCost, Params, true );
			}

			updatePop( NewCost, Params, sH );

			ScutHist.incr( SelScut );
			MethodHist.incr( SelMethod );
			DrawHist.incr( SelDraw );
			D3Hist.incr( SelM3 );
			Gen1Hist.incr( SelGen1 );

			StallCount = 0;
		}

		return( StallCount );
	}

protected:
	static const int MantSize = 54; ///< Mantissa size of the bitmask
		///< operations.
		///<
	static const uint64_t MantSizeMask = ( 1ULL << MantSize ) - 1; ///< Mask
		///< that corresponds to mantissa.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double MantMultI; ///< =1/MantMult.
		///<
	double ScutCntr; ///< Short-cut probability counter.
		///<
	double AllpCntr; ///< All-parameter randomization probability counter.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double AllpProbDamp; ///< AllpProbs damping. Applied for higher dimensions
		///< as the "all parameter" randomization is ineffective in higher
		///< dimensions.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<
	uint64_t* IntParams; ///< Temporary integer value parameter buffer.
		///<
	CBiteOptHist< 4, 8, 2 > ScutHist; ///< Short-cut method's histogram.
		///<
	CBiteOptHist< 4, 4, 2 > MethodHist; ///< Population generator method
		///< histogram.
		///<
	CBiteOptHist< 3, 3, 1 > DrawHist; ///< Method draw histogram.
		///<
	CBiteOptHist< 3, 4, 2 > D3Hist; ///< Draw method 3's histogram.
		///<
	CBiteOptHist< 8, 16, 1 > Gen1Hist; ///< Generator method 1's histogram.
		///<
	int PrevSelMethod; ///< Previously successfully used method; contains
		///< random method index if optimization was not successful.
		///<
	int SelGen1; ///< Generator method 1's selector (temporary).
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
		const double* const MinParams = getParamsOrdered( mpi );
		memcpy( Params, MinParams, ParamCount * sizeof( MinParams[ 0 ]));

		// Select a single random parameter or all parameters for further
		// operations.

		int i;
		int a;
		int b;

		SelGen1 = Gen1Hist.select( rnd );

		static const double AllpProbs[ 2 ] = { 0.42, 0.81 };
		AllpCntr += AllpProbDamp * AllpProbs[ SelGen1 & 1 ];

		if( AllpCntr >= 1.0 )
		{
			AllpCntr -= 1.0;
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

		static const double Probs[ 2 ] = { 0.61, 0.88 };
		CentCntr += Probs[( SelGen1 >> 1 ) & 1 ];

		if( CentCntr >= 1.0 )
		{
			CentCntr -= 1.0;

			// Random move around random previous solution vector.

			static const double SpanMults[ 2 ] = { 1.5, 2.5 };
			const double m = SpanMults[ SelGen1 >> 2 ] * rnd.getRawScaleInv();
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
		const double* const MinParams = getParamsOrdered( mpi );
		const int si1 = (int) ( mp2 * PopSize );
		const double* const rp1 = getParamsOrdered( si1 );
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const double m1 = getBit( rnd );
			const double m2 = 1.0 - m1;

			Params[ i ] = CentParams[ i ] * m1 +
				( MinParams[ i ] + ( MinParams[ i ] - rp1[ i ])) * m2;
		}
	}

	/**
	 * "Entropy bit mixing"-based solution generator. Performs crossing-over
	 * of an odd (important) number of random solutions via XOR operation.
	 * Slightly less effective than the DE-based mixing, but makes the
	 * optimization method more diverse overall.
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
	 */

	void generateSol5( CBiteRnd& rnd )
	{
		const double* const CrossParams1 =
			getParamsOrdered( (int) ( mp2 * PopSize ));

		const double p2 = rnd.getRndValue();
		const double* const CrossParams2 =
			getParamsOrdered( (int) ( p2 * p2 * PopSize ));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// Produce a random bit mixing mask.

			const uint64_t crpl = ( rnd.getUniformRaw() |
				( (uint64_t) rnd.getUniformRaw() << 30 ));

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
		return( Opts[ BestOpt ] -> getBestParams() );
	}

	virtual double getBestCost() const
	{
		return( Opts[ BestOpt ] -> getBestCost() );
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

	void updateDims( const int aParamCount, const int M = 16,
		const int PopSize0 = 0 )
	{
		if( aParamCount == ParamCount && M == OptCount )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		OptCount = M;
		OptCountRnd = (double) OptCount / CBiteRnd :: getRawScale();
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

			if( i == 0 || Opts[ i ] -> getBestCost() <
				Opts[ BestOpt ] -> getBestCost() )
			{
				BestOpt = i;
			}
		}

		CurOpt = 0;
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
		}
		else
		{
			int PushOpt;

			if( OptCount == 2 )
			{
				PushOpt = ( CurOpt + 1 ) & 1;
			}
			else
			{
				while( true )
				{
					PushOpt = (int) ( rnd.getUniformRaw() * OptCountRnd );

					if( PushOpt != CurOpt )
					{
						break;
					}
				}
			}

			const int sc = Opts[ CurOpt ] -> optimize( rnd, Opts[ PushOpt ]);

			if( Opts[ CurOpt ] -> getBestCost() <
				Opts[ BestOpt ] -> getBestCost() )
			{
				BestOpt = CurOpt;
			}

			StallCount = ( sc == 0 ? 0 : StallCount + 1 );

			CurOpt++;

			if( CurOpt == OptCount )
			{
				CurOpt = 0;
			}
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
	double OptCountRnd; ///< Multiplier used to scale "raw" random value to
		///< obtain random Opts index.
		///<
	CBiteOptWrap** Opts; ///< Optimization objects.
		///<
	int BestOpt; ///< Optimizer that contains the best solution.
		///<
	int CurOpt; ///< Current optimization object index.
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

		if( k == 0 || opt.getBestCost() < *minf )
		{
			for( i = 0; i < N; i++ )
			{
				x[ i ] = opt.getBestParams()[ i ];
			}

			*minf = opt.getBestCost();
		}
	}
}

#endif // BITEOPT_INCLUDED
