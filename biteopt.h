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
 * @version 2021.3
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
	double AllpProb; ///< All parameters randomization probability.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double CentSpan; ///< Centroid move range multiplier.
		///<

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: MantMult( 1ULL << MantSize )
		, MantMultI( 1.0 / ( 1ULL << MantSize ))
		, IntParams( NULL )
	{
		AllpProb = 0.5;
		CentProb = 0.8;
		CentSpan = 2.5;
		ParOpt.CentPow = 6.0;
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

		InitEvals = aPopSize;
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();
		resetCentroid();
		updateDiffValues();

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		int i;
		int j;

		for( j = 0; j < PopSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double v = wrapParam( rnd,
					( j == 0 && InitParams != NULL ?
					( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ] : getGaussian( rnd ) * 0.25 + 0.5 ));

				CurParams[ j ][ i ] = v;
				CentParams[ i ] += v / PopSize;
				NewParams[ i ] = getRealValue( CurParams[ j ], i );
			}

			insertPopOrder( optcost( NewParams ), j, j );
			updateBestCost( CurCosts[ j ], NewParams );
		}

		AllpCntr = rnd.getRndValue();
		CentCntr = rnd.getRndValue();
		ParamCountRnd = (double) ParamCount / rnd.getRawScale();
		ParamCntr = (int) ( rnd.getUniformRaw() * ParamCountRnd );

		CentSpanRnd = CentSpan / rnd.getRawScale();
		AllpProbDamp = AllpProb * 2.0 / ParamCount;

		ParOpt.init( rnd );

		PrevUseMethod = (int) ( rnd.getRndValue() * MethodCount );
		mp = 0.0;
		mp2 = 0.0;
		mpi = 0;

		memset( MethodHist, 0, sizeof( MethodHist ));
		updateMethodProbs();
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
	 * threshold value is getInitEvals() * 8. When this value was reached the
	 * probability of plateau is high. This value however should not be solely
	 * relied upon when considering a stopping criteria: a hard iteration
	 * limit should be always used as in some cases convergence time may be
	 * very high with small but frequent improving steps. This value is best
	 * used to allocate iteration budget between optimization attempts more
	 * efficiently.
	 */

	int optimize( CBiteRnd& rnd, CBiteOpt* const PushOpt = NULL )
	{
		double NewCost;
		bool DoEval = true;
		int UseMethod = -1;
		int i;

		if( rnd.getRndValue() < 0.06 )
		{
			// Parameter value short-cuts, they considerably reduce
			// convergence time for some functions while not severely
			// impacting performance for other functions.
			//
			// Reuses any previously generated "mp" value.

			i = (int) ( rnd.getUniformRaw() * ParamCountRnd );

			if( getBit( rnd ))
			{
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
			if( getBit( rnd ))
			{
				UseMethod = (int) ( rnd.getRndValue() * MethodCount );
			}
			else
			if( getBit( rnd ))
			{
				UseMethod = PrevUseMethod;
			}
			else
			{
				const double mv = rnd.getRndValue() * MethodProbSum;
				UseMethod = 0;

				for( i = 1; i < MethodCount; i++ )
				{
					if( mv >= MethodProbs[ i - 1 ] && mv < MethodProbs[ i ])
					{
						UseMethod = i;
						break;
					}
				}
			}

			mp = rnd.getRndValue();
			mp2 = mp * mp;
			mpi = (int) ( mp * mp2 * 4 );

			if( UseMethod == 0 )
			{
				ParOpt.optimize( rnd, &NewCost, Params );
				DoEval = false;
			}
			else
			if( UseMethod == 1 )
			{
				generateSol1( rnd );
			}
			else
			if( UseMethod == 2 )
			{
				generateSol2( rnd );
			}
			else
			{
				if( getBit( rnd ))
				{
					generateSol3( rnd );
				}
				else
				{
					generateSol4( rnd );
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

			if( PushOpt != NULL && PushOpt != this )
			{
				PushOpt -> updatePop( NewCost, Params );
			}
		}

		const int sH = PopOrder[ PopSize1 ];
		int hincr = 0;

		if( NewCost >= CurCosts[ sH ])
		{
			// Upper bound cost constraint check failed, reject this solution.

			if( UseMethod >= 0 )
			{
				hincr = -1;
				PrevUseMethod = (int) ( rnd.getRndValue() * MethodCount );
			}

			StallCount++;
		}
		else
		{
			if( UseMethod >= 0 )
			{
				hincr = 1;
				PrevUseMethod = UseMethod;
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

			StallCount = 0;
		}

		if( hincr != 0 )
		{
			MethodHist[ UseMethod ] += hincr * 2;

			for( i = 0; i < MethodCount; i++ )
			{
				MethodHist[ i ] += hincr;
			}

			updateMethodProbs();
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
	double AllpCntr; ///< All-parameter randomization probability counter.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double AllpProbDamp; ///< Damped AllpProb value that depends on
		///< Damping is applied for higher dimensions as the "all parameter"
		///< randomization is ineffective in higher dimensions.
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	double CentSpanRnd; ///< CentSpan multiplier converted into "raw"
		///< random value scale.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<
	uint64_t* IntParams; ///< Temporary integer value parameter buffer.
		///<
	double mp, mp2; ///< Temporary variables.
		///<
	int mpi; ///< Temporary variable.
		///<
	static const int MethodCount = 4; ///< The number of solution generation
		///< methods in use.
		///< 
	int MethodHist[ MethodCount ]; ///< Methods' histogram.
		///<
	double MethodProbs[ MethodCount ]; ///< Method's probabilities,
		///< cumulative.
		///<
	double MethodProbSum; ///< Sum of method probabilities.
		///<
	int PrevUseMethod; ///< Previously successfully used method; contains
		///< random method index if optimization was not successful.
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
	 * Function updates solution generation method probabilities, based on
	 * histogram.
	 */

	void updateMethodProbs()
	{
		int MinHist = MethodHist[ 0 ];
		int i;

		for( i = 1; i < MethodCount; i++ )
		{
			if( MethodHist[ i ] < MinHist )
			{
				MinHist = MethodHist[ i ];
			}
		}

		MinHist--;
		double HistSum = 0.0;

		for( i = 0; i < MethodCount; i++ )
		{
			MethodProbs[ i ] = MethodHist[ i ] - MinHist;
			HistSum += MethodProbs[ i ];
		}

		HistSum /= MethodCount;
		MethodProbSum = 0.0;

		for( i = 0; i < MethodCount; i++ )
		{
			MethodProbs[ i ] = ( MethodProbs[ i ] < HistSum ? HistSum :
				MethodProbs[ i ]) + MethodProbSum;

			MethodProbSum = MethodProbs[ i ];
		}
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

		AllpCntr += AllpProbDamp;

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

		CentCntr += CentProb;

		if( CentCntr >= 1.0 )
		{
			CentCntr -= 1.0;

			// Random move around random previous solution vector.

			const double m1 = rnd.getTPDFRaw() * CentSpanRnd;
			const double m2 = rnd.getTPDFRaw() * CentSpanRnd;
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
	 *
	 * Not currently used.
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
			const uint64_t crpl = rnd.getUniformRaw() |
				( (uint64_t) rnd.getUniformRaw() << 30 );

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

	virtual int getInitEvals() const
	{
		int ec = 0;
		int i;

		for( i = 0; i < OptCount; i++ )
		{
			ec += Opts[ i ] -> getInitEvals();
		}

		return( ec );
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
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		int i;

		for( i = 0; i < OptCount; i++ )
		{
			Opts[ i ] -> init( rnd, InitParams );

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
	 * threshold value is getInitEvals() / M * 8.
	 */

	int optimize( CBiteRnd& rnd )
	{
		if( OptCount == 1 )
		{
			StallCount = Opts[ 0 ] -> optimize( rnd );
		}
		else
		{
			int NextOpt;

			if( OptCount < 3 )
			{
				NextOpt = ( CurOpt + 1 ) % 2;
			}
			else
			{
				while( true )
				{
					NextOpt = (int) ( rnd.getRndValue() * OptCount );

					if( NextOpt != CurOpt )
					{
						break;
					}
				}
			}

			const int sc = Opts[ CurOpt ] -> optimize( rnd, Opts[ NextOpt ]);

			if( Opts[ CurOpt ] -> getBestCost() <
				Opts[ BestOpt ] -> getBestCost() )
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
				CurOpt = NextOpt;
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

	const int useiter = (int) ( iter * sqrt( (double) M )) -
		opt.getInitEvals();

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
