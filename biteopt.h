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
 * @version 2021.2
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
	double RandProb[ 2 ]; ///< Parameter value randomization probability.
		///<
	double RandProb2[ 2 ]; ///< Alt parameter value randomization probability.
		///<
	double AllpProb[ 2 ]; ///< All parameters randomization probability.
		///<
	double CentProb[ 2 ]; ///< Centroid move probability.
		///<
	double CentSpan[ 2 ]; ///< Centroid move range multiplier.
		///<
	double ScutProb; ///< Short-cut probability.
		///<
	double MantSizeSh; ///< MantSize used in "bitmask inversion" shift, used
		///< to shrink or extend the range.
		///<
	double MantSizeSh2; ///< MantSize used in "bitmask inversion" shift, used
		///< to shrink or extend the range, for second operand.
		///<
	double PopSizeBase; ///< Minimal population size.
		///<
	double PopSizeMult; ///< Dimensional population size multiplier.
		///<
	double ParOptProb[ 2 ]; ///< Parallel optimizer's engagement probablity.
		///<
	double EntmProb[ 2 ]; ///< "Entropy bit mixing" method probability.
		///<

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: MantMult( 1LL << MantSize )
		, MantMultI( 1.0 / ( 1LL << MantSize ))
		, TmpParams( NULL )
	{
		// Cost=2.340515
		RandProb[ 0 ] = 0.31183247;
		RandProb[ 1 ] = 0.96389921;
		RandProb2[ 0 ] = 0.37889711;
		RandProb2[ 1 ] = 0.32030958;
		AllpProb[ 0 ] = 0.45959446;
		AllpProb[ 1 ] = 0.98226970;
		CentProb[ 0 ] = 0.95205787;
		CentProb[ 1 ] = 0.19277299;
		CentSpan[ 0 ] = 2.96375815;
		CentSpan[ 1 ] = 0.18697579;
		ScutProb = 0.06575653;
		MantSizeSh = 20.24787030;
		MantSizeSh2 = 76.14566404;
		PopSizeBase = 11.45111332;
		PopSizeMult = 1.98384173;
		ParOptProb[ 0 ] = 0.06342471;
		ParOptProb[ 1 ] = 0.31881036;
		EntmProb[ 0 ] = 0.16989104;
		EntmProb[ 1 ] = 0.56884090;
		ParOpt.CentPow = 9.94313896;
		ParOpt.RadPow = 29.85162390;
	}

	/**
	 * Function returns pointer to the parallel optimizer, used to optimize
	 * its hyper-parameters.
	 */

	CSpherOpt* getParOpt()
	{
		return( &ParOpt );
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
			(int) ( PopSizeBase + aParamCount * PopSizeMult ));

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();
		initBaseBuffers( aParamCount, aPopSize );

		Params = CurParams[ aPopSize ];
		TmpParams = new uint64_t[ aParamCount ];

		ParOpt.Owner = this;
		ParOpt.updateDims( ParamCount );

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

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}

		resetCommonVars();
		resetCentroid();

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		int j;

		for( j = 0; j < PopSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double r = pow( fabs( rnd.getRndValue() -
					rnd.getRndValue() ), 0.125 ) *
					( getBit( rnd ) ? -0.5 : 0.5 ) + 0.5;

				const double v = wrapParam( rnd,
					( j == 0 && InitParams != NULL ?
					( InitParams[ i ] - MinValues[ i ]) / DiffValues[ i ] :
					r ));

				CurParams[ j ][ i ] = v;
				CentParams[ i ] += v / PopSize;
				NewParams[ i ] = getRealValue( v, i );
			}

			insertPopOrder( optcost( NewParams ), j, j );

			if( j == 0 || CurCosts[ j ] < BestCost )
			{
				BestCost = CurCosts[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					BestParams[ i ] = NewParams[ i ];
				}
			}
		}

		RandCntr = rnd.getRndValue();
		RandCntr2 = rnd.getRndValue();
		AllpCntr = rnd.getRndValue();
		CentCntr = rnd.getRndValue();
		EntmCntr = rnd.getRndValue();
		ParamCountRnd = (double) ParamCount / rnd.getRawScale();
		ParamCntr = (int) ( rnd.getUniformRaw() * ParamCountRnd );
		RandSwitch = 0;

		CentSpanRnd[ 0 ] = CentSpan[ 0 ] / rnd.getRawScale();
		CentSpanRnd[ 1 ] = CentSpan[ 1 ] / rnd.getRawScale();
		AllpProbDamp = 2.0 / ParamCount;

		ParOpt.init( rnd );
		ParOptProbM = 1;
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
		int i;

		if( rnd.getRndValue() <
			ParOptProbM * ParOptProb[( RandSwitch >> 6 ) & 1 ])
		{
			double NewCost;
			ParOpt.optimize( rnd, &NewCost, Params );

			const int sH = PopOrder[ PopSize1 ];

			if( NewCost >= CurCosts[ sH ])
			{
				RandSwitch &= ~64;

				if( ParOptProbM > 1 )
				{
					ParOptProbM--;
				}

				StallCount++;
			}
			else
			{
				RandSwitch |= 64;

				if( NewCost < BestCost )
				{
					// Record the best solution.

					for( i = 0; i < ParamCount; i++ )
					{
						BestParams[ i ] = getRealValue( Params[ i ], i );
					}

					BestCost = NewCost;
				}

				updatePop( NewCost, Params, sH );
				StallCount = 0;
			}

			return( StallCount );
		}

		// Random selection between best solutions, reduces sensitivity to
		// noise.

		const double mp = rnd.getRndValue(); // Also reused later.
		const double mp2 = mp * mp; // Used later.
		const double mp3 = mp2 * mp;
		const int mpi = (int) ( mp3 * 4 );
		const double* const MinParams = getParamsOrdered( mpi );

		int UseRandSwitch = RandSwitch; // RandSwitch to use next.
		int RaiseFlags = 0; // RandSwitch flags to raise on optimization
			// improvement.

		RandCntr += RandProb[ RandSwitch & 1 ];

		if( RandCntr >= 1.0 )
		{
			RaiseFlags |= 1;
			RandCntr -= 1.0;

			RandCntr2 += RandProb2[( RandSwitch >> 1 ) & 1 ];

			if( RandCntr2 >= 1.0 )
			{
				RaiseFlags |= 2;
				RandCntr2 -= 1.0;

				// Alternative randomization method, works well for convex
				// functions. Use the very best solution and a random previous
				// solution. "mp*mp" is equivalent of giving more weight to
				// better solutions.

				const int si1 = (int) ( mp2 * PopSize );
				const double* const rp1 = getParamsOrdered( si1 );

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = ( getBit( rnd ) ? CentParams[ i ] :
						MinParams[ i ] + ( MinParams[ i ] - rp1[ i ]));
				}
			}
			else
			{
				memcpy( Params, MinParams, ParamCount * sizeof( Params[ 0 ]));

				// Select a single random parameter or all parameters for
				// further operations.

				int a;
				int b;

				AllpCntr += AllpProb[( RandSwitch >> 2 ) & 1 ] * AllpProbDamp;
					// Apply probability damping for higher dimensions as
					// "all parameter" randomization is ineffective in higher
					// dimensions.

				if( AllpCntr >= 1.0 )
				{
					RaiseFlags |= 4;
					AllpCntr -= 1.0;
					a = 0;
					b = ParamCount - 1;
				}
				else
				{
					a = ParamCntr;
					b = ParamCntr;
					ParamCntr = ( ParamCntr == 0 ?
						ParamCount : ParamCntr ) - 1;
				}

				// Bitmask inversion operation, works as a "driver" of
				// optimization process.

				const int imasks = (int) ( mp2 * mp2 * MantSizeSh );
				const uint64_t imask =
					( imasks > 63 ? 0 : MantSizeMask >> imasks );

				const double rr = rnd.getRndValue();
				const int imask2s = (int) ( rr * rr * MantSizeSh2 );
				const uint64_t imask2 =
					( imask2s > 63 ? 0 : MantSizeMask >> imask2s );

				const int si0 = (int) ( mp3 * PopSize );
				const double* const rp0 = getParamsOrdered( si0 );

				for( i = a; i <= b; i++ )
				{
					const uint64_t v1 = (uint64_t) ( Params[ i ] * MantMult );
					const uint64_t v2 = (uint64_t) ( rp0[ i ] * MantMult );
					uint64_t v0 = (( v1 ^ imask ) + ( v2 ^ imask2 )) >> 1;
					Params[ i ] = v0 * MantMultI;
				}

				const int ci = ( RandSwitch >> 3 ) & 1;
				CentCntr += CentProb[ ci ];

				if( CentCntr >= 1.0 )
				{
					RaiseFlags |= 8;
					CentCntr -= 1.0;

					// Random move around random previous solution vector.

					const double m1 = rnd.getTPDFRaw() * CentSpanRnd[ ci ];
					const double m2 = rnd.getTPDFRaw() * CentSpanRnd[ ci ];
					const int si = (int) ( mp2 * PopSize );
					const double* const rp1 = getParamsOrdered( si );

					for( i = a; i <= b; i++ )
					{
						Params[ i ] -= ( Params[ i ] - rp1[ i ]) * m1;
						Params[ i ] -= ( Params[ i ] - rp1[ i ]) * m2;
					}
				}
			}
		}
		else
		{
			EntmCntr += EntmProb[( RandSwitch >> 5 ) & 1 ];

			if( EntmCntr >= 1.0 )
			{
				RaiseFlags |= 32;
				EntmCntr -= 1.0;

				// "Entropy bit mixing"-based population cross-over. Slightly
				// less effective than the DE-based mixing below, but makes
				// the optimization method more diverse overall.

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
							TmpParams[ i ] =
								(uint64_t) ( rp[ i ] * MantMult );
						}
					}
					else
					{
						for( i = 0; i < ParamCount; i++ )
						{
							TmpParams[ i ] ^=
								(uint64_t) ( rp[ i ] * MantMult );
						}
					}
				}

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = TmpParams[ i ] * MantMultI;
				}
			}
			else
			{
				// Select worst and a random previous solution from the
				// ordered list, apply offsets to reduce sensitivity to noise.

				const int si0 = mpi + (int) ( mp * ( PopSize1 - mpi ));
				const double* const OrigParams = getParamsOrdered( si0 );
				const double* const MaxParams = getParamsOrdered(
					PopSize1 - mpi );

				// Select two more previous solutions to be used in the mix.

				const double r = rnd.getRndValue();
				const int si1 = (int) ( r * r * PopSize );
				const double* const rp1 = getParamsOrdered( si1 );
				const double* const rp2 = getParamsOrdered( PopSize1 - si1 );

				for( i = 0; i < ParamCount; i++ )
				{
					// The "step in the right direction" (Differential
					// Evolution "mutation") operation towards the best
					// (minimal) and away from the worst (maximal) parameter
					// vector, plus a difference of two random vectors.

					Params[ i ] = MinParams[ i ] -
						(( MaxParams[ i ] - OrigParams[ i ]) -
						( rp1[ i ] - rp2[ i ])) * 0.5;
				}
			}
		}

		if( rnd.getRndValue() < ScutProb )
		{
			UseRandSwitch = 0;
			RaiseFlags = 16;

			// Low-probability parameter value short-cuts, they considerably
			// reduce convergence time for some functions while not severely
			// impacting performance for other functions.

			i = (int) ( rnd.getUniformRaw() * ParamCountRnd );
			const double v = getRealValue( Params[ i ], i );

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = ( v - MinValues[ i ]) / DiffValues[ i ];
			}
		}

		// Wrap parameter values so that they stay in the [0; 1] range.

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = wrapParam( rnd, Params[ i ]);
			NewParams[ i ] = getRealValue( Params[ i ], i );
		}

		// Evaluate objective function with new parameters.

		const double NewCost = optcost( NewParams );

		if( PushOpt != NULL && PushOpt != this )
		{
			PushOpt -> pushParams( NewCost, Params, 0 );
		}

		const int sH = PopOrder[ PopSize1 ];

		if( NewCost >= CurCosts[ sH ])
		{
			// Upper bound cost constraint check failed, reject this solution.

			if( RaiseFlags != 0 )
			{
				RandSwitch = 0;
			}
			else
			{
				RandSwitch = UseRandSwitch | 1; // Raise RandProb flag. This
					// is not critically important, but introduces a kind of
					// "order" useful when optimizing method's
					// hyper-parameters.
			}

			if( ParOptProbM < 3 )
			{
				ParOptProbM++;
			}

			StallCount++;

			return( StallCount );
		}

		RandSwitch = UseRandSwitch | RaiseFlags;

		if( NewCost < BestCost )
		{
			// Record the best solution.

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;
		}

		updatePop( NewCost, Params, sH );

		StallCount = 0;

		return( StallCount );
	}

protected:
	static const int MantSize = 54; ///< Mantissa size of the bitmask
		///< operations.
		///<
	static const uint64_t MantSizeMask = ( 1LL << MantSize ) - 1; ///< Mask
		///< that corresponds to mantissa.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double MantMultI; ///< =1/MantMult.
		///<
	double RandCntr; ///< Parameter randomization probability counter.
		///<
	double RandCntr2; ///< Alt parameter randomization probability counter.
		///<
	double AllpCntr; ///< All-parameter randomization probability counter.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double EntmCntr; ///< Entropy mixing method probability counter.
		///<
	double AllpProbDamp; ///< AllpProb damping that depends on ParamCount.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	int RandSwitch; ///< Index flags for probability values switching.
		///< State automata-like.
		///<
	double CentSpanRnd[ 2 ]; ///< CentSpan multiplier converted into "raw"
		///< random value scale.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	uint64_t* TmpParams; ///< Temporary integer value parameter buffer.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] TmpParams;
	}

	/**
	 * Function "pushes" externally-provided parameters to *this object.
	 *
	 * @param NewCost Cost of the solution being pushed.
	 * @param PushParams Parameter vector being "pushed".
	 * @param NewRandSwitch New "RandSwitch" value.
	 */

	void pushParams( const double NewCost, const double* const PushParams,
		const int NewRandSwitch )
	{
		const int sH = PopOrder[ PopSize1 ];

		if( NewCost < CurCosts[ sH ])
		{
			double* const rp = CurParams[ sH ];
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] += ( PushParams[ i ] - rp[ i ]) * PopSizeI;
				rp[ i ] = PushParams[ i ];
			}

			insertPopOrder( NewCost, sH, PopSize1 );
			RandSwitch = NewRandSwitch;
		}
	}

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
	int ParOptProbM; ///< Parallel optimizer's engagement probablity
		///< multiplier.
		///<
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
