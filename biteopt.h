//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBiteOpt class.
 *
 * @section license License
 *
 * Copyright (c) 2016-2018 Aleksey Vaneev
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

#include <math.h>
#include "biternd.h"

/**
 * BiteOpt optimization class. Implements a stochastic non-linear
 * bound-constrained derivative-free optimization method.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CBiteOpt
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

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: MantMult( 1LL << MantSize )
		, MantMultI( 1.0 / ( 1LL << MantSize ))
		, ParamCount( 0 )
		, PopSize( 0 )
		, PopOrder( NULL )
		, CurParamsBuf( NULL )
		, CurParams( NULL )
		, CurCosts( NULL )
		, CentParams( NULL )
		, MinValues( NULL )
		, MaxValues( NULL )
		, DiffValues( NULL )
		, BestParams( NULL )
		, Params( NULL )
		, NewParams( NULL )
	{
		// Cost=2.101823
		RandProb[ 0 ] = 0.40891214;
		RandProb[ 1 ] = 0.98786512;
		RandProb2[ 0 ] = 0.45861012;
		RandProb2[ 1 ] = 0.32664134;
		AllpProb[ 0 ] = 0.52827167;
		AllpProb[ 1 ] = 0.98183325;
		CentProb[ 0 ] = 0.94613304;
		CentProb[ 1 ] = 0.01094351;
		CentSpan[ 0 ] = 2.87005732;
		CentSpan[ 1 ] = 1.50000000;
		ScutProb = 0.09000000;
		MantSizeSh = 21.97054644;
		MantSizeSh2 = 81.69065658;
		PopSizeBase = 11.10471664;
		PopSizeMult = 1.92913120;
	}

	~CBiteOpt()
	{
		deleteBuffers();
	}

	/**
	 * Function updates dimensionality of *this object. Function does nothing
	 * if dimensionality has not changed since the last call. This function
	 * should be called at least once before calling the init() function.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0, the default formula will be used.
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

		ParamCount = aParamCount;
		PopSize = aPopSize;
		PopSize1 = aPopSize - 1;
		PopSizeI = 1.0 / aPopSize;
		PopOrder = new int[ PopSize ];
		CurParamsBuf = new double[ PopSize * ParamCount ];
		CurParams = new double*[ PopSize ];
		CurCosts = new double[ PopSize ];
		CentParams = new double[ ParamCount ];
		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		DiffValues = new double[ ParamCount ];
		BestParams = new double[ ParamCount ];
		Params = new double[ ParamCount ];
		NewParams = new double[ ParamCount ];

		int i;

		for( i = 0; i < PopSize; i++ )
		{
			CurParams[ i ] = CurParamsBuf + i * ParamCount;
		}
	}

	/**
	 * @return The number of initial objective function evaluations.
	 * Corresponds to the population size.
	 */

	int getInitEvals() const
	{
		return( PopSize );
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
			CentParams[ i ] = 0.0;
		}

		// Initialize solution vectors randomly, calculate objective function
		// values of these solutions.

		int j;

		for( j = 0; j < PopSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double r = pow( fabs( rnd.getRndValue() -
					rnd.getRndValue() ), 0.125 ) *
					( rnd.getRndValue() < 0.5 ? -0.5 : 0.5 ) + 0.5;

				const double v = ( j == 0 && InitParams != NULL ?
					wrapParam( rnd, ( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ]) : r );

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
		ParamCntr = (int) ( rnd.getRndValue() * ParamCount );
		RandSwitch = 0;
		StallCount = 0;

		CentSpanRnd[ 0 ] = CentSpan[ 0 ] / rnd.getRawScale();
		CentSpanRnd[ 1 ] = CentSpan[ 1 ] / rnd.getRawScale();
		PopSizeRnd = (double) PopSize / rnd.getRawScale();
		ParamCountRnd = (double) ParamCount / rnd.getRawScale();
		AllpProbDamp = 2.0 / ParamCount;
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
		// Random selection between best solutions, reduces sensitivity to
		// noise.

		const double mp = rnd.getRndValue(); // Also reused later.
		const double mp2 = mp * mp; // Used later.
		const int mpi = (int) ( mp * mp2 * 4 );
		const double* const MinParams = CurParams[ PopOrder[ mpi ]];

		int i;

		int RaiseFlags = 0; // Which RandSwitch flags to raise on
			// optimization improvement.

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

				const int si = (int) ( mp2 * PopSize );
				const double* const rp1 = CurParams[ PopOrder[ si ]];

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = ( rnd.getRndValue() < 0.5 ?
						CentParams[ i ] :
						MinParams[ i ] + ( MinParams[ i ] - rp1[ i ]));
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = MinParams[ i ];
				}

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

				const int64_t imask =
					MantSizeMask >> (int64_t) ( mp2 * mp2 * MantSizeSh );

				const double rr = rnd.getRndValue();
				const int64_t imask2 =
					MantSizeMask >> (int64_t) ( rr * rr * MantSizeSh2 );

				const int si = (int) ( mp * mp2 * PopSize );
				const double* const rp0 = CurParams[ PopOrder[ si ]];

				for( i = a; i <= b; i++ )
				{
					const int64_t v1 = (int64_t) ( Params[ i ] * MantMult );
					const int64_t v2 = (int64_t) ( rp0[ i ] * MantMult );
					int64_t v0 = (( v1 ^ imask ) + ( v2 ^ imask2 )) >> 1;
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
					const double* const rp1 = CurParams[ PopOrder[ si ]];

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
			// Select worst and a random previous solution from the ordered
			// list, apply offsets to reduce sensitivity to noise.

			const int op = (int) ( mp * 3 );
			const int si = mpi + (int) ( mp * ( PopSize1 - mpi ));
			const double* const OrigParams = CurParams[ PopOrder[ si ]];
			const double* const MaxParams = CurParams[ PopOrder[
				PopSize1 - op ]];

			// Select two more previous solutions to be used in the mix.

			const int si2 = (int) ( rnd.getUniformRaw() * PopSizeRnd );
			const double* const rp1 = CurParams[ PopOrder[ si2 ]];
			const double* const rp2 = CurParams[ PopOrder[ PopSize1 - si2 ]];

			for( i = 0; i < ParamCount; i++ )
			{
				// The "step in the right direction" (Differential Evolution
				// "mutation") operation towards the best (minimal) and away
				// from the worst (maximal) parameter vector, plus a
				// difference of two random vectors.

				Params[ i ] = MinParams[ i ] -
					(( MaxParams[ i ] - OrigParams[ i ]) -
					( rp1[ i ] - rp2[ i ])) * 0.5;
			}
		}

		if( mp < ScutProb )
		{
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

		if( PushOpt != NULL )
		{
			const int sH = PushOpt -> PopOrder[ PushOpt -> PopSize1 ];

			if( NewCost <= PushOpt -> CurCosts[ sH ])
			{
				double* const rp = PushOpt -> CurParams[ sH ];

				for( i = 0; i < ParamCount; i++ )
				{
					PushOpt -> CentParams[ i ] +=
						( Params[ i ] - rp[ i ]) * PushOpt -> PopSizeI;

					rp[ i ] = Params[ i ];
				}

				PushOpt -> insertPopOrder( NewCost, sH, PushOpt -> PopSize1 );
				PushOpt -> RandSwitch = RandSwitch;
			}
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
				RandSwitch |= 1; // Raise RandProb flag. This is not
					// critically important, but introduces a kind of "order"
					// useful when optimizing method's hyper-parameters.
			}

			StallCount++;

			return( StallCount );
		}

		RandSwitch |= RaiseFlags;

		if( NewCost < BestCost )
		{
			// Record the best solution.

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;
		}

		// Replace the highest-cost previous solution, update centroid.

		double* const rp = CurParams[ sH ];

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( Params[ i ] - rp[ i ]) * PopSizeI;
			rp[ i ] = Params[ i ];
		}

		insertPopOrder( NewCost, sH, PopSize1 );
		StallCount = 0;

		return( StallCount );
	}

	/**
	 * @return Best parameter vector.
	 */

	const double* getBestParams() const
	{
		return( BestParams );
	}

	/**
	 * @return Cost of the best parameter vector.
	 */

	double getBestCost() const
	{
		return( BestCost );
	}

	/**
	 * Virtual function that should fill minimal parameter value vector.
	 *
	 * @param[out] p Minimal value vector.
	 */

	virtual void getMinValues( double* const p ) const = 0;

	/**
	 * Virtual function that should fill maximal parameter value vector.
	 *
	 * @param[out] p Maximal value vector.
	 */

	virtual void getMaxValues( double* const p ) const = 0;

	/**
	 * Virtual function (objective function) that should calculate parameter
	 * vector's optimization cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) = 0;

protected:
	static const int64_t MantSize = 54; ///< Mantissa size of the "bitmask
		///< inversion" operation.
		///<
	static const int64_t MantSizeMask = ( 1LL << MantSize ) - 1; ///< Mask
		///< that corresponds to mantissa.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double MantMultI; ///< =1/MantMult.
		///<
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int PopSize; ///< The size of population in use.
		///<
	int PopSize1; ///< = PopSize - 1.
		///<
	double PopSizeI; ///< = 1/PopSize.
		///<
	double RandCntr; ///< Parameter randomization probability counter.
		///<
	double RandCntr2; ///< Alt parameter randomization probability counter.
		///<
	double AllpCntr; ///< All-parameter randomization probability counter.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double AllpProbDamp; ///< AllpProb damping that depends on ParamCount.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	int RandSwitch; ///< Index flags for probability values switching.
		///< State automata-like.
		///<
	int StallCount; ///< The number of iterations without improvement.
		///<
	double CentSpanRnd[ 2 ]; ///< CentSpan multiplier converted into "raw"
		///< random value scale.
		///<
	double PopSizeRnd; ///< PopSize converted into "raw" random value
		///< scale.
		///<
	double ParamCountRnd; ///< ParamCount converted into "raw" random value
		///< scale.
		///<
	int* PopOrder; ///< The current solution vectors ordering,
		///< ascending-sorted by cost.
		///<
	double* CurParamsBuf; ///< CurParams buffer.
		///<
	double** CurParams; ///< Current working parameter vectors.
		///<
	double* CurCosts; ///< Best costs of current working parameter vectors.
		///<
	double* CentParams; ///< Centroid of the current parameter vectors.
		///<
	double* MinValues; ///< Minimal parameter values.
		///<
	double* MaxValues; ///< Maximal parameter values.
		///<
	double* DiffValues; ///< Difference between maximal and minimal parameter
		///< values.
		///<
	double* BestParams; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<
	double* NewParams; ///< Temporary new parameter buffer.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	void deleteBuffers()
	{
		delete[] PopOrder;
		delete[] CurParamsBuf;
		delete[] CurParams;
		delete[] CurCosts;
		delete[] CentParams;
		delete[] MinValues;
		delete[] MaxValues;
		delete[] DiffValues;
		delete[] BestParams;
		delete[] Params;
		delete[] NewParams;
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0] range, by wrapping it over the boundaries using random
	 * operator. This operation improves convergence in comparison to
	 * clamping.
	 *
	 * @param v Parameter value to wrap.
	 * @return Wrapped parameter value.
	 */

	static double wrapParam( CBiteRnd& rnd, double v )
	{
		if( v < 0.0 )
		{
			if( v > -1.0 )
			{
				return( rnd.getRndValue() * -v );
			}

			return( rnd.getRndValue() );
		}

		if( v > 1.0 )
		{
			if( v < 2.0 )
			{
				return( 1.0 - rnd.getRndValue() * ( v - 1.0 ));
			}

			return( rnd.getRndValue() );
		}

		return( v );
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param Params Parameter vector of interest.
	 * @param i Parameter index.
	 */

	double getRealValue( const double v, const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] * v );
	}

	/**
	 * Function inserts the specified solution into the PopOrder
	 * array at the appropriate offset, increasing the number of items by 1.
	 *
	 * @param Cost Solution's cost.
	 * @param f Solution's index.
	 * @param ItemCount The current number of items in the array.
	 */

	void insertPopOrder( const double Cost, const int f, const int ItemCount )
	{
		CurCosts[ f ] = Cost;

		// Perform binary search.

		int z = 0;
		int hi = ItemCount;

		while( z < hi )
		{
			int mid = ( z + hi ) >> 1;

			if( CurCosts[ PopOrder[ mid ]] >= Cost )
			{
				hi = mid;
			}
			else
			{
				z = mid + 1;
			}
		}

		// Insert element at the correct sorted position.

		int i;

		for( i = ItemCount; i > z; i-- )
		{
			PopOrder[ i ] = PopOrder[ i - 1 ];
		}

		PopOrder[ z ] = f;
	}
};

/**
 * Deep optimization class. Based on an array of M CBiteOpt objects. This
 * "deep" method pushes the newly-obtained solution to the next CBiteOpt
 * object which is then optimized.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CBiteOptDeep
{
public:
	CBiteOptDeep()
		: ParamCount( 0 )
		, BiteCount( 0 )
		, Opts( NULL )
	{
	}

	~CBiteOptDeep()
	{
		deleteBuffers();
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
		if( aParamCount == ParamCount && M == BiteCount )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		BiteCount = M;
		Opts = new CBiteOptWrap*[ BiteCount ];

		int i;

		for( i = 0; i < BiteCount; i++ )
		{
			Opts[ i ] = new CBiteOptWrap( this );
			Opts[ i ] -> updateDims( aParamCount, PopSize0 );
		}
	}

	/**
	 * @return The number of initial objective function evaluations.
	 */

	int getInitEvals() const
	{
		int ec = 0;
		int i;

		for( i = 0; i < BiteCount; i++ )
		{
			ec += Opts[ i ] -> getInitEvals();
		}

		return( ec );
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

		for( i = 0; i < BiteCount; i++ )
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
		if( BiteCount == 1 )
		{
			StallCount = Opts[ 0 ] -> optimize( rnd );
		}
		else
		{
			int NextOpt = ( CurOpt == BiteCount - 1 ? 0 : CurOpt + 1 );

			if( Opts[ CurOpt ] -> optimize( rnd, Opts[ NextOpt ]) == 0 )
			{
				StallCount = 0;
			}
			else
			{
				StallCount++;
				NextOpt = ( CurOpt - 1 + (int) ( rnd.getRndValue() * 3 ) +
					BiteCount ) % BiteCount;
			}

			if( Opts[ CurOpt ] -> getBestCost() <
				Opts[ BestOpt ] -> getBestCost() )
			{
				BestOpt = CurOpt;
			}

			CurOpt = NextOpt;
		}

		return( StallCount );
	}

	/**
	 * @return Best parameter vector.
	 */

	const double* getBestParams() const
	{
		return( Opts[ BestOpt ] -> getBestParams() );
	}

	/**
	 * @return Cost of the best parameter vector.
	 */

	double getBestCost() const
	{
		return( Opts[ BestOpt ] -> getBestCost() );
	}

	/**
	 * Virtual function that should fill minimal parameter value vector.
	 *
	 * @param[out] p Minimal value vector.
	 */

	virtual void getMinValues( double* const p ) const = 0;

	/**
	 * Virtual function that should fill maximal parameter value vector.
	 *
	 * @param[out] p Maximal value vector.
	 */

	virtual void getMaxValues( double* const p ) const = 0;

	/**
	 * Virtual function (objective function) that should calculate parameter
	 * vector's optimization cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) = 0;

protected:
	/**
	 * Wrapper class for CBiteOpt class.
	 */

	class CBiteOptWrap : public CBiteOpt
	{
	public:
		CBiteOptDeep* Owner; ///< Owning object.
			///<

		CBiteOptWrap( CBiteOptDeep* const aOwner = NULL )
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
	int BiteCount; ///< The total number of optimization objects in use.
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

			for( i = 0; i < BiteCount; i++ )
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
