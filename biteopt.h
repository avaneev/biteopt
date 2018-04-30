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
 * BiteOpt stochastic optimization class. Implements a stochastic non-linear
 * bound-constrained derivative-free optimization strategy. It maintains a
 * cost-ordered population list of previously evaluated solutions that are
 * evolved towards a lower cost. On every iteration, the highest-cost solution
 * in the list can be replaced with a new solution, and the list reordered. A
 * population of solutions allows the strategy to space solution vectors apart
 * from each other thus making them cover a larger parameter search space
 * collectively. Beside that, parameter randomization and the "step in the
 * right direction" operation are used that move the solutions into position
 * with a probabilistically lower objective function value.
 *
 * The benefit of this strategy is a relatively high robustness: it can
 * successfully optimize a wide range of multi-dimensional test functions.
 * Another benefit is a low convergence time which depends on the complexity
 * of the objective function. Like many stochastic optimization strategies
 * with fast convergence, this strategy can't solve problems with narrow or
 * rogue optimums. Hard (multi-modal) problems may require many optimization
 * attempts to reach optimum. The name "BiteOpt" is an acronym for "BITmask
 * Evolution OPTimization".
 *
 * Algorithm's description is available at https://github.com/avaneev/biteopt
 */

class CBiteOpt
{
public:
	double MinxMult; ///< Minimal/maximal solution move range multiplier.
		///<
	double RandProb; ///< Parameter value randomization probability.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double CentSpan; ///< Centroid move range multiplier.
		///<
	double AllpProb; ///< All parameters randomization probability.
		///<
	double ScutProb; ///< Short-cut probability.
		///<
	double RandProb2; ///< Alt parameter value randomization probability.
		///<

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: MantMult( 1 << MantSize )
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
		, BestAuxData( NULL )
		, CurAuxData( NULL )
		, Params( NULL )
		, NewParams( NULL )
	{
		MinxMult = 0.50000000;
		RandProb = 0.60981671;
		CentProb = 0.87995568;
		CentSpan = 2.94153018;
		AllpProb = 0.20000000;
		ScutProb = 0.11000000;
		RandProb2 = 0.25000000;
	}

	~CBiteOpt()
	{
		deleteAuxData( CurAuxData );
		deleteAuxData( BestAuxData );
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
			12 + aParamCount * 2 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		PopSize = aPopSize;
		PopSize1 = aPopSize - 1;
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
				const double v = ( j == 0 && InitParams != NULL ?
					wrapParam( rnd, ( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ]) : rnd.getRndValue() );

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

				if( BestAuxData != NULL )
				{
					deleteAuxData( BestAuxData );
				}

				BestAuxData = CurAuxData;
				CurAuxData = NULL;
			}
		}

		CentCntr = rnd.getRndValue();
		RandCntr = rnd.getRndValue();
		ParamCntr = (int) ( rnd.getRndValue() * ParamCount );
		RandCntr2 = rnd.getRndValue();
		StallCount = 0;
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
		const int mpi = (int) ( mp * mp * mp * 4 );
		const double* const MinParams = CurParams[ PopOrder[ mpi ]];

		int i;

		RandCntr += RandProb;

		if( RandCntr >= 1.0 )
		{
			RandCntr -= 1.0;

			RandCntr2 += RandProb2;

			if( RandCntr2 >= 1.0 )
			{
				RandCntr2 -= 1.0;

				// Alternative randomization method, works well for convex
				// functions.

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = MinParams[ i ] +
						( CentParams[ i ] - MinParams[ i ]) *
						( rnd.getRndValue() < 0.5 ? 1.0 : -1.0 );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = MinParams[ i ];
				}

				// Select a single random parameter or all parameters.

				int a;
				int b;

				if( rnd.getRndValue() < AllpProb )
				{
					a = 0;
					b = ParamCount - 1;
				}
				else
				{
					a = ParamCntr;
					b = a;
					ParamCntr = ( ParamCntr == 0 ?
						ParamCount : ParamCntr ) - 1;
				}

				// Bitmask inversion operation, works as a "driver" of
				// optimization process, applied to 1 random parameter at a
				// time.

				const double r = mp;
				const int imask = ( 2 <<
					(int) (( 0.999999997 - r * r * r * r ) * MantSize )) - 1;

				for( i = a; i <= b; i++ )
				{
					Params[ i ] = ( (int) ( Params[ i ] * MantMult ) ^
						imask ) / MantMult;
				}

				CentCntr += CentProb;

				if( CentCntr >= 1.0 )
				{
					CentCntr -= 1.0;

					// Move towards centroid vector or beyond it, randomly.

					const double m = rnd.getRndValue() * CentSpan;
					const double m2 = rnd.getRndValue() * CentSpan;

					for( i = a; i <= b; i++ )
					{
						Params[ i ] -= ( Params[ i ] - CentParams[ i ]) * m;
						Params[ i ] -= ( Params[ i ] - CentParams[ i ]) * m2;
					}
				}
			}
		}
		else
		{
			// Select a random previous solution from the ordered list,
			// apply offset to reduce sensitivity to noise.

			const double r = mp;
			const int op = (int) ( r * r * 3 );
			const int si = mpi + (int) ( mp * ( PopSize1 - op - mpi ));
			const double* const OrigParams = CurParams[ PopOrder[ si ]];
			const double* const MaxParams = CurParams[ PopOrder[
				PopSize1 - op ]];

			for( i = 0; i < ParamCount; i++ )
			{
				// The "step in the right direction" operation towards the
				// best (minimal) and away from the worst (maximal) parameter
				// vector.

				Params[ i ] = MinParams[ i ] -
					( MaxParams[ i ] - OrigParams[ i ]) * MinxMult;
			}
		}

		if( mp < ScutProb )
		{
			// Low-probability parameter value short-cuts, they considerably
			// reduce convergence time for some functions while not severely
			// impacting performance for other functions.

			i = (int) ( rnd.getRndValue() * ParamCount );
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

		if( PushOpt != NULL && PushOpt != this &&
			PushOpt -> ParamCount == ParamCount )
		{
			const int sH = PushOpt -> PopOrder[ PushOpt -> PopSize1 ];

			if( NewCost <= PushOpt -> CurCosts[ sH ])
			{
				double* const rp = PushOpt -> CurParams[ sH ];

				for( i = 0; i < ParamCount; i++ )
				{
					PushOpt -> CentParams[ i ] +=
						( Params[ i ] - rp[ i ]) / PushOpt -> PopSize;

					rp[ i ] = Params[ i ];
				}

				PushOpt -> insertPopOrder( NewCost, sH, PushOpt -> PopSize1 );
			}
		}

		const int sH = PopOrder[ PopSize1 ];

		if( NewCost >= CurCosts[ sH ])
		{
			// Upper bound cost constraint check failed, reject this solution.

			StallCount++;

			return( StallCount );
		}

		if( NewCost < BestCost )
		{
			// Record the best solution.

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;

			if( BestAuxData != NULL )
			{
				deleteAuxData( BestAuxData );
			}

			BestAuxData = CurAuxData;
			CurAuxData = NULL;
		}

		// Replace the highest-cost previous solution, update centroid.

		double* const rp = CurParams[ sH ];

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( Params[ i ] - rp[ i ]) / PopSize;
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
	 * Function returns pointer to aux data provided with the best solution.
	 * This pointer should not be deleted.
	 */

	void* getBestAuxData() const
	{
		return( BestAuxData );
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

	/**
	 * Function deletes auxiliary data previously provided with the solution.
	 *
	 * @param p Pointer to aux data, can be NULL. Should be type-casted to the
	 * actual data type.
	 */

	virtual void deleteAuxData( void* const p )
	{
	}

protected:
	static const int MantSize = 29; ///< Mantissa size of bitmask inversion
		///< operation. Must be lower than the random number generator's
		///< precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int PopSize; ///< The size of population in use.
		///<
	int PopSize1; ///< = PopSize - 1.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double RandCntr; ///< Parameter randomization probability counter.
		///<
	int ParamCntr; ///< Parameter randomization index counter.
		///<
	double RandCntr2; ///< Alt parameter randomization probability counter.
		///<
	int StallCount; ///< The number of iterations without improvement.
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
	void* BestAuxData; ///< Auxiliary data of the best parameter vector.
		///<
	void* CurAuxData; ///< Aux data provided with the latest solution, will be
		///< deleted if new solution is not the best solution.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<
	double* NewParams; ///< Temporary new parameter buffer.
		///<

	/**
	 * Function provides auxiliary data, this function should be called in the
	 * optcost() function. The deleteAuxData() should be implemented to delete
	 * unused aux data.
	 *
	 * @param AuxData Pointer to aux data. Can be NULL.
	 */

	void provideAuxData( void* const AuxData )
	{
		deleteAuxData( CurAuxData );
		CurAuxData = AuxData;
	}

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
 * Deep stochastic optimization class. Based on an array of M CBiteOpt
 * objects. This "deep" strategy pushes the newly-obtained solution to the
 * next CBiteOpt object which is then optimized. This strategy while
 * increasing the convergence time by a factor of about sqrt(M) is able to
 * solve even the most noisy non-linear functions.
 *
 * This strategy is most effective on stochastic functions or functions with
 * huge fluctuations near the global solution that are not very expensive to
 * calculate and that have a large iteration budget. Tests have shown that on
 * smooth functions that have many strongly competing minima this strategy
 * increases the chance to find a global solution by a factor of sqrt(M)
 * relative to the CBiteOpt class, but still requires several runs at
 * different random seeds. However, the number of required runs in most cases
 * is lower by about sqrt(M). So, in many cases it's more efficient to
 * increase the iteration budget by a factor of sqrt(M) which will in turn
 * increase the chance to find a global optimum by a factor of sqrt(M) and
 * also reduce the number of optimization attempts by the same number. In
 * practice, the chance to find a global optimum is increased even more
 * considerably with this "deep" strategy.
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
	 */

	void updateDims( const int aParamCount, const int M = 16 )
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
			Opts[ i ] -> updateDims( aParamCount );
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
	 * threshold value is getInitEvals() * 8.
	 */

	int optimize( CBiteRnd& rnd )
	{
		const int NextOpt = ( CurOpt == BiteCount - 1 ? 0 : CurOpt + 1 );

		if( Opts[ CurOpt ] -> optimize( rnd, Opts[ NextOpt ]) == 0 )
		{
			StallCount = 0;
		}
		else
		{
			StallCount++;
		}

		if( Opts[ CurOpt ] -> getBestCost() <
			Opts[ BestOpt ] -> getBestCost() )
		{
			BestOpt = CurOpt;
		}

		CurOpt = NextOpt;

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
