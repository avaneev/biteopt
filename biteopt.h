//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBiteOpt class.
 *
 * @section license License
 *
 * Copyright (c) 2016-2017 Aleksey Vaneev
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
 * cost-ordered population of previously evaluated solutions that are evolved
 * towards a lower cost. On every iteration, the highest-cost solution in the
 * list is unconditionally replaced with a new solution, and the list is
 * reordered. A population of solutions allows the strategy to space solution
 * vectors apart from each other thus making them cover a larger parameter
 * search space collectively. Beside that, parameter randomization and the
 * "step in the right direction" operation are used that move the solutions
 * into position with a probabilistically lower objective function value.
 *
 * The benefit of this strategy is a relatively high robustness: it can
 * successfully optimize a wide range of 1-10 dimensional test functions.
 * Another benefit is a low convergence time which depends on the complexity
 * of the objective function. Like many stochastic optimization strategies
 * with fast convergence, this strategy can't solve problems with narrow or
 * rogue optimums. Harder problems may require dozens of optimization attempts
 * to reach optimum.
 *
 * The algorithm consists of the following elements:
 *
 * 1. A population of previous solutions is maintained. A solution is an
 * independent parameter vector which is evolved towards a better solution.
 * Also a cost-ordered list of solutions is maintaned. On every iteration
 * a single randomly-selected previous solution is evolved.
 *
 * 2. On every iteration the "step in the right direction" operation is
 * performed using the current best and worst solutions.
 *
 * 3. Depending on the RandProb probability, also an individual parameter
 * value randomization is performed using "bitmask inversion" operation.
 *
 * 4. With CentProb probability the "step in the right direction" operation is
 * performed using the centroid vector.
 *
 * 5. After each objective function evaluation, the highest-cost previous
 * solution is replaced using the cost constraint.
 */

class CBiteOpt
{
public:
	double CostMult; ///< Solution rejection cost threshold multiplier.
		///<
	double MinxMult; ///< Minimal/maximal solution move range multiplier.
		///<
	double RandProb; ///< Parameter value randomization probability.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double CentOffs; ///< Centroid move range multiplier offset.
		///<
	double CentSpan; ///< Centroid move range random multiplier.
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
		, Params( NULL )
		, NewParams( NULL )
	{
		// Cost=14.337680
		CostMult = 1.38433223;
		MinxMult = 0.50005377;
		RandProb = 0.50360804;
		CentProb = 0.29195602;
		CentOffs = 0.45708476;
		CentSpan = 1.78321594;
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
			(int) ( 18.0 + aParamCount * aParamCount / 6.0 ));

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
		// values of these solutions and find the best cost.

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
			}
		}

		CentCntr = rnd.getRndValue();
		RandCntr = rnd.getRndValue();
		ParamCntr = (int) ( rnd.getRndValue() * ParamCount );
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return "True" if optimizer's state was improved on this iteration.
	 * Many successive "false" results means optimizer has reached a plateau.
	 */

	bool optimize( CBiteRnd& rnd )
	{
		// Select a random previous solution from the ordered list.

		const int si = (int) ( rnd.getRndValue() * PopSize );
		const double* const OrigParams = CurParams[ PopOrder[ si ]];
		const double* const MinParams = CurParams[ PopOrder[ 0 ]];
		const double* const MaxParams = CurParams[ PopOrder[ PopSize1 ]];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			// The "step in the right direction" operation towards the best
			// (minimal) and away from the worst (maximal) parameter vector.

			Params[ i ] = MinParams[ i ] -
				( MaxParams[ i ] - OrigParams[ i ]) * MinxMult;
		}

		RandCntr += RandProb;

		if( RandCntr >= 1.0 )
		{
			RandCntr -= 1.0;

			// Bitmask inversion operation, works as a "driver" of
			// optimization process, applied to 1 random parameter at a time.

			const int imask =
				( 2 << (int) ( rnd.getRndValue() * MantSize )) - 1;

			Params[ ParamCntr ] = ( (int) ( Params[ ParamCntr ] * MantMult ) ^
				imask ) / MantMult;

			ParamCntr = ( ParamCntr == 0 ? ParamCount : ParamCntr ) - 1;
		}

		CentCntr += CentProb;

		if( CentCntr >= 1.0 )
		{
			CentCntr -= 1.0;

			// Move towards centroid vector or beyond it, randomly.

			for( i = 0; i < ParamCount; i++ )
			{
				const double m = CentOffs + rnd.getRndValue() * CentSpan;
				Params[ i ] -= ( Params[ i ] - CentParams[ i ]) * m;
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

		const int sH = PopOrder[ PopSize1 ];
		const double cT = CurCosts[ sH ] +
			( CurCosts[ sH ] - CurCosts[ PopOrder[ 0 ]]) * CostMult;

		if( NewCost > cT )
		{
			return( false );
		}

		if( NewCost < BestCost )
		{
			// Record the best solution.

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;
		}

		// Replace the highest-cost previous solution.

		double* const rp = CurParams[ sH ];

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( Params[ i ] - rp[ i ]) / PopSize;
			rp[ i ] = Params[ i ];
		}

		insertPopOrder( NewCost, sH, PopSize1 );

		return( true );
	}

	/**
	 * Function performs iterative optimization function until the required
	 * minimum cost is reached or the maximum number of iterations is
	 * exceeded.
	 *
	 * @param rnd Random number generator.
	 * @param MinCost Target minimum cost.
	 * @param PlateauIters The number of successive iterations which have not
	 * produced a useful change, to treat as reaching a plateau. When plateau
	 * is reached, *this object will be reinitialized (if MaxIters>0).
	 * Effective values depend on the PopSize value, lower PopSize values
	 * require higher PlateauIters values.
	 * @param MaxIters The maximal number of iterations to perform. Set to 0
	 * to run any number of iterations until plateau is reached.
	 * @return The number of iterations performed, including reinitialization
	 * calls. May be greater than MaxIters.
	 */

	int optimizePlateau( CBiteRnd& rnd, const double MinCost,
		const int PlateauIters, const int MaxIters )
	{
		double* TmpBestParams = new double[ ParamCount ];
		double TmpBestCost;
		int Iters = PopSize;
		int i;

		init( rnd );
		bool IsFirstInit = true;
		int PlateauCount = 0;

		while( true )
		{
			if( BestCost <= MinCost )
			{
				delete[] TmpBestParams;
				return( Iters );
			}

			if( Iters >= MaxIters && MaxIters > 0 )
			{
				break;
			}

			const bool WasImproved = optimize( rnd );
			Iters++;

			if( WasImproved )
			{
				PlateauCount = 0;
			}
			else
			{
				PlateauCount++;

				if( PlateauCount >= PlateauIters )
				{
					if( IsFirstInit || BestCost < TmpBestCost )
					{
						TmpBestCost = BestCost;

						for( i = 0; i < ParamCount; i++ )
						{
							TmpBestParams[ i ] = BestParams[ i ];
						}
					}

					IsFirstInit = false;

					if( MaxIters <= 0 )
					{
						break;
					}

					init( rnd );
					Iters += PopSize;
					PlateauCount = 0;
				}
			}
		}

		if( !IsFirstInit && TmpBestCost < BestCost )
		{
			BestCost = TmpBestCost;

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = TmpBestParams[ i ];
			}
		}

		delete[] TmpBestParams;
		return( Iters );
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

	virtual double optcost( const double* const p ) const = 0;

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
	 * operator. This operation increases convergence in comparison to
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

#endif // BITEOPT_INCLUDED
