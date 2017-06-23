//$ nocpp

/**
 * @file bitefan.h
 *
 * @brief The inclusion file for the CBEOOptimizerFan class.
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
 * BiteOpt stochastic function optimization class. This optimization strategy
 * is based on now outdated CBEOOptimizer and CBEOOptimizer2, CBEOOptimizerFan
 * stochastic derivative-less strategies. It uses an ordered list of several
 * current parameter vectors (called "fan elements") that are evolved towards
 * a lower cost. Highest-cost "fan element" in the list can be replaced with a
 * new solution if "fan element's" cost (plus some margin) is higher than that
 * of the new solution's. Having several "fan elements" allows the strategy
 * to space them apart from each other thus making them cover a larger
 * parameter search space collectively. Beside that, parameter randomization
 * and several "step in the right direction" operations are used that move the
 * parameter vector into position with a probabilistically lower objective
 * function value.
 *
 * The benefit of this strategy is a relatively high robustness: it can
 * successfully optimize a wide range of test functions. Another benefit is a
 * low convergence time which depends on the complexity of the objective
 * function. Like many stochastic optimization strategies with fast
 * convergence, this strategy can't solve problems with narrow or rogue
 * optimums. Harder problems may require dozens of optimization attempts to
 * reach optimum.
 *
 * The algorithm consists of the following elements:
 *
 * 1. A complex of "fan elements" is maintained. A "fan element" is an
 * independent parameter vector which is evolved towards a better solution.
 * Also a cost-ordered list of "fan elements" is maintaned. On every iteration
 * a single random "fan element" is evolved.
 *
 * 2. A centroid vector of all "fan elements" is maintained.
 *
 * 3. On every iteration (not together with N.4) the "step in the right
 * direction" operation is performed using the current best solution. Also
 * the "step in the right direction" operation is performed using the current
 * worst solution.
 *
 * 4. With CentProb probability the "step in the right direction" operation is
 * performed using the centroid vector.
 *
 * 5. With RandProb probability the parameter value randomization operation is
 * performed, which is the main driver of the evolutionary process.
 *
 * 6. After each objective function evaluation, an attempt to replace the
 * highest-cost "fan element" is performed using the cost constraint. This
 * approach is based on an assumption that the later solutions tend to be
 * statistically better than the earlier solutions.
 */

class CBiteOpt
{
public:
	double CostMult; ///< "Fan element" cost threshold multiplier.
		///<
	double MinpMult; ///< Best solution's move range multiplier.
		///<
	double MaxpMult; ///< Worst solution's move range multiplier.
		///<
	double CentMult; ///< Centroid move range multiplier.
		///<
	double CentOffs; ///< Centroid move range shift.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double RandProb; ///< Parameter value randomization probability.
		///<
	double RandMult; ///< Parameter value randomization multiplier.
		///<

	/**
	 * Constructor.
	 */

	CBiteOpt()
		: CentProb( 0.25 )
		, RandProb( 0.65 )
		, ParamCount( 0 )
		, FanSize( 0 )
		, FanOrder( NULL )
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
		// Cost=6.861364 ItAvg=348.792489 ItRtAvg=0.249729
		CostMult = 0.69643934;
		MinpMult = 0.59410786;
		MaxpMult = 0.71090436;
		CentMult = 0.61069564;
		CentOffs = 0.73554560;
		RandMult = 9.39813923;
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
	 * @param FanSize0 The number of "fan elements" to use. If set to 0, the
	 * default formula will be used.
	 */

	void updateDims( const int aParamCount, const int FanSize0 = 0 )
	{
		const int aFanSize = ( FanSize0 > 0 ? FanSize0 :
			16 + aParamCount * aParamCount / 5 );

		if( aParamCount == ParamCount && aFanSize == FanSize )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		FanSize = aFanSize;
		FanSize1 = aFanSize - 1;
		FanOrder = new int[ FanSize ];
		CurParamsBuf = new double[ FanSize * ParamCount ];
		CurParams = new double*[ FanSize ];
		CurCosts = new double[ FanSize ];
		CentParams = new double[ ParamCount ];
		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		DiffValues = new double[ ParamCount ];
		BestParams = new double[ ParamCount ];
		Params = new double[ ParamCount ];
		NewParams = new double[ ParamCount ];

		int i;

		for( i = 0; i < FanSize; i++ )
		{
			CurParams[ i ] = CurParamsBuf + i * ParamCount;
		}
	}

	/**
	 * Function initializes *this optimizer. Performs N=FanSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		CentCntr = 1.0;
		RandCntr = 1.0;
		ParamCntr = 0;

		getMinValues( MinValues );
		getMaxValues( MaxValues );
		int i;
		int j;

		for( i = 0; i < ParamCount; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}

		// Initialize "fan element" parameter vectors.

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] = 0.0;
		}

		if( InitParams != NULL )
		{
			for( j = 0; j < FanSize; j++ )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const double v = ( j == 0 ?
						wrapParam( rnd, ( InitParams[ i ] - MinValues[ i ]) /
						DiffValues[ i ]) : rnd.getRndValue() );

					CurParams[ j ][ i ] = v;
					CentParams[ i ] += CurParams[ j ][ i ];
				}
			}
		}
		else
		{
			for( j = 0; j < FanSize; j++ )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					CurParams[ j ][ i ] = rnd.getRndValue();
					CentParams[ i ] += CurParams[ j ][ i ];
				}
			}
		}

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] /= FanSize;
		}

		// Calculate costs of "fan elements" and find the best cost.

		for( j = 0; j < FanSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				NewParams[ i ] = getRealValue( CurParams[ j ][ i ], i );
			}

			insertFanOrder( optcost( NewParams ), j, j );

			if( j == 0 || CurCosts[ j ] < BestCost )
			{
				BestCost = CurCosts[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					BestParams[ i ] = NewParams[ i ];
				}
			}
		}
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
		// Select a random "fan element" from the ordered list.

		const int s = FanOrder[ (int) ( rnd.getRndValue() * FanSize )];
		int i;

		const double* const OrigParams = CurParams[ s ];
		const double* const MinParams = CurParams[ FanOrder[ 0 ]];
		const double* const MaxParams = CurParams[ FanOrder[ FanSize1 ]];

		if( CentCntr >= 1.0 )
		{
			CentCntr -= 1.0;

			// Move towards centroid vector or beyond it, randomly.

			for( i = 0; i < ParamCount; i++ )
			{
				const double m = CentOffs + rnd.getRndValue() * CentMult;
				Params[ i ] = OrigParams[ i ] -
					( OrigParams[ i ] - CentParams[ i ]) * m;
			}
		}
		else
		for( i = 0; i < ParamCount; i++ )
		{
			// The "step in the right direction" operation towards the best
			// (minimal) parameter vector.

			Params[ i ] = OrigParams[ i ] -
				( OrigParams[ i ] - MinParams[ i ]) * MinpMult;

			// The "step in the right direction" operation away from the worst
			// (maximal) parameter vector.

			Params[ i ] -= ( MaxParams[ i ] - Params[ i ]) * MaxpMult;
		}

		CentCntr += CentProb;

		if( RandCntr >= 1.0 )
		{
			RandCntr -= 1.0;

			// Parameter randomization operation, works as a "driver" of
			// optimization process, applied to a random parameter.

			const double p = Params[ ParamCntr ] +
				( rnd.getRndValue() - 0.5 ) * fabs( CentParams[ ParamCntr ] -
				MinParams[ ParamCntr ]) * RandMult;

			// A very interesting approach: mix (partially replace) another
			// random parameter with previously randomized parameter. Such
			// approach probably works due to possible mutual correlation
			// between parameters, especially in multi-dimensional functions.
			// Mix constant is randomized, with adjusted probability
			// distribution.

			const int rp = (int) ( rnd.getRndValue() * ParamCount );
			double rr = rnd.getRndValue();
			Params[ rp ] -= ( Params[ rp ] - p ) * ( 1.0 - rr * rr );

			Params[ ParamCntr ] = Params[ rp ];
			ParamCntr = ( ParamCntr == 0 ? ParamCount : ParamCntr ) - 1;
		}

		RandCntr += RandProb;

		// Wrap parameter values so that they stay in the [0; 1] range.

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = wrapParam( rnd, Params[ i ]);
		}

		// Evaluate objective function with new parameters.

		for( i = 0; i < ParamCount; i++ )
		{
			NewParams[ i ] = getRealValue( Params[ i ], i );
		}

		const double NewCost = optcost( NewParams );

		// Try to replace the highest cost "fan element".

		const int sH = FanOrder[ FanSize1 ];
		const double cT = CurCosts[ sH ] +
			( CurCosts[ sH ] - CurCosts[ FanOrder[ 0 ]]) * CostMult;

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

		// Replace highest cost "fan element".

		double* const rp = CurParams[ sH ];

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( Params[ i ] - rp[ i ]) / FanSize;
			rp[ i ] = Params[ i ];
		}

		insertFanOrder( NewCost, sH, FanSize1 );

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
	 * Effective values depend on the FanSize value, lower FanSize values
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
		int Iters = FanSize;
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
					Iters += FanSize;
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
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int FanSize; ///< The number of "fan elements" in use.
		///<
	int FanSize1; ///< = FanSize - 1.
		///<
	double CentCntr; ///< Centroid move probability counter.
		///<
	double RandCntr; ///< Randomization operation probability counter.
		///<
	int ParamCntr; ///< Parameter index counter.
		///<
	int* FanOrder; ///< The current "fan element" ordering, ascending-sorted
		///< by cost.
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
		delete[] FanOrder;
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
	 * Function inserts the specified "fan element" index into the FanOrder
	 * array at the appropriate offset, increasing the number of items by 1.
	 *
	 * @param Cost "Fan element's" cost.
	 * @param f "Fan element's" index.
	 * @param ItemCount The current number of items in the array.
	 */

	void insertFanOrder( const double Cost, const int f, const int ItemCount )
	{
		CurCosts[ f ] = Cost;

		// Perform binary search.

		int z = 0;
		int hi = ItemCount;

		while( z < hi )
		{
			int mid = ( z + hi ) >> 1;

			if( CurCosts[ FanOrder[ mid ]] >= Cost )
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
			FanOrder[ i ] = FanOrder[ i - 1 ];
		}

		FanOrder[ z ] = f;
	}
};

#endif // BITEOPT_INCLUDED
