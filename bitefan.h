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

#ifndef BITEFAN_INCLUDED
#define BITEFAN_INCLUDED

#include <math.h>
#include "biternd.h"

/**
 * "Bitmask evolution" version "fan" optimization class. This strategy is
 * based on now outdated CBEOOptimizer and CBEOOptimizer2 stochastic
 * derivative-less strategies, and uses several current parameter vectors
 * ("fan elements"). Highest cost "fan element" can be replaced with a new
 * solution if "fan element's" cost (plus some margin) is higher than that of
 * the new solution's. Having several "fan elements" allows parameter vectors
 * to be spaced apart from each other thus making them cover a larger
 * parameter search space collectively. The strategy was named as "bitmask
 * evolution", because at its core an operation of inversion of a random
 * segment of parameter value's lowest bits is used. Beside that, several
 * "step in the right direction" operations are used that move the solution
 * vector into position with a probably lower objective function value.
 *
 * The benefit of this strategy is a relatively high robustness: it can
 * successfully optimize a wide range of test functions. Another benefit is a
 * low convergence time which depends on the complexity of the objective
 * function. This strategy does not solve all global optimization problems
 * successfully, but strives to provide the "minimum among minima" solution.
 * Like many stochastic optimization strategies, this strategy can't solve
 * problems with narrow or rogue optimums.
 *
 * The strategy consists of the following elements:
 *
 * 1. A set of "fan elements" is maintained. A "fan element" is an independent
 * parameter vector which is evolved towards a better solution. Also an
 * ordered list of "fan elements" is maintaned.
 *
 * 2. The previous attempted/rejected solution parameter vector for each
 * "fan element" is maintained.
 *
 * 3. A single previous (outdated) historic solution is maintained, which is
 * shared among all "fan elements".
 *
 * 4. A centroid vector of all "fan elements" is maintained.
 *
 * 5. On every step the "best move" operation is performed which involves the
 * current best solution.
 *
 * 6. On every step an "away from history move" operation is performed which
 * involves an outdated historic solution.
 *
 * 7. With 50% probability the "bitmask evolution" (inversion of a random
 * range of the lowest bits of a single random parameter) operation is
 * performed, which is the main driver of the evolutionary process, followed
 * by the "step in the right direction" operation using a previous (rejected)
 * solution.
 *
 * 8. Additionally, with 33% probability the "step in the right direction"
 * operation is performed using the centroid vector.
 *
 * 9. After each objective function evaluation, an attempt to replace the
 * highest cost "fan element" is performed using the cost constraint. This
 * method is based on an assumption that the later solutions tend to be
 * statistically better than the earlier solutions. History is updated with a
 * previous (replaced) solution whenever a "fan element" is replaced.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal parameter values assigned to
 * each optimization parameter. Set to 2, 3 or 4 to better solve more complex
 * objective functions. Not all functions will benefit from an increased
 * value. Note that the overhead is increased proportionally to this value.
 * @tparam FanSize The number of "fan elements" to use. Higher values increase
 * convergence time, but improve solution quality.
 */

template< int ParamCount0, int ValuesPerParam = 1,
	int FanSize = 6 + ParamCount0 * ParamCount0 / 2 >
class CBEOOptimizerFan
{
public:
	double CostMult; ///< "Fan element" cost threshold multiplier.
		///<
	double BestMult; ///< Best move range multiplier.
		///<
	double HistMult; ///< History move range multiplier.
		///<
	double PrevMult; ///< Previous move range multiplier.
		///<
	double CentMult; ///< Centroid move range multiplier.
		///<
	double CentOffs; ///< Centroid move range shift.
		///<

	/**
	 * Constructor.
	 */

	CBEOOptimizerFan()
		: MantMult( 1 << MantSize )
	{
		CostMult = 1.44176476;
		BestMult = 0.65528962;
		HistMult = 0.53627432;
		PrevMult = 0.43633471;
		CentMult = 1.26765785;
		CentOffs = 0.67884726;
	}

	/**
	 * Function initializes *this optimizer. Performs N=FanSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBEORnd& rnd, const double* const InitParams = NULL )
	{
		PrevCnt = 0;
		CentCnt = 0;
		ParamCnt = 0;

		getMinValues( MinValues );
		getMaxValues( MaxValues );
		int i;
		int j;

		for( i = 0; i < ParamCount0; i++ )
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
					const int k = i / ValuesPerParam;
					const double v = ( j == 0 ?
						wrapParam(( InitParams[ k ] - MinValues[ k ]) /
						DiffValues[ k ]) : rnd.getRndValue() );

					CurParams[ j ][ i ] = v;
					PrevParams[ j ][ i ] = CurParams[ j ][ i ];
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
					PrevParams[ j ][ i ] = CurParams[ j ][ i ];
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
			double Params[ ParamCount0 ];

			for( i = 0; i < ParamCount0; i++ )
			{
				Params[ i ] = getParamValue( CurParams[ j ], i );
			}

			insertFanOrder( optcost( Params ), j, j );

			if( j == 0 || CurCosts[ j ] < BestCost )
			{
				BestCost = CurCosts[ j ];

				for( i = 0; i < ParamCount0; i++ )
				{
					BestParams[ i ] = Params[ i ];
				}
			}
		}

		// Initialize history with random values.

		for( i = 0; i < ParamCount; i++ )
		{
			HistParams[ i ] = rnd.getRndValue();
		}
	}

	/**
	 * Function performs the parameter optimization step that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return "True" if optimizer's state was improved on this step. Many
	 * successive "false" results means optimizer has reached a plateau.
	 */

	bool optimize( CBEORnd& rnd )
	{
		const int s = FanOrder[ (int) ( rnd.getRndValue() * FanSize )];
		double Params[ ParamCount ];
		int i;

		// The "step in the right direction" operation towards the best
		// (minimal) parameter vector.

		const double* const OrigParams = CurParams[ s ];
		const double* const MinParams = CurParams[ FanOrder[ 0 ]];

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = OrigParams[ i ] -
				( OrigParams[ i ] - MinParams[ i ]) * BestMult;
		}

		// Move away from a previous historic solution.

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] -= ( HistParams[ i ] - Params[ i ]) * HistMult;
		}

		PrevCnt ^= 1;

		if( PrevCnt == 1 )
		{
			// Bitmask inversion operation, works as a "driver" of
			// optimization process, applied to 1 random parameter at a time.

			const int imask =
				( 2 << (int) ( rnd.getRndValue() * MantSize )) - 1;

			Params[ ParamCnt ] = ( (int) ( Params[ ParamCnt ] * MantMult ) ^
				imask ) / MantMult;

			ParamCnt = ( ParamCnt == 0 ? ParamCount : ParamCnt ) - 1;

			// The "step in the right direction" operation, away from the
			// previously rejected solution.

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] -=
					( PrevParams[ s ][ i ] - Params[ i ]) * PrevMult;
			}
		}

		CentCnt = ( CentCnt + 1 ) % 3;

		if( CentCnt == 1 )
		{
			// Move towards centroid vector randomly.

			for( i = 0; i < ParamCount; i++ )
			{
				const double m = CentOffs + rnd.getRndValue() * CentMult;
				Params[ i ] -= ( Params[ i ] - CentParams[ i ]) * m;
			}
		}

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = wrapParam( Params[ i ]);
		}

		// Evaluate objective function with new parameters.

		double NewParams[ ParamCount0 ];

		for( i = 0; i < ParamCount0; i++ )
		{
			NewParams[ i ] = getParamValue( Params, i );
		}

		const double NewCost = optcost( NewParams );

		// Try to replace the highest cost "fan element".

		const int sH = FanOrder[ FanSize1 ];
		const double cT = CurCosts[ sH ] +
			( CurCosts[ sH ] - CurCosts[ FanOrder[ 0 ]]) * CostMult;

		if( NewCost > cT )
		{
			if( PrevCnt == 1 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					PrevParams[ s ][ i ] = Params[ i ];
				}
			}

			return( false );
		}

		if( NewCost < BestCost )
		{
			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;
		}

		// Replace highest cost "fan element".

		double* const rp = CurParams[ sH ];
		double* const pp = PrevParams[ sH ];

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( Params[ i ] - rp[ i ]) / FanSize;
			HistParams[ i ] = rp[ i ];
			pp[ i ] = rp[ i ];
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

	int optimizePlateau( CBEORnd& rnd, const double MinCost,
		const int PlateauIters, const int MaxIters )
	{
		double TmpBestParams[ ParamCount0 ];
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

						for( i = 0; i < ParamCount0; i++ )
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

			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = TmpBestParams[ i ];
			}
		}

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
	static const int ParamCount = ParamCount0 * ValuesPerParam; ///< The total
		///< number of internal parameter values in use.
		///<
	static const int FanSize1 = FanSize - 1; ///< = FanSize - 1.
		///<
	static const int MantSize = 29; ///< Mantissa size of bitmask inversion
		///< operation. Must be lower than the random number generator's
		///< precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	int PrevCnt; ///< Previous move counter.
		///<
	int CentCnt; ///< Centroid move counter.
		///<
	int ParamCnt; ///< Parameter index counter.
		///<
	int FanOrder[ FanSize ]; ///< The current "fan element" ordering,
		///< ascending-sorted by cost.
		///<
	double CurParams[ FanSize ][ ParamCount ]; ///< Current working parameter
		///< vectors.
		///<
	double CurCosts[ FanSize ]; ///< Best costs of current working parameter
		///< vectors.
		///<
	double CentParams[ ParamCount ]; ///< Centroid of the current parameter
		///< vectors.
		///<
	double PrevParams[ FanSize ][ ParamCount ]; ///< Previously evaluated
		///< (and rejected) parameters.
		///<
	double HistParams[ ParamCount ]; ///< Last better parameter values.
		///<
	double MinValues[ ParamCount0 ]; ///< Minimal parameter values.
		///<
	double MaxValues[ ParamCount0 ]; ///< Maximal parameter values.
		///<
	double DiffValues[ ParamCount0 ]; ///< Difference between maximal and
		///< minimal parameter values.
		///<
	double BestParams[ ParamCount0 ]; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0] range, by wrapping it over the boundaries. This operation
	 * increases convergence in comparison to clamping.
	 *
	 * @param v Parameter value to wrap.
	 * @return Wrapped parameter value.
	 */

	static double wrapParam( double v )
	{
		while( true )
		{
			if( v < 0.0 )
			{
				v = -v;
			}

			if( v > 1.0 )
			{
				v = 2.0 - v;
				continue;
			}

			return( v );
		}
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param Params Parameter vector of interest.
	 * @param i Parameter index.
	 */

	double getParamValue( const double* const Params, const int i ) const
	{
		const double* const p = Params + i * ValuesPerParam;
		double v = p[ 0 ];
		int k;

		for( k = 1; k < ValuesPerParam; k++ )
		{
			v += p[ k ];
		}

		return( MinValues[ i ] + DiffValues[ i ] * v / ValuesPerParam );
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

#endif // BITEFAN_INCLUDED
