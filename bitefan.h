//$ nocpp

/**
 * @file bitefan.h
 *
 * @brief The inclusion file for the CBEOOptimizerFan class.
 *
 * @section license License
 * 
 * Copyright (c) 2016 Aleksey Vaneev
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
 * based on the CBEOOptimizer2 strategy, but uses several current parameter
 * vectors ("fan elements"). Highest cost "fan element" can be replaced with a
 * new solution if "fan element's" cost (plus some margin) is higher than that
 * of the new solution's. Having several "fan elements" allows parameter
 * vectors to be spaced apart from each other thus making them cover a larger
 * parameter search space collectively. The "fan elements" are used unevenly:
 * the lower cost ones are evolved more frequently than the higher cost ones.
 *
 * The benefit of this strategy is increased robustness: it can successfully
 * optimize a wider range of functions. Another benefit is a considerably
 * decreased convergence time in deeper optimizations. This strategy does not
 * solve all global optimization problems successfully, but strives to provide
 * the "minimum among minima" solution.
 *
 * This strategy is associated with a high overhead per function evaluation.
 * Due to this fact, for simple functions and not deep optimization it may be
 * more beneficial to use the CBEOOptimizer2 class.
 *
 * The strategy consists of the following elements. Most operations utilize a
 * log-like distributed random number values. The parameter space is itself
 * log-like transformed on each function evaluation.
 *
 * 1. A set of "fan elements" is maintained. A "fan element" is an independent
 * parameter vector which is randomly evolved towards a better solution. The
 * "fan element" with the lowest cost is evolved more frequently than
 * "fan elements" with the higher costs.
 *
 * 2. The previous attempted solution parameter vector for each "fan element"
 * is maintained.
 *
 * 3. A set of 14 better historic solutions is maintained. The history is
 * shared among all "fan elements". History is at first initialized with
 * random values that work as an additional source of randomization on initial
 * optimization steps.
 *
 * 4. A running-average centroid vector of better solutions is maintained.
 *
 * 5. With 48% probability a "history move" operation is performed which
 * involves a random historic solution. This operation consists of the
 * "step in the right direction" operation.
 *
 * 6. With the remaining probability the "bitmask evolution" (inversion of a
 * random range of the lowest bits) operation is performed, which is the main
 * driver of the evolutionary process, followed by the "step in the right
 * direction" operation using the previous solution.
 *
 * 7. Additionally, with 33% probability the "step in the right direction"
 * operation is performed using the centroid vector.
 *
 * 8. After each function evaluation, an attempt to replace the highest cost
 * "fan element" is performed using cost constraint. This method is based on
 * an assumption that the later solutions tend to be statistically better than
 * the earlier solutions. History is updated with a previous (replaced)
 * solution whenever a "fan element" is replaced.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal parameter values assigned to
 * each optimization parameter. Set to 2, 3 or 4 to better solve more complex
 * functions. Not all functions will benefit from an increased value. Note
 * that the overhead is increased proportionally to this value.
 */

template< int ParamCount0, int ValuesPerParam = 1 >
class CBEOOptimizerFan
{
public:
	double HistProb; ///< History move probability.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double CentTime; ///< Centroid averaging time (samples).
		///<
	double CostMult; ///< "Fan element" cost threshold multiplier.
		///<
	double HistMult; ///< History move range multiplier.
		///<
	double CentMult; ///< Centroid move range multiplier.
		///<
	double PrevMult; ///< Previous move range multiplier.
		///<
	double BestMult; ///< Best move range multiplier.
		///<
	double SpaceMult; ///< Parameter space adjustment multiplier.
		///<
	double HistRMult; ///< History move adjustment multiplier.
		///<
	double PrevRMult; ///< Previous move adjustment multiplier.
		///<

	/**
	 * Constructor.
	 */

	CBEOOptimizerFan()
		: MantMult( 1 << MantSize )
	{
		// Machine-optimized values.

		SpaceMult = 0.268346;
		HistRMult = 0.433420;
		PrevRMult = 0.573783;
		HistProb = 0.480000;
		CentProb = 0.333333;
		CentTime = 8.144686;
		CostMult = 1.370876;
		HistMult = 0.840073;
		CentMult = 0.845718;
		PrevMult = 1.001365;
		BestMult = 0.744041;
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBEORnd& rnd, const double* const InitParams = NULL )
	{
		CentCoeff = calcAvgCoeff( CentTime );

		HistM1 = HistRMult * sqr( sqr( HistMult ));
		HistM2 = ( 1.0 - HistRMult ) * sqr( sqr( HistMult ));

		PrevM1 = PrevRMult * sqr( sqr( PrevMult ));
		PrevM2 = ( 1.0 - PrevRMult ) * sqr( sqr( PrevMult ));
		SpaceMult2 = 1.0 - SpaceMult;

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
						clampParam(( InitParams[ k ] - MinValues[ k ]) /
						DiffValues[ k ]) : rnd.getRndValue() );

					CurParams[ j ][ i ] = getParamInv( v );
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
					CurParams[ j ][ i ] = getParamInv( rnd.getRndValue() );
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

			CurCosts[ j ] = optcost( Params );
			insertFanOrder( CurCosts[ j ], j, j );

			if( j == 0 || CurCosts[ j ] < BestCost )
			{
				BestCost = CurCosts[ j ];

				for( i = 0; i < ParamCount0; i++ )
				{
					BestParams[ i ] = Params[ i ];
				}
			}
		}

		// Initialize history with random values. This works as an additional
		// source of initial randomization.

		HistPos = 0;

		for( j = 0; j < HistSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				HistParams[ j ][ i ] = getParamInv( rnd.getRndValue() );
			}

			HistCosts[ j ] = BestCost; // Not entirely correct, but works.
		}
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 * @return "True" if optimizer's state was improved on this step. Many
	 * successive "false" results means optimizer has reached a plateau.
	 */

	bool optimize( CBEORnd& rnd )
	{
		const int s = FanOrder[ (int) ( sqr( rnd.getRndValue() ) * FanSize )];
		double Params[ ParamCount ];
		int i;

		if( true )
		{
			const double* const OrigParams = CurParams[ s ];
			const double* const UseParams = CurParams[ FanOrder[ 0 ]];
			const double m = rnd.getRndValue() * BestMult;

			// The "step in the right direction" operation towards the
			// best parameter vector.

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = OrigParams[ i ] -
					( OrigParams[ i ] - UseParams[ i ]) * m;
			}
		}

		if( rnd.getRndValue() < HistProb )
		{
			// Move towards one of the historic solutions.

			const int Pos = (int) ( rnd.getRndValue() * HistSize );
			const double* const UseParams = HistParams[ Pos ];
			const double r = rnd.getRndValue();
			const double m = sqrt( sqrt( r * ( HistM1 + r * r * HistM2 )));

			if( CurCosts[ s ] < HistCosts[ Pos ])
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( Params[ i ] -
						( UseParams[ i ] - Params[ i ]) * m );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( Params[ i ] -
						( Params[ i ] - UseParams[ i ]) * m );
				}
			}
		}
		else
		{
			const double r = rnd.getRndValue();
			const double m = sqrt( sqrt( r * ( PrevM1 + r * r * PrevM2 )));

			for( i = 0; i < ParamCount; i++ )
			{
				// Bitmask inversion operation, works as a "driver" of
				// optimization process.

				const int imask = ( 2 <<
					(int) ( sqrt( rnd.getRndValue() ) * MantSize )) - 1;

				double np = ( (int) ( Params[ i ] * MantMult ) ^ imask ) /
					MantMult;

				// The "step in the right direction" operation.

				Params[ i ] = wrapParam( np -
					( PrevParams[ s ][ i ] - np ) * m );
			}
		}

		if( rnd.getRndValue() < CentProb )
		{
			// Move towards centroid vector.

			const double m = rnd.getRndValue() * CentMult;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = wrapParam( Params[ i ] -
					( Params[ i ] - CentParams[ i ]) * m );
			}
		}

		// Evaluate function with new parameters.

		double NewParams[ ParamCount0 ];

		for( i = 0; i < ParamCount0; i++ )
		{
			NewParams[ i ] = getParamValue( Params, i );
		}

		const double NewCost = optcost( NewParams );

		if( NewCost < BestCost )
		{
			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;
		}

		// Try to replace the highest cost "fan element".

		const int sH = FanOrder[ FanSize1 ];
		const double cT = CurCosts[ sH ] +
			( CurCosts[ sH ] - CurCosts[ FanOrder[ 0 ]]) * CostMult;

		if( NewCost > cT )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				PrevParams[ s ][ i ] = Params[ i ];
			}

			return( false );
		}
		else
		{
			HistPos = ( HistPos + 1 ) % HistSize;
			HistCosts[ HistPos ] = CurCosts[ sH ];
			double* const hp = HistParams[ HistPos ];
			double* const rp = CurParams[ sH ];

			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] += ( rp[ i ] - CentParams[ i ]) * CentCoeff;
				hp[ i ] = rp[ i ];
				PrevParams[ sH ][ i ] = rp[ i ];
				rp[ i ] = Params[ i ];
			}

			CurCosts[ sH ] = NewCost;
			insertFanOrder( NewCost, sH, FanSize1 );

			return( true );
		}
	}

	/**
	 * Function performs iterative optimization function until the required
	 * minimum cost is reached or the maximum number of iterations was
	 * exceeded.
	 *
	 * @param rnd Random number generator.
	 * @param MinCost Target minimum cost.
	 * @param PlateauIters The number of successive iterations which have not
	 * produced a useful change, to treat as reaching a plateau. When plateau
	 * is reached, *this object will be reinitialized.
	 * @param MaxIters The maximal number of iterations to perform.
	 * @return The number of iterations performed, excluding reinitialization
	 * calls.
	 */

	int optimizePlateau( CBEORnd& rnd, const double MinCost,
		const int PlateauIters, const int MaxIters )
	{
		double TmpBestParams[ ParamCount0 ];
		double TmpBestCost;
		int Iters = 0;
		int i;

		init( rnd );
		bool IsFirstInit = true;
		int PlateauCount = 0;

		while( true )
		{
			const bool WasImproved = optimize( rnd );
			Iters++;

			if( BestCost <= MinCost )
			{
				return( Iters );
			}

			if( Iters >= MaxIters )
			{
				break;
			}

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

					init( rnd );
					IsFirstInit = false;
					PlateauCount = 0;
				}
			}
		}

		if( !IsFirstInit )
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
	 * Virtual function that should calculate parameter vector's optimization
	 * cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) const = 0;

protected:
	static const int ParamCount = ParamCount0 * ValuesPerParam; ///< The total
		///< number of internal parameter values in use.
		///<
	static const int FanSize = 4; ///< The number of "fan elements" to use.
		///<
	static const int FanSize1 = FanSize - 1; ///< = FanSize - 1.
		///<
	static const int HistSize = 14; ///< The size of the history.
		///<
	static const int MantSize = 29; ///< Mantissa size of bitmask inversion
		///< operation. Must be lower than the random number generator's
		///< precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double HistM1; ///< History move multiplier 1.
		///<
	double HistM2; ///< History move multiplier 1.
		///<
	double PrevM1; ///< Previous move multiplier 1.
		///<
	double PrevM2; ///< Previous move multiplier 2.
		///<
	double SpaceMult2; ///< Parameter space adjustment multiplier 2.
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
	double CentParams[ ParamCount ]; ///< Centroid of the best parameter
		///< vectors.
		///<
	double CentCoeff; ///< Averaging coefficient for update of CentParams
		///< values.
		///<
	double PrevParams[ FanSize ][ ParamCount ]; ///< Previously evaluated
		///< parameters.
		///<
	double HistParams[ HistSize ][ ParamCount ]; ///< Best historic parameter
		///< values.
		///<
	double HistCosts[ HistSize ]; ///< The costs of the best historic
		///< parameter values.
		///<
	int HistPos; ///< Best parameter value history position.
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
	 * Function calculates the averaging coefficient for 1st order low-pass
	 * filtering.
	 *
	 * @param Count The approximate number of values to average.
	 */

	static double calcAvgCoeff( const double Count )
	{
		const double theta = 2.79507498389883904 / Count;
		const double costheta2 = 2.0 - cos( theta );
		return( 1.0 - ( costheta2 - sqrt( costheta2 * costheta2 - 1.0 )));
	}

	/**
	 * @param x Value to square.
	 * @return Square of the argument.
	 */

	static double sqr( const double x )
	{
		return( x * x );
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0) range, by wrapping it over the boundaries. This operation
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

			if( v > 0.9999999999 )
			{
				v = 0.9999999999 - v;
				continue;
			}

			return( v );
		}
	}

	/**
	 * Function clamps the specified parameter value so that it stays in the
	 * [0.0; 1.0) range.
	 *
	 * @param v Parameter value to clamp.
	 * @return Clamped parameter value.
	 */

	static double clampParam( double v )
	{
		if( v < 0.0 )
		{
			return( 0.0 );
		}

		if( v > 0.9999999999 )
		{
			return( 0.9999999999 );
		}

		return( v );
	}

	/**
	 * @param v Value to inverse.
	 * @return Inverse of the parameter value space function.
	 */

	double getParamInv( const double v ) const
	{
		return(( sqrt( 4.0 * sqr( sqr( v )) * SpaceMult2 + sqr( SpaceMult )) -
			SpaceMult ) / ( 2.0 * SpaceMult2 ));
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

		v /= ValuesPerParam;
		v = sqrt( sqrt( v * ( SpaceMult + v * SpaceMult2 )));

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
		int z;

		for( z = 0; z < ItemCount; z++ )
		{
			if( Cost <= CurCosts[ FanOrder[ z ]])
			{
				break;
			}
		}

		int i;

		for( i = ItemCount; i > z; i-- )
		{
			FanOrder[ i ] = FanOrder[ i - 1 ];
		}

		FanOrder[ z ] = f;
	}
};

#endif // BITEFAN_INCLUDED
