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
#include <stdlib.h>
#include <string.h>
#include "biternd.h"

/**
 * "Bitmask evolution" version "fan" optimization class. This strategy is
 * based on the CBEOOptimizer2 strategy, but uses several current parameter
 * vectors ("fan elements"). Any parameter vector can be replaced with a new
 * solution if parameter vector's cost (with some margin) is higher than that
 * of the new solution's. Having several "fan elements" allows parameter
 * vectors to be spaced apart from each other thus making them cover a larger
 * parameter search space collectively. The "fan elements" are used unevenly:
 * lower cost ones are evolved more frequently than the others.
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
 * square root distributed random number values. The parameter space is itself
 * square rooted on each function evaluation.
 *
 * 1. A set of "fan elements" is maintained. A "fan element" is an independent
 * parameter vector which is randomly evolved towards a better solution. The
 * "fan element" with the lowest cost is evolved more frequently than
 * "fan elements" with the higher costs.
 *
 * 2. A set of 14 best historic solutions is maintained. The history is shared
 * among all "fan elements". History is at first initialized with random
 * values that work as an additional source of randomization on initial
 * optimization steps.
 *
 * 3. A running-average centroid vector of best solutions is maintained.
 *
 * 4. The previous attempted solution parameter vector for each "fan element"
 * is maintained.
 *
 * 5. With 48% probability a crossing-over operation is performed which
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
 * 8. After each function evaluation, an attempt to replace one of the "fan
 * elements" is performed using cost constraints. This method is based on an
 * assumption that the later solutions tend to be better than the earlier
 * solutions.
 *
 * 9. History is updated with a previous solution whenever a better solution
 * is found or when a "fan element" is replaced.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal parameter values assigned to
 * each optimization parameter. Set to 2 or 3 to better solve more complex
 * functions. Not all functions will benefit from an increased value. Note
 * that the overhead is increased proportionally to this value.
 */

template< int ParamCount0, int ValuesPerParam = 1 >
class CBEOOptimizerFan
{
public:
	double CrossProb; ///< Crossing-over probability.
		///<
	double CentProb; ///< Centroid move probability.
		///<
	double CentTime; ///< Centroid averaging time (samples).
		///<
	double AvgCostMult; ///< Average "fan element" cost threshold multiplier.
		///<
	double CrossMults[ 3 ]; ///< Crossing-over range multipliers for each
		///< "fan element".
		///<
	double CentMult; ///< Centroid move range multiplier.
		///<

	/**
	 * Constructor.
	 */

	CBEOOptimizerFan()
		: MantMult( 1 << MantSize )
	{
		// Machine-optimized values.

		CrossProb = 0.480000;
		CentProb = 0.333333;
		CentTime = 9.570437;
		AvgCostMult = 2.721252;
		CrossMults[ 0 ] = 0.969624;
		CrossMults[ 1 ] = 0.977170;
		CrossMults[ 2 ] = 0.869048;
		CentMult = 0.942157;
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBEORnd& rnd, const double* const InitParams = NULL )
	{
		AvgCoeff = calcAvgCoeff( CentTime );

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

					CurParams[ j ][ i ] = v * v;
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
					const double v = rnd.getRndValue();
					CurParams[ j ][ i ] = v * v;
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
				const double v = rnd.getRndValue();
				HistParams[ j ][ i ] = v * v;
			}

			HistCosts[ j ] = BestCost; // Not entirely correct, but works.
		}

		updateAvgCost();
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 */

	void optimize( CBEORnd& rnd )
	{
		int s = (int) ( sqr( rnd.getRndValue() ) * FanSize );

		// "Fan element" sorting, for 3 items.

		if( CurCosts[ 1 ] < CurCosts[ 0 ])
		{
			if( CurCosts[ 2 ] < CurCosts[ 1 ])
			{
				static const int sorder[ 3 ] = { 2, 1, 0 };
				s = sorder[ s ];
			}
			else
			if( CurCosts[ 2 ] < CurCosts[ 0 ])
			{
				static const int sorder[ 3 ] = { 1, 2, 0 };
				s = sorder[ s ];
			}
			else
			{
				static const int sorder[ 3 ] = { 1, 0, 2 };
				s = sorder[ s ];
			}
		}
		else
		{
			if( CurCosts[ 2 ] < CurCosts[ 1 ])
			{
				if( CurCosts[ 2 ] < CurCosts[ 0 ])
				{
					static const int sorder[ 3 ] = { 2, 0, 1 };
					s = sorder[ s ];
				}
				else
				{
					static const int sorder[ 3 ] = { 0, 2, 1 };
					s = sorder[ s ];
				}
			}
		}

		const double cm = CrossMults[ s ];
		double* const Params = CurParams[ s ];
		double SaveParams[ ParamCount ];
		int i;

		if( rnd.getRndValue() < CrossProb )
		{
			// Crossing-over with one of the historic solutions.

			const int CrossHistPos = (int) ( rnd.getRndValue() * HistSize );
			const double* const UseParams = HistParams[ CrossHistPos ];
			const double m = sqrt( rnd.getRndValue() ) * cm;

			if( CurCosts[ s ] < HistCosts[ CrossHistPos ])
			{
				for( i = 0; i < ParamCount; i++ )
				{
					SaveParams[ i ] = Params[ i ];

					// The "step in the right direction" operation.

					Params[ i ] = wrapParam( Params[ i ] -
						( UseParams[ i ] - Params[ i ]) * m );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					SaveParams[ i ] = Params[ i ];

					// The "step in the right direction" operation.

					Params[ i ] = wrapParam( Params[ i ] -
						( Params[ i ] - UseParams[ i ]) * m );
				}
			}
		}
		else
		{
			const double m = sqrt( rnd.getRndValue() );

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

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
			const double m = rnd.getRndValue() * cm * CentMult;

			// The "step in the right direction" operation.

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

		double CopyParams[ ParamCount ];

		for( i = 0; i < ParamCount; i++ )
		{
			CopyParams[ i ] = Params[ i ];
			Params[ i ] = SaveParams[ i ];
		}

		// Possibly replace another least-performing "fan element".

		int f = -1;
		double MaxCost;

		for( i = 0; i < FanSize; i++ )
		{
			double d = ( CurCosts[ i ] - AvgCost ) * AvgCostMult;

			if( d < 0.0 )
			{
				d = 0.0;
			}

			if( NewCost < CurCosts[ i ] + d )
			{
				if( f == -1 || CurCosts[ i ] > MaxCost )
				{
					MaxCost = CurCosts[ i ];
					f = i;
				}
			}
		}

		if( f == -1 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				PrevParams[ s ][ i ] = CopyParams[ i ];
			}
		}
		else
		{
			double* const hp = advanceHist();

			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] +=
					( CurParams[ f ][ i ] - CentParams[ i ]) * AvgCoeff;

				hp[ i ] = CurParams[ f ][ i ];
				PrevParams[ f ][ i ] = CurParams[ f ][ i ];
				CurParams[ f ][ i ] = CopyParams[ i ];
			}

			HistCosts[ HistPos ] = CurCosts[ f ];
			CurCosts[ f ] = NewCost;
			updateAvgCost();
		}
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
	static const int FanSize = 3; ///< The number of "fan elements" to use.
		///<
	static const int HistSize = 14; ///< The size of the history.
		///<
	static const int MantSize = 29; ///< Mantissa size of bitmask inversion
		///< operation. Must be lower than the random number generator's
		///< precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
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
	double AvgCoeff; ///< Averaging coefficient for update of CentParams
		///< values.
		///<
	double AvgCost; ///< Average cost of all working parameter vectors.
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

		return( MinValues[ i ] +
			DiffValues[ i ] * sqrt( v / ValuesPerParam ));
	}

	/**
	 * @return Function advances the history position and returns pointer to
	 * the history vector.
	 */

	double* advanceHist()
	{
		HistPos = ( HistPos + 1 ) % HistSize;
		return( HistParams[ HistPos ]);
	}

	/**
	 * Function updates distances of all "fan elements".
	 */

	void updateAvgCost()
	{
		AvgCost = CurCosts[ 0 ];
		int i;

		for( i = 1; i < FanSize; i++ )
		{
			AvgCost += CurCosts[ i ];
		}

		AvgCost /= FanSize;
	}
};

#endif // BITEFAN_INCLUDED
