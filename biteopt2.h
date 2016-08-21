//$ nocpp

/**
 * @file biteopt2.h
 *
 * @brief The inclusion file for the CBEOOptimizer2 class.
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

#ifndef BITEOPT2_INCLUDED
#define BITEOPT2_INCLUDED

#include <math.h>
#include "biternd.h"

/**
 * "Bitmask evolution" version 2 optimization class. Implements a very simple
 * stochastic evolutionary optimization method (strategy) which involves
 * inversion of a random segment of parameter value's lowest bits at each
 * step. Additionally includes the "step in the right direction" operation.
 *
 * This version provides a faster convergence time.
 *
 * This strategy was tested on several classic 2-parameter optimization
 * problems and it performed fairly well. Global problems (with multiple local
 * minima) may not be handled well by this strategy.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal parameter values assigned to
 * each optimization parameter. Set to 2 or 3 to better solve more complex
 * functions. Not all functions will benefit from an increased value. Note
 * that the overhead is increased proportionally to this value.
 */

template< int ParamCount0, int ValuesPerParam = 1 >
class CBEOOptimizer2
{
public:
	/**
	 * Constructor.
	 */

	CBEOOptimizer2()
		: MantMult( 1 << MantSize )
	{
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBEORnd& rnd, const double* const InitParams = NULL )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );
		int i;

		for( i = 0; i < ParamCount0; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}

		if( InitParams != NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const int k = i / ValuesPerParam;
				const double v = ( InitParams[ k ] - MinValues[ k ]) /
					DiffValues[ k ];

				if( v < 0.0 )
				{
					Params[ i ] = 0.0;
				}
				else
				if( v > 0.9999999999 )
				{
					Params[ i ] = 0.9999999999;
				}
				else
				{
					Params[ i ] = v * v;
				}

				PrevParams[ i ] = Params[ i ];
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double v = rnd.getRndValue();
				Params[ i ] = v * v;
				PrevParams[ i ] = Params[ i ];
			}
		}

		HistPos = 0;
		int j;

		for( j = 0; j < HistSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double v = rnd.getRndValue();
				HistParams[ j ][ i ] = v * v;
			}
		}

		for( i = 0; i < ParamCount0; i++ )
		{
			BestParams[ i ] = getParamValue( i );
		}

		BestCost = optcost( BestParams );
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 */

	void optimize( CBEORnd& rnd )
	{
		double SaveParams[ ParamCount ];
		int i;

		if( rnd.getRndValue() < 0.27 )
		{
			const int Pos = (int) ( rnd.getRndValue() * HistSize );
			const double* UseParams = HistParams[ Pos ];
			const double m = rnd.getRndValue();

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
			const double m = rnd.getRndValue();

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// Bitmask inversion operation, works as a "driver" of
				// optimization process.

				const int imask = ( 2 <<
					(int) ( rnd.getRndValue() * MantSize )) - 1;

				double np = ( (int) ( Params[ i ] * MantMult ) ^ imask ) /
					MantMult;

				// The "step in the right direction" operation.

				Params[ i ] = wrapParam( np - ( PrevParams[ i ] - np ) * m );
			}
		}

		double NewParams[ ParamCount0 ];

		for( i = 0; i < ParamCount0; i++ )
		{
			NewParams[ i ] = getParamValue( i );
		}

		const double NewCost = optcost( NewParams );

		if( NewCost >= BestCost )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				PrevParams[ i ] = Params[ i ];
				Params[ i ] = SaveParams[ i ];
			}
		}
		else
		{
			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;

			HistPos = ( HistPos + 1 ) % HistSize;
			double* const hp = HistParams[ HistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				hp[ i ] = SaveParams[ i ];
			}
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
	static const int HistSize = 14; ///< The size of the history.
		///<
	static const int MantSize = 29; ///< Mantissa size of bitmask inversion
		///< operation. Must be lower than the random number generator's
		///< precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double Params[ ParamCount ]; ///< Current working parameter vector.
		///<
	double PrevParams[ ParamCount ]; ///< Previously evaluated parameter
		///< vector.
		///<
	double HistParams[ HistSize ][ ParamCount ]; ///< Best historic parameter
		///< vectors.
		///<
	int HistPos; ///< Best parameter value history position.
		///<
	double BestParams[ ParamCount0 ]; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	double MinValues[ ParamCount0 ]; ///< Minimal parameter values.
		///<
	double MaxValues[ ParamCount0 ]; ///< Maximal parameter values.
		///<
	double DiffValues[ ParamCount0 ]; ///< Difference between maximal and
		///< minimal parameter values.
		///<

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
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param i Parameter index.
	 */

	double getParamValue( const int i ) const
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
};

#endif // BITEOPT2_INCLUDED
