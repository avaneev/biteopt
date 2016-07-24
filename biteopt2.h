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
 * evolutionary optimization method (strategy) which involves inversion of a
 * random segment of parameter value's lowest bits at each step. Additionally
 * includes a mysterious "previous attempt intermix" operation. Does not
 * require ancestors nor history of previous best solutions.
 *
 * This version provides a quite fast convergence time, a very small code
 * size and minimal memory requirement. The only drawback is that this
 * strategy requires 3 instead of 2 random number generator calls per
 * parameter on each step.
 *
 * This strategy was tested on several classic 2-parameter optimization
 * problems and it performed fairly well. Global problems (with multiple local
 * minima) may not be handled well by this strategy.
 *
 * @tparam ParamCount The number of parameters being optimized.
 */

template< int ParamCount, int HistSize = 16 >
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

		for( i = 0; i < ParamCount; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}

		HistPos = 0;
		HistCount = 0;

		if( InitParams != NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double v = ( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ];

				Params[ i ] = (int) ( v * v * MantMult );
				PrevParams[ i ] = Params[ i ];
				HistParams[ 0 ][ i ] = Params[ i ];
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = (int) ( rnd.getRndValue() * MantMult );
				PrevParams[ i ] = Params[ i ];
				HistParams[ 0 ][ i ] = Params[ i ];
			}
		}

		for( i = 0; i < ParamCount; i++ )
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
		int SaveParams[ ParamCount ];
		int i;

		if( rnd.getRndValue() < 0.1 )
		{
			// Crossing-over with the historic best solutions.

			const int CrossHistPos = (int) ( rnd.getRndValue() * HistCount );
			const int* UseParams =
				HistParams[( HistPos + CrossHistPos ) % HistSize ];

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// Replace lower bits with the historic best solution's ones.

				const int icmask = ( 2 <<
					(int) ( rnd.getRndValue() * MantSize )) - 1;

				Params[ i ] &= ~icmask;
				Params[ i ] |= UseParams[ i ] & icmask;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// The mysterious "previous attempt intermix" operation.

				const double r = rnd.getRndValue();
				Params[ i ] = (int) ( Params[ i ] * r +
					( Params[ i ] * 2.0 - PrevParams[ i ]) * ( 1.0 - r ));

				// Bitmask inversion operation with value clamping.

				const int imask = ( 2 <<
					(int) ( rnd.getRndValue() * MantSize )) - 1;

				if( Params[ i ] < 0 )
				{
					Params[ i ] = 0 ^ imask;
				}
				else
				if( Params[ i ] > MantSize1 )
				{
					Params[ i ] = MantSize1 ^ imask;
				}
				else
				{
					Params[ i ] ^= imask;
				}

				// Reduce swing of randomization by 20%.

				Params[ i ] = (int) ( SaveParams[ i ] * 0.2 +
					Params[ i ] * 0.8 );
			}
		}

		double NewParams[ ParamCount ];

		for( i = 0; i < ParamCount; i++ )
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
			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			BestCost = NewCost;

			HistPos = ( HistPos == 0 ? HistSize : HistPos ) - 1;
			int* const hp = HistParams[ HistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				hp[ i ] = SaveParams[ i ];
			}

			if( HistCount < HistSize )
			{
				HistCount++;
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
	static const int MantSize = 30; ///< Mantissa size of values.
		///<
	static const int MantSize1 = ( 1 << MantSize ) - 1; ///<
		///< Equals to ( 1 << MantSize ) - 1.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	int Params[ ParamCount ]; ///< Current working parameter states.
		///<
	int PrevParams[ ParamCount ]; ///< Previously evaluated parameters.
		///<
	int HistParams[ HistSize ][ ParamCount ]; ///< Best historic parameter
		///< values.
		///<
	int HistPos; ///< Best parameter value history position.
		///<
	int HistCount; ///< The total number of history additions performed.
		///< Always <= HistSize.
		///<
	double BestParams[ ParamCount ]; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	double MinValues[ ParamCount ]; ///< Minimal parameter values.
		///<
	double MaxValues[ ParamCount ]; ///< Maximal parameter values.
		///<
	double DiffValues[ ParamCount ]; ///< Difference between maximal and
		///< minimal parameter values.
		///<

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param i Parameter index.
	 */

	double getParamValue( const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] *
			sqrt( Params[ i ] / MantMult ));
	}
};

#endif // BITEOPT2_INCLUDED
