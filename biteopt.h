//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The inclusion file for the CBEOOptimizer class.
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

#ifndef BITEOPT_INCLUDED
#define BITEOPT_INCLUDED

#include "biternd.h"

/**
 * "Bitmask evolution" optimization class. Implements a very simple stochastic
 * evolutionary optimization method (strategy) which involves inversion of a
 * random segment of parameter value's lowest bits at each step. Additionally
 * includes a crossing-over operation which in some cases improves convergence
 * considerably. For more robustness it is possible to assign several internal
 * values to each optimization parameter. This strategy is associated with
 * a very small code size and a minimal memory requirement.
 *
 * This strategy was tested on several classic 2-parameter optimization
 * problems and it performed fairly well. Global problems (with multiple local
 * minima) may not be handled well by this strategy, but in practice this
 * strategy strives to provide "minimum among minima" nevertheless.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal parameter values assigned to
 * each optimization parameter. Set to 2 or 3 to better solve more complex
 * functions. Not all functions will benefit from an increased value. Note
 * that the overhead is increased proportionally to this value.
 * @tparam HistSize Best parameter values history size. Affects convergence
 * time. Setting too low or too high values increases convergence time.
 */

template< int ParamCount0, int ValuesPerParam = 1, int HistSize = 64 >
class CBEOOptimizer
{
public:
	/**
	 * Constructor.
	 *
	 * @param aCrossProb Crossing-over probability [0; 1].
	 */

	CBEOOptimizer( const double aCrossProb = 0.4 )
		: CrossProb( aCrossProb )
		, MantMult( 1 << MantSize )
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

		HistPos = 0;
		int i;

		if( InitParams != NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const int k = i / ValuesPerParam;
				const double v = ( InitParams[ k ] - MinValues[ k ]) /
					( MaxValues[ k ] - MinValues[ k ]);

				Params[ i ] = (int) ( v * MantMult );

				if( Params[ i ] < 0 )
				{
					Params[ i ] = 0;
				}
				else
				if( Params[ i ] >= ( 1 << MantSize ))
				{
					Params[ i ] = ( 1 << MantSize ) - 1;
				}
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = (int) ( rnd.getRndValue() * MantMult );
			}
		}

		int j;

		for( j = 0; j < HistSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				HistParams[ j ][ i ] = (int) ( rnd.getRndValue() * MantMult );
			}
		}

		for( i = 0; i < ParamCount0; i++ )
		{
			DiffValues[ i ] = ( MaxValues[ i ] - MinValues[ i ]) /
				ValuesPerParam / MantMult;

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

		if( rnd.getRndValue() < CrossProb )
		{
			// Crossing-over with the historic best solutions.

			const int CrossHistPos = (int) ( rnd.getRndValue() * HistSize );
			const int* UseParams = HistParams[ CrossHistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// Crossing-over operation: copies lower bits of the historic
				// best solution.

				const int icmask =
					( 2 << (int) ( rnd.getRndValue() * MantSize )) - 1;

				Params[ i ] &= ~icmask;
				Params[ i ] |= UseParams[ i ] & icmask;

				// Bitmask inversion operation, works as a "driver" of
				// optimization process.

				const int imask =
					( 2 << (int) ( rnd.getRndValue() * MantSize )) - 1;

				Params[ i ] ^= imask;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// Bitmask inversion operation, works as a "driver" of
				// optimization process.

				const int imask =
					( 2 << (int) ( rnd.getRndValue() * MantSize )) - 1;

				Params[ i ] ^= imask;
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
				Params[ i ] = SaveParams[ i ];
			}
		}
		else
		{
			BestCost = NewCost;

			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			HistPos = ( HistPos == 0 ? HistSize : HistPos ) - 1;
			int* const hp = HistParams[ HistPos ];

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
	static const int ParamCount = ParamCount0 * ValuesPerParam; ///< Overall
		///< parameter count.
		///<
	static const int MantSize = 30; ///< Mantissa size of values.
		///<
	double CrossProb; ///< Crossing-over probability [0; 1].
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	int Params[ ParamCount ]; ///< Current working parameter states.
		///<
	int HistParams[ HistSize ][ ParamCount ]; ///< Best historic parameter
		///< states.
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
		///< minimal parameter values, multiplied by value's scaling factor.
		///<

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param i Parameter index.
	 */

	double getParamValue( const int i ) const
	{
		const int* p = Params + i * ValuesPerParam;
		double v = p[ 0 ];
		int j;

		for( j = 1; j < ValuesPerParam; j++ )
		{
			v += p[ j ];
		}

		return( MinValues[ i ] + DiffValues[ i ] * v );
	}
};

#endif // BITEOPT_INCLUDED
