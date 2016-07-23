//$ nocpp

/**
 * @file biteopt.h
 *
 * @brief The "main" inclusion file.
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
 * "Bitmask evolution" optimization class. Implements a very simple
 * evolutionary optimization method (strategy) which involves inversion of a
 * random segment of parameter value's lowest bits at each step. Additionally
 * includes a crossing-over operation which in some cases improves convergence
 * considerably. In some cases crossing-over reduces convergence, but only
 * slightly. For more robustness it is possible to assign several internal
 * values to each optimization parameter.
 *
 * This strategy was tested on several classic 2-parameter optimization
 * problems and it performed fairly well. Global problems (with multiple local
 * minima) may not be handled well by this strategy, but in practice this
 * strategy strives to provide "minimum among minima" nevertheless.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal optimization values to use
 * for each parameter (1 or more). Increasing this value increases
 * optimization robustness, but may increase the time of convergence.
 * @tparam HistSize Best parameter values history size. Affects convergence
 * time. Setting too low or too high values increases convergence time.
 */

template< int ParamCount0, int ValuesPerParam = 1, int HistSize = 64 >
class CBEOOptimizer
{
public:
	bool WasBestCost; ///< "True" if the best cost was found on the last
		///< optimize() function call.
		///<

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
	 */

	void init( CBEORnd& rnd )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		HistPos = 0;
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = (int) ( rnd.getRndValue() * MantMult );
			int j;

			for( j = 0; j < HistSize; j++ )
			{
				HistParams[ j ][ i ] = Params[ i ];
			}
		}

		for( i = 0; i < ParamCount0; i++ )
		{
			DiffValues[ i ] = ( MaxValues[ i ] - MinValues[ i ]) /
				ValuesPerParam / MantMult;

			BestParams[ i ] = getParamValue( i );
		}

		BestCost = optcost( BestParams );
		WasBestCost = true;
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 */

	void optimize( CBEORnd& rnd )
	{
		double NewParams[ ParamCount0 ];
		double NewCost;
		const bool DoCrossover = ( rnd.getRndValue() < CrossProb );
		int SaveParams[ ParamCount ];
		int i;

		if( DoCrossover )
		{
			const int CrossHistPos = (int) ( rnd.getRndValue() * HistSize );
			const int* UseParams =
				HistParams[( HistPos + CrossHistPos ) % HistSize ];

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];
				int icmask = ( 1 << ( MantSize -
					(int) ( rnd.getRndValue() * MantSize ))) - 1;

				Params[ i ] &= ~icmask;
				Params[ i ] |= UseParams[ i ] & icmask;

				const int imask = ( 1 << ( MantSize -
					(int) ( rnd.getRndValue() * MantSize ))) - 1;

				Params[ i ] ^= imask;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];
				const int imask = ( 1 << ( MantSize -
					(int) ( rnd.getRndValue() * MantSize ))) - 1;

				Params[ i ] ^= imask;
			}
		}

		for( i = 0; i < ParamCount0; i++ )
		{
			NewParams[ i ] = getParamValue( i );
		}

		NewCost = optcost( NewParams );

		if( NewCost >= BestCost )
		{
			WasBestCost = false;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = SaveParams[ i ];
			}
		}
		else
		{
			WasBestCost = true;
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
	 * Function adds history parameters derived from another optimizer.
	 *
	 * @param s Source optimizer.
	 * @param SetCurrent "True" if current parameters should be replaced
	 * instead of adding them into history.
	 */

	void addHistParams( const CBEOOptimizer& s, const bool SetCurrent )
	{
		int i;

		if( SetCurrent )
		{
			if( s.getBestCost() < BestCost )
			{
				BestCost = s.getBestCost();

				for( i = 0; i < ParamCount0; i++ )
				{
					BestParams[ i ] = s.BestParams[ i ];
				}
			}

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = s.Params[ i ];
			}
		}
		else
		{
			HistPos = ( HistPos == 0 ? HistSize : HistPos ) - 1;
			int* const hp = HistParams[ HistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				hp[ i ] = s.Params[ i ];
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
