//$ nocpp

/**
 * @file bitehive.h
 *
 * @brief The "hive" strategy inclusion file.
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

#ifndef BITEHIVE_INCLUDED
#define BITEHIVE_INCLUDED

#include <stdlib.h>

/**
 * "Bitmask Evolution Hive" optimization strategy class.
 *
 * The CBEOHive class implements "hive" optimization strategy which utilizes
 * several CBEOOptimizer objects in parallel exchanging solutions between
 * them. This can offer a benefit for some "hard" functions by reducing the
 * number of required function evaluations by a factor of 2. For other
 * functions there may be a negative benefit. CBEOHive is best used with
 * low (0.1-0.2) crossing-over probabilities and ValuesPerParam=3 or 4.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal optimization values to use
 * for each parameter (1 or more). Increasing this value increases
 * optimization robustness, but may increase the time of convergence.
 * @tparam HiveSize The number of optimizers to use in the "hive". Must be
 * above 2.
 * @tparam HistSize Best parameter values history size. Affects convergence
 * time. Setting too low or too high values increases convergence time.
 */

template< int ParamCount0, int ValuesPerParam = 4, int HiveSize = 5,
	int HistSize = 64 >
class CBEOHive
{
public:
	/**
	 * Constructor.
	 *
	 * @param aCrossProb Crossing-over probability [0; 1].
	 */

	CBEOHive( const double aCrossProb = 0.15 )
	{
		int i;

		for( i = 0; i < HiveSize; i++ )
		{
			Opts[ i ] = new COpt( this, aCrossProb );
		}
	}

	~CBEOHive()
	{
		int i;

		for( i = 0; i < HiveSize; i++ )
		{
			delete Opts[ i ];
		}
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 */

	void init( CBEORnd& rnd )
	{
		BestCost = 1e100;
		int i;

		for( i = 0; i < HiveSize; i++ )
		{
			Opts[ i ] -> init( rnd );

			if( Opts[ i ] -> getBestCost() < BestCost )
			{
				BestCost = Opts[ i ] -> getBestCost();
			}
		}
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 */

	void optimize( CBEORnd& rnd )
	{
		int i;

		for( i = 0; i < HiveSize; i++ )
		{
			if( !Opts[ i ] -> WasBestCost )
			{
				continue;
			}

			int j;

			if( Opts[ i ] -> getBestCost() <= BestCost )
			{
				BestCost = Opts[ i ] -> getBestCost();

				for( j = 0; j < ParamCount0; j++ )
				{
					BestParams[ j ] = Opts[ i ] -> getBestParams()[ j ];
				}

				const int t = HiveSize / 2 - 1;

				if( i >= t )
				{
					for( j = 0; j < t; j++ )
					{
						Opts[ j ] -> addHistParams( *Opts[ i ], true );
					}
				}
				else
				{
					for( j = t; j < HiveSize; j++ )
					{
						Opts[ j ] -> addHistParams( *Opts[ i ], false );
					}
				}
			}
		}

		for( i = 0; i < HiveSize; i++ )
		{
			Opts[ i ] -> optimize( rnd );
		}
	}

	/**
	 * @return Best parameter vector.
	 */

	double getBestCost() const
	{
		return( BestCost );
	}

	/**
	 * @return Cost of the best parameter vector.
	 */

	const double* getBestParams() const
	{
		return( BestParams );
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
	/**
	 * Internal "overloaded" optimizer class.
	 */

	class COpt : public CBEOOptimizer< ParamCount0, ValuesPerParam, HistSize >
	{
	public:
		COpt( CBEOHive* const aOwner, const double aCrossProb )
			: CBEOOptimizer< ParamCount0, ValuesPerParam, HistSize >(
				aCrossProb )
			, Owner( aOwner )
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

		virtual double optcost( const double* const p ) const
		{
			return( Owner -> optcost( p ));
		}

	protected:
		CBEOHive* Owner;
	};

	COpt* Opts[ HiveSize ]; ///< Optimizers.
		///<
	double BestCost; ///< The global best cost found so far.
		///<
	double BestParams[ ParamCount0 ]; ///< The global best parameters found
		///< so far.
		///<
};

#endif // BITEHIVE_INCLUDED
