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

#include <stdint.h>

/**
 * Class that implements random number generation. The default implementation
 * includes a relatively fast pseudo-random number generator (PRNG) using a
 * classic formula "seed = ( a * seed + c ) % m" (LCG). This implementation
 * uses bits 32-62 (30 bits) of the state variable which ensures at least 2^32
 * period in the lowest significant bit of the resulting pseudo-random
 * sequence. See https://en.wikipedia.org/wiki/Linear_congruential_generator
 * for more details.
 */

class CBEORnd
{
public:
	/**
	 * Function initializes *this PRNG object.
	 *
	 * @param NewSeed New random seed value. Lower 10 bits are used to select
	 * pseudo-random sequence.
	 */

	void init( const int NewSeed )
	{
		seed = NewSeed;

		// Skip first values to make PRNG "settle down".

		seed = 500009 * seed + 300119;
		seed = 500009 * seed + 300119;
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	virtual double getRndValue()
	{
		seed = 500009 * seed + 300119;

		return( (double) (int) (( seed >> 32 ) & 0x3FFFFFFF ) / 0x40000000 );
	}

private:
	uint64_t seed; ///< The current random seed value.
		///<
};

/**
 * Bitmask evolution optimization class. Implements a very simple evolution
 * strategy which involves inversion of a segment of parameter value's bits
 * at each step. Additionally includes crossing over operation which in some
 * cases improves convergence considerably. In some cases crossing over
 * reduces convergence, but only slightly. For more robustness it is possible
 * to assign several internal values to each optimization parameter.
 *
 * This strategy was tested on several classic 2-parameter optimization
 * problems and it performed fairly well. Global (multiple local minima)
 * problems may not be handled well by this strategy, but in practice this
 * strategy strives to provide "minimum among minimums" nevertheless.
 *
 * @tparam ParamCount0 The number of parameters being optimized.
 * @tparam ValuesPerParam The number of internal optimization values to use
 * for each parameter (1 or more). Increasing this value increases
 * optimization robustness, but may considerably increase the time of
 * convergence.
 * @tparam HistSize Best parameter values history size. Affects convergence
 * time setting to too low or too high values increases convergence time.
 */

template< int ParamCount0, int ValuesPerParam = 1, int HistSize = 32 >
class CBEOOptimizer
{
public:
	/**
	 * Constructor.
	 *
	 * @param aCrossProb Crossing-over probability [0; 1]. Should be usually
	 * set to a low value like 0.05 or 0.1.
	 */

	CBEOOptimizer( const double aCrossProb = 0.1 )
		: CrossProb( aCrossProb )
		, MantMult( 1 << MantSize )
		, BestCost( 1e10 )
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
				PrevBestParams[ j ][ i ] = Params[ i ];
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
		double NewParams[ ParamCount0 ];
		double NewCost;
		const bool DoCrossover = ( rnd.getRndValue() < CrossProb );
		int ipm[ ParamCount ]; // Bit inversion masks.
		int SaveParams[ ParamCount ];
		int i;

		if( DoCrossover )
		{
			const int CrossHistPos = (int) ( rnd.getRndValue() * HistSize );
			const int* UseParams = PrevBestParams[ CrossHistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];
				const int rmask = (int) ( rnd.getRndValue() * MantMult );
				Params[ i ] &= ( UseParams[ i ] | rmask );
			}

			for( i = 0; i < ParamCount0; i++ )
			{
				NewParams[ i ] = getParamValue( i );
			}

			NewCost = optcost( NewParams );

			for( i = 0; i < ParamCount; i++ )
			{
				const int t = Params[ i ];
				Params[ i ] = SaveParams[ i ];
				SaveParams[ i ] = t;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				ipm[ i ] = ( 1 << ( MantSize -
					(int) ( rnd.getRndValue() * MantSize ))) - 1;

				Params[ i ] ^= ipm[ i ];
			}

			for( i = 0; i < ParamCount0; i++ )
			{
				NewParams[ i ] = getParamValue( i );
			}

			NewCost = optcost( NewParams );

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] ^= ipm[ i ];
			}
		}

		if( NewCost < BestCost )
		{
			BestCost = NewCost;

			for( i = 0; i < ParamCount0; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}

			int* const hp = PrevBestParams[ HistPos ];
			HistPos = ( HistPos + 1 ) % HistSize;

			for( i = 0; i < ParamCount; i++ )
			{
				hp[ i ] = Params[ i ];
			}

			// Revert to the best parameter setting.

			if( DoCrossover )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const int t = Params[ i ];
					Params[ i ] = SaveParams[ i ];
					SaveParams[ i ] = t;
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] ^= ipm[ i ];
				}
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
	int PrevBestParams[ HistSize ][ ParamCount ]; ///< Best previous parameter
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
