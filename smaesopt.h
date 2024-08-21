//$ nocpp

/**
 * @file smaesopt.h
 *
 * @version 2024.6
 *
 * @brief The inclusion file for the CSMAESOpt class.
 *
 * @section license License
 *
 * Copyright (c) 2016-2024 Aleksey Vaneev
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

#ifndef SMAESOPT_INCLUDED
#define SMAESOPT_INCLUDED

#include "biteort.h"

/**
 * Sigma Adaptation Evolution Strategy class. Fundamentally similar to CMA-ES,
 * but mainly focuses on sigma adaptation.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CSMAESOpt : public CBiteOptBase< double >
{
public:
	/**
	 * Function updates dimensionality of *this object.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0 or negative, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 : 13 + aParamCount );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		Ort.updateDims( aParamCount, aPopSize );
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values, only used as centroid,
	 * not evaluated.
	 * @param InitRadius Initial radius, multiplier relative to the default
	 * sigma value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		initCommonVars( rnd );

		cure = 0;
		curem = (int) ceil( CurPopSize * Ort.EvalFac );

		// Provide initial centroid and sigma (PopParams is used here
		// temporarily, otherwise initially undefined).

		const double sd = 0.25 * InitRadius;
		int i;

		if( InitParams == NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				PopParams[ 0 ][ i ] = MinValues[ i ] + DiffValues[ i ] * 0.5;
				PopParams[ 1 ][ i ] = fabs( DiffValues[ i ]) * sd;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				PopParams[ 0 ][ i ] = InitParams[ i ];
				PopParams[ 1 ][ i ] = fabs( DiffValues[ i ]) * sd;
			}
		}

		Ort.init( PopParams[ 0 ], PopParams[ 1 ]);
	}

	/**
	 * Function samples a random population vector based on the current
	 * distribution, with feasibility guarantee.
	 *
	 * @param rnd Random number generator.
	 * @param[out] op Resulting parameter vector.
	 */

	void sample( CBiteRnd& rnd, double* const op ) const
	{
		// Generate vector, check its feasibility, and resample it up to 10
		// times.

		int infcount = 0;
		int i;

		while( true )
		{
			Ort.sample( rnd, op );

			if( isFeasible( op ))
			{
				break;
			}

			infcount++;

			if( infcount == 10 )
			{
				// Force bound constraints.

				for( i = 0; i < ParamCount; i++ )
				{
					op[ i ] = wrapParamReal( rnd, op[ i ], i );
				}

				break;
			}
		}
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return The number of non-improving iterations so far.
	 */

	int optimize( CBiteRnd& rnd )
	{
		double* const Params = getCurParams();

		sample( rnd, Params );

		const double NewCost = fixCostNaN( optcost( Params ));
		NewCosts[ 0 ] = NewCost;
		LastValues = Params;

		updatePop( NewCost, Params );
		updateBestCost( NewCost, Params );

		AvgCost += NewCost;
		cure++;

		if( cure >= curem )
		{
			AvgCost /= cure;

			if( AvgCost < HiBound )
			{
				HiBound = AvgCost;
			}

			resetCurPopPos();
			AvgCost = 0.0;
			cure = 0;

			Ort.update( *this );
		}

		StallCount = ( NewCost < HiBound ? 0 : StallCount + 1 );

		return( StallCount );
	}

protected:
	CBiteOrt Ort; ///< Rotation vector and orthogonalization calculator.
	int cure; ///< Current evaluation index, greater or equal to
		///< "curem" if population distribution needs to be updated.
	int curem; ///< "cure" value threshold.

	/**
	 * Function returns "true" if the supplied vector is feasible.
	 *
	 * @param x Vector to check.
	 */

	bool isFeasible( const double* const x ) const
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			if( x[ i ] < MinValues[ i ] || x[ i ] > MaxValues[ i ])
			{
				return( false );
			}
		}

		return( true );
	}
};

#endif // SMAESOPT_INCLUDED
