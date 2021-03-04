//$ nocpp

/**
 * @file smaesopt.h
 *
 * @brief The inclusion file for the CSMAESOpt class.
 *
 * @section license License
 *
 * Copyright (c) 2016-2021 Aleksey Vaneev
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
 *
 * @version 2021.8
 */

#ifndef SMAESOPT_INCLUDED
#define SMAESOPT_INCLUDED

#include "biteoptort.h"

/**
 * Sigma Adaptation Evolution Strategy class. Fundamentally similar to CMA-ES,
 * but mainly focuses on sigma adaptation.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CSMAESOpt : public CBiteOptBase
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

		deleteBuffers();
		initBaseBuffers( aParamCount, aPopSize );

		EvalFac = 2.0;

		Ort.updateDims( aParamCount, aPopSize, EvalFac );
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
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();

		curpi = 0;
		cure = 0;

		// Provide initial centroid and sigma (CurParams is used here
		// temporarily, otherwise initially undefined).

		const double sd = 0.25 * InitRadius;
		int i;

		if( InitParams == NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CurParams[ 0 ][ i ] =
					( MinValues[ i ] + MaxValues[ i ]) * 0.5;

				const double d = MaxValues[ i ] - MinValues[ i ];
				CurParams[ 1 ][ i ] = fabs( d ) * sd;
				DiffValues[ i ] = 1.0 / d;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CurParams[ 0 ][ i ] = InitParams[ i ];
				const double d = MaxValues[ i ] - MinValues[ i ];
				CurParams[ 1 ][ i ] = fabs( d ) * sd;
				DiffValues[ i ] = 1.0 / d;
			}
		}

		UsePopSize = Ort.init( CurParams[ 0 ], CurParams[ 1 ]);
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
					if( op[ i ] < MinValues[ i ])
					{
						op[ i ] = MinValues[ i ];
					}
					else
					if( op[ i ] > MaxValues[ i ])
					{
						op[ i ] = MaxValues[ i ];
					}
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
	 * @param OutCost If not NULL, pointer to variable that receives cost
	 * of the newly-evaluated solution.
	 * @param OutParams If not NULL, pointer to array that receives
	 * newly-evaluated parameter vector, in normalized scale.
	 * @return The number of non-improving iterations so far.
	 */

	int optimize( CBiteRnd& rnd, double* const OutCost = NULL,
		double* const OutParams = NULL )
	{
		double* const Params = CurParams[ curpi ];
		int i;

		sample( rnd, Params );

		const double NewCost = optcost( Params );

		if( OutCost != NULL )
		{
			*OutCost = NewCost;
		}

		if( OutParams != NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				OutParams[ i ] = ( Params[ i ] - MinValues[ i ]) *
					DiffValues[ i ];
			}
		}

		updateBestCost( NewCost, Params );

		if( curpi < UsePopSize )
		{
			sortPop( NewCost, curpi );
			curpi++;
		}
		else
		{
			const int ps1 = UsePopSize - 1;

			if( NewCost <= CurCosts[ ps1 ])
			{
				memcpy( CurParams[ ps1 ], Params,
					ParamCount * sizeof( CurParams[ 0 ]));

				sortPop( NewCost, ps1 );
			}
		}

		AvgCost += NewCost;
		cure++;

		if( cure >= UsePopSize * EvalFac )
		{
			AvgCost /= cure;

			if( AvgCost < HiBound )
			{
				HiBound = AvgCost;
				StallCount = 0;
			}
			else
			{
				StallCount += cure;
			}

			AvgCost = 0.0;
			curpi = 0;
			cure = 0;
			UsePopSize = Ort.update( CurParams );
		}

		return( StallCount );
	}

protected:
	int UsePopSize; ///< Current population size.
		///<
	double EvalFac; ///< Function evalutions factor.
		///<
	CBiteOptOrt Ort; ///< Rotation vector and orthogonalization calculator.
		///<
	int curpi; ///< Current parameter index.
		///<
	int cure; ///< Current evaluation index, equals UsePopSize if population
		///< distribution needs to be updated.
		///<

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
