//$ nocpp

/**
 * @file spheropt.h
 *
 * @brief The inclusion file for the CSpherOpt class.
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
 * @version 2021.2
 */

#ifndef SPHEROPT_INCLUDED
#define SPHEROPT_INCLUDED

#include "biteaux.h"

/**
 * "Converging hyper-spheroid" optimizer class. Converges quite fast, but is
 * not effective for dimensions below 4.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CSpherOpt : public CBiteOptBase
{
public:
	double CentPow; ///< Centroid power factor.
		///<
	double RadPow; ///< Radius power factor.
		///<
	double EvalFac; ///< Evaluations factor.
		///<

	CSpherOpt()
		: WPopCent( NULL )
		, WPopRad( NULL )
	{
		CentPow = 6.0;
		RadPow = 18.0;
		EvalFac = 2.0;
	}

	/**
	 * Function updates dimensionality of *this object.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0 or negative, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 : 14 + aParamCount );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();
		initBaseBuffers( aParamCount, aPopSize );

		WPopCent = new double[ PopSize ];
		WPopRad = new double[ PopSize ];

		double s = 0.0;
		double s2 = 0.0;
		int i;

		for( i = 0; i < PopSize; i++ )
		{
			const double l = 1.0 - (double) i / ( PopSize * EvalFac );
			const double v = pow( l, CentPow );
			WPopCent[ i ] = v;
			s += v;
			const double v2 = pow( l, RadPow );
			WPopRad[ i ] = v2;
			s2 += v2;
		}

		for( i = 0; i < PopSize; i++ )
		{
			WPopCent[ i ] /= s;
			WPopRad[ i ] /= s2;
		}
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();

		Radius = 0.5;
		curpi = 0;
		cure = 0;

		// Provide initial centroid and sigma (CurParams is used temporarily,
		// otherwise initially undefined).

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
			CentParams[ i ] = 0.5;
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
		double s2 = 1e-300;
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = rnd.getRndValue() - 0.5;
			s2 += Params[ i ] * Params[ i ];
		}

		const double d = Radius / sqrt( s2 );

		for( i = 0; i < ParamCount; i++ )
		{
			Params[ i ] = wrapParam( rnd, CentParams[ i ] + Params[ i ] * d );
			NewParams[ i ] = getRealValue( Params[ i ], i );
		}

		const double NewCost = optcost( NewParams );

		if( OutCost != NULL )
		{
			*OutCost = NewCost;
		}

		if( OutParams != NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				OutParams[ i ] = Params[ i ];
			}
		}

		if( NewCost < BestCost )
		{
			BestCost = NewCost;

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = NewParams[ i ];
			}
		}

		if( curpi < PopSize )
		{
			insertPopOrder( NewCost, curpi, curpi );
			curpi++;
		}
		else
		{
			const int sH = PopOrder[ PopSize1 ];

			if( NewCost < CurCosts[ sH ])
			{
				for( i = 0; i < ParamCount; i++ )
				{
					CurParams[ sH ][ i ] = Params[ i ];
				}

				insertPopOrder( NewCost, sH, PopSize1 );
			}
		}

		AvgCost += NewCost;
		cure++;

		if( cure >= PopSize * EvalFac )
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
			update();
		}

		return( StallCount );
	}

protected:
	double* WPopCent; ///< Weighting coefficients for centroid.
		///<
	double* WPopRad; ///< Weighting coefficients for radius.
		///<
	double Radius; ///< Current radius.
		///<
	int curpi; ///< Current parameter index.
		///<
	int cure; ///< Current evaluation index.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] WPopCent;
		delete[] WPopRad;
	}

	/**
	 * Function updates centroid and radius.
	 */

	void update()
	{
		const double* ip = CurParams[ PopOrder[ 0 ]];
		double* const cp = CentParams;
		double w = WPopCent[ 0 ];
		int i;
		int j;

		for( i = 0; i < ParamCount; i++ )
		{
			cp[ i ] = ip[ i ] * w;
		}

		for( j = 1; j < PopSize; j++ )
		{
			ip = CurParams[ PopOrder[ j ]];
			w = WPopCent[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += ip[ i ] * w;
			}
		}

		Radius = 0.0;

		for( j = 0; j < PopSize; j++ )
		{
			ip = CurParams[ PopOrder[ j ]];
			double s = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const double d = ip[ i ] - cp[ i ];
				s += d * d;
			}

			Radius += s * WPopRad[ j ];
		}

		Radius = sqrt( Radius );
	}
};

#endif // SPHEROPT_INCLUDED
