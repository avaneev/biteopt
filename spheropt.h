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
 * @version 2021.6
 */

#ifndef SPHEROPT_INCLUDED
#define SPHEROPT_INCLUDED

#include "biteaux.h"

/**
 * "Converging hyper-spheroid" optimizer class. Simple, converges quite fast.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CSpherOpt : public CBiteOptBase
{
public:
	double Jitter; ///< Solution sampling random jitter, improves convergence
		///< at low dimensions. Usually, a fixed value.
		///<
	double EvalFac; ///< Evaluations factor. Usually, a fixed value.
		///<

	CSpherOpt()
		: WBuf( NULL )
	{
		Jitter = 2.5;
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

		WBuf = new double[ PopSize * 6 ];
		WPopCent[ 0 ] = &WBuf[ 0 ];
		WPopCent[ 1 ] = &WBuf[ PopSize ];
		WPopCent[ 2 ] = &WBuf[ PopSize * 2 ];
		WPopRad[ 0 ] = &WBuf[ PopSize * 3 ];
		WPopRad[ 1 ] = &WBuf[ PopSize * 4 ];
		WPopRad[ 2 ] = &WBuf[ PopSize * 5 ];

		double s1[ 3 ] = { 0.0 };
		double s2[ 3 ] = { 0.0 };
		int i;
		int k;

		for( i = 0; i < PopSize; i++ )
		{
			const double l = 1.0 - (double) i / ( PopSize * EvalFac );

			static const double WCent[ 3 ] = { 5.0, 7.5, 10.0 };
			static const double WRad[ 3 ] = { 16.0, 18.0, 20.0 };

			for( k = 0; k < 3; k++ )
			{
				const double v1 = pow( l, WCent[ k ]);
				WPopCent[ k ][ i ] = v1;
				s1[ k ] += v1;

				const double v2 = pow( l, WRad[ k ]);
				WPopRad[ k ][ i ] = v2;
				s2[ k ] += v2;
			}
		}

		for( i = 0; i < PopSize; i++ )
		{
			for( k = 0; k < 3; k++ )
			{
				WPopCent[ k ][ i ] /= s1[ k ];
				WPopRad[ k ][ i ] /= s2[ k ];
			}
		}

		JitMult = 2.0 * Jitter / aParamCount;
		JitOffs = 1.0 - JitMult * 0.5;
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid.
	 * @param InitRadius Initial radius, multiplier relative to the default
	 * sigma value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();
		updateDiffValues();

		Radius = 0.5 * InitRadius;
		curpi = 0;
		cure = 0;

		// Provide initial centroid and sigma (CurParams is used temporarily,
		// otherwise initially undefined).

		int i;

		if( InitParams == NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] = 0.5;
			}

			DoCentEval = false;
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] = wrapParam( rnd,
					( InitParams[ i ] - MinValues[ i ]) / DiffValues[ i ]);
			}

			DoCentEval = true;
		}

		CentPowHist.reset();
		RadPowHist.reset();
		cpm = 0;
		rpm = 0;
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

		if( DoCentEval )
		{
			DoCentEval = false;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = CentParams[ i ];
				NewParams[ i ] = getRealValue( CentParams, i );
			}
		}
		else
		{
			double s2 = 1e-300;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rnd.getRndValue() - 0.5;
				s2 += Params[ i ] * Params[ i ];
			}

			const double d = Radius / sqrt( s2 );

			if( ParamCount > 4 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( rnd,
						CentParams[ i ] + Params[ i ] * d );

					NewParams[ i ] = getRealValue( Params, i );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const double m = JitOffs + rnd.getRndValue() * JitMult;

					Params[ i ] = wrapParam( rnd,
						CentParams[ i ] + Params[ i ] * d * m );

					NewParams[ i ] = getRealValue( Params, i );
				}
			}
		}

		const double NewCost = optcost( NewParams );

		if( OutCost != NULL )
		{
			*OutCost = NewCost;
		}

		if( OutParams != NULL )
		{
			memcpy( OutParams, Params, ParamCount * sizeof( Params[ 0 ]));
		}

		updateBestCost( NewCost, NewParams );

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
				memcpy( CurParams[ sH ], Params,
					ParamCount * sizeof( Params[ 0 ]));

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

				CentPowHist.incr( cpm );
				RadPowHist.incr( rpm );
			}
			else
			{
				StallCount += cure;

				CentPowHist.decr( cpm );
				RadPowHist.decr( rpm );
			}

			cpm = CentPowHist.select( rnd );
			rpm = RadPowHist.select( rnd );

			AvgCost = 0.0;
			curpi = 0;
			cure = 0;

			update();
		}

		return( StallCount );
	}

protected:
	double* WBuf; ///< Buffer for weighting coefficients.
		///<
	double* WPopCent[ 3 ]; ///< Weighting coefficients for centroid.
		///<
	double* WPopRad[ 3 ]; ///< Weighting coefficients for radius.
		///<
	double JitMult; ///< Jitter multiplier.
		///<
	double JitOffs; ///< Jitter multiplier offset.
		///<
	double Radius; ///< Current radius.
		///<
	int curpi; ///< Current parameter index.
		///<
	int cure; ///< Current evaluation index.
		///<
	bool DoCentEval; ///< "True" if an initial objective function evaluation
		///< at centroid point is required.
	CBiteOptHist< 3, 3, 2 > CentPowHist; ///< Centroid power factor histogram.
		///<
	CBiteOptHist< 3, 3, 2 > RadPowHist; ///< Radius power factor histogram.
		///<
	int cpm; ///< Centroid power factor selector.
		///<
	int rpm; ///< Radius power factor selector.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] WBuf;
	}

	/**
	 * Function updates centroid and radius.
	 */

	void update()
	{
		const double* ip = CurParams[ PopOrder[ 0 ]];
		double* const cp = CentParams;
		const double* const wc = WPopCent[ cpm ];
		double w = wc[ 0 ];
		int i;
		int j;

		for( i = 0; i < ParamCount; i++ )
		{
			cp[ i ] = ip[ i ] * w;
		}

		for( j = 1; j < PopSize; j++ )
		{
			ip = CurParams[ PopOrder[ j ]];
			w = wc[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += ip[ i ] * w;
			}
		}

		const double* const rc = WPopRad[ rpm ];
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

			Radius += s * rc[ j ];
		}

		Radius = sqrt( Radius );
	}
};

#endif // SPHEROPT_INCLUDED
