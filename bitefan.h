//$ nocpp

/**
 * @file bitefan.h
 *
 * @brief The inclusion file for the CBEOOptimizerFan class.
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

#ifndef BITEFAN_INCLUDED
#define BITEFAN_INCLUDED

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "biternd.h"

/**
 * "Bitmask evolution" version "fan" optimization class. This strategy is
 * based on the CBEOOptimizer2 strategy, but uses several current parameter
 * vectors ("fan elements"). Any parameter vector can be replaced with a new
 * solution if parameter vector's cost is higher than that of the new
 * solution's and the "distance" of the new solution is not considerably low.
 * The "distance" constraint allows parameter vectors to be spaced apart from
 * each other thus making them cover a larger parameter search space
 * collectively. The "fan elements" are used unevenly: some are used more
 * frequently than the others.
 *
 * The benefit of this strategy is increased robustness: it can successfully
 * optimize a wider range of functions. Another benefit is a considerably
 * decreased convergence time in deeper optimizations. This strategy does not
 * solve all global optimization problems successfully, but strives to provide
 * the "minimum among minima" solution.
 *
 * This strategy is associated with a high overhead per function evaluation.
 * Due to this fact, for simple functions and not deep optimization it may be
 * more beneficial to use the CBEOOptimizer2 class.
 *
 * The strategy consists of the following elements. Most operations utilize a
 * square root distributed random number values. The parameter space is itself
 * square rooted on each function evaluation.
 *
 * 1. A set of "fan elements" is maintained. A "fan element" is an independent
 * parameter vector which is randomly evolved towards a better solution. The
 * "fan element" with the highest cost is evolved more frequently than the
 * "fan elements" with lower costs.
 *
 * 2. A set of 14 best historic solutions is maintained. The history is shared
 * among all "fan elements". History is at first initialized with random
 * values that work as an additional source of randomization on initial
 * optimization steps.
 *
 * 3. The previous attempted solution parameter vector for each "fan element"
 * is maintained.
 *
 * 4. A running average (centroid) and average centroid distance of all
 * attempted solutions for each "fan element" is maintained.
 *
 * 5. With 57% probability a crossing-over operation is performed which
 * involves a random historic best solution. This operation consists of the
 * "step in the right direction" operation.
 *
 * 6. With the remaining probability the "step in the right direction"
 * operation is performed using the previous solution, followed by the
 * "bitmask evolution" operation which is the main driver of the evolutionary
 * process.
 *
 * 7. If a better solution was not found on the current step, an attempt to
 * replace one of the "fan elements" is performed using cost and parameter
 * distance constraints.
 *
 * @tparam ParamCount The number of parameters being optimized.
 */

template< int ParamCount >
class CBEOOptimizerFan
{
public:
	/**
	 * Constructor.
	 */

	CBEOOptimizerFan()
		: MantMult( 1 << MantSize )
		, MantDiv08( 0.8 / ( 1 << MantSize ))
		, AvgCoeff( calcAvgCoeff( 28 ))
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

		double Centr[ ParamCount ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			Centr[ i ] = 0.0;
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}

		HistPos = 0;
		int j;

		if( InitParams != NULL )
		{
			for( j = 0; j < FanSize; j++ )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const double v = ( j == 0 ?
						( InitParams[ i ] - MinValues[ i ]) /
						DiffValues[ i ] : rnd.getRndValue() );

					CurParams[ j ][ i ] = v * v;
					PrevParams[ j ][ i ] = CurParams[ j ][ i ];
					Centr[ i ] += CurParams[ j ][ i ];
				}
			}
		}
		else
		{
			for( j = 0; j < FanSize; j++ )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const double v = rnd.getRndValue();
					CurParams[ j ][ i ] = v * v;
					PrevParams[ j ][ i ] = CurParams[ j ][ i ];
					Centr[ i ] += CurParams[ j ][ i ];
				}
			}
		}

		// Calculate costs of "fan elements" and find the best cost.

		for( j = 0; j < FanSize; j++ )
		{
			double Params[ ParamCount ];

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = getParamValue( CurParams[ j ], i );
			}

			CurCosts[ j ] = optcost( Params );
			PrevCosts[ j ] = CurCosts[ j ];

			if( j == 0 || CurCosts[ j ] < BestCost )
			{
				BestCost = CurCosts[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					BestParams[ i ] = Params[ i ];
				}
			}
		}

		// Initialize history with random values. This works as an additional
		// source of initial randomization.

		for( j = 0; j < HistSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double v = rnd.getRndValue();
				HistParams[ j ][ i ] = v * v;
			}

			HistCosts[ j ] = BestCost; // Not entirely correct, but works.
		}

		// Set the same centroid in all "fan elements".

		for( i = 0; i < ParamCount; i++ )
		{
			Centr[ i ] /= FanSize;

			for( j = 0; j < FanSize; j++ )
			{
				AvgParams[ j ][ i ] = Centr[ i ];
			}
		}

		// Calculate initial average centroid distance.

		double Dist = 0.0;

		for( j = 0; j < FanSize; j++ )
		{
			double s = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const double d = Centr[ i ] - CurParams[ j ][ i ];
				s += d * d;
			}

			Dist += s;
		}

		Dist /= FanSize;

		for( j = 0; j < FanSize; j++ )
		{
			AvgCentrDists[ j ] = Dist;
		}

		updateDistances();
	}

	/**
	 * Function performs 1 parameter optimization step.
	 *
	 * @param rnd Random number generator.
	 */

	void optimize( CBEORnd& rnd )
	{
		const int s = (int) ( isqr( rnd.getRndValue() ) * FanSize );
		double* const Params = CurParams[ s ];
		double SaveParams[ ParamCount ];
		const double rp = rnd.getRndValue();
		int i;

		if( rp < 0.57 )
		{
			// Crossing-over with one of the historic best solutions.

			const int CrossHistPos = (int) ( rnd.getRndValue() * HistSize );
			const double* const UseParams = HistParams[ CrossHistPos ];

			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// The "step in the right direction" operation, with reduction
				// of swing by 42%, and with value clamping.

				double d;

				if( CurCosts[ s ] < HistCosts[ CrossHistPos ])
				{
					d = UseParams[ i ] - Params[ i ];
				}
				else
				{
					d = Params[ i ] - UseParams[ i ];
				}

				Params[ i ] = clampParam( Params[ i ] - d *
					sqrt( rnd.getRndValue() ) * 0.58 );
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				SaveParams[ i ] = Params[ i ];

				// The "step in the right direction" operation.

				const double r = sqrt( rnd.getRndValue() );
				double d;

				if( CurCosts[ s ] < PrevCosts[ s ])
				{
					d = PrevParams[ s ][ i ] - Params[ i ];
				}
				else
				{
					d = Params[ i ] - PrevParams[ s ][ i ];
				}

				double np = clampParam( Params[ i ] - d * r );

				// Bitmask inversion operation with value clamping, works as
				// a "driver" of optimization process.

				const int imask = ( 2 << (int) ( r * MantSize )) - 1;
				np = (int) ( np * MantMult ) ^ imask;

				// Reduce swing of randomization by 20%.

				Params[ i ] = SaveParams[ i ] * 0.2 + np * MantDiv08;
			}
		}

		// Keep average evaluated parameter values (centroid), and average
		// centroid distance (squared).

		double CentrDist = 0.0;

		for( i = 0; i < ParamCount; i++ )
		{
			AvgParams[ s ][ i ] +=
				( Params[ i ] - AvgParams[ s ][ i ]) * AvgCoeff;

			const double d = AvgParams[ s ][ i ] - Params[ i ];
			CentrDist += d * d;
		}

		AvgCentrDists[ s ] += ( CentrDist - AvgCentrDists[ s ]) * AvgCoeff;

		if( CentrDist > AvgCentrDists[ s ])
		{
			// If new parameter values are located far from the centroid,
			// bring them closer to the centroid.

			const double m =
				sqrt( sqrt( AvgCentrDists[ s ] / CentrDist )) * 0.73;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = AvgParams[ s ][ i ] +
					( Params[ i ] - AvgParams[ s ][ i ]) * m;
			}
		}

		// Evaluate function with new parameters.

		double NewParams[ ParamCount ];

		for( i = 0; i < ParamCount; i++ )
		{
			NewParams[ i ] = getParamValue( Params, i );
		}

		const double NewCost = optcost( NewParams );

		if( NewCost >= CurCosts[ s ])
		{
			double CopyParams[ ParamCount ];
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				PrevParams[ s ][ i ] = Params[ i ];
				CopyParams[ i ] = Params[ i ];
				Params[ i ] = SaveParams[ i ];
			}

			PrevCosts[ s ] = NewCost;

			// Possibly replace another least-performing "fan element".

			int f = -1;
			double MaxDist = 0.0;

			for( i = 0; i < FanSize; i++ )
			{
				double d = ( CurCosts[ i ] - AvgCost ) * 2.25;

				if( d < 0.0 )
				{
					d = 0.0;
				}

				if( NewCost < CurCosts[ i ] + d )
				{
					const double NewDist = calcDistance( CopyParams, i );

					if( NewDist > MaxDist && NewDist > AvgDist * 0.21 )
					{
						MaxDist = NewDist;
						f = i;
					}
				}
			}

			if( f != -1 )
			{
				double* const hp = advanceHist();

				for( i = 0; i < ParamCount; i++ )
				{
					hp[ i ] = CurParams[ f ][ i ];
					CurParams[ f ][ i ] = CopyParams[ i ];
				}

				HistCosts[ HistPos ] = CurCosts[ f ];
				CurCosts[ f ] = NewCost;
				updateDistances();
			}
		}
		else
		{
			if( NewCost < BestCost )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					BestParams[ i ] = NewParams[ i ];
				}

				BestCost = NewCost;
			}

			double* const hp = advanceHist();

			for( i = 0; i < ParamCount; i++ )
			{
				hp[ i ] = SaveParams[ i ];
			}

			HistCosts[ HistPos ] = CurCosts[ s ];
			CurCosts[ s ] = NewCost;
			updateDistances();
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
	static const int FanSize = 3; ///< The number of "fan elements" to use.
		///<
	static const int HistSize = 14; ///< The size of the history.
		///<
	static const int MantSize = 30; ///< Mantissa size of values. Must be
		///< synchronized with the random number generator's precision.
		///<
	double MantMult; ///< Mantissa multiplier (1 << MantSize).
		///<
	double MantDiv08; ///< Mantissa divisor (0.8 / MantMult).
		///<
	double CurParams[ FanSize ][ ParamCount ]; ///< Current working parameter
		///< vectors.
		///<
	double CurCosts[ FanSize ]; ///< Best costs of current working parameter
		///< vectors.
		///<
	double CurDists[ FanSize ]; ///< Average distances to other working
		///< parameters vectors.
		///<
	double AvgParams[ FanSize ][ ParamCount ]; ///< Averaged parameter vectors
		///< used. Resembles a centroid.
		///<
	double AvgCentrDists[ FanSize ]; ///< Average distances of parameter
		///< vectors to centroid, squared.
		///<
	double AvgCoeff; ///< Averaging coefficient for update of AvgParams and
		///< AvgCentrDists values.
		///<
	double AvgDist; ///< Average distance of all working parameter vectors.
		///<
	double AvgCost; ///< Average cost of all working parameter vectors.
		///<
	double PrevParams[ FanSize ][ ParamCount ]; ///< Previously evaluated
		///< parameters.
		///<
	double PrevCosts[ FanSize ]; ///< Costs of previously evaluated
		///< parameters.
		///<
	double HistParams[ HistSize ][ ParamCount ]; ///< Best historic parameter
		///< values.
		///<
	double HistCosts[ HistSize ]; ///< The costs of the best historic
		///< parameter values.
		///<
	int HistPos; ///< Best parameter value history position.
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
	 * Structure used for "fan element" sorting.
	 */

	struct CFanSortStruct
	{
		int i; ///< The index of the element.
			///<
		double Cost; ///< The cost of the element.
			///<
	};

	/**
	 * Function calculates the averaging coefficient for 1st order low-pass
	 * filtering.
	 *
	 * @param Count The approximate number of values to average.
	 */

	static double calcAvgCoeff( const int Count )
	{
		const double theta = 2.79507498389883904 / Count;
		const double costheta2 = 2.0 - cos( theta );
		return( 1.0 - ( costheta2 - sqrt( costheta2 * costheta2 - 1.0 )));
	}

	/**
	 * @param x Value to invert square, in the range 0 to 1.
	 * @return Inverted square of the argument.
	 */

	static double isqr( const double x )
	{
		const double x1 = 1.0 - x;
		return( 1.0 - x1 * x1 );
	}

	/**
	 * Function clamps the specified parameter value so that it stays in the
	 * 0.0 to 1.0 range, inclusive.
	 *
	 * @param v Parameter value to clamp.
	 * @return Clamped parameter value.
	 */

	static double clampParam( const double v )
	{
		if( v < 0.0 )
		{
			return( 0.0 );
		}

		if( v > 1.0 )
		{
			return( 1.0 );
		}

		return( v );
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param Params Parameter vector of interest.
	 * @param i Parameter index.
	 */

	double getParamValue( const double* const Params, const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] * sqrt( Params[ i ]));
	}

	/**
	 * @return Function advances the history position and returns pointer to
	 * the history vector.
	 */

	double* advanceHist()
	{
		HistPos = ( HistPos + 1 ) % HistSize;
		return( HistParams[ HistPos ]);
	}

	/**
	 * "Fan element" sorting function, used in the qsort() function call.
	 *
	 * @param p1 Element 1.
	 * @param p2 Element 2.
	 */

	static int FanElementSortFn( const void* p1, const void* p2 )
	{
		const double c1 = ( (CFanSortStruct*) p1 ) -> Cost;
		const double c2 = ( (CFanSortStruct*) p2 ) -> Cost;

		if( c1 < c2 )
		{
			return( 1 );
		}

		if( c1 > c2 )
		{
			return( -1 );
		}

		return( 0 );
	}

	/**
	 * Function sorts "fan elements" by cost.
	 */

	void sortFanElements()
	{
		CFanSortStruct fe[ FanSize ];
		int i;

		for( i = 0; i < FanSize; i++ )
		{
			fe[ i ].i = i;
			fe[ i ].Cost = CurCosts[ i ];
		}

		qsort( fe, FanSize, sizeof( fe[ 0 ]), FanElementSortFn );

		double CurParamsS[ FanSize ][ ParamCount ];
		double CurCostsS[ FanSize ];
		double CurDistsS[ FanSize ];
		double AvgParamsS[ FanSize ][ ParamCount ];
		double AvgCentrDistsS[ FanSize ];
		double PrevParamsS[ FanSize ][ ParamCount ];
		double PrevCostsS[ FanSize ];

		memcpy( CurParamsS, CurParams, sizeof( CurParamsS ));
		memcpy( CurCostsS, CurCosts, sizeof( CurCostsS ));
		memcpy( CurDistsS, CurDists, sizeof( CurDistsS ));
		memcpy( AvgParamsS, AvgParams, sizeof( AvgParamsS ));
		memcpy( AvgCentrDistsS, AvgCentrDists, sizeof( AvgCentrDistsS ));
		memcpy( PrevParamsS, PrevParams, sizeof( PrevParamsS ));
		memcpy( PrevCostsS, PrevCosts, sizeof( PrevCostsS ));

		for( i = 0; i < FanSize; i++ )
		{
			const int s = fe[ i ].i;

			memcpy( CurParams[ i ], CurParamsS[ s ], sizeof( CurParams[ i ]));
			CurCosts[ i ] = CurCostsS[ s ];
			CurDists[ i ] = CurDistsS[ s ];
			memcpy( AvgParams[ i ], AvgParamsS[ s ], sizeof( AvgParams[ i ]));
			AvgCentrDists[ i ] = AvgCentrDistsS[ s ];
			memcpy( PrevParams[ i ], PrevParamsS[ s ],
				sizeof( PrevParams[ i ]));

			PrevCosts[ i ] = PrevCostsS[ s ];
		}
	}

	/**
	 * Function calculates distance of the specified parameter values to all
	 * other parameter vectors ("fan elements").
	 *
	 * @param Params Parameters whose distance to calculate.
	 * @param Skip The index of "fan element" to skip from calculation.
	 * @return Distance to "fan elements".
	 */

	double calcDistance( const double* const Params, const int Skip ) const
	{
		double Dist = 0.0;
		int j;
		int i;

		for( j = 0; j < Skip; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Dist += fabs( CurParams[ j ][ i ] - Params[ i ]);
			}
		}

		for( j = Skip + 1; j < FanSize; j++ )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Dist += fabs( CurParams[ j ][ i ] - Params[ i ]);
			}
		}

		return( Dist / ( ParamCount * ( FanSize - 1 )));
	}

	/**
	 * Function updates distances of all "fan elements".
	 */

	void updateDistances()
	{
		AvgDist = 0.0;
		int i;

		for( i = 0; i < FanSize; i++ )
		{
			CurDists[ i ] = calcDistance( CurParams[ i ], i );
			AvgDist += CurDists[ i ];
		}

		AvgDist /= FanSize;
		AvgCost = CurCosts[ 0 ];
		double MinCost = CurCosts[ 0 ];

		for( i = 1; i < FanSize; i++ )
		{
			AvgCost += CurCosts[ i ];

			if( CurCosts[ i ] < MinCost )
			{
				MinCost = CurCosts[ i ];
			}
		}

		AvgCost /= FanSize;

		sortFanElements();
	}
};

#endif // BITEFAN_INCLUDED
