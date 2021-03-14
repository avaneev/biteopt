//$ nocpp

/**
 * @file nmsopt.h
 *
 * @brief The inclusion file for the CNMSeqOpt class.
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
 * @version 2021.17
 */

#ifndef NMSOPT_INCLUDED
#define NMSOPT_INCLUDED

#include "biteaux.h"

/**
 * Sequential Nelder-Mead simplex method.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CNMSeqOpt : public CBiteOptBase< double >
{
public:
	CNMSeqOpt()
		: x2( NULL )
	{
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
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 :
			( aParamCount + 1 ) * 4 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		N = aParamCount;
		M = aPopSize;
		M1m = 1.0 / ( M - 1 );

		initBuffers( N, M );
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid.
	 * @param InitRadius Initial radius, relative to the default value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars( rnd );

		// Initialize parameter vectors, costs and centroid.

		double* const xx = x[ 0 ];
		int i;

		if( InitParams != NULL )
		{
			memcpy( xx, InitParams, N * sizeof( x[ 0 ]));
		}
		else
		{
			for( i = 0; i < N; i++ )
			{
				xx[ i ] = MinValues[ i ] + DiffValues[ i ] * 0.5;
			}
		}

		xlo = 0;

		const double sd = 0.25 * InitRadius;
		int j;

		for( j = 1; j < M; j++ )
		{
			double* const xj = x[ j ];

			for( i = 0; i < N; i++ )
			{
				xj[ i ] = xx[ i ] + DiffValues[ i ] * getGaussian( rnd ) * sd;
			}
		}

		State = stReflection;
		DoInitEvals = true;
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @param OutCost If not NULL, pointer to variable that receives cost
	 * of the newly-evaluated solution.
	 * @param OutValues If not NULL, pointer to array that receives
	 * newly-evaluated parameter vector, in real scale.
	 * @return The number of non-improving iterations so far.
	 */

	int optimize( CBiteRnd& rnd, double* const OutCost = NULL,
		double* const OutValues = NULL )
	{
		int i;

		if( DoInitEvals )
		{
			y[ CurPopPos ] = eval( rnd, x[ CurPopPos ], OutCost, OutValues );

			if( y[ CurPopPos ] < y[ xlo ])
			{
				xlo = CurPopPos;
			}

			CurPopPos++;

			if( CurPopPos == M )
			{
				DoInitEvals = false;
				calccent();
			}

			return( 0 );
		}

		StallCount++;

		static const double alpha = 1.0; // Reflection coeff.
		static const double gamma = 2.0; // Expansion coeff.
		static const double rho = -0.5; // Contraction coeff.
		static const double sigma = 0.5; // Reduction coeff.

		double* const xH = x[ xhi ]; // Highest cost parameter vector.

		switch( State )
		{
			case stReflection:
			{
				for( i = 0; i < N; i++ )
				{
					x1[ i ] = x0[ i ] + alpha * ( x0[ i ] - xH[ i ]);
				}

				y1 = eval( rnd, x1, OutCost, OutValues );

				if( y1 >= y[ xlo ] && y1 < y[ xhi2 ])
				{
					copy( x1, y1 );
					StallCount = 0;
				}
				else
				{
					if( y1 < y[ xlo ])
					{
						State = stExpansion;
						StallCount = 0;
					}
					else
					{
						State = stContraction;
					}
				}

				break;
			}

			case stExpansion:
			{
				for( i = 0; i < N; i++ )
				{
					x2[ i ] = x0[ i ] + gamma * ( x0[ i ] - xH[ i ]);
				}

				const double y2 = eval( rnd, x2, OutCost, OutValues );
				xlo = xhi;

				if( y2 < y1 )
				{
					copy( x2, y2 );
				}
				else
				{
					copy( x1, y1 );
				}

				State = stReflection;
				break;
			}

			case stContraction:
			{
				for( i = 0; i < N; i++ )
				{
					x2[ i ] = x0[ i ] + rho * ( x0[ i ] - xH[ i ]);
				}

				const double y2 = eval( rnd, x2, OutCost, OutValues );

				if( y2 < y[ xhi ])
				{
					if( y2 < y[ xlo ])
					{
						xlo = xhi;
					}

					copy( x2, y2 );
					State = stReflection;
					StallCount = 0;
				}
				else
				{
					rx = x[ xlo ];
					rj = 0;
					State = stReduction;
				}

				break;
			}

			case stReduction:
			{
				double* xx = x[ rj ];

				if( xx == rx )
				{
					rj++;
					xx = x[ rj ];
				}

				for( i = 0; i < N; i++ )
				{
					xx[ i ] = rx[ i ] + sigma * ( xx[ i ] - rx[ i ]);
				}

				y[ rj ] = eval( rnd, xx, OutCost, OutValues );

				if( y[ rj ] < y[ xlo ])
				{
					xlo = rj;
					StallCount = 0;
				}

				rj++;

				if( rj == M || ( rj == M - 1 && x[ rj ] == rx ))
				{
					calccent();
					State = stReflection;
				}

				break;
			}
		}

		return( StallCount );
	}

private:
	int N; ///< The total number of internal parameter values in use.
		///<
	int M; ///< The number of points in a simplex.
		///<
	double M1m; ///< = 1 / ( M - 1 ).
		///<
	int xlo; ///< Current lowest cost parameter vector.
		///<
	int xhi; ///< Current highest cost parameter vector.
		///<
	int xhi2; ///< Current second highest cost parameter vector.
		///<
	double** x; ///< Parameter vectors for all points.
		///<
	double* y; ///< Parameter vector costs.
		///<
	double* x0; // Centroid parameter vector.
		///<
	double* x1; ///< Temporary parameter vector 1.
		///<
	double y1; ///< Cost of temporary parameter vector 1.
		///<
	double* x2; ///< Temporary parameter vector 2.
		///<
	double* rx; ///< Lowest cost parameter vector used during reduction.
		///<
	int rj; ///< Current vector during reduction.
		///<
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
		///<

	/**
	 * Algorithm's state automata states.
	 */

	enum EState
	{
		stReflection, // Reflection.
		stExpansion, // Expansion.
		stContraction, // Contraction.
		stReduction // Reduction.
	};

	EState State; ///< Current optimization state.
		///<

	virtual void initBuffers( const int aParamCount, const int aPopSize )
	{
		CBiteOptBase :: initBuffers( aParamCount, aPopSize );

		x = PopParams;
		y = PopCosts;
		x0 = CentParams;
		x1 = TmpParams;
		x2 = new double[ N ];
	}

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] x2;
	}

	/**
	 * Function finds the highest-cost vector.
	 */

	void findhi()
	{
		if( y[ 0 ] > y[ 1 ])
		{
			xhi = 0;
			xhi2 = 1;
		}
		else
		{
			xhi = 1;
			xhi2 = 0;
		}

		int j;

		for( j = 2; j < M; j++ )
		{
			if( y[ j ] > y[ xhi ])
			{
				xhi2 = xhi;
				xhi = j;
			}
			else
			if( y[ j ] > y[ xhi2 ])
			{
				xhi2 = j;
			}
		}
	}

	/**
	 * Function calculates the centroid vector.
	 */

	void calccent()
	{
		findhi();
		int i;

		for( i = 0; i < N; i++ )
		{
			x0[ i ] = 0.0;
		}

		int j;

		for( j = 0; j < M; j++ )
		{
			if( j != xhi )
			{
				const double* const xx = x[ j ];

				for( i = 0; i < N; i++ )
				{
					x0[ i ] += xx[ i ];
				}
			}
		}

		for( i = 0; i < N; i++ )
		{
			x0[ i ] *= M1m;
		}
	}

	/**
	 * Function replaces the highest-cost vector with a new vector.
	 *
	 * @param ip Input vector.
	 * @param cost Input vector's cost.
	 */

	void copy( const double* const ip, const double cost )
	{
		y[ xhi ] = cost;
		double* const xH = x[ xhi ];
		int i;

		memcpy( xH, ip, N * sizeof( xH[ 0 ]));
		findhi();

		const double* const nxH = x[ xhi ];

		if( xH != nxH )
		{
			for( i = 0; i < N; i++ )
			{
				x0[ i ] += ( xH[ i ] - nxH[ i ]) * M1m;
			}
		}
	}

	/**
	 * Function evaluates parameter vector and applies value range wrapping,
	 * also records a new best solution.
	 *
	 * @param rnd Random number generator.
	 * @param p Parameter vector to evaluate.
	 * @param OutCost If not NULL, pointer to variable that receives cost
	 * of the newly-evaluated solution.
	 * @param OutValues If not NULL, pointer to array that receives
	 * newly-evaluated parameter vector, in real scale.
	 */

	double eval( CBiteRnd& rnd, const double* p,
		double* const OutCost = NULL, double* const OutValues = NULL )
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			NewValues[ i ] = wrapParamReal( rnd, p[ i ], i );
		}

		const double cost = optcost( NewValues );

		if( OutCost != NULL )
		{
			*OutCost = cost;
		}

		if( OutValues != NULL )
		{
			memcpy( OutValues, NewValues, N * sizeof( OutValues[ 0 ]));
		}

		updateBestCost( cost, NewValues );

		return( cost );
	}
};

#endif // NMSOPT_INCLUDED
