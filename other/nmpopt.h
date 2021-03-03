//$ nocpp

/**
 * @file nmpopt.h
 *
 * @brief The inclusion file for the CNelderMeadPlusOpt class.
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
 * @version 2021.7
 */

#ifndef NMPOPT_INCLUDED
#define NMPOPT_INCLUDED

#include "../biteaux.h"

/**
 * Sequential Nelder-Mead simplex method with "Z" parameter improvement.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CNelderMeadPlusOpt : public CBiteOptBase
{
public:
	CNelderMeadPlusOpt()
		: Z( 0 )
		, x2( NULL )
	{
	}

	/**
	 * Function updates dimensionality of *this object.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param aZ Z-parameter, >= 1. Improves convergence and increases
	 * convergence time.
	 */

	void updateDims( const int aParamCount, const int aZ = 1 )
	{
		if( aParamCount * aZ == ParamCount && aZ == Z )
		{
			return;
		}

		deleteBuffers();

		Z = aZ;
		N = aParamCount * Z;
		M = N + 1;

		initBaseBuffers( N, M );

		ParamCount0 = aParamCount;
		M1m = 1.0 / ( M - 1 );
		x = CurParams;
		y = CurCosts;
		x0 = CentParams;
		x1 = CurParams[ M ];
		x2 = new double[ N ];
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 * @param InitRadius Initial radius, relative to the default value
	 * (<= 1.0).
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		resetCommonVars();
		updateDiffValues();

		// Initialize parameter vectors, costs and centroid.

		double* const xx = x[ 0 ];
		const double w = 1.0 / Z;
		int i;
		int k;

		if( InitParams != NULL )
		{
			for( i = 0; i < ParamCount0; i++ )
			{
				const double v = ( InitParams[ i ] - MinValues[ i ]) /
					DiffValues[ i ];

				for( k = 0; k < Z; k++ )
				{
					xx[ i * Z + k ] = v * w;
				}
			}
		}
		else
		{
			for( i = 0; i < ParamCount0; i++ )
			{
				for( k = 0; k < Z; k++ )
				{
					xx[ i * Z + k ] = 0.5 * w;
				}
			}
		}

		xlo = 0;

		const double sd = InitRadius * 0.25 * w;
		int j;

		for( j = 1; j < M; j++ )
		{
			double* const xj = x[ j ];

			for( i = 0; i < ParamCount0; i++ )
			{
				for( k = 0; k < Z; k++ )
				{
					const double v = xx[ i * Z + k ] +
						getGaussian( rnd ) * sd;

					xj[ i * Z + k ] = v;
				}
			}
		}

		State = stReflection;
		DoInitEvals = true;
		InitEvalIndex = 0;
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
		int i;

		if( DoInitEvals )
		{
			y[ InitEvalIndex ] = eval( rnd, x[ InitEvalIndex ], OutCost,
				OutParams );

			if( y[ InitEvalIndex ] < y[ xlo ])
			{
				xlo = InitEvalIndex;
			}

			InitEvalIndex++;

			if( InitEvalIndex == M )
			{
				DoInitEvals = false;
				calccent();
			}

			return( 0 );
		}

		StallCount++;

		const double alpha = 1.0; // Reflection coeff.
		const double gamma = 2.0; // Expansion coeff.
		const double rho = -0.5; // Contraction coeff.
		const double sigma = 0.5; // Reduction coeff.

		double* const xH = x[ xhi ]; // Highest cost parameter vector.

		switch( State )
		{
			case stReflection:
			{
				for( i = 0; i < N; i++ )
				{
					x1[ i ] = x0[ i ] + alpha * ( x0[ i ] - xH[ i ]);
				}

				y1 = eval( rnd, x1, OutCost, OutParams );

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

				const double y2 = eval( rnd, x2, OutCost, OutParams );
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

				const double y2 = eval( rnd, x2, OutCost, OutParams );

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

				y[ rj ] = eval( rnd, xx, OutCost, OutParams );

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
	int ParamCount0; ///< Original parameter count.
		///<
	int Z; ///< Z-parameter.
		///<
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
	int InitEvalIndex; ///< Current initial population index.
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

	/**
	 * Function deletes previously allocated buffers.
	 */

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

		for( i = 0; i < N; i++ )
		{
			xH[ i ] = ip[ i ];
		}

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
	 * @param OutParams If not NULL, pointer to array that receives
	 * newly-evaluated parameter vector, in normalized scale.
	 */

	double eval( CBiteRnd& rnd, const double* p,
		double* const OutCost = NULL, double* const OutParams = NULL )
	{
		int i;

		for( i = 0; i < ParamCount0; i++ )
		{
			double v = 0.0;
			int k;

			for( k = 0; k < Z; k++ )
			{
				v += *p;
				p++;
			}

			NewParams[ i ] = wrapParam( rnd, v );
		}

		if( OutParams != NULL )
		{
			memcpy( OutParams, NewParams,
				ParamCount0 * sizeof( NewParams[ 0 ]));
		}

		for( i = 0; i < ParamCount0; i++ )
		{
			NewParams[ i ] = getRealValue( NewParams, i );
		}

		const double cost = optcost( NewParams );

		if( OutCost != NULL )
		{
			*OutCost = cost;
		}

		updateBestCost( cost, NewParams );

		return( cost );
	}
};

#endif // NMPOPT_INCLUDED
