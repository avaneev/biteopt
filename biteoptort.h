//$ nocpp

/**
 * @file biteoptort.h
 *
 * @brief The inclusion file for the CBiteOptOrt class.
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

#ifndef BITEOPTORT_INCLUDED
#define BITEOPTORT_INCLUDED

#include "biteaux.h"

/**
 * Rotation matrix calculation class, based on the Eigen decomposition of the
 * covariance matrix. Used to "orthogonalize" sample population. Applies
 * weighting to population.
 */

class CBiteOptOrt
{
public:
	double CentPow; ///< Centroid weighting power coefficient.
		///<
	double CovUpdSlow; ///< Covariance matrix update rate, slow.
		///<
	double CovUpdFast; ///< Covariance matrix update rate, fast.
		///<
	double SigmaMulBase1; ///< Sigma damping multiplier, for sphere.
		///<
	double SigmaMulBase2; ///< Sigma damping multiplier, for needle.
		///<
	double SigmaMulExp2; ///< Sigma expansion multiplier, for needle.
		///<

	CBiteOptOrt()
		: ParamCount( 0 )
		, PopSize( 0 )
		, CovParamsBuf( NULL )
		, CovParams( NULL )
		, BParamsBuf( NULL )
		, BParams( NULL )
		, DParams( NULL )
		, DParamsN( NULL )
		, CentParams( NULL )
		, PrevCentParams( NULL )
		, WPopCent( NULL )
		, WPopCov( NULL )
		, PopParamsBuf( NULL )
		, PopParams( NULL )
		, SortedParams( NULL )
		, TmpParams( NULL )
	{
		CentPow = 6.0;
		CovUpdSlow = 5.0;
		CovUpdFast = 18.0;
		SigmaMulBase1 = 0.4;
		SigmaMulBase2 = 0.9;
		SigmaMulExp2 = 1.15;
	}

	/**
	 * Function updates dimensions of *this object.
	 *
	 * @param aParamCount The number of elements in the parameter vector.
	 * @param aPopSize Size of population used to calculate the rotation
	 * matrix.
	 * @param aEvalFac The multiplier of the population size - the actual
	 * number of function evaluations performed, >=1.
	 */

	void updateDims( const int aParamCount, const int aPopSize,
		const double aEvalFac )
	{
		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		PopSize = aPopSize;
		EvalFac = aEvalFac;
		CovParamsBuf = new double[ ParamCount * ParamCount ];
		CovParams = new double*[ ParamCount ];
		BParamsBuf = new double[ ParamCount * ParamCount ];
		BParams = new double*[ ParamCount ];
		DParams = new double[ ParamCount ];
		DParamsN = new double[ ParamCount ];
		CentParams = new double[ ParamCount ];
		PrevCentParams = new double[ ParamCount ];
		WPopCent = new double[ PopSize ];
		WPopCov = new double[ PopSize ];
		PopParamsBuf = new double[ ParamCount * PopSize ];
		PopParams = new double*[ ParamCount ];
		SortedParams = new double*[ PopSize ];
		TmpParams = new double[ ParamCount ];

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CovParams[ i ] = CovParamsBuf + i * ParamCount;
			BParams[ i ] = BParamsBuf + i * ParamCount;
			PopParams[ i ] = PopParamsBuf + i * PopSize;
		}
	}

	/**
	 * Function updates centroid and covariance estimation weights.
	 */

	void updateWeights()
	{
		// Calculate weights for centroid and covariance calculation.

		double s = 0.0;
		int i;

		for( i = 0; i < UsePopSize; i++ )
		{
			const double l = 1.0 - (double) i / ( UsePopSize * EvalFac );
			WPopCent[ i ] = pow( l, CentPow );
			s += WPopCent[ i ];
		}

		for( i = 0; i < UsePopSize; i++ )
		{
			WPopCent[ i ] /= s;
			WPopCov[ i ] = sqrt( WPopCent[ i ]);
		}
	}

	/**
	 * Function initializes *this object to "no rotation" state, with the
	 * specified initial centroid and sigma.
	 *
	 * @param InitCent Initial centroids per parameter (NULL for origin).
	 * @param InitSigma Initial sigmas per parameter (NULL for 0.5).
	 * @return Population size to use.
	 */

	int init( const double* const InitCent = NULL,
		const double* const InitSigma = NULL )
	{
		memset( CovParamsBuf, 0, ParamCount * ParamCount *
			sizeof( CovParamsBuf[ 0 ]));

		memset( BParamsBuf, 0, ParamCount * ParamCount *
			sizeof( BParamsBuf[ 0 ]));

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CovParams[ i ][ i ] = 1.0;
			BParams[ i ][ i ] = 1.0;
			CentParams[ i ] = ( InitCent == NULL ? 0.0 : InitCent[ i ]);
			DParams[ i ] = ( InitSigma == NULL ? 0.25 : InitSigma[ i ]);
			DParamsN[ i ] = DParams[ i ];
		}

		spc = 0.0;
		UsePopSize = PopSize;
		updateWeights();

		return( UsePopSize );
	}

	/**
	 * Function performs rotation matrix update.
	 *
	 * @param CurParams Population parameter vectors.
	 * @param PopOrder Population ordering, used to obtain unordered
	 * population indices.
	 * @return Current population size.
	 */

	int update( double** const CurParams, const int* const PopOrder )
	{
		// Prepare PopParams (vector of per-parameter population deviations),
		// later used to calculate weighted covariances, use current centroid
		// vector.

		int i;
		int j;

		for( i = 0; i < UsePopSize; i++ )
		{
			SortedParams[ i ] = CurParams[ PopOrder[ i ]]; // For fast access.
		}

		for( i = 0; i < ParamCount; i++ )
		{
			double* const op = PopParams[ i ];
			const double c = CentParams[ i ];

			for( j = 0; j < UsePopSize; j++ )
			{
				op[ j ] = ( SortedParams[ j ][ i ] - c ) * WPopCov[ j ];
			}
		}

		// Store older centroid and calculate new centroid.

		memcpy( PrevCentParams, CentParams,
			ParamCount * sizeof( PrevCentParams[ 0 ]));

		calcCent();

		// Update covariance matrix, the left-handed triangle only. Uses leaky
		// integrator averaging filter. "avgc" selects corner frequency
		// of the filter.

		// Select covariance update rate based on geometry parameter,
		// depends on sphericity.

		const double spc2 = spc * spc;
		const double spc3 = spc2 * spc;
		const double spc4 = spc2 * spc2;

		const double avgc = 1.0 - pow( 0.01,
			( CovUpdFast * spc4 + CovUpdSlow * ( 1.0 - spc4 )) /
			( UsePopSize * EvalFac ));

		// Perform covariance matrix update using a leaky integrator filter.

		for( j = 0; j < ParamCount; j++ )
		{
			double* const cp = CovParams[ j ];

			for( i = 0; i <= j; i++ )
			{
				const double cov = dotp( PopParams[ j ], PopParams[ i ]);
				cp[ i ] += ( cov - cp[ i ]) * avgc;
			}
		}

		// Calculate rotation matrix and standard deviations.

		eigen();

		// Calculate distribution's geometry parameter.

		double maxd = 1e-100; // Maximal deviation among all parameters.

		for( i = 0; i < ParamCount; i++ )
		{
			if( DParams[ i ] > maxd )
			{
				maxd = DParams[ i ];
			}
		}

		maxd = 1.0 / maxd;
		spc = 0.0;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = DParams[ i ] * maxd;
			spc += TmpParams[ i ];
		}

		spc /= ParamCount;

		// Apply extension or contraction to positive and negative halfs of
		// standard deviations, depending on centroid step size.

		for( i = 0; i < ParamCount; i++ )
		{
			PrevCentParams[ i ] = CentParams[ i ] - PrevCentParams[ i ];
		}

		ortnc( PrevCentParams, DParamsN ); // Orthogonalize centroid step.

		for( i = 0; i < ParamCount; i++ )
		{
			const double d = DParamsN[ i ];

			DParamsN[ i ] = DParams[ i ] - d;
			DParams[ i ] += d;
		}

		// Modify standard deviations, apply reduction, or extension to
		// overly shrunk dimensions. Depends on sphericity.

		for( i = 0; i < ParamCount; i++ )
		{
			const double c = 1.0 - TmpParams[ i ];
			const double m = ( SigmaMulBase2 +
				( SigmaMulExp2 - SigmaMulBase2 ) * c * c * c ) *
				( 1.0 - spc3 ) + SigmaMulBase1 * spc3;

			DParams[ i ] *= m;
			DParamsN[ i ] *= m;
		}

		return( UsePopSize );
	}

	/**
	 * Function applies orthogonalization (rotation) to the parameter vector,
	 * subtracting centroid vector.
	 *
	 * @param ip Input unorthogonalized vector.
	 * @param[out] op Output orthogonalized vector, can be equal to "ip".
	 */

	void ort( const double* const ip, double* const op )
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = ip[ i ] - CentParams[ i ];
		}

		for( i = 0; i < ParamCount; i++ )
		{
			double s = 0.0;
			int j;

			for( j = 0; j < ParamCount; j++ )
			{
				s += BParams[ j ][ i ] * TmpParams[ j ];
			}

			op[ i ] = s;
		}
	}

	/**
	 * Function applies orthogonalization (rotation) to the parameter vector,
	 * without subtracting centroid vector.
	 *
	 * @param ip Input unorthogonalized vector.
	 * @param[out] op Output orthogonalized vector, should not be equal to
	 * "ip".
	 */

	void ortnc( const double* const ip, double* const op )
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			double s = 0.0;
			int j;

			for( j = 0; j < ParamCount; j++ )
			{
				s += BParams[ j ][ i ] * ip[ j ];
			}

			op[ i ] = s;
		}
	}

	/**
	 * Function undos orthogonalization (rotation) of the parameter vector,
	 * adding centroid vector.
	 *
	 * @param ip Input orthogonalized vector.
	 * @param[out] op Output unorthogonalized vector, can be equal to "ip".
	 */

	void unort( const double* const ip, double* const op )
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = ip[ i ];
		}

		for( i = 0; i < ParamCount; i++ )
		{
			double s = 0.0;
			int j;

			for( j = 0; j < ParamCount; j++ )
			{
				s += BParams[ i ][ j ] * TmpParams[ j ];
			}

			op[ i ] = s + CentParams[ i ];
		}
	}

	/**
	 * Function "samples" new random population vector making a random draw
	 * from the current distribution.
	 *
	 * @param rnd Random number generator.
	 * @param op Resulting vector.
	 */

	void sample( CBiteRnd& rnd, double* const op ) const
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = getGaussian( rnd );

			TmpParams[ i ] *= ( TmpParams[ i ] < 0.0 ?
				DParamsN[ i ] : DParams[ i ]);
		}

		for( i = 0; i < ParamCount; i++ )
		{
			double s = 0.0;
			int j;

			for( j = 0; j < ParamCount; j++ )
			{
				s += BParams[ i ][ j ] * TmpParams[ j ];
			}

			op[ i ] = s + CentParams[ i ];
		}
	}

protected:
	int ParamCount; ///< The number of parameters in parameter vector.
		///<
	int PopSize; ///< Population size (max).
		///<
	double EvalFac; ///< Function evaluations factor.
		///<
	int UsePopSize; ///< Current population size.
		///<
	double* CovParamsBuf; ///< CovParams buffer.
		///<
	double** CovParams; ///< Covariance matrix.
		///<
	double* BParamsBuf; ///< BParams buffer.
		///<
	double** BParams; ///< Rotation matrix.
		///<
	double* DParams; ///< Std. deviations vector.
		///<
	double* DParamsN; ///< Std. deviations vector, for negative values.
		///<
	double* CentParams; ///< Centroid vector, weighted.
		///<
	double* PrevCentParams; ///< Previous centroid vector.
		///<
	double* WPopCent; ///< Weighting coefficients for ordered population, for
		///< centroid calculation.
		///<
	double* WPopCov; ///< Weighting coefficients for covariance calculation,
		///< squared.
		///<
	double* PopParamsBuf; ///< PopParams buffer.
		///<
	double** PopParams; ///< Population vectors per parameter (deviations,
		///< weighted).
		///<
	double** SortedParams; ///< Temporary sorted parameter vectors.
		///<
	double* TmpParams; ///< Temporary parameter vector.
		///<
	double spc; ///< Distribution's sphericity coefficient. 1 - fully
		///< spherical.

	/**
	 * Function deletes previously allocated buffers.
	 */

	void deleteBuffers()
	{
		delete[] CovParamsBuf;
		delete[] CovParams;
		delete[] BParamsBuf;
		delete[] BParams;
		delete[] DParams;
		delete[] DParamsN;
		delete[] CentParams;
		delete[] PrevCentParams;
		delete[] WPopCent;
		delete[] WPopCov;
		delete[] PopParamsBuf;
		delete[] PopParams;
		delete[] SortedParams;
		delete[] TmpParams;
	}

	/**
	 * Householder triagonalization routine from JAMA package.
	 */

	static void tred2( const int n, double* d, double** V, double* e )
	{
		//  This is derived from the Algol procedures tred2 by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int j = 0; j < n; j++) {
			d[j] = V[n-1][j];
		}

		// Householder reduction to tridiagonal form.

		for (int i = n-1; i > 0; i--) {

			// Scale to avoid under/overflow.

			double scale = 0.0;
			double h = 0.0;
			for (int k = 0; k < i; k++) {
				scale = scale + abs(d[k]);
			}
			if (scale == 0.0) {
				e[i] = d[i-1];
				for (int j = 0; j < i; j++) {
					d[j] = V[i-1][j];
					V[i][j] = 0.0;
					V[j][i] = 0.0;
				}
			} else {

				// Generate Householder vector.

				for (int k = 0; k < i; k++) {
					d[k] /= scale;
					h += d[k] * d[k];
				}
				double f = d[i-1];
				double g = sqrt(h);
				if (f > 0) {
					g = -g;
				}
				e[i] = scale * g;
				h = h - f * g;
				d[i-1] = f - g;
				for (int j = 0; j < i; j++) {
					e[j] = 0.0;
				}

				// Apply similarity transformation to remaining columns.

				for (int j = 0; j < i; j++) {
					f = d[j];
					V[j][i] = f;
					g = e[j] + V[j][j] * f;
					for (int k = j+1; k <= i-1; k++) {
						g += V[k][j] * d[k];
						e[k] += V[k][j] * f;
					}
					e[j] = g;
				}
				f = 0.0;
				for (int j = 0; j < i; j++) {
					e[j] /= h;
					f += e[j] * d[j];
				}
				double hh = f / (h + h);
				for (int j = 0; j < i; j++) {
					e[j] -= hh * d[j];
				}
				for (int j = 0; j < i; j++) {
					f = d[j];
					g = e[j];
					for (int k = j; k <= i-1; k++) {
						V[k][j] -= (f * e[k] + g * d[k]);
					}
					d[j] = V[i-1][j];
					V[i][j] = 0.0;
				}
			}
			d[i] = h;
		}

		// Accumulate transformations.

		for (int i = 0; i < n-1; i++) {
			V[n-1][i] = V[i][i];
			V[i][i] = 1.0;
			double h = d[i+1];
			if (h != 0.0) {
				for (int k = 0; k <= i; k++) {
					d[k] = V[k][i+1] / h;
				}
				for (int j = 0; j <= i; j++) {
					double g = 0.0;
					for (int k = 0; k <= i; k++) {
						g += V[k][i+1] * V[k][j];
					}
					for (int k = 0; k <= i; k++) {
						V[k][j] -= g * d[k];
					}
				}
			}
			for (int k = 0; k <= i; k++) {
				V[k][i+1] = 0.0;
			}
		}
		for (int j = 0; j < n; j++) {
			d[j] = V[n-1][j];
			V[n-1][j] = 0.0;
		}
		V[n-1][n-1] = 1.0;
		e[0] = 0.0;
	}

	/**
	 * sqrt(a^2 + b^2) with possible under/overflow. Adequate for this
	 * method's purposes.
	 */

	static double hypot_( double a, double b )
	{
		return( sqrt( a * a + b * b ));
	}

	/**
	 * Function returns maximum of 2 values.
	 */

	static double max( const double a, const double b )
	{
		return( a > b ? a : b );
	}

	/**
	 * Symmetric tridiagonal QL algorithm from JAMA package. Std. deviation
	 * vector sorting removed.
	 */

	static void tql2( const int n, double* d, double** V, double* e )
	{
		//  This is derived from the Algol procedures tql2, by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int i = 1; i < n; i++) {
			e[i-1] = e[i];
		}
		e[n-1] = 0.0;

		double f = 0.0;
		double tst1 = 0.0;
		double eps = pow(2.0,-52.0);
		for (int l = 0; l < n; l++) {

			// Find small subdiagonal element

			tst1 = max(tst1,abs(d[l]) + abs(e[l]));
			int m = l;
			while (m < n) {
				if (abs(e[m]) <= eps*tst1) {
					break;
				}
				m++;
			}

			// If m == l, d[l] is an eigenvalue,
			// otherwise, iterate.

			if (m > l) {
				do {
					// Compute implicit shift

					double g = d[l];
					double p = (d[l+1] - g) / (2.0 * e[l]);
					double r = hypot_(p,1.0);
					if (p < 0) {
						r = -r;
					}
					d[l] = e[l] / (p + r);
					d[l+1] = e[l] * (p + r);
					double dl1 = d[l+1];
					double h = g - d[l];
					for (int i = l+2; i < n; i++) {
						d[i] -= h;
					}
					f = f + h;

					// Implicit QL transformation.

					p = d[m];
					double c = 1.0;
					double c2 = c;
					double c3 = c;
					double el1 = e[l+1];
					double s = 0.0;
					double s2 = 0.0;
					for (int i = m-1; i >= l; i--) {
						c3 = c2;
						c2 = c;
						s2 = s;
						g = c * e[i];
						h = c * p;
						r = hypot_(p,e[i]);
						e[i+1] = s * r;
						s = e[i] / r;
						c = p / r;
						p = c * d[i] - s * g;
						d[i+1] = h + s * (c * g + s * d[i]);

						// Accumulate transformation.

						for (int k = 0; k < n; k++) {
							h = V[k][i+1];
							V[k][i+1] = s * V[k][i] + c * h;
							V[k][i] = c * V[k][i] - s * h;
						}
					}
					p = -s * s2 * c3 * el1 * e[l] / dl1;
					e[l] = s * p;
					d[l] = c * p;

					// Check for convergence.

				} while (abs(e[l]) > eps*tst1);
			}
			d[l] = d[l] + f;
			e[l] = 0.0;
		}
	}

	/**
	 * Function calculates eigenpairs of the current covariance matrix and
	 * updates DParams standard deviations vector.
	 */

	void eigen()
	{
		// Replicate left-handed triangular matrix.

		int i;
		int j;

		for( j = 0; j < ParamCount; j++ )
		{
			for( i = 0; i <= j; i++ )
			{
				BParams[ j ][ i ] = BParams[ i ][ j ] = CovParams[ j ][ i ];
			}
		}

		// Calculate eigen-pairs.

		tred2( ParamCount, DParams, BParams, TmpParams );
		tql2( ParamCount, DParams, BParams, TmpParams );

		for( i = 0; i < ParamCount; i++ )
		{
			if( DParams[ i ] < 1e-30 )
			{
				DParams[ i ] = 1e-15;
			}
			else
			{
				DParams[ i ] = sqrt( DParams[ i ]);
			}
		}
	}

	/**
	 * Function calculates centroid vector of population, with weighting.
	 * Requires SortedParams filled with pointers to ordered population
	 * vectors.
	 */

	void calcCent()
	{
		const double* ip = SortedParams[ 0 ];
		double w = WPopCent[ 0 ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] = ip[ i ] * w;
		}

		int j;

		for( j = 1; j < UsePopSize; j++ )
		{
			ip = SortedParams[ j ];
			w = WPopCent[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] += ip[ i ] * w;
			}
		}
	}

	/**
	 * Function calculates dot product of two per-parameter population
	 * vectors.
	 *
	 * @param p1 Parameter A's population vector.
	 * @param p2 Parameter B's population vector.
	 */

	double dotp( const double* const p1, const double* const p2 ) const
	{
		double s = 0.0;
		int i;

		for( i = 0; i < UsePopSize; i++ )
		{
			s += p1[ i ] * p2[ i ];
		}

		return( s );
	}

	/**
	 * Function generates a Gaussian-distributed pseudo-random number with
	 * mean=0 and std.dev=1.
	 *
	 * @param rnd Uniform PRNG.
	 */

	static double getGaussian( CBiteRnd& rnd )
	{
		double q, u, v;

		do
		{
			u = rnd.getRndValue();
			v = rnd.getRndValue();

			if( u <= 0.0 || v <= 0.0 )
			{
				u = 1.0;
				v = 1.0;
			}

			v = 1.7156 * ( v - 0.5 );
			const double x = u - 0.449871;
			const double y = fabs( v ) + 0.386595;
			q = x * x + y * ( 0.19600 * y - 0.25472 * x );

			if( q < 0.27597 )
			{
				break;
			}
		} while(( q > 0.27846 ) || ( v * v > -4.0 * log( u ) * u * u ));

		return( v / u );
	}
};

#endif // BITEOPTORT_INCLUDED
