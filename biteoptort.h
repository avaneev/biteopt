//$ nocpp

#ifndef BITEOPTORT_INCLUDED
#define BITEOPTORT_INCLUDED

#include <math.h>
#include <string.h>
#include "biternd.h"

/**
 * Rotation matrix calculation class, based on the Eigen decomposition of the
 * covariance matrix. Used to "orthogonalize" sample population.
 */

class CBiteOptOrt
{
public:
	CBiteOptOrt()
		: ParamCount( 0 )
		, PopSize( 0 )
		, CovParamsBuf( NULL )
		, CovParams( NULL )
		, BParamsBuf( NULL )
		, BParams( NULL )
		, DParams( NULL )
		, CentParams( NULL )
		, WPop( NULL )
		, WPopSqrt( NULL )
		, PopParamsBuf( NULL )
		, PopParams( NULL )
		, TmpParams( NULL )
	{
	}

	/**
	 * Function updates dimensions of *this object.
	 *
	 * @param aParamCount The number of elements in the parameter vector.
	 * @param aPopSize Size of population used to calculate the rotation
	 * matrix.
	 */

	void updateDims( const int aParamCount, const int aPopSize )
	{
		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		PopSize = aPopSize;
		CovParamsBuf = new double[ ParamCount * ParamCount ];
		CovParams = new double*[ ParamCount ];
		BParamsBuf = new double[ ParamCount * ParamCount ];
		BParams = new double*[ ParamCount ];
		DParams = new double[ ParamCount ];
		CentParams = new double[ ParamCount ];
		WPop = new double[ PopSize ];
		WPopSqrt = new double[ PopSize ];
		PopParamsBuf = new double[ ParamCount * PopSize ];
		PopParams = new double*[ ParamCount ];
		TmpParams = new double[ ParamCount ];

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CovParams[ i ] = CovParamsBuf + i * ParamCount;
			BParams[ i ] = BParamsBuf + i * ParamCount;
			PopParams[ i ] = PopParamsBuf + i * PopSize;
		}

		double s = 0.0;

		for( i = 0; i < PopSize; i++ )
		{
			const double x = 1.0 - 0.5 * i / ( PopSize - 1 );
			WPop[ i ] = x * x;
			s += WPop[ i ];
		}

		for( i = 0; i < PopSize; i++ )
		{
			WPop[ i ] /= s;
			WPopSqrt[ i ] = sqrt( WPop[ i ]);
		}
	}

	/**
	 * Function initializes *this object to "no rotation" state.
	 */

	void init()
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
			CentParams[ i ] = 0.0;
			DParams[ i ] = 0.5;
		}

		WAccum = 0.0;
	}

	/**
	 * Function updates weight accumulator and performs rotation matrix
	 * update.
	 *
	 * @param ordi Ordered population vector index which was updated.
	 * @param CurParams Population parameter vectors.
	 * @param[out] OrtParams Orthogonalized population parameter vectors.
	 * @param PopOrder Population ordering, used to obtain unordered
	 * population indices.
	 */

	void update( const int ordi, double** CurParams, double** OrtParams,
		const int* PopOrder )
	{
		WAccum += WPop[ ordi ];

		if( WAccum < 2.0 )
		{
			const int k = PopOrder[ ordi ];
			ort( CurParams[ k ], OrtParams[ k ]);
			return;
		}

		WAccum = 0.0;

		calcAvg( CurParams, PopOrder );

		// Prepare PopParams, vector of per-parameter population deviations,
		// later used to calculate covariances.

		int i;
		int j;

		for( i = 0; i < ParamCount; i++ )
		{
			const double avg = CentParams[ i ];
			double* const op = PopParams[ i ];

			for( j = 0; j < PopSize; j++ )
			{
				const int k = PopOrder[ j ];
				op[ j ] = ( CurParams[ k ][ i ] - avg ) * WPopSqrt[ j ];
			}
		}

		// Calculate covariance matrix, left-handed triangle only.

		for( j = 0; j < ParamCount; j++ )
		{
			for( i = 0; i <= j; i++ )
			{
				CovParams[ j ][ i ] = calcCov( PopParams[ j ], PopParams[ i ]);
//				CovParams[ j ][ i ] +=
//					( calcCov( PopParams[ j ], PopParams[ i ]) -
//					CovParams[ j ][ i ]) * 0.5;
			}
		}

		// Calculate rotation matrix.

		eigen();

		// Apply orthogonalization.

		for( j = 0; j < PopSize; j++ )
		{
			ort( CurParams[ j ], OrtParams[ j ]);
		}

		// Calculate standard deviations of each orthogonalized parameter.

		calcStdDev( OrtParams, PopOrder );
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
				s += BParams[ i ][ j ] * TmpParams[ j ];
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
				s += BParams[ j ][ i ] * TmpParams[ j ];
			}

			op[ i ] = s + CentParams[ i ];
		}
	}

	/**
	 * Function adds "noise" to the specified parameter vector. The noise
	 * approximately follows the probability distribution of population.
	 */

	void addNoise( CBiteRnd& rnd, double* const op,const double*mp )
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			TmpParams[ i ] = DParams[ i ] * rnd.getTPDFRaw() /
				rnd.getRawScale() * 2;
		}

		for( i = 0; i < ParamCount; i++ )
		{
			double s = 0.0;
			int j;

			for( j = 0; j < ParamCount; j++ )
			{
				s += BParams[ j ][ i ] * TmpParams[ j ];
			}

			op[ i ] = s + CentParams[ i ];
		}
	}

protected:
	int ParamCount; ///< The number of parameters in parameter vector.
		///<
	int PopSize; ///< Population size.
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
	double* CentParams; ///< Centroid vector.
		///<
	double* WPop; ///< Weighting coefficients for ordered population.
		///<
	double* WPopSqrt; ///< sqrt(WPop) for faster covariance calculation.
		///<
	double* PopParamsBuf; ///< PopParams buffer.
		///<
	double** PopParams; ///< Population vectors per parameter (deviations,
		///< weighted).
		///<
	double* TmpParams; ///< Temporary parameter vector.
		///<
	double WAccum; ///< Update weight accumulator. When this value reaches a
		///< certain threshold, the rotation matrix is recalculated.
		///<

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
		delete[] CentParams;
		delete[] WPop;
		delete[] WPopSqrt;
		delete[] PopParamsBuf;
		delete[] PopParams;
		delete[] TmpParams;
	}

	/**
	 * Householder triagonalization routine from JAMA package. Slightly
	 * optimized (and reduced precision) by replacing divisions by inverse
	 * multiplications.
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

				const double scalei = 1.0 / scale;
				for (int k = 0; k < i; k++) {
//					d[k] /= scale;
					d[k] *= scalei;
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
				const double hi = 1.0 / h;
				for (int j = 0; j < i; j++) {
//					e[j] /= h;
					e[j] *= hi;
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
				const double hi=1.0/h;
				for (int k = 0; k <= i; k++) {
//					d[k] = V[k][i+1] / h;
					d[k] = V[k][i+1] * hi;
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
	 * sqrt(a^2 + b^2) without under/overflow, commented out and used a simple
	 * more efficient version.
	 */

	static double hypot_( double a, double b )
	{
		return(sqrt(a*a+b*b));
/*		double r;
		if (abs(a) > abs(b)) {
			r = b/a;
			r = abs(a)*sqrt(1+r*r);
		} else if (b != 0) {
			r = a/b;
			r = abs(b)*sqrt(1+r*r);
		} else {
			r = 0.0;
		}
		return r;*/
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
	 * Function calculates eigenpairs of the current covariance matrix.
	 */

	void eigen()
	{
		// Replicate left-handed triangular matrix.

		int j;

		for( j = 0; j < ParamCount; j++ )
		{
			int i;

			for( i = 0; i <= j; i++ )
			{
				BParams[ j ][ i ] = BParams[ i ][ j ] = CovParams[ j ][ i ];
			}
		}

		tred2( ParamCount, DParams, BParams, TmpParams );
		tql2( ParamCount, DParams, BParams, TmpParams );
	}

	/**
	 * Function calculates centroid (average) vector of population, without
	 * weighting.
	 *
	 * @param CurParams Population parameter vectors.
	 * @param PopOrder Population ordering.
	 */

	void calcAvg( const double** const CurParams, const int* const PopOrder )
	{
		const double* ip = CurParams[ PopOrder[ 0 ]];
		double w = WPop[ 0 ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] = ip[ i ] * w;
		}

		int j;

		for( j = 1; j < PopSize; j++ )
		{
			ip = CurParams[ PopOrder[ j ]];
			w = WPop[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] += ip[ i ] * w;
			}
		}

/*		const double m = 1.0 / PopSize;

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] *= m;
		}*/
	}

	/**
	 * Function calculates weighted covariance of two per-parameter population
	 * vectors.
	 *
	 * @param p1 Parameter A's population vector.
	 * @param p2 Parameter B's population vector.
	 */

	double calcCov( const double* const p1, const double* const p2 ) const
	{
		double s = 0.0;
		int i;

		for( i = 0; i < PopSize; i++ )
		{
			s += p1[ i ] * p2[ i ];
		}

		return( s );
	}

	/**
	 * Function calculates standard deviations of parameters.
	 *
	 * @param OrtParams Orthogonalized population parameter vectors.
	 * @param PopOrder Population ordering.
	 */

	void calcStdDev( double** const OrtParams,
		const int* const PopOrder )
	{
		const double* ip = OrtParams[ PopOrder[ 0 ]];
		double w = WPop[ 0 ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			DParams[ i ] = ip[ i ] * ip[ i ] * w;
		}

		int j;

		for( j = 1; j < PopSize; j++ )
		{
			ip = OrtParams[ PopOrder[ j ]];
			w = WPop[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				DParams[ i ] += ip[ i ] * ip[ i ] * w;
			}
		}

		for( i = 0; i < ParamCount; i++ )
		{
			DParams[ i ] = sqrt( DParams[ i ]);
		}
	}
};

#endif // BITEOPTORT_INCLUDED
