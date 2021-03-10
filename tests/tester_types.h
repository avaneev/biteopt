#ifndef TESTER_TYPES_INCLUDED
#define TESTER_TYPES_INCLUDED

#include "../biteaux.h"

#if !defined( M_PI )
	#define M_PI 3.14159265358979324
#endif // !defined( M_PI )

CBiteRnd fnrnd( 1000000000 );

inline double sqr( const double x )
{
	return( x * x );
}

template< class T >
inline T roundf( const T d )
{
	return( d < 0.0 ? -floor( (T) 0.5 - d ) : floor( d + (T) 0.5 ));
}

/**
 * Structure holds details about test function.
 */

struct CTestFn
{
	const char* Name; ///< Function's name.
	int Dims; ///< The number of dimensions, 0 - unlimited.
	double RangeMin; ///< Dimensions range min.
	double RangeMax; ///< Dimensions range max.
	double OptValue; ///< Optimal value.

	double (*CalcFunc)( const double* const x, const int N ); ///< Calculation
		///< function.

	void (*ParamFunc)( double* const minv, double* const maxv,
		const int N, double* const optv ); ///< Range and optimal value
		///< function, can be NULL. "optv" value is initalized by OptValue
		///< before ParamFunc call, may not need update in ParamFunc.
		///<
};

#endif // TESTER_TYPES_INCLUDED
