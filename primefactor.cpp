// Biprime-number factoring method. Uses BiteOptDeep with a huge depth factor.
// This is purely a proof-of-concept example, may not work for factoring of
// huge numbers, so far works for 19-digit numbers, does not work for all
// numbers. Since this method is probabilistic, various variables should be
// changed to obtain solutions: depth, attempt count, MaxFactor/MaxDiff
// values. The resulting "cost" is equal to 0 when a precise solution was
// found.
//
// This method is highly-parallelizable, each attempt can be executed in
// parallel in its own thread. Tests indicated that this method finds factors
// extremely fast when factors have similar magnitude.

#include <stdio.h>
#include <stdint.h>
#include "biteopt.h"

const uint64_t val = 1012322327LL * 1115382761LL; //322532699LL * 512532767LL;//7152347LL * 3152399LL;//n=pq
const int Depth = 100; // BiteOptDeep depth factor, affects lower bound of
	// achieved cost function. Possibly should be increased if a larger number
	// have to be factored.
const int Attempts = 200000; // Factorization attempts to perform.
const int MaxIters = (int) ( 2000000.0 * sqrt( (double) Depth )); // Maximal
	// number of iterations to perform.
const int StallFactor = 50; // Factor to stop attempt at when optimizer stalls.
const int N = 2; // The number of dimensions in optimizer.
const double MaxFactor = 2000000000.0 / 2; // Max value of factors.
const double MaxDiff = 300000000.0 / 2; // Max difference between factors.
bool DoPrint = false;

/**
 * Optimizer class.
 */

class CPrimeFactorizerOpt : public CBiteOptDeep
{
public:
	CPrimeFactorizerOpt()
	{
		updateDims( N, Depth );
	}

	virtual void getMinValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = 0.0;
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			p[ i ] = 1.0;
		}
	}

	virtual double optcost( const double* const p )
	{
		uint64_t a1 = (uint64_t) ( p[ 0 ] * MaxFactor );
		uint64_t a2 = a1 + (uint64_t) ( p[ 1 ] * MaxDiff );

		a1 = ( a1 << 1 ) | 1; // Produce odd factor.
		a2 = ( a2 << 1 ) | 1; // Produce odd factor.

		const uint64_t nval = a1 * a2;
		const uint64_t d = ( val > nval ? val - nval : nval - val );

		if( DoPrint )
		{
			printf( "    f1 = %llu\n", a1 );
			printf( "    f2 = %llu\n", a2 );
			printf( "n = pq = %llu\n", val );
			printf( " f1*f2 = %llu\n", nval );
		}

		return( (double) d );
	}
};

int main()
{
	int j;

	for( j = 0; j < Attempts; j++ )
	{
		CBiteRnd rnd;
		rnd.init( j + 1 );

		CPrimeFactorizerOpt opt;
		opt.init( rnd );
		const int MaxStallCount = opt.getInitEvals() / Depth * StallFactor;

		int i;

		for( i = 0; i < MaxIters; i++ )
		{
			if( opt.optimize( rnd ) > MaxStallCount )
			{
				// Optimization stalled.

				break;
			}
		}

		printf( "attempt\t%i\tcost\t%.8g\titers\t%i\n", j + 1,
			opt.getBestCost(), i );

		if( opt.getBestCost() == 0.0 )
		{
			printf( "Solution found\n" );

			DoPrint = true;
			opt.optcost( opt.getBestParams() );
			break;
		}
	}
}
