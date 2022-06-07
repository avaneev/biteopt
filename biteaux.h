//$ nocpp

/**
 * @file biteaux.h
 *
 * @brief The inclusion file for the CBiteRnd, CBiteOptPop, CBiteOptParPops,
 * CBiteOptInterface, and CBiteOptBase classes.
 *
 * @section license License
 * 
 * Copyright (c) 2016-2022 Aleksey Vaneev
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
 * @version 2022.16
 */

#ifndef BITEAUX_INCLUDED
#define BITEAUX_INCLUDED

#include <stdint.h>
#include <math.h>
#include <string.h>

/**
 * Type for an externally-provided random number generator, to be used instead
 * of the default PRNG. Note that if the external produces 64-bit random
 * values, they can be safely truncated/typecasted to the "uint32_t" type. If
 * the external PRNG produces floating-point values, they should be scaled to
 * the 32-bit unsigned integer range. 32-bit PRNG output is required for
 * compatibility with older 32-bit PRNGs (but CBiteRnd is 64-bit PRNG).
 */

typedef uint32_t( *biteopt_rng )( void* rng_data );

/**
 * Class that implements a pseudo-random number generator (PRNG). The default
 * implementation includes a fast high-quality PRNG (2^159 period). See
 * https://github.com/avaneev/prvhash for more details.
 */

class CBiteRnd
{
public:
	/**
	 * Default constructor, calls the init() function.
	 */

	CBiteRnd()
	{
		init( 1 );
	}

	/**
	 * Constructor, calls the init() function.
	 *
	 * @param NewSeed Random seed value.
	 */

	CBiteRnd( const int NewSeed )
	{
		init( NewSeed );
	}

	/**
	 * Function initializes *this PRNG object.
	 *
	 * @param NewSeed New random seed value. Ignored, if "arf" is non-NULL.
     * @param arf External random number generator to use; NULL: use the
	 * default PRNG. Note that the external RNG should be seeded externally.
	 * @param ardata Data pointer to pass to the "arf" function.
	 */

	void init( const int NewSeed, biteopt_rng const arf = NULL,
		void* const ardata = NULL )
	{
		rf = arf;
		rdata = ardata;

		BitsLeft = 0;
		Seed = (uint64_t) NewSeed;
		lcg = 0;
		Hash = 0;

		// Skip first values to make PRNG "settle down".

		int i;

		for( i = 0; i < 5; i++ )
		{
			advance();
		}
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	double get()
	{
		return(( advance() >> ( 64 - 53 )) * 0x1p-53 );
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1). N1 denotes the
	 * number of bins, not the maximal returned value.
	 */

	int getInt( const int N1 )
	{
		return( (int) ( get() * N1 ));
	}

	/**
	 * @return Random number in the range [0; 1). Beta distribution with
	 * Alpha=0.5, Beta=1 (squared).
	 */

	double getSqr()
	{
		const double v = ( advance() >> ( 64 - 53 )) * 0x1p-53;

		return( v * v );
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1). Beta distribution
	 * with Alpha=0.5, Beta=1 (squared). N1 denotes the number of bins, not
	 * the maximal returned value.
	 */

	int getSqrInt( const int N1 )
	{
		return( (int) ( getSqr() * N1 ));
	}

	/**
	 * @return Uniformly-distributed random number in the "raw" scale.
	 */

	uint64_t getRaw()
	{
		return( advance() );
	}

	/**
	 * @return TPDF random number in the range (-1; 1).
	 */

	double getTPDF()
	{
		const int64_t v1 = (int64_t) ( advance() >> ( 64 - 53 ));
		const int64_t v2 = (int64_t) ( advance() >> ( 64 - 53 ));

		return(( v1 - v2 ) * 0x1p-53 );
	}

	/**
	 * Function generates a Gaussian-distributed pseudo-random number with
	 * mean=0 and std.dev=1.
	 *
	 * Algorithm is adopted from "Leva, J. L. 1992. "A Fast Normal Random
	 * Number Generator", ACM Transactions on Mathematical Software, vol. 18,
	 * no. 4, pp. 449-453".
	 */

	double getGaussian()
	{
		double q, u, v;

		do
		{
			u = get();
			v = get();

			if( u == 0.0 || v == 0.0 )
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

		} while( q > 0.27846 || v * v > -4.0 * log( u ) * u * u );

		return( v / u );
	}

	/**
	 * Function returns the next random bit, usually used for 50% probability
	 * evaluations efficiently.
	 */

	int getBit()
	{
		if( BitsLeft == 0 )
		{
			BitPool = advance();

			const int b = (int) ( BitPool & 1 );

			BitsLeft = 63;
			BitPool >>= 1;

			return( b );
		}

		const int b = (int) ( BitPool & 1 );

		BitsLeft--;
		BitPool >>= 1;

		return( b );
	}

protected:
	biteopt_rng rf; ///< External random number generator to use; NULL: use
		///< the default PRNG.
		///<
	void* rdata; ///< Data pointer to pass to the "rf" function.
		///<
	uint64_t Seed, lcg, Hash; ///< PRNG state variables.
		///<
	uint64_t BitPool; ///< Bit pool.
		///<
	int BitsLeft; ///< The number of bits left in the bit pool.
		///<

	/**
	 * Function advances the PRNG and returns the next PRNG value.
	 */

	uint64_t advance()
	{
		if( rf != NULL )
		{
			uint64_t r = ( *rf )( rdata );
			r |= (uint64_t) ( *rf )( rdata ) << 32;

			return( r );
		}

		Seed *= lcg * 2 + 1;
		const uint64_t rs = Seed >> 32 | Seed << 32;
		Hash += rs + 0xAAAAAAAAAAAAAAAA;
		lcg += Seed + 0x5555555555555555;
		Seed ^= Hash;

		return( lcg ^ rs );
	}
};

/**
 * Histogram class. Used to keep track of success of various choices. Updates
 * probabilities of future choices based on the selection outcome.
 *
 * Called "histogram" for historic reasons. The current implementation uses
 * bubble-sort-alike method to update a vector of possible choices. The
 * selection is made as a weighted-random draw of a value from this vector.
 */

class CBiteOptHistBase
{
public:
	/**
	 * Constructor.
	 *
	 * @param Count The number of possible choices, greater than 1.
	 */

	CBiteOptHistBase( const int aCount )
		: Count( aCount )
		, CountSp( Count * SparseMul )
		, CountSp1( CountSp - 1 )
	{
	}

	/**
	 * This function resets *this object, should be called before calling
	 * other functions, including after object's construction.
	 *
	 * @param rnd PRNG object.
	 */

	void reset( CBiteRnd& rnd )
	{
		int j;

		for( j = 0; j < SlotCount; j++ )
		{
			int i;

			for( i = 0; i < Count; i++ )
			{
				int k;

				for( k = 0; k < SparseMul; k++ )
				{
					Sels[ j ][ i * SparseMul + k ] = i;
				}
			}

			// Randomized swap-mixing.

			for( i = 0; i < CountSp * SparseMul; i++ )
			{
				const int i1 = rnd.getInt( CountSp );
				const int i2 = rnd.getInt( CountSp );

				const int t = Sels[ j ][ i1 ];
				Sels[ j ][ i1 ] = Sels[ j ][ i2 ];
				Sels[ j ][ i2 ] = t;
			}
		}

		Slot = 0;

		select( rnd );
	}

	/**
	 * An auxiliary function that returns choice count.
	 */

	int getChoiceCount() const
	{
		return( Count );
	}

	/**
	 * This function should be called when a certain choice is successful.
	 * This function should only be called after a prior select() calls.
	 *
	 * @param rnd PRNG object. May not be used.
	 * @param v Histogram increment value, [0; 1].
	 */

	void incr( CBiteRnd& rnd, const double v = 1.0 )
	{
		if( Selp > 0 && rnd.get() < v * v ) // Boost choice with a good cost.
		{
			Sels[ Slot ][ Selp ] = Sels[ Slot ][ Selp - 1 ];
			Sels[ Slot ][ Selp - 1 ] = Sel;
		}

		if( Selp > CountSp / 2 && Slot + 1 < SlotCount )
		{
			Slot++;
		}
	}

	/**
	 * This function should be called when a certain choice is a failure.
	 * This function should only be called after a prior select() calls.
	 *
	 * @param rnd PRNG object. May not be used.
	 */

	void decr( CBiteRnd& rnd )
	{
		if( Selp < CountSp1 ) // Demote inefficient choice.
		{
			Sels[ Slot ][ Selp ] = Sels[ Slot ][ Selp + 1 ];
			Sels[ Slot ][ Selp + 1 ] = Sel;
		}

		if( Selp < CountSp / 2 && Slot > 0 )
		{
			Slot--;
		}
	}

	/**
	 * Function produces a random choice index based on the current *this
	 * object's state. Note that "select" functions can only be called once
	 * for a given *this object during the optimize() function call.
	 *
	 * @param rnd PRNG object.
	 */

	int select( CBiteRnd& rnd )
	{
		const double r = rnd.get();
		Selp = (int) ( r * sqrt( r ) * CountSp );
		Sel = Sels[ Slot ][ Selp ];

		return( Sel );
	}

	/**
	 * Function returns the latest made choice index.
	 */

	int getSel() const
	{
		return( Sel );
	}

protected:
	static const int MaxCount = 4; ///< The maximal number of choices
		///< supported.
		///<
	static const int SparseMul = 5; ///< Multiplier used to obtain an actual
		///< length of the choice vector.
		///<
	static const int SlotCount = 4; ///< The number of choice vectors in use.
		///<
	int Count; ///< The number of choices in use.
		///<
	int CountSp; ///< = Count * SparseMul. The actual length of the choice
		///< vector.
		///<
	int CountSp1; ///< = CountSp - 1.
		///<
	int Sels[ SlotCount ][ MaxCount * SparseMul ]; ///< Choice vectors.
		///<
	int Sel; ///< The latest selected choice. Available only after the
		///< select() function calls.
		///<
	int Selp; ///< The index of the choice in the Sels vector.
		///<
	int Slot; ///< The current Sels vector, depending on incr/decr.
		///<
};

/**
 * Templated histogram class, for convenient constructor's Count parameter
 * specification.
 *
 * @tparam tCount The number of possible choices, greater than 1.
 */

template< int tCount >
class CBiteOptHist : public CBiteOptHistBase
{
public:
	CBiteOptHist()
		: CBiteOptHistBase( tCount )
	{
	}
};

/**
 * Class implements storage of population parameter vectors, costs, centroid,
 * and ordering.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBiteOptPop
{
public:
	CBiteOptPop()
		: ParamCount( 0 )
		, PopSize( 0 )
		, PopParamsBuf( NULL )
		, PopParams( NULL )
		, PopCosts( NULL )
		, CentParams( NULL )
		, NeedCentUpdate( false )
	{
	}

	CBiteOptPop( const CBiteOptPop& s )
		: PopParamsBuf( NULL )
		, PopParams( NULL )
		, PopCosts( NULL )
		, CentParams( NULL )
	{
		initBuffers( s.ParamCount, s.PopSize );
		copy( s );
	}

	virtual ~CBiteOptPop()
	{
		deleteBuffers();
	}

	CBiteOptPop& operator = ( const CBiteOptPop& s )
	{
		copy( s );
		return( *this );
	}

	/**
	 * Function initializes all common buffers, and "PopSize" variables. This
	 * function should be called when population's dimensions were changed.
	 * This function calls the deleteBuffers() function to release any
	 * derived classes' allocated buffers. Allocates an additional vector for
	 * temporary use, which is at the same the last vector in the PopParams
	 * array. Derived classes should call this function of the base class.
	 *
	 * @param aParamCount New parameter count.
	 * @param aPopSize New population size. If <= 0, population buffers will
	 * not be allocated.
	 */

	virtual void initBuffers( const int aParamCount, const int aPopSize )
	{
		deleteBuffers();

		ParamCount = aParamCount;
		ParamCountI = 1.0 / ParamCount;
		PopSize = aPopSize;
		CurPopSize = aPopSize;
		CurPopSizeI = 1.0 / CurPopSize;
		CurPopSize1 = aPopSize - 1;

		PopParamsBuf = new ptype[( aPopSize + 1 ) * aParamCount ];
		PopParams = new ptype*[ aPopSize + 1 ]; // Last element is temporary.
		PopCosts = new double[ aPopSize ];
		CentParams = new ptype[ aParamCount ];

		int i;

		for( i = 0; i <= aPopSize; i++ )
		{
			PopParams[ i ] = PopParamsBuf + i * aParamCount;
		}

		TmpParams = PopParams[ aPopSize ];
	}

	/**
	 * Function copies population from the specified source population. If
	 * *this population has a different size, or is uninitialized, it will
	 * be initialized to source's population size.
	 *
	 * @param s Source population to copy. Should be initalized.
	 */

	void copy( const CBiteOptPop& s )
	{
		if( ParamCount != s.ParamCount || PopSize != s.PopSize )
		{
			initBuffers( s.ParamCount, s.PopSize );
		}

		CurPopSize = s.CurPopSize;
		CurPopSizeI = s.CurPopSizeI;
		CurPopSize1 = s.CurPopSize1;
		CurPopPos = s.CurPopPos;
		NeedCentUpdate = s.NeedCentUpdate;

		int i;

		for( i = 0; i < CurPopSize; i++ )
		{
			memcpy( PopParams[ i ], s.PopParams[ i ],
				ParamCount * sizeof( PopParams[ i ][ 0 ]));
		}

		memcpy( PopCosts, s.PopCosts, CurPopSize * sizeof( PopCosts[ 0 ]));

		if( !NeedCentUpdate )
		{
			memcpy( CentParams, s.CentParams, ParamCount *
				sizeof( CentParams[ 0 ]));
		}
	}

	/**
	 * Function recalculates centroid based on the current population size.
	 * The NeedCentUpdate variable can be checked if centroid update is
	 * needed. This function resets the NeedCentUpdate to "false".
	 */

	void updateCentroid()
	{
		NeedCentUpdate = false;

		const int BatchCount = ( 1 << IntOverBits ) - 1;
		const double m = CurPopSizeI;
		ptype* const cp = CentParams;
		int i;
		int j;

		if( CurPopSize <= BatchCount )
		{
			memcpy( cp, PopParams[ 0 ], ParamCount * sizeof( cp[ 0 ]));

			for( j = 1; j < CurPopSize; j++ )
			{
				const ptype* const p = PopParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					cp[ i ] += p[ i ];
				}
			}

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] = (ptype) ( cp[ i ] * m );
			}
		}
		else
		{
			// Batched centroid calculation, for more precision and no integer
			// overflows.

			ptype* const tp = TmpParams;
			int pl = CurPopSize;
			j = 0;
			bool DoCopy = true;

			while( pl > 0 )
			{
				int c = ( pl > BatchCount ? BatchCount : pl );
				pl -= c;
				c--;

				memcpy( tp, PopParams[ j ], ParamCount * sizeof( tp[ 0 ]));

				while( c > 0 )
				{
					j++;
					const ptype* const p = PopParams[ j ];

					for( i = 0; i < ParamCount; i++ )
					{
						tp[ i ] += p[ i ];
					}

					c--;
				}

				if( DoCopy )
				{
					DoCopy = false;

					for( i = 0; i < ParamCount; i++ )
					{
						cp[ i ] = (ptype) ( tp[ i ] * m );
					}
				}
				else
				{
					for( i = 0; i < ParamCount; i++ )
					{
						cp[ i ] += (ptype) ( tp[ i ] * m );
					}
				}
			}
		}
	}

	/**
	 * Function returns pointer to the centroid vector. The NeedUpdateCent
	 * should be checked and and if it is equal to "true", the
	 * updateCentroid() function called.
	 */

	const ptype* getCentroid() const
	{
		return( CentParams );
	}

	/**
	 * Function returns population's parameter vector by the specified index
	 * from the ordered list.
	 *
	 * @param i Parameter vector index.
	 */

	const ptype* getParamsOrdered( const int i ) const
	{
		return( PopParams[ i ]);
	}

	/**
	 * Function returns a pointer to array of population vector pointers,
	 * which are sorted in the ascending cost order.
	 */

	const ptype** getPopParams() const
	{
		return( (const ptype**) PopParams );
	}

	/**
	 * Function returns current population size.
	 */

	int getCurPopSize() const
	{
		return( CurPopSize );
	}

	/**
	 * Function returns current population position.
	 */

	int getCurPopPos() const
	{
		return( CurPopPos );
	}

	/**
	 * Function resets the current population position to zero. This function
	 * is usually called when the population needs to be completely changed.
	 * This function should be called before any updates to *this population
	 * (usually during optimizer's initialization).
	 */

	void resetCurPopPos()
	{
		CurPopPos = 0;
		NeedCentUpdate = false;
	}

	/**
	 * Function returns "true" if the specified cost meets population's
	 * cost constraint. The check is synchronized with the sortPop() function.
	 *
	 * @param Cost Cost value to evaluate.
	 */

	bool isAcceptedCost( const double Cost ) const
	{
		return( Cost <= PopCosts[ CurPopSize1 ]);
	}

	/**
	 * Function replaces the highest-cost previous solution, updates centroid.
	 * This function considers the value of the CurPopPos variable - if it is
	 * smaller than the CurPopSize, the new solution will be added to
	 * population without any checks.
	 *
	 * @param NewCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 * @param DoUpdateCentroid "True" if centroid should be updated using
	 * running sum. This update is done for parallel populations.
	 * @param DoCostCheck "True" if the cost contraint should be checked.
	 * Function returns "CurPopSize" if the cost constraint was not met;
	 * insertion position otherwise.
	 */

	int updatePop( const double NewCost, const ptype* const UpdParams,
		const bool DoUpdateCentroid, const bool DoCostCheck )
	{
		if( CurPopPos < CurPopSize )
		{
			memcpy( PopParams[ CurPopPos ], UpdParams,
				ParamCount * sizeof( PopParams[ CurPopPos ][ 0 ]));

			const int p = sortPop( NewCost, CurPopPos );
			CurPopPos++;

			return( p );
		}

		if( DoCostCheck )
		{
			if( !isAcceptedCost( NewCost ))
			{
				return( CurPopSize );
			}
		}

		ptype* const rp = PopParams[ CurPopSize1 ];

		if( DoUpdateCentroid )
		{
			ptype* const cp = CentParams;
			const double m = CurPopSizeI;
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += (ptype) (( UpdParams[ i ] - rp[ i ]) * m );
				rp[ i ] = UpdParams[ i ];
			}
		}
		else
		{
			memcpy( rp, UpdParams, ParamCount * sizeof( rp[ 0 ]));
			NeedCentUpdate = true;
		}

		return( sortPop( NewCost, CurPopSize1 ));
	}

	/**
	 * Function increases current population size, and updates the required
	 * variables. This function can only be called if CurPopSize is less than
	 * PopSize.
	 *
	 * @param CopyVec If >=, a parameter vector with this index will be copied
	 * to the newly-added vector. The maximal population cost will be copied
	 * as well.
	 */

	void incrCurPopSize( const int CopyVec = -1 )
	{
		if( CopyVec >= 0 )
		{
			PopCosts[ CurPopSize ] = PopCosts[ CurPopSize1 ];
			memcpy( PopParams[ CurPopSize ], PopParams[ CopyVec ],
				ParamCount * sizeof( PopParams[ CurPopSize ][ 0 ]));
		}

		CurPopSize++;
		CurPopSizeI = 1.0 / CurPopSize;
		CurPopSize1++;
		NeedCentUpdate = true;
	}

	/**
	 * Function decreases current population size, and updates the required
	 * variables.
	 */

	void decrCurPopSize()
	{
		CurPopSize--;
		CurPopSizeI = 1.0 / CurPopSize;
		CurPopSize1--;
		NeedCentUpdate = true;
	}

protected:
	static const int IntOverBits = ( sizeof( ptype ) > 4 ? 5 : 3 ); ///< The
		///< number of bits of precision required for centroid calculation and
		///< overflows.
		///<
	static const int IntMantBits = sizeof( ptype ) * 8 - 1 - IntOverBits; ///<
		///< Mantissa size of the integer parameter values (higher by 1 bit in
		///< practice for real value 1.0). Accounts for a sign bit, and
		///< possible accumulation overflows.
		///<
	static const int64_t IntMantMult = 1LL << IntMantBits; ///< Mantissa
		///< multiplier.
		///<
	static const int64_t IntMantMultM = -IntMantMult; ///< Negative
		///< IntMantMult.
		///<
	static const int64_t IntMantMult2 = ( IntMantMult << 1 ); ///< =
		///< IntMantMult * 2.
		///<
	static const int64_t IntMantMask = IntMantMult - 1; ///< Mask that
		///< corresponds to mantissa.
		///<

	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	double ParamCountI; ///< = 1.0 / ParamCount.
		///<
	int PopSize; ///< The size of population in use (maximal).
		///<
	int CurPopSize; ///< Current population size.
		///<
	int CurPopSize1; ///< = CurPopSize - 1.
		///<
	double CurPopSizeI; ///< = 1.0 / CurPopSize.
		///<
	int CurPopPos; ///< Current population position, for initial population
		///< update. This variable should be initialized by the optimizer.
		///<
	ptype* PopParamsBuf; ///< Buffer for all PopParams vectors.
		///<
	ptype** PopParams; ///< Population parameter vectors. Always kept sorted
		///< in ascending cost order.
		///<
	double* PopCosts; ///< Costs of population parameter vectors, sorting
		///< order corresponds to PopParams.
		///<
	ptype* CentParams; ///< Centroid of the current parameter vectors.
		///<
	bool NeedCentUpdate; ///< "True" if centroid update is needed.
		///<
	ptype* TmpParams; ///< Temporary parameter vector, points to the last
		///< element of the PopParams array.
		///<

	/**
	 * Function deletes buffers previously allocated via the initBuffers()
	 * function. Derived classes should call this function of the base class.
	 */

	virtual void deleteBuffers()
	{
		delete[] PopParamsBuf;
		delete[] PopParams;
		delete[] PopCosts;
		delete[] CentParams;
	}

	/**
	 * Function performs re-sorting of the population based on the cost of a
	 * newly-added solution, and stores new cost.
	 *
	 * @param Cost Solution's cost.
	 * @param i Solution's index (usually, CurPopSize1).
	 * @return Ordered insertion index.
	 */

	int sortPop( const double Cost, int i )
	{
		ptype* const InsertParams = PopParams[ i ];

		while( i > 0 )
		{
			const double c1 = PopCosts[ i - 1 ];

			if( c1 < Cost )
			{
				break;
			}

			PopCosts[ i ] = c1;
			PopParams[ i ] = PopParams[ i - 1 ];
			i--;
		}

		PopCosts[ i ] = Cost;
		PopParams[ i ] = InsertParams;

		return( i );
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0] range (including in integer range), by wrapping it over the
	 * boundaries using random operator. This operation improves convergence
	 * in comparison to clamping.
	 *
	 * @param v Parameter value to wrap.
	 * @return Wrapped parameter value.
	 */

	static ptype wrapParam( CBiteRnd& rnd, const ptype v )
	{
		if( (ptype) 0.25 == 0 )
		{
			if( v < 0 )
			{
				if( v > IntMantMultM )
				{
					return( (ptype) ( rnd.get() * -v ));
				}

				return( (ptype) ( rnd.getRaw() & IntMantMask ));
			}

			if( v > IntMantMult )
			{
				if( v < IntMantMult2 )
				{
					return( (ptype) ( IntMantMult -
						rnd.get() * ( v - IntMantMult )));
				}

				return( (ptype) ( rnd.getRaw() & IntMantMask ));
			}

			return( v );
		}
		else
		{
			if( v < 0.0 )
			{
				if( v > -1.0 )
				{
					return( (ptype) ( rnd.get() * -v ));
				}

				return( (ptype) rnd.get() );
			}

			if( v > 1.0 )
			{
				if( v < 2.0 )
				{
					return( (ptype) ( 1.0 - rnd.get() * ( v - 1.0 )));
				}

				return( (ptype) rnd.get() );
			}

			return( v );
		}
	}

	/**
	 * Function generates a Gaussian-distributed pseudo-random number, in
	 * integer scale, with the specified mean and std.dev.
	 *
	 * @param rnd Uniform PRNG.
	 * @param sd Standard deviation multiplier.
	 * @param meanInt Mean value, in integer scale.
	 */

	static ptype getGaussianInt( CBiteRnd& rnd, const double sd,
		const ptype meanInt )
	{
		while( true )
		{
			const double r = rnd.getGaussian() * sd;

			if( r > -8.0 && r < 8.0 )
			{
				return( (ptype) ( r * IntMantMult + meanInt ));
			}
		}
	}
};

/**
 * Population class that embeds a dynamically-allocated parallel population
 * objects.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBiteOptParPops : virtual public CBiteOptPop< ptype >
{
public:
	CBiteOptParPops()
		: ParPopCount( 0 )
	{
	}

	virtual ~CBiteOptParPops()
	{
		int i;

		for( i = 0; i < ParPopCount; i++ )
		{
			delete ParPops[ i ];
		}
	}

protected:
	using CBiteOptPop< ptype > :: ParamCount;

	static const int MaxParPopCount = 8; ///< The maximal number of parallel
		///< population supported.
		///<
	CBiteOptPop< ptype >* ParPops[ MaxParPopCount ]; ///< Parallel
		///< population orbiting *this population.
		///<
	int ParPopCount; ///< Parallel population count. This variable should only
		///< be changed via the setParPopCount() function.
		///<

	/**
	 * Function changes the parallel population count, and reallocates
	 * buffers.
	 *
	 * @param NewCount New parallel population count, >= 0.
	 */

	void setParPopCount( const int NewCount )
	{
		while( ParPopCount > NewCount )
		{
			ParPopCount--;
			delete ParPops[ ParPopCount ];
		}

		while( ParPopCount < NewCount )
		{
			ParPops[ ParPopCount ] = new CBiteOptPop< ptype >();
			ParPopCount++;
		}
	}

	/**
	 * Function returns index of the parallel population that is most close
	 * to the specified parameter vector. Function returns -1 if the cost
	 * constraint is not met in all parallel populations.
	 *
	 * @param Cost Cost of parameter vector, used to filter considered
	 * parallel population pool.
	 * @param p Parameter vector.
	 */

	int getMinDistParPop( const double Cost, const ptype* const p ) const
	{
		int ppi[ MaxParPopCount ];
		int ppc = 0;
		int i;

		for( i = 0; i < ParPopCount; i++ )
		{
			if( ParPops[ i ] -> isAcceptedCost( Cost ))
			{
				ppi[ ppc ] = i;
				ppc++;
			}
		}

		if( ppc == 0 )
		{
			return( -1 );
		}

		if( ppc == 1 )
		{
			return( ppi[ 0 ]);
		}

		double s[ MaxParPopCount ];

		if( ppc == 4 )
		{
			const ptype* const c0 = ParPops[ ppi[ 0 ]] -> getCentroid();
			const ptype* const c1 = ParPops[ ppi[ 1 ]] -> getCentroid();
			const ptype* const c2 = ParPops[ ppi[ 2 ]] -> getCentroid();
			const ptype* const c3 = ParPops[ ppi[ 3 ]] -> getCentroid();
			double s0 = 0.0;
			double s1 = 0.0;
			double s2 = 0.0;
			double s3 = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const ptype v = p[ i ];
				const double d0 = (double) ( v - c0[ i ]);
				const double d1 = (double) ( v - c1[ i ]);
				const double d2 = (double) ( v - c2[ i ]);
				const double d3 = (double) ( v - c3[ i ]);
				s0 += d0 * d0;
				s1 += d1 * d1;
				s2 += d2 * d2;
				s3 += d3 * d3;
			}

			s[ 0 ] = s0;
			s[ 1 ] = s1;
			s[ 2 ] = s2;
			s[ 3 ] = s3;
		}
		else
		if( ppc == 3 )
		{
			const ptype* const c0 = ParPops[ ppi[ 0 ]] -> getCentroid();
			const ptype* const c1 = ParPops[ ppi[ 1 ]] -> getCentroid();
			const ptype* const c2 = ParPops[ ppi[ 2 ]] -> getCentroid();
			double s0 = 0.0;
			double s1 = 0.0;
			double s2 = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const ptype v = p[ i ];
				const double d0 = (double) ( v - c0[ i ]);
				const double d1 = (double) ( v - c1[ i ]);
				const double d2 = (double) ( v - c2[ i ]);
				s0 += d0 * d0;
				s1 += d1 * d1;
				s2 += d2 * d2;
			}

			s[ 0 ] = s0;
			s[ 1 ] = s1;
			s[ 2 ] = s2;
		}
		else
		if( ppc == 2 )
		{
			const ptype* const c0 = ParPops[ ppi[ 0 ]] -> getCentroid();
			const ptype* const c1 = ParPops[ ppi[ 1 ]] -> getCentroid();
			double s0 = 0.0;
			double s1 = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const ptype v = p[ i ];
				const double d0 = (double) ( v - c0[ i ]);
				const double d1 = (double) ( v - c1[ i ]);
				s0 += d0 * d0;
				s1 += d1 * d1;
			}

			s[ 0 ] = s0;
			s[ 1 ] = s1;
		}

		int pp = 0;
		double d = s[ pp ];

		for( i = 1; i < ppc; i++ )
		{
			if( s[ i ] <= d )
			{
				pp = i;
				d = s[ i ];
			}
		}

		return( ppi[ pp ]);
	}
};

/**
 * Base virtual abstract class that defines common optimizer interfacing
 * functions.
 */

class CBiteOptInterface
{
public:
	CBiteOptInterface()
	{
	}

	virtual ~CBiteOptInterface()
	{
	}

	/**
	 * @return Best parameter vector.
	 */

	virtual const double* getBestParams() const = 0;

	/**
	 * @return Cost of the best parameter vector.
	 */

	virtual double getBestCost() const = 0;

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
	 * Virtual function (objective function) that should calculate parameter
	 * vector's optimization cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) = 0;
};

/**
 * The base class for optimizers of the "biteopt" project.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBiteOptBase : public CBiteOptInterface,
	virtual protected CBiteOptParPops< ptype >
{
private:
	CBiteOptBase( const CBiteOptBase& )
	{
		// Copy-construction unsupported.
	}

	CBiteOptBase& operator = ( const CBiteOptBase& )
	{
		// Copying unsupported.
		return( *this );
	}

public:
	CBiteOptBase()
		: MinValues( NULL )
		, MaxValues( NULL )
		, DiffValues( NULL )
		, DiffValuesI( NULL )
		, BestValues( NULL )
		, NewValues( NULL )
		, HistCount( 0 )
	{
	}

	virtual const double* getBestParams() const
	{
		return( BestValues );
	}

	virtual double getBestCost() const
	{
		return( BestCost );
	}

	static const int MaxHistCount = 64; ///< The maximal number of histograms
		///< that can be added (for static arrays).
		///<

	/**
	 * Function returns a pointer to an array of histograms in use.
	 */

	CBiteOptHistBase** getHists()
	{
		return( Hists );
	}

	/**
	 * Function returns a pointer to an array of histogram names.
	 */

	const char** getHistNames() const
	{
		return( (const char**) HistNames );
	}

	/**
	 * Function returns the number of histograms in use.
	 */

	int getHistCount() const
	{
		return( HistCount );
	}

protected:
	using CBiteOptParPops< ptype > :: IntMantMult;
	using CBiteOptParPops< ptype > :: ParamCount;
	using CBiteOptParPops< ptype > :: PopSize;
	using CBiteOptParPops< ptype > :: CurPopSize;
	using CBiteOptParPops< ptype > :: CurPopSizeI;
	using CBiteOptParPops< ptype > :: CurPopSize1;
	using CBiteOptParPops< ptype > :: CurPopPos;
	using CBiteOptParPops< ptype > :: NeedCentUpdate;
	using CBiteOptParPops< ptype > :: resetCurPopPos;

	double* MinValues; ///< Minimal parameter values.
		///<
	double* MaxValues; ///< Maximal parameter values.
		///<
	double* DiffValues; ///< Difference between maximal and minimal parameter
		///< values.
		///<
	double* DiffValuesI; ///< Inverse DiffValues.
		///<
	double* BestValues; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	double* NewValues; ///< Temporary new parameter buffer, with real values.
		///<
	int StallCount; ///< The number of iterations without improvement.
		///<
	double HiBound; ///< Higher cost bound, for StallCount estimation. May not
		///< be used by the optimizer.
		///<
	double AvgCost; ///< Average cost in the latest batch. May not be used by
		///< the optimizer.
		///<
	CBiteOptHistBase* Hists[ MaxHistCount ]; ///< Pointers to histogram
		///< objects, for indexed access in some cases.
		///<
	const char* HistNames[ MaxHistCount ]; ///< Histogram names.
		///<
	int HistCount; ///< The number of histograms in use.
		///<
	static const int MaxApplyHists = 32; /// The maximal number of histograms
		///< that can be used during the optimize() function call.
		///<
	CBiteOptHistBase* ApplyHists[ MaxApplyHists ]; ///< Histograms used in
		///< "selects" during the optimize() function call.
		///<
	int ApplyHistsCount; ///< The number of "selects" used during the
		///< optimize() function call.
		///<

	virtual void initBuffers( const int aParamCount, const int aPopSize )
	{
		CBiteOptParPops< ptype > :: initBuffers( aParamCount, aPopSize );

		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		DiffValues = new double[ ParamCount ];
		DiffValuesI = new double[ ParamCount ];
		BestValues = new double[ ParamCount ];
		NewValues = new double[ ParamCount ];
	}

	virtual void deleteBuffers()
	{
		CBiteOptParPops< ptype > :: deleteBuffers();

		delete[] MinValues;
		delete[] MaxValues;
		delete[] DiffValues;
		delete[] DiffValuesI;
		delete[] BestValues;
		delete[] NewValues;
	}

	/**
	 * Function resets common variables used by optimizers to their default
	 * values, including registered histograms, calls the resetCurPopPos()
	 * and updateDiffValues() functions. This function is usually called in
	 * the init() function of the optimizer.
	 */

	void resetCommonVars( CBiteRnd& rnd )
	{
		updateDiffValues();
		resetCurPopPos();

		CurPopSize = PopSize;
		CurPopSizeI = 1.0 / PopSize;
		CurPopSize1 = PopSize - 1;
		BestCost = 1e300;
		StallCount = 0;
		HiBound = 1e300;
		AvgCost = 0.0;
		ApplyHistsCount = 0;

		int i;

		for( i = 0; i < HistCount; i++ )
		{
			Hists[ i ] -> reset( rnd );
		}
	}

	/**
	 * Function updates values in the DiffValues array, based on values in the
	 * MinValues and MaxValues arrays.
	 */

	void updateDiffValues()
	{
		int i;

		if( (ptype) 0.25 == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double d = MaxValues[ i ] - MinValues[ i ];

				DiffValues[ i ] = d / IntMantMult;
				DiffValuesI[ i ] = IntMantMult / d;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
				DiffValuesI[ i ] = 1.0 / DiffValues[ i ];
			}
		}
	}

	/**
	 * Function updates BestCost value and BestValues array, if the specified
	 * NewCost is better.
	 *
	 * @param NewCost New solution's cost.
	 * @param UpdValues New solution's values. The values should be in the
	 * "real" value range.
	 */

	void updateBestCost( const double NewCost, const double* const UpdValues )
	{
		if( NewCost <= BestCost )
		{
			BestCost = NewCost;

			memcpy( BestValues, UpdValues,
				ParamCount * sizeof( BestValues[ 0 ]));
		}
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param NormParams Parameter vector of interest, in normalized scale.
	 * @param i Parameter index.
	 */

	double getRealValue( const ptype* const NormParams, const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] * NormParams[ i ]);
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [MinValue; MaxValue] real range, by wrapping it over the boundaries
	 * using random operator. This operation improves convergence in
	 * comparison to clamping.
	 *
	 * @param v Parameter value to wrap.
	 * @param i Parameter index.
	 * @return Wrapped parameter value.
	 */

	double wrapParamReal( CBiteRnd& rnd, const double v, const int i ) const
	{
		const double minv = MinValues[ i ];

		if( v < minv )
		{
			const double dv = DiffValues[ i ];

			if( v > minv - dv )
			{
				return( minv + rnd.get() * ( minv - v ));
			}

			return( minv + rnd.get() * dv );
		}

		const double maxv = MaxValues[ i ];

		if( v > maxv )
		{
			const double dv = DiffValues[ i ];

			if( v < maxv + dv )
			{
				return( maxv - rnd.get() * ( v - dv ));
			}

			return( maxv - rnd.get() * dv );
		}

		return( v );
	}

	/**
	 * Function adds a histogram to the Hists list.
	 *
	 * @param h Histogram object to add.
	 * @param hname Histogram's name, should be a static constant.
	 */

	void addHist( CBiteOptHistBase& h, const char* const hname )
	{
		Hists[ HistCount ] = &h;
		HistNames[ HistCount ] = hname;
		HistCount++;
	}

	/**
	 * Function performs choice selection based on the specified histogram,
	 * and adds the histogram to apply list.
	 *
	 * @param Hist Histogram.
	 * @param rnd PRNG object.
	 */

	template< class T >
	int select( T& Hist, CBiteRnd& rnd )
	{
		ApplyHists[ ApplyHistsCount ] = &Hist;
		ApplyHistsCount++;

		return( Hist.select( rnd ));
	}

	/**
	 * Function applies histogram increments on optimization success.
	 *
	 * @param rnd PRNG object.
	 * @param v Increment value, [0; 1].
	 */

	void applyHistsIncr( CBiteRnd& rnd, const double v = 1.0 )
	{
		const int c = ApplyHistsCount;
		ApplyHistsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplyHists[ i ] -> incr( rnd, v );
		}
	}

	/**
	 * Function applies histogram decrements on optimization fail.
	 *
	 * @param rnd PRNG object.
	 */

	void applyHistsDecr( CBiteRnd& rnd )
	{
		const int c = ApplyHistsCount;
		ApplyHistsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplyHists[ i ] -> decr( rnd );
		}
	}
};

#endif // BITEAUX_INCLUDED
