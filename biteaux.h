//$ nocpp

/**
 * @file biteaux.h
 *
 * @brief The inclusion file for the CBiteRnd, CBiteOptPop, CBiteOptInterface,
 * and CBiteOptBase classes.
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
 * @version 2021.3
 */

#ifndef BITEAUX_INCLUDED
#define BITEAUX_INCLUDED

#include <stdint.h>
#include <math.h>
#include <string.h>

/**
 * Class that implements random number generation. The default implementation
 * includes a relatively fast pseudo-random number generator (PRNG) using a
 * classic formula "seed = ( a * seed + c ) % m" (LCG). This implementation
 * uses bits 32-61 (30 bits) of the state variable which ensures at least 2^32
 * period in the lowest significant bit of the resulting pseudo-random
 * sequence. See https://en.wikipedia.org/wiki/Linear_congruential_generator
 * for more details.
 */

class CBiteRnd
{
public:
	/**
	 * Function initializes *this PRNG object.
	 *
	 * @param NewSeed New random seed value.
	 */

	void init( const int NewSeed )
	{
		seed = NewSeed;

		// Skip first values to make PRNG "settle down".

		seed = 500009 * seed + 300119;
		seed = 500009 * seed + 300119;
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	double getRndValue()
	{
		seed = 500009 * seed + 300119;

		return((( seed >> 32 ) & 0x3FFFFFFF ) * 9.31322574615478516e-10 );
	}

	/**
	 * @return Scale of "raw" random values returned by functions with the
	 * "raw" suffix.
	 */

	static int32_t getRawScale()
	{
		return( 0x40000000 );
	}

	/**
	 * @return Uniformly-distributed random number in the range [0; 1), in the
	 * "raw" scale.
	 */

	int32_t getUniformRaw()
	{
		seed = 500009 * seed + 300119;

		return( (int32_t) (( seed >> 32 ) & 0x3FFFFFFF ));
	}

	/**
	 * @return TPDF random number in the range (-1; 1), in the "raw" scale.
	 */

	int32_t getTPDFRaw()
	{
		seed = 500009 * seed + 300119;
		const int32_t v1 = (int32_t) (( seed >> 32 ) & 0x3FFFFFFF );
		seed = 500009 * seed + 300119;
		const int32_t v2 = (int32_t) (( seed >> 32 ) & 0x3FFFFFFF );

		return( v1 - v2 );
	}

private:
	uint64_t seed; ///< The current random seed value.
		///<
};

/**
 * Class implements storage of population parameter vectors, costs, centroid,
 * and ordering.
 */

class CBiteOptPop
{
public:
	CBiteOptPop()
		: ParamCount( 0 )
		, PopSize( 0 )
		, PopOrder( NULL )
		, CurParamsBuf( NULL )
		, CurParams( NULL )
		, CurCosts( NULL )
		, CentParams( NULL )
	{
	}

	~CBiteOptPop()
	{
		deletePopBuffers();
	}

	/**
	 * Function initializes population storage buffers, include 1 vector for
	 * temporary use. This function can only be called after the construction
	 * of *this object, or after the deletePopBuffers() function call.
	 *
	 * @param aParamCount New parameter count.
	 * @param aPopSize New population size.
	 */

	void initPopBuffers( const int aParamCount, const int aPopSize )
	{
		deletePopBuffers();

		ParamCount = aParamCount;
		PopSize = aPopSize;
		PopSize1 = aPopSize - 1;
		PopSizeI = 1.0 / aPopSize;

		PopOrder = new int[ PopSize ];
		CurParamsBuf = new double[( PopSize + 1 ) * ParamCount ];
		CurParams = new double*[ PopSize + 1 ]; // Last element is temporary.
		CurCosts = new double[ PopSize ];
		CentParams = new double[ ParamCount ];

		int i;

		for( i = 0; i <= PopSize; i++ )
		{
			CurParams[ i ] = CurParamsBuf + i * ParamCount;
		}
	}

	/**
	 * Function deletes buffers previously allocated via the initPopBuffers()
	 * function.
	 */

	void deletePopBuffers()
	{
		delete[] PopOrder;
		delete[] CurParamsBuf;
		delete[] CurParams;
		delete[] CurCosts;
		delete[] CentParams;
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
			initPopBuffers( s.ParamCount, s.PopSize );
		}

		memcpy( PopOrder, s.PopOrder, PopSize * sizeof( s.PopOrder[ 0 ]));
		memcpy( CurParamsBuf, s.CurParamsBuf, PopSize * ParamCount *
			sizeof( s.CurParamsBuf[ 0 ]));

		memcpy( CurCosts, s.CurCosts, PopSize * sizeof( s.CurCosts[ 0 ]));
		memcpy( CentParams, s.CentParams, ParamCount *
			sizeof( s.CentParams[ 0 ]));
	}

	/**
	 * Function resets centroid vector.
	 */

	void resetCentroid()
	{
		memset( CentParams, 0, ParamCount * sizeof( CentParams[ 0 ]));
	}

	/**
	 * Function returns pointer to the centroid vector.
	 */

	const double* getCentroid() const
	{
		return( CentParams );
	}

	/**
	 * Function returns population's parameter vector by the specified index
	 * from the ordered list.
	 *
	 * @param i Parameter vector index.
	 */

	const double* getParamsOrdered( const int i ) const
	{
		return( CurParams[ PopOrder[ i ]]);
	}

	/**
	 * Function calculates Euclidean distance of the specifed vector to *this
	 * population's centroid. Function returns the square of the distance.
	 *
	 * @param p Parameter vector.
	 */

	double getDistanceSqr( const double* const p ) const
	{
		double s = 0.0;
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const double d = CentParams[ i ] - p[ i ];
			s += d * d;
		}

		return( s );
	}

	/**
	 * Function replaces the highest-cost previous solution (sH), updates
	 * centroid.
	 *
	 * @param NewCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 * @param sH Index of vector to update. If equal to -1, the NewCost value
	 * will be first compared to worst cost solution present in the
	 * population.
	 * @return If sH is not specified or negative the function returns "false"
	 * if cost constraint is not met. "True" otherwise.
	 */

	bool updatePop( const double NewCost, const double* const UpdParams,
		int sH = -1 )
	{
		if( sH < 0 )
		{
			sH = PopOrder[ PopSize1 ];

			if( NewCost >= CurCosts[ sH ])
			{
				return( false );
			}
		}

		double* const rp = CurParams[ sH ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CentParams[ i ] += ( UpdParams[ i ] - rp[ i ]) * PopSizeI;
			rp[ i ] = UpdParams[ i ];
		}

		insertPopOrder( NewCost, sH, PopSize1 );

		return( true );
	}

protected:
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int PopSize; ///< The size of population in use.
		///<
	int PopSize1; ///< = PopSize - 1.
		///<
	double PopSizeI; ///< = 1/PopSize.
		///<
	int* PopOrder; ///< The current solution vectors ordering,
		///< ascending-sorted by cost.
		///<
	double* CurParamsBuf; ///< CurParams buffer.
		///<
	double** CurParams; ///< Current working parameter vectors.
		///<
	double* CurCosts; ///< Best costs of current working parameter vectors.
		///<
	double* CentParams; ///< Centroid of the current parameter vectors.
		///<

	/**
	 * Function inserts the specified solution into the PopOrder
	 * array at the appropriate offset, increasing the number of items by 1.
	 *
	 * @param Cost Solution's cost.
	 * @param f Solution's index.
	 * @param ItemCount The current number of items in the array.
	 */

	void insertPopOrder( const double Cost, const int f, const int ItemCount )
	{
		CurCosts[ f ] = Cost;

		// Perform binary search.

		int z = 0;
		int hi = ItemCount;

		while( z < hi )
		{
			int mid = ( z + hi ) >> 1;

			if( CurCosts[ PopOrder[ mid ]] >= Cost )
			{
				hi = mid;
			}
			else
			{
				z = mid + 1;
			}
		}

		// Insert element at the correct sorted position.

		int i;

		for( i = ItemCount; i > z; i-- )
		{
			PopOrder[ i ] = PopOrder[ i - 1 ];
		}

		PopOrder[ z ] = f;
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
	 * @return The number of initial objective function evaluations.
	 * May correspond to the population size.
	 */

	virtual int getInitEvals() const = 0;

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
 */

class CBiteOptBase : public CBiteOptInterface, virtual protected CBiteOptPop
{
public:
	CBiteOptBase()
		: MinValues( NULL )
		, MaxValues( NULL )
		, DiffValues( NULL )
		, NewParams( NULL )
		, BestParams( NULL )
		, BestCost( 0.0 )
		, InitEvals( 0 )
		, BitsLeft( 0 )
	{
	}

	virtual ~CBiteOptBase()
	{
		deleteBuffers();
	}

	virtual int getInitEvals() const
	{
		return( InitEvals );
	}

	virtual const double* getBestParams() const
	{
		return( BestParams );
	}

	virtual double getBestCost() const
	{
		return( BestCost );
	}

protected:
	double* MinValues; ///< Minimal parameter values.
		///<
	double* MaxValues; ///< Maximal parameter values.
		///<
	double* DiffValues; ///< Difference between maximal and minimal parameter
		///< values.
		///<
	double* NewParams; ///< Temporary new parameter buffer, with real values.
		///<
	double* BestParams; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	int InitEvals; ///< Initial number of function evaluations performed by
		///< the opimzizer.
		///<
	int StallCount; ///< The number of iterations without improvement.
		///<
	int BitPool; ///< Bit pool.
		///<
	int BitsLeft; ///< The number of bits left in the bit pool. This variable
		///< should be reset to 0 on each optimizer's init() function call.
		///<
	double HiBound; ///< Higher cost bound, for StallCount estimation. May not
		///< be used by the optimizer.
		///<
	double AvgCost; ///< Average cost in the latest batch. May not be used by
		///< the optimizer.
		///<

	/**
	 * Function initializes all common buffers, PopSize1, and PopSizeI
	 * variables.
	 *
	 * @param aParamCount New parameter count.
	 * @param aPopSize New population size. If <= 0, population buffers will
	 * not be allocated.
	 */

	void initBaseBuffers( const int aParamCount, const int aPopSize )
	{
		if( aPopSize > 0 )
		{
			initPopBuffers( aParamCount, aPopSize );
		}
		else
		{
			ParamCount = aParamCount;
			PopSize = 0;
			PopSize1 = 0;
			PopSizeI = 0.0;
		}

		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		DiffValues = new double[ ParamCount ];
		NewParams = new double[ ParamCount ];
		BestParams = new double[ ParamCount ];
	}

	/**
	 * Function deletes previously allocated buffers.
	 */

	virtual void deleteBuffers()
	{
		delete[] MinValues;
		delete[] MaxValues;
		delete[] DiffValues;
		delete[] NewParams;
		delete[] BestParams;
	}

	/**
	 * Function resets common variables used by optimizers to their default
	 * values. This function is usually called in the init() function of the
	 * optimizer.
	 */

	void resetCommonVars()
	{
		BestCost = 1e300;
		StallCount = 0;
		BitsLeft = 0;
		HiBound = 1e300;
		AvgCost = 0.0;
	}

	/**
	 * Function updates values in the DiffValues array, based on values in the
	 * MinValues and MaxValues arrays.
	 */

	void updateDiffValues()
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
		}
	}

	/**
	 * Function updates BestCost and BestParams values, if the specified
	 * NewCost is better.
	 *
	 * @param NewCost New solution's cost.
	 * @param UpdParams New solution's values.
	 * @param IsNormalized "True" if values are normalized, and should be
	 * converted with the getRealValue() function.
	 */

	void updateBestCost( const double NewCost, const double* const UpdParams,
		const bool IsNormalized = false )
	{
		if( NewCost < BestCost )
		{
			BestCost = NewCost;

			if( IsNormalized )
			{
				int i;

				for( i = 0; i < ParamCount; i++ )
				{
					BestParams[ i ] = getRealValue( UpdParams, i );
				}
			}
			else
			{
				memcpy( BestParams, UpdParams,
					ParamCount * sizeof( UpdParams[ 0 ]));
			}
		}
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param NormParams Parameter vector of interest, in normalized scale.
	 * @param i Parameter index.
	 */

	double getRealValue( const double* const NormParams, const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] * NormParams[ i ]);
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0] range, by wrapping it over the boundaries using random
	 * operator. This operation improves convergence in comparison to
	 * clamping.
	 *
	 * @param v Parameter value to wrap.
	 * @return Wrapped parameter value.
	 */

	static double wrapParam( CBiteRnd& rnd, double v )
	{
		if( v < 0.0 )
		{
			if( v > -1.0 )
			{
				return( rnd.getRndValue() * -v );
			}

			return( rnd.getRndValue() );
		}

		if( v > 1.0 )
		{
			if( v < 2.0 )
			{
				return( 1.0 - rnd.getRndValue() * ( v - 1.0 ));
			}

			return( rnd.getRndValue() );
		}

		return( v );
	}

	/**
	 * Function returns the next random bit, for use in 50% probability
	 * evaluations efficiently.
	 */

	int getBit( CBiteRnd& rnd )
	{
		if( BitsLeft == 0 )
		{
			BitPool = rnd.getUniformRaw();
			BitsLeft = 29;

			const int b = ( BitPool & 1 );
			BitPool >>= 1;

			return( b );
		}
		else
		{
			const int b = ( BitPool & 1 );
			BitPool >>= 1;
			BitsLeft--;

			return( b );
		}
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

#endif // BITEAUX_INCLUDED
