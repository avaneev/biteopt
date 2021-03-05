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
 * @version 2021.9
 */

#ifndef BITEAUX_INCLUDED
#define BITEAUX_INCLUDED

#include <stdint.h>
#include <math.h>
#include <string.h>

/**
 * Class that implements a pseudo-random number generator (PRNG). The default
 * implementation includes a relatively fast PRNG that uses a classic formula
 * "seed = ( a * seed + c ) % m" (LCG). This implementation uses bits 34-63
 * (30 bits) of the state variable which ensures at least 2^34 period in the
 * lowest significant bit of the resulting pseudo-random sequence. See
 * https://en.wikipedia.org/wiki/Linear_congruential_generator for more
 * details.
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
	 * @param NewSeed New random seed value.
	 */

	void init( const int NewSeed )
	{
		BitsLeft = 0;
		seed = (uint64_t) NewSeed;

		// Skip first values to make PRNG "settle down".

		advance();
		advance();
	}

	/**
	 * @return Scale of "raw" random values returned by functions with the
	 * "raw" suffix.
	 */

	inline static int getRawScale()
	{
		return( 0x40000000 );
	}

	/**
	 * Function returns the number of significant bits in the "raw" random
	 * value.
	 */

	inline static int getRawBitCount()
	{
		return( 30 );
	}

	/**
	 * @return Inverse scale of "raw" random values returned by functions with
	 * the "raw" suffix.
	 */

	inline static double getRawScaleInv()
	{
		return( 9.31322574615478516e-10 );
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	double getRndValue()
	{
		advance();

		return(( seed >> 34 ) * getRawScaleInv() );
	}

	/**
	 * @return Uniformly-distributed random number in the range [0; 1), in the
	 * "raw" scale.
	 */

	int getUniformRaw()
	{
		advance();

		return( (int) ( seed >> 34 ));
	}

	/**
	 * @return TPDF random number in the range (-1; 1), in the "raw" scale.
	 */

	int getTPDFRaw()
	{
		advance();
		const int v1 = (int) ( seed >> 34 );

		advance();
		const int v2 = (int) ( seed >> 34 );

		return( v1 - v2 );
	}

	/**
	 * Function returns the next random bit, usually used for 50% probability
	 * evaluations efficiently.
	 */

	int getBit()
	{
		if( BitsLeft == 0 )
		{
			BitPool = getUniformRaw();
			BitsLeft = getRawBitCount() - 1;

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

private:
	uint64_t seed; ///< The current random seed value.
		///<
	int BitPool; ///< Bit pool.
		///<
	int BitsLeft; ///< The number of bits left in the bit pool.
		///<

	/**
	 * Function advances the PRNG.
	 */

	void advance()
	{
		seed = 500009 * seed + 300119;
	}
};

/**
 * Base histogram class.
 */

class CBiteOptHistBase
{
public:
	/**
	 * This function should be called when a certain choice is successful.
	 * This function should only be called after a prior select() calls.
	 */

	virtual void incr() = 0;

	/**
	 * This function should be called when a certain choice is a failure.
	 * This function should only be called after a prior select() calls.
	 */

	virtual void decr() = 0;

protected:
	int Sel; ///< The latest selected choice. Available only after the
		///< select() function calls.
		///<
};

/**
 * Histogram class. Used to keep track of success of various choices. Updates
 * probabilities of future choices based on the current histogram state.
 *
 * If the choice is used to select probability value, with a later probability
 * evalutaion, the choice should be set to a random choice if the probability
 * does not realise.
 *
 * @tparam Count The number of possible choices, greater than 1.
 * @tparam Divisor Divisor used to obtain a minimal histogram value. Controls
 * the ratio between minimal and maximal possible probabilities of all
 * choices. Should be usually equal to Count, but can be set up to Count * 2
 * if a certain choice is likely to be effective most of the time.
 * @tparam IncrDecr Histogram increment (on success) or decrement
 * (on failure). Usually equals to 1, but a higher value can be used if a
 * steep momentary probability change of a certain choice is effective.
 */

template< int Count, int Divisor, int IncrDecr >
class CBiteOptHist : virtual public CBiteOptHistBase
{
public:
	CBiteOptHist()
		: m( 1.0 / Divisor )
		, rcm( (double) Count / CBiteRnd :: getRawScale() )
	{
	}

	/**
	 * This function resets histogram, should be called before calling other
	 * functions, including after object's construction.
	 *
	 * @param rnd PRNG object. Not used.
	 */

	void reset( CBiteRnd& rnd )
	{
		memset( Hist, 0, sizeof( Hist ));
		updateProbs();
	}

	virtual void incr()
	{
		Hist[ Sel ] += IncrDecr;
		updateProbs();
	}

	virtual void decr()
	{
		Hist[ Sel ] -= IncrDecr;
		updateProbs();
	}

	/**
	 * Function produces a random choice index based on the current histogram
	 * state. Note that "select" functions can only be called once for a given
	 * histogram during the optimize() function call.
	 *
	 * @param rnd PRNG object.
	 */

	int select( CBiteRnd& rnd )
	{
		const double rv = rnd.getRndValue() * ProbSum;

		if( Count == 2 )
		{
			Sel = ( rv >= Probs[ 0 ]);
			return( Sel );
		}
		else
		{
			int i;

			for( i = 1; i < Count; i++ )
			{
				if( rv >= Probs[ i - 1 ] && rv <= Probs[ i ])
				{
					Sel = i;
					return( i );
				}
			}

			Sel = 0;
			return( 0 );
		}
	}

	/**
	 * Function makes a uniformly-distributed choice.
	 *
	 * @param rnd PRNG object.
	 */

	int selectRandom( CBiteRnd& rnd )
	{
		Sel = (int) ( rnd.getUniformRaw() * rcm );
		return( Sel );
	}

	/**
	 * Function forces a specified choice.
	 *
	 * @param s Choice index, must be in Count range.
	 */

	int selectForce( const int s )
	{
		Sel = s;
		return( s );
	}

	/**
	 * Function "unselects" a previously selected choice so that the choice on
	 * either the incr() or decr() call is randomized.
	 */

	void unselect( CBiteRnd& rnd )
	{
		if( Count == 2 )
		{
			Sel = rnd.getBit();
		}
		else
		{
			Sel = (int) ( rnd.getUniformRaw() * rcm );
		}
	}

protected:
	double m; ///< Multiplier (depends on Divisor)
		///<
	double rcm; ///< Raw random value multiplier that depends on Count.
		///<
	int Hist[ Count ]; ///< Histogram.
		///<
	double Probs[ Count ]; ///< Probabilities, cumulative.
		///<
	double ProbSum; ///< Sum of probabilities, for random variable scaling.
		///<

	/**
	 * Function updates probabilities of choices based on the histogram state.
	 */

	void updateProbs()
	{
		int MinHist = Hist[ 0 ];
		int i;

		for( i = 1; i < Count; i++ )
		{
			if( Hist[ i ] < MinHist )
			{
				MinHist = Hist[ i ];
			}
		}

		MinHist--;
		double HistSum = 0.0;

		for( i = 0; i < Count; i++ )
		{
			Probs[ i ] = Hist[ i ] - MinHist;
			HistSum += Probs[ i ];
		}

		HistSum *= m;
		ProbSum = 0.0;

		for( i = 0; i < Count; i++ )
		{
			const double v = ( Probs[ i ] < HistSum ? HistSum : Probs[ i ]) +
				ProbSum;

			Probs[ i ] = v;
			ProbSum = v;
		}
	}
};

/**
 * Histogram class for binary variables. A lot more computationally-efficient,
 * but functionally similar to the CBiteOptHistogram class in performance.
 *
 * The choice should be set to a random choice if the probability does not
 * realise.
 */

class CBiteOptHistBinary : virtual public CBiteOptHistBase
{
public:
	/**
	 * This function resets histogram to a random state, should be called
	 * before calling other functions, including after object's construction.
	 *
	 * @param rnd PRNG object.
	 */

	void reset( CBiteRnd& rnd )
	{
		const int b = rnd.getBit();

		Hist[ 0 ] = b;
		Hist[ 1 ] = 1 - b;
	}

	/**
	 * This function should be called when a certain choice is successful.
	 *
	 * @param Index Choice index. Should be equal to 0 or 1.
	 */

	virtual void incr()
	{
		Hist[ Sel ]++;
	}

	/**
	 * This function should be called when a certain choice is a failure.
	 *
	 * @param Index Choice index. Should be equal to 0 or 1.
	 */

	virtual void decr()
	{
		Hist[ Sel ]--;
	}

	/**
	 * Function produces a fixed binary choice index based on the current
	 * histogram state.
	 *
	 * @param rnd PRNG object. Not used.
	 */

	int select( CBiteRnd& rnd )
	{
		Sel = ( Hist[ 1 ] > Hist[ 0 ]);
		return( Sel );
	}

	/**
	 * Function "unselects" a previously selected choice so that the choice on
	 * either incr() or decr() call is randomized.
	 */

	void unselect( CBiteRnd& rnd )
	{
		Sel = rnd.getBit();
	}

protected:
	int Hist[ 2 ]; ///< Histogram.
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
		, CurParamsBuf( NULL )
		, CurParams( NULL )
		, CurCosts( NULL )
		, CentParams( NULL )
	{
	}

	CBiteOptPop( const CBiteOptPop& s )
	{
		initPopBuffers( s.ParamCount, s.PopSize );
		copy( s );
	}

	~CBiteOptPop()
	{
		deletePopBuffers();
	}

	/**
	 * Function initializes population storage buffers, include 1 vector for
	 * temporary use.
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
		CurPopSize = aPopSize;
		CurPopSize1 = aPopSize - 1;
		NeedCentUpdate = false;

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

		CurPopSize = s.CurPopSize;
		CurPopSize1 = s.CurPopSize1;
		NeedCentUpdate = s.NeedCentUpdate;

		int i;

		for( i = 0; i < CurPopSize; i++ )
		{
			memcpy( CurParams[ i ], s.CurParams[ i ],
				ParamCount * sizeof( CurParams[ i ][ 0 ]));
		}

		memcpy( CurCosts, s.CurCosts, CurPopSize * sizeof( CurCosts[ 0 ]));

		if( !NeedCentUpdate )
		{
			memcpy( CentParams, s.CentParams, ParamCount *
				sizeof( CentParams[ 0 ]));
		}
	}

	/**
	 * Function resets centroid vector.
	 */

	void resetCentroid()
	{
		memset( CentParams, 0, ParamCount * sizeof( CentParams[ 0 ]));
	}

	/**
	 * Function recalculates centroid based on the current population size.
	 * The NeedCentUpdate variable can be checked if centroid update is
	 * needed. This function resets the NeedCentUpdate to "false".
	 */

	void updateCentroid()
	{
		resetCentroid();

		double* const cp = CentParams;
		int i;
		int j;

		for( j = 0; j < CurPopSize; j++ )
		{
			const double* const p = CurParams[ j ];

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += p[ i ];
			}
		}

		const double m = 1.0 / CurPopSize;

		for( i = 0; i < ParamCount; i++ )
		{
			cp[ i ] *= m;
		}

		NeedCentUpdate = false;
	}

	/**
	 * Function returns pointer to the centroid vector. The NeedUpdateCent
	 * should be checked and and if it is equal to "true", the
	 * updateCentroid() function called.
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
		return( CurParams[ i ]);
	}

	/**
	 * Function returns current population size.
	 */

	int getCurPopSize() const
	{
		return( CurPopSize );
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
	 * Function replaces the highest-cost previous solution, updates centroid.
	 *
	 * @param NewCost Cost of the new solution.
	 * @param UpdParams New parameter values.
	 * @param DoUpdateCentroid "True" if centroid should be updated using
	 * running sum. This update is done for parallel populations.
	 * @param DoCostCheck "True" if the cost contraint should be checked.
	 * Function returns "false" if the cost constraint was not met, "true"
	 * otherwise.
	 */

	bool updatePop( const double NewCost, const double* const UpdParams,
		const bool DoUpdateCentroid, const bool DoCostCheck )
	{
		if( DoCostCheck )
		{
			if( NewCost > CurCosts[ CurPopSize1 ])
			{
				return( false );
			}
		}

		double* const rp = CurParams[ CurPopSize1 ];

		if( DoUpdateCentroid )
		{
			double* const cp = CentParams;
			const double m = 1.0 / CurPopSize;
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += ( UpdParams[ i ] - rp[ i ]) * m;
				rp[ i ] = UpdParams[ i ];
			}
		}
		else
		{
			memcpy( rp, UpdParams, ParamCount * sizeof( rp[ 0 ]));
			NeedCentUpdate = true;
		}

		sortPop( NewCost, CurPopSize1 );

		return( true );
	}

protected:
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int PopSize; ///< The size of population in use (maximal).
		///<
	int PopSize1; ///< = PopSize - 1.
		///<
	int CurPopSize; ///< Current population size.
		///<
	int CurPopSize1; ///< = CurPopSize - 1.
		///<
	bool NeedCentUpdate; ///< "True" if centroid update is needed.
		///<
	double* CurParamsBuf; ///< CurParams buffer.
		///<
	double** CurParams; ///< Current working parameter vectors. Always kept
		///< sorted in ascending cost order.
		///<
	double* CurCosts; ///< Best costs of current working parameter vectors.
		///<
	double* CentParams; ///< Centroid of the current parameter vectors.
		///<

	/**
	 * Function performs re-sorting of the population based on the cost of a
	 * newly-added solution, and stores new cost.
	 *
	 * @param Cost Solution's cost.
	 * @param i Solution's index (usually, CurPopSize1).
	 */

	void sortPop( const double Cost, int i )
	{
		double* const InsertParams = CurParams[ i ];

		while( i > 0 )
		{
			const double c1 = CurCosts[ i - 1 ];

			if( c1 < Cost )
			{
				break;
			}

			CurCosts[ i ] = c1;
			CurParams[ i ] = CurParams[ i - 1 ];
			i--;
		}

		CurCosts[ i ] = Cost;
		CurParams[ i ] = InsertParams;
	}

private:
	/**
	 * Function deletes buffers previously allocated via the initPopBuffers()
	 * function.
	 */

	void deletePopBuffers()
	{
		delete[] CurParamsBuf;
		delete[] CurParams;
		delete[] CurCosts;
		delete[] CentParams;
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
	{
	}

	virtual ~CBiteOptBase()
	{
		deleteBuffers();
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
	int StallCount; ///< The number of iterations without improvement.
		///<
	double HiBound; ///< Higher cost bound, for StallCount estimation. May not
		///< be used by the optimizer.
		///<
	double AvgCost; ///< Average cost in the latest batch. May not be used by
		///< the optimizer.
		///<
	static const int MaxApplyHists = 16; /// The maximal number of histograms
		///< that can be used during the optimize() function call.
		///<
	CBiteOptHistBase* ApplyHists[ MaxApplyHists ]; ///< Histograms used in
		///< "selects" during the optimize() function call.
		///<
	int ApplyHistsCount; ///< The number of "selects" used during the
		///< optimize() function call.
		///<

	/**
	 * Function initializes all common buffers, and "PopSize" variables.
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
			CurPopSize = 0;
			CurPopSize1 = 0;
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
		CurPopSize = PopSize;
		CurPopSize1 = PopSize1;
		NeedCentUpdate = false;
		BestCost = 1e300;
		StallCount = 0;
		HiBound = 1e300;
		AvgCost = 0.0;
		ApplyHistsCount = 0;
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
		if( NewCost <= BestCost )
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
					ParamCount * sizeof( BestParams[ 0 ]));
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

	static double wrapParam( CBiteRnd& rnd, const double v )
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
	 * Function performs random choice selection based on histogram's choice
	 * count, and adds the histogram to apply list.
	 *
	 * @param Hist Histogram.
	 * @param rnd PRNG object.
	 */

	template< class T >
	int selectRandom( T& Hist, CBiteRnd& rnd )
	{
		ApplyHists[ ApplyHistsCount ] = &Hist;
		ApplyHistsCount++;

		return( Hist.selectRandom( rnd ));
	}

	/**
	 * Function forces a specified choice selection to the histogram, and adds
	 * the histogram to apply list.
	 *
	 * @param Hist Histogram.
	 * @param s Choice index. Will be returned unchanged.
	 */

	template< class T >
	int selectForce( T& Hist, const int s )
	{
		ApplyHists[ ApplyHistsCount ] = &Hist;
		ApplyHistsCount++;

		return( Hist.selectForce( s ));
	}

	/**
	 * Function performs choice selection based on the specified histogram,
	 * and adds the histogram to apply list. Specialized for binary
	 * histograms.
	 *
	 * @param Hist Histogram.
	 * @param rnd PRNG object.
	 */

	int select( CBiteOptHistBinary& Hist, CBiteRnd& rnd )
	{
		ApplyHists[ ApplyHistsCount ] = &Hist;
		ApplyHistsCount++;

		return( Hist.select( rnd ));
	}

	/**
	 * Function "unselects" a previously selected choice in the specified
	 * histogram so that the choice on either the applyHistsIncr() or
	 * applyHistsDecr() call is randomized.
	 *
	 * @param Hist Histogram.
	 * @param rnd PRNG object.
	 */

	template< class T >
	void unselect( T& Hist, CBiteRnd& rnd )
	{
		Hist.unselect( rnd );
	}

	/**
	 * Function applies histogram increments on optimization success.
	 */

	void applyHistsIncr()
	{
		const int c = ApplyHistsCount;
		ApplyHistsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplyHists[ i ] -> incr();
		}
	}

	/**
	 * Function applies histogram decrements on optimization fail.
	 */

	void applyHistsDecr()
	{
		const int c = ApplyHistsCount;
		ApplyHistsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplyHists[ i ] -> decr();
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
