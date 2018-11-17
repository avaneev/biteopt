//$ nocpp

#ifndef SAESOPT_INCLUDED
#define SAESOPT_INCLUDED

#include "biteoptort.h"

/**
 * Sigma Adaptation Evolution Strategy class. Fundamentally similar to CMA-ES,
 * but performs sigma adaptation.
 */

class CSAESOpt
{
public:
	CSAESOpt()
		: ParamCount( 0 )
		, PopSize( 0 )
		, PopOrder( NULL )
		, CurParamsBuf( NULL )
		, CurParams( NULL )
		, CurCosts( NULL )
		, MinValues( NULL )
		, MaxValues( NULL )
		, BestParams( NULL )
	{
	}

	~CSAESOpt()
	{
		deleteBuffers();
	}

	/**
	 * Function updates dimensionality of *this object.

	 * @param aParamCount The number of parameters being optimized.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		if( aParamCount == ParamCount )
		{
			return;
		}

		deleteBuffers();

		ParamCount = aParamCount;
		PopSize = 25 + ParamCount * 2.0;
		PopSize1 = PopSize - 1;

		PopOrder = new int[ PopSize ];
		CurParamsBuf = new double[ PopSize * ParamCount ];
		CurParams = new double*[ PopSize ];
		CurCosts = new double[ PopSize ];
		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		BestParams = new double[ ParamCount ];

		int i;

		for( i = 0; i < PopSize; i++ )
		{
			CurParams[ i ] = CurParamsBuf + i * ParamCount;
		}

		Ort.updateDims( ParamCount, PopSize );
	}

	/**
	 * @return The number of initial objective function evaluations.
	 */

	int getInitEvals() const
	{
		return( 0 );
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );

		BestCost = 1e100;
		cure = 0;

		// Provide initial centroid and sigma (CurParams is used temporarily,
		// otherwise initially undefined).

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			CurParams[ 0 ][ i ] = ( MinValues[ i ] + MaxValues[ i ]) * 0.5;
			CurParams[ 1 ][ i ] = fabs( MaxValues[ i ] - MinValues[ i ]) / 6.0;
		}

		Ort.init( CurParams[ 0 ], CurParams[ 1 ]);
	}

	/**
	 * Function samples a random population vector based on the current
	 * distribution, with feasibility guarantee.
	 *
	 * @param rnd Random number generator.
	 * @param[out] op Resulting parameter vector.
	 */

	void sample( CBiteRnd& rnd, double* const op ) const
	{
		// Generate vector, check its feasibility, and resample it up to 10
		// times.

		int infcount = 0;
		int i;

		while( true )
		{
			Ort.sample( rnd, op );

			if( isFeasible( op ))
			{
				break;
			}

			infcount++;

			if( infcount == 10 )
			{
				// Force bound constraints.

				for( i = 0; i < ParamCount; i++ )
				{
					if( op[ i ] < MinValues[ i ])
					{
						op[ i ] = MinValues[ i ];
					}
					else
					if( op[ i ] > MaxValues[ i ])
					{
						op[ i ] = MaxValues[ i ];
					}
				}

				break;
			}
		}
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return The number of non-improving iterations so far. Always 0.
	 */

	int optimize( CBiteRnd& rnd )
	{
		double* const Params = CurParams[ cure ];
		sample( rnd, Params );

		const double NewCost = optcost( Params );

		if( NewCost < BestCost )
		{
			BestCost = NewCost;
			int i;

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = Params[ i ];
			}
		}

		insertPopOrder( NewCost, cure, cure );

		if( cure == PopSize1 )
		{
			cure = 0;
			Ort.update( CurParams, PopOrder );
		}
		else
		{
			cure++;
		}

		return( 0 );
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
	 * Virtual function (objective function) that should calculate parameter
	 * vector's optimization cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) = 0;

protected:
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int PopSize; ///< The size of population in use.
		///<
	int PopSize1; ///< = PopSize - 1.
		///<
	CBiteOptOrt Ort; ///< Rotation vector and orthogonalization calculator.
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
	double* MinValues; ///< Minimal parameter values.
		///<
	double* MaxValues; ///< Maximal parameter values.
		///<
	double* BestParams; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	int cure; ///< Current evaluation index, equals PopSize if population
		///< distribution needs to be updated.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	void deleteBuffers()
	{
		delete[] PopOrder;
		delete[] CurParamsBuf;
		delete[] CurParams;
		delete[] CurCosts;
		delete[] MinValues;
		delete[] MaxValues;
		delete[] BestParams;
	}

	/**
	 * Function returns "true" if the supplied vector is feasible.
	 *
	 * @param x Vector to check.
	 */

	bool isFeasible( const double* const x ) const
	{
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			if( x[ i ] < MinValues[ i ] || x[ i ] > MaxValues[ i ])
			{
				return( false );
			}
		}

		return( true );
	}

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
			const int mid = ( z + hi ) >> 1;

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

#endif // SAESOPT_INCLUDED
