//$ nocpp

#ifndef CCMAES_INCLUDED
#define CCMAES_INCLUDED

#include "../biteaux.h"
#include "../cmaes/cmaes_interface.h"

/**
 * Interface to CMA-ES optimization method.
 */

class CCMAESOpt : public CBiteOptInterface
{
public:
	CCMAESOpt()
		: ParamCount( 0 )
		, IsInit( false )
		, MinValues( NULL )
		, MaxValues( NULL )
		, BestParams( NULL )
		, Params( NULL )
	{
	}

	virtual ~CCMAESOpt()
	{
		if( IsInit )
		{
			cmaes_exit( &cma );
		}

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

	/**
	 * Function updates dimensionality of *this object. Function does nothing

	 * @param aParamCount The number of parameters being optimized.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		if( aParamCount == ParamCount )
		{
			return;
		}

		lambda = 20 + ParamCount * 2;

		deleteBuffers();

		ParamCount = aParamCount;
		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		BestParams = new double[ ParamCount ];
		Params = new double[ ParamCount ];
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams Initial parameter values.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL )
	{
		if( IsInit )
		{
			cmaes_exit( &cma );
		}

		getMinValues( MinValues );
		getMaxValues( MaxValues );

		double* const initial_solution = BestParams;
		double* const initial_sigma = Params;
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			initial_solution[ i ] = MinValues[ i ] +
				( MaxValues[ i ] - MinValues[ i ]) * rnd.getRndValue();

			initial_sigma[ i ] = ( MaxValues[ i ] - MinValues[ i ]) / 6.0;
		}

		cmaes_init_para( &cma, ParamCount, initial_solution, initial_sigma,
			1 + (int) ( rnd.getRndValue() * 1000000.0 ), lambda, "no" );

		cma.sp.filename = strdup( "no" );
		y = cmaes_init_final( &cma );
		IsInit = true;
		curx = lambda;
		WasBestCostSet = false;
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return The number of non-improving iterations so far. The plateau
	 * threshold value is ParamCount * 16.
	 */

	int optimize( CBiteRnd& rnd )
	{
		if( curx == lambda )
		{
			X = cmaes_SamplePopulation( &cma );
			curx = 0;
		}

		int infcount = 0;

		while( !is_feasible( X[ curx ], MinValues, MaxValues, ParamCount ) &&
			infcount < 10 )
		{
			cmaes_ReSampleSingle( &cma, curx );
			infcount++;
		}

		const double* const p = X[ curx ];
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			if( p[ i ] < MinValues[ i ])
			{
				Params[ i ] = MinValues[ i ];
			}
			else
			if( p[ i ] > MaxValues[ i ])
			{
				Params[ i ] = MaxValues[ i ];
			}
			else
			{
				Params[ i ] = p[ i ];
			}
		}

		y[ curx ] = optcost( Params );

		if( !WasBestCostSet || y[ curx ] <= BestCost )
		{
			WasBestCostSet = true;
			BestCost = y[ curx ];

			for( i = 0; i < ParamCount; i++ )
			{
				BestParams[ i ] = Params[ i ];
			}
		}

		curx++;

		if( curx == lambda )
		{
			cmaes_UpdateDistribution( &cma, y );
		}

		return( 0 );
	}

protected:
	int ParamCount; ///< The total number of internal parameter values in use.
		///<
	int lambda; ///< CMA-ES "lambda" parameter, specifies population size.
		///<
	cmaes_t cma; ///< CMA-ES object.
		///<
	bool IsInit; ///< "True" if "cma" object was initialized.
		///<
	double* y; ///< Array of objective function values (size = lambda).
		///<
	double* const* X; ///< Array of parameter vectors (size = lambda).
		///<
	int curx; ///< Current parameter vector being evaluated, equals lambda
		///< if population distribution needed to be updated.
		///<
	double* MinValues; ///< Minimal parameter values.
		///<
	double* MaxValues; ///< Maximal parameter values.
		///<
	double* BestParams; ///< Best parameter vector.
		///<
	double BestCost; ///< Cost of the best parameter vector.
		///<
	bool WasBestCostSet; ///< "True" if initial BestCost was set.
		///<
	double* Params; ///< Temporary parameter buffer.
		///<

	/**
	 * Function deletes previously allocated buffers.
	 */

	void deleteBuffers()
	{
		delete[] MinValues;
		delete[] MaxValues;
		delete[] BestParams;
		delete[] Params;
	}

	static bool is_feasible( const double *x, const double *lower,
		const double *upper, size_t number_of_variables )
	{
		for( size_t i = 0; i < number_of_variables; ++i )
		{
			if( x[i] < lower[i] || x[i] > upper[i])
			{
				return( false );
			}
		}

		return( true );
	}
};

#endif // CCMAES_INCLUDED
