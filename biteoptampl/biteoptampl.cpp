// AMPL interface for BiteOpt derivative-free optimization method.
//
// AMPL NL parser sources and compiled library "amplsolv" should be put into
// the "solvers" directory. AMPL parser can be acquired at
// https://ampl.com/netlib/ampl/
//
// The model `.mod` file should first be converted to the `.nl` format via
// the `ampl -og` command.

//$ lib "solvers/amplsolv"
//$ skip_include "z|solvers/getstub.h"

#define USE_SOLDB 0 // For internal use: 1 = use solution database.

#if USE_SOLDB
	#include "updsol.h"
#endif // USE_SOLDB

#include "solvers/getstub.h"
#include "../biteopt.h"
#include <math.h>

static fint depth = 8;
static fint attc = 10;
static int maxdim = 60;
static int nprob = -1;
static int progr = 0;
static double itmult = 1.0;
static double tol = 1e-5; // Constraint tolerance.

real* tmpx; // Temporary holder for solution and rounding.
real* tmpcon; // Temprorary holder for constraint bodies.
real* fc; // Temporary buffer of constraint penalties.
int con_notmet; // no. contraints not met.
double last_ov = 0.0; // Last evaluated obj value.
int negate; // negate objective.
ASL *asl;

static keyword keywds[] = {	/* must be in alphabetical order */
	KW("attc", L_val, &attc, "no. of attempts (default 10)"),
	KW("depth", L_val, &depth, "solver's depth (default 8), expected value 1 to 32"),
	KW("itmult", D_val, &itmult, "iteration number multiplier (default 1.0)"),
	KW("maxdim", I_val, &maxdim, "maximum number of dimensions to accept (default 60)"),
	KW("nprob", I_val, &nprob, "objective choice: 1 (default) = 1st"),
	KW("progr", I_val, &progr, "1 - print progress (default 0)"),
	KW("tol", D_val, &tol, "constraint tolerance (default 1e-5)"),
	KW("version", Ver_val, 0, "report version"),
	KW("wantsol", WS_val, 0, WSu_desc_ASL+5)
};

static char biteoptvers[] =
	"AMPL/BITEOPT\0\nAMPL/BITEOPT Driver Version 2022.19\n";

static Option_Info Oinfo = {
	"biteoptampl", "BITEOPT-2022.19", "biteopt_options", keywds, nkeywds, 1.,
	biteoptvers, 0,0,0,0,0, 202219
};

int xround( real* x, int n )
{
	int i;

	for( i = 0; i < n; i++ )
	{
		x[ i ] = ( x[ i ] < 0.0 ?
			-floor( 0.5 - x[ i ]) : floor( x[ i ] + 0.5 ));
	}

	return( n );
}

int solround( real* x )
{
	int nround = 0;
	int nint = niv + nbv;

	if( nint > 0 )
	{
		nround = xround( x + n_var - nint, nint );
	}

	if( nlvbi > 0 )
	{
		nround += xround( x + ( nlvb - nlvbi ), nlvbi );
	}

	if( nlvci > 0 )
	{
		nround += xround( x + ( nlvc - nlvci ), nlvci );
	}

	if( nlvoi > 0 )
	{
		nround += xround( x + ( nlvo - nlvoi ), nlvoi );
	}

	return( nround );
}

static double objfn( int N, const double* const x )
{
	memcpy( tmpx, x, N * sizeof( tmpx[ 0 ]));
	solround( tmpx );

	fint nerr = 0;
	double f = objval( nprob, tmpx, &nerr );
	f = ( nerr > 0 ? 1e300 : ( negate ? -f : f ));

	con_notmet = 0;

	if( n_con > 0 )
	{
		nerr = 0;
		conval( tmpx, tmpcon, &nerr );

		int i;

		if( nerr > 0 )
		{
			for( i = 0; i < n_con; i++ )
			{
				fc[ i ] = 0.0;
			}

			f += 1e300;
			con_notmet = n_con;
		}
		else
		{
			for( i = 0; i < n_con; i++ )
			{
				fc[ i ] = 0.0;

				if( LUrhs[ i ] > negInfinity && tmpcon[ i ] < LUrhs[ i ])
				{
					double a = LUrhs[ i ] - tmpcon[ i ];

					if( a > tol )
					{
						fc[ i ] = a - tol;
						con_notmet++;
					}
				}

				if( Urhsx[ i ] < Infinity && tmpcon[ i ] > Urhsx[ i ])
				{
					double a = tmpcon[ i ] - Urhsx[ i ];

					if( a > tol )
					{
						fc[ i ] = a - tol;
						con_notmet++;
					}
				}
			}
		}
	}

	return( f );
}

class CBiteOptAMPL : public CBiteOptDeep
{
public:
	virtual void getMinValues( double* const p ) const
	{
		int i;

		for( i = 0; i < n_var; i++ )
		{
			p[ i ] = LUv[ i ];
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		int i;

		for( i = 0; i < n_var; i++ )
		{
			p[ i ] = Uvx[ i ];
		}
	}

	virtual double optcost( const double* const p )
	{
		last_ov = objfn( n_var, p );

		if( con_notmet > 0 )
		{
			const double ps = pow( 3.0, 1.0 / n_con );
			const double pnsi = 1.0 / sqrt( (double) n_con );
			double pns = 0.0;
			double pnsm = 0.0;
			int i;

			for( i = 0; i < n_con; i++ )
			{
				pns = pns * ps + pnsi + fc[ i ] + fc[ i ] * fc[ i ] * fc[ i ];
				pnsm = pnsm * ps + pnsi;
			}

			return( last_ov + 1e10 * ( 1.0 + ( pns - pnsm )));
		}

		return( last_ov );
	}
};

int main( int argc, char* argv[])
{
	FILE *nl;
	char buf[ 2048 ], *stub;
	char* msg = "";
	int msgo = 0;
	int infc = 0;

	asl = ASL_alloc( ASL_read_fg );
	want_derivs = 0;

	stub = getstops( argv, &Oinfo );
	nl = jac0dim( stub, (fint) strlen( stub ));

	if( nprob < 0 )
	{
		nprob = 0;
	}
	else
	{
		nprob--;
	}

	if( nprob < 0 || nprob >= n_obj )
	{
		msg = "No objective.";
		goto done;
	}

	if( n_obj > 1 )
	{
		msgo += sprintf( buf + msgo,
			"Model contains %d objectives. This solver supports \n"
			"single-objective optimization only. You can specify the \n",
			"objective index via the \"nprob\" parameter.", n_obj );
	}

	if( n_con > 0 )
	{
		msgo += sprintf( buf + msgo, "Model contains %d constraints.\n",
			n_con );
	}

	if( n_var > maxdim )
	{
		msg = "Too many dimensions (more than \"maxdim\").";
		goto done;
	}

	X0 = (real*) Malloc(( 4 * n_var + 4 * n_con ) * sizeof( real ));
	LUv = X0 + n_var;
	Uvx = LUv + n_var;
	tmpx = Uvx + n_var;
	fc = tmpx + n_var;
	tmpcon = fc + n_con;
	LUrhs = tmpcon + n_con;
	Urhsx = LUrhs + n_con;
	fg_read( nl, 0 );

	negate = ( objtype[ nprob ] != 0 );
	int i;

	for( i = 0; i < n_var; i++ )
	{
		if( LUv[ i ] <= negInfinity )
		{
			infc++;
			LUv[ i ] = -1e9;
		}

		if( Uvx[ i ] >= Infinity )
		{
			infc++;
			Uvx[ i ] = 1e9;
		}
	}

	if( infc > 0 )
	{
		msgo += sprintf( buf + msgo,
			"Infinity var ranges were limited to [-1e9; 1e9] range.\n" );
	}

	goto start;

done:
	printf( "model: %s (%s)\n", stub,
		( negate ? "maximization" : "minimization" ));

	msgo += sprintf( buf + msgo, msg );

	write_sol( buf, X0, 0, &Oinfo );
	fclose( nl );

	return( 0 );

start:
	CBiteOptAMPL opt;
	opt.updateDims( n_var, depth );

	CBiteRnd rnd;
	rnd.init( 1 );

	int fnevals = 0;
	const double p = ( n_var >= 60 ? 1.6 : 4.0 * exp( -0.015 * n_var ));
	const int hardlim = (int) ( 1000.0 * itmult * pow( (double) n_var, p ) *
		sqrt( (double) depth ));

	const int sc_thresh = n_var * 512;

	#if USE_SOLDB
		bool UseBestSol = false;
		double MinObj = 0.0;

		VOXERRSKIP( loadBestSols() );

		CMap< CString, double > :: iterator it2 = BestSols.find( stub );

		if( it2 != BestSols.end() )
		{
			MinObj = it2.value();
			UseBestSol = true;
		}
	#endif // USE_SOLDB

	double f;
	int f_notmet;
	int f_iters;
	int khl = 0;
	int k;

	for( k = 0; k < attc; k++ )
	{
		opt.init( rnd );

		double thr = 200.0;
		i = 0;

		while( i < hardlim )
		{
			int sc = opt.optimize( rnd );
			fnevals++;
			i++;

			#if USE_SOLDB
			if( UseBestSol && con_notmet == 0 )
			{
				const double d = ( MinObj == 0.0 ?
					SolTol0 : fabs( MinObj ) * SolTol );

				if( negate )
				{
					if( -last_ov >= MinObj - d )
					{
						k = attc;
						break;
					}
				}
				else
				{
					if( last_ov <= MinObj + d )
					{
						k = attc;
						break;
					}
				}
			}
			#endif // USE_SOLDB

			if( progr && i > thr )
			{
				const double nf = objfn( n_var, opt.getBestParams() );

				printf( "Attempt %i/%i n_var=%i n_con=%i iter %i, "
					"Objective = %.15g\n", k + 1, attc, n_var, n_con, i,
					( negate ? -nf : nf ));

				thr *= 1.4;
			}

			if( sc > sc_thresh )
			{
				break;
			}
		}

		if( i >= hardlim )
		{
			khl++;
		}

		const double nf = objfn( n_var, opt.getBestParams() );

		if( k == 0 || con_notmet < f_notmet ||
			( con_notmet == f_notmet && nf < f ))
		{
			f_iters = i;

			for( i = 0; i < n_var; i++ )
			{
				X0[ i ] = opt.getBestParams()[ i ];
			}

			f = nf;
			f_notmet = con_notmet;
		}
	}

	msgo += sprintf( buf + msgo, "%s:\n%ld function evaluations "
		"(%ld attempts, depth=%i)\n", Oinfo.bsname, fnevals, attc, depth );

	msgo += sprintf( buf + msgo, "Hard iteration limit is %i per attempt.\n",
		hardlim );

	msgo += sprintf( buf + msgo, "Hard iteration limit achieved in %i of %i "
		"attempts.\n", khl, attc );

	f = ( negate ? -f : f );

	msgo += sprintf( buf + msgo, "Objective = %.*g\n", obj_prec(), f );

	if( f_notmet > 0 )
	{
		msgo += sprintf( buf + msgo,
			"!!! %i constraint(s) not met, infeasible solution\n", f_notmet );
	}

	#if USE_SOLDB
		VOXERRSKIP( updateSol( stub, n_var, n_con, negate, f, f_notmet,
			f_iters, khl ));
	#endif // USE_SOLDB

	solround( X0 );
	goto done;
}
