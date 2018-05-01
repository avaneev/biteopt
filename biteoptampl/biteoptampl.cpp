// AMPL interface for BiteOpt derivative-free optimization method.
//
// AMPL NL parser sources and compiled library "amplsolv" should be put into
// the "solvers" directory. AMPL parser can be acquired at
// http://www.netlib.org/ampl/

//$ lib "solvers/amplsolv"
//$ skip_include "z|solvers/getstub.h"

#include "solvers/getstub.h"
#include "../biteopt.h"

static fint depth = 9;
static fint attc = 10;
static int nprob = -1;
static double itmult = 1.0;
real* tmpx; // Temporary holder for solution and rounding.
real* tmpcon; // Temprorary holder for constraint bodies.
int tmpcon_notmet; // no. contraints not met.
int negate;
ASL *asl;

static keyword keywds[] = {	/* must be in alphabetical order */
	KW("attcnt", L_val, &attc, "no. of attempts (default 10)"),
	KW("depth", L_val, &depth, "solver's depth (default 9), expected value 1 to 32"),
	KW("itmult", D_val, &itmult, "iteration number multiplier (default 1.0)"),
	KW("nprob", I_val, &nprob, "objective choice: 1 (default) = 1st"),
	KW("version", Ver_val, 0, "report version"),
	KW("wantsol", WS_val, 0, WSu_desc_ASL+5)
};

static char biteoptvers[] =
	"AMPL/BITEOPT\0\nAMPL/BITEOPT Driver Version 20180501\n";

static Option_Info Oinfo = {
	"biteoptampl", "BITEOPT", "biteopt_options", keywds, nkeywds, 1.,
	biteoptvers, 0,0,0,0,0, 20180501
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
	int nint;
	int nround = 0;
	real* x2;

	if( nint = niv + nbv )
	{
		x2 = x + n_var - nint;
		nround = xround( x2, nint );
	}

	if( nint = nlvbi )
	{
		x2 = x + ( nlvb - nint );
		nround += xround( x2, nint );
	}

	if( nint = nlvci )
	{
		x2 = x + ( nlvc - nint );
		nround += xround( x2, nint );
	}

	if( nint = nlvoi )
	{
		x2 = x + (nlvo - nint);
		nround += xround( x2, nint );
	}

	return( nround );
}

static double objfn( int N, const double* x )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		tmpx[ i ] = x[ i ];
	}

	solround( tmpx );

	double f = objval( nprob, tmpx, 0 );
	f = ( negate ? -f : f );
	tmpcon_notmet = 0;

	if( n_con > 0 )
	{
		conval( tmpx, tmpcon, 0 );

		for( i = 0; i < n_con; i++ )
		{
			if( fabs( LUrhs[ i ] - Urhsx[ i ]) <= 1e-11 )
			{
				double a = fabs( tmpcon[ i ] - LUrhs[ i ]);

				if( a > 1e-6 )
				{
					f += ( a + a * a ) * 10000.0;
					tmpcon_notmet++;
				}
			}
			else
			{
				if( LUrhs[ i ] > negInfinity &&
					tmpcon[ i ] < LUrhs[ i ])
				{
					double a = LUrhs[ i ] - tmpcon[ i ];
					f += ( a + a * a + a * a * a ) * 10000.0;
					tmpcon_notmet++;
				}

				if( Urhsx[ i ] < Infinity &&
					tmpcon[ i ] > Urhsx[ i ])
				{
					double a = tmpcon[ i ] - Urhsx[ i ];
					f += ( a + a * a + a * a * a ) * 10000.0;
					tmpcon_notmet++;
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
		return( objfn( n_var, p ));
	}
};

int main( int argc, char* argv[])
{
	FILE *nl;
	char buf[ 2048 ], *stub;
	char* msg = "";
	int msgo = 0;

	asl = ASL_alloc( ASL_read_fg );
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

	if( n_con > 0 )
	{
		msgo += sprintf( buf + msgo,
			"NOTE: model contains %d constraints. Constraint support \n"
			"of this solver is experimental, constraints are applied as \n"
			"penalties to the objective function.\n", n_con );
	}

	X0 = (real*) Malloc(( 4 * n_var + 3 * n_con ) * sizeof( real ));
	LUv = X0 + n_var;
	Uvx = LUv + n_var;
	tmpx = Uvx + n_var;
	tmpcon = tmpx + n_var;
	LUrhs = tmpcon + n_con;
	Urhsx = LUrhs + n_con;
	fg_read( nl, 0 );

	negate = ( objtype[ nprob ] != 0 );
	int i;
	int infc = 0;

	for( i = 0; i < n_var; i++ )
	{
		if( LUv[ i ] <= negInfinity )
		{
			infc++;
			LUv[ i ] = -1e8;
		}

		if( Uvx[ i ] >= Infinity )
		{
			infc++;
			Uvx[ i ] = 1e8;
		}
	}

	if( infc > 0 )
	{
		msgo += sprintf( buf + msgo,
			"Infinity var ranges were limited to [-1e8; 1e8] range.\n" );
	}

	double f;
	CBiteOptAMPL opt;
	opt.updateDims( n_var, depth );

	CBiteRnd rnd;
	rnd.init( 1 );

	int fnevals = 0;
	const int hardlim = (int) ( itmult * 2000.0 * pow( (double) n_var, 1.75 ) *
		sqrt( (double) depth ));

	const int sc_thresh = (int) ( opt.getInitEvals() * 10.0 / depth );

	int kmet = ( n_con > 0 ? 0 : 1 );
	int k;

	for( k = 0; k < attc; k++ )
	{
		opt.init( rnd );
		fnevals += opt.getInitEvals();
		i = 0;
		double thr = 200.0;

		while( i < hardlim )
		{
			int sc = opt.optimize( rnd );
			fnevals++;
			i++;

/*			if( i > (int) thr )
			{
				printf( "Attempt %i/%i n_vars=%i n_cosntr=%i iter %i, "
					"obj.value = %.15g\n", k + 1, attc, n_var, n_con, i,
					opt.getBestCost() );

				thr *= 1.4;
			}
*/
			if( sc > sc_thresh )
			{
				break;
			}
		}

		objfn( n_var, opt.getBestParams() );

		if( k == 0 || ( kmet == 0 && tmpcon_notmet == 0 ) ||
			( opt.getBestCost() < f && tmpcon_notmet == 0 ))
		{
			for( i = 0; i < n_var; i++ )
			{
				X0[ i ] = opt.getBestParams()[ i ];
			}

			f = opt.getBestCost();
			kmet += ( tmpcon_notmet == 0 ? 1 : 0 );
		}
	}

	msgo += sprintf( buf + msgo, "%s:\n%ld function evaluations "
		"(%ld attempts, depth=%i)\n", Oinfo.bsname, fnevals, attc, depth );

	solround( X0 );
	f = objfn( n_var, X0 );

	msgo += sprintf( buf + msgo, "Objective = %.*g\n", obj_prec(),
		negate ? -f : f );

	if( n_con > 0 && tmpcon_notmet > 0 )
	{
		msgo += sprintf( buf + msgo,
			"(%i constraint(s) not met, infeasible solution)\n",
			tmpcon_notmet );
	}

done:
	printf( "model: %s\n", stub );
	msgo += sprintf( buf + msgo, msg );
	write_sol( buf, X0, 0, &Oinfo );
	fclose( nl );
}
