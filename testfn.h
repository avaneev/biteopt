//$ nocpp

// Global optimization test functions.
//
// Sources:
// http://infinity77.net/global_optimization/test_functions.html
// https://www.sfu.ca/~ssurjano/optimization.html
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page364.htm
// http://benchmarkfcns.xyz/fcns
// http://al-roomi.org/benchmarks/unconstrained/2-dimensions/120-damavandi-s-function
// https://github.com/numbbo/coco/tree/master/code-experiments/src
// https://github.com/alusov/mathexplib/tree/master/testfuncs

#include <math.h>

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

#if !defined( M_PI )
	#define M_PI 3.14159265358979324
#endif // !defined( M_PI )

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

	double (*ParamFunc)( double* const minv, double* const maxv,
		const int N ); ///< Range and optimal value function, can be NULL.
		///<
};

static double calcThreeHumpCamel( const double* const x, const int N )
{
	return( 2*sqr(x[0])-1.05*sqr(sqr(x[0]))+sqr(sqr(sqr(x[0])))/6+
		x[0]*x[1]+sqr(x[1]) );
}

static const CTestFn TestFnThreeHumpCamel = { "ThreeHumpCamel", 2, -10.0,
	10.0, 0.0, &calcThreeHumpCamel };

static double calcBooth( const double* const x, const int N )
{
	return( sqr(x[0]+2*x[1]-7)+sqr(2*x[0]+x[1]-5));
}

static const CTestFn TestFnBooth = { "Booth", 2, -10.0, 10.0, 0.0,
	&calcBooth };

static double calcMatyas( const double* const x, const int N )
{
	return( 0.26*(sqr(x[0])+sqr(x[1]))-0.48*x[0]*x[1] );
}

static const CTestFn TestFnMatyas = { "Matyas", 2, -10.0, 10.0, 0.0,
	&calcMatyas };

static double calcSphere( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( s );
}

static const CTestFn TestFnSphere = { "Sphere", 0, -10.0, 10.0, 0.0,
	&calcSphere };

static double calcLevy13( const double* const x, const int N )
{
	return( sqr(x[0]-1.0)*(sqr(sin(3.0*M_PI*x[1]))+1.0)+
		sqr(x[1]-1.0)*(sqr(sin(2.0*M_PI*x[1]))+1.0)+sqr(sin(3.0*M_PI*x[0])));
}

static const CTestFn TestFnLevy13 = { "Levy13", 2, -10.0, 10.0, 0.0,
	&calcLevy13 };

static double calcSchaffer01( const double* const x, const int N )
{
	return( 0.5+(sqr(sin(sqr(sqr(x[0])+sqr(x[1]))))-0.5)/
		sqr(1.0+0.001*sqr(sqr(x[0])+sqr(x[1]))));
}

static const CTestFn TestFnSchaffer01 = { "Schaffer01", 2, -100.0, 100.0,
	0.0, &calcSchaffer01 };

static double calcSchaffer02( const double* const x, const int N )
{
	return( 0.5+(sqr(sin(sqr(x[0])-sqr(x[1])))-0.5)/
		sqr(1.0+0.001*sqr(sqr(x[0])+sqr(x[1]))));
}

static const CTestFn TestFnSchaffer02 = { "Schaffer02", 2, -100.0, 100.0,
	0.0, &calcSchaffer02 };

static double calcSchaffer03( const double* const x, const int N )
{
	return( 0.5+(sqr(sin(cos(fabs(sqr(x[0])-sqr(x[1])))))-0.5)/
		sqr(1.0+0.001*sqr(sqr(x[0])+sqr(x[1]))));
}

static const CTestFn TestFnSchaffer03 = { "Schaffer03", 2, -100.0, 100.0,
	0.0024558581695, &calcSchaffer03 };

static double calcSchaffer04( const double* const x, const int N )
{
	return( 0.5+(sqr(cos(sin(fabs(sqr(x[0])-sqr(x[1])))))-0.5)/
		sqr(1.0+0.001*sqr(sqr(x[0])+sqr(x[1]))));
}

static const CTestFn TestFnSchaffer04 = { "Schaffer04", 2, -100.0, 100.0,
	0.2929486652748, &calcSchaffer04 };

static double calcSchaffer06( const double* const x, const int N )
{
	return( 0.5+(sqr(sin(sqrt(sqr(x[0])+sqr(x[1]))))-0.5)/
		sqr(1.0+0.001*(sqr(x[0])+sqr(x[1]))) );
}

static const CTestFn TestFnSchaffer06 = { "Schaffer06", 2, -100.0, 100.0, 0.0,
	&calcSchaffer06 };

static double calcAckley( const double* const x, const int N )
{
	const double a = 20.0;
	const double b = 0.2;
	const double c = 2.0 * M_PI;
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr(x[i]);
		s2 += cos(c*x[i]);
	}

	return( -a*exp(-b*sqrt(s1/N))-exp(s2/N)+a+exp(1.0));
}

static const CTestFn TestFnAckley = { "Ackley", 0, -32.0, 32.0, 0.0,
	&calcAckley };

static double calcAckley2( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( -200.0*exp(-0.02*sqrt(s)) );
}

static const CTestFn TestFnAckley2 = { "Ackley2", 0, -32.0, 32.0, -200.0,
	&calcAckley2 };

static double calcAckley3( const double* const x, const int N )
{
	return( -200.0*exp(-0.02*sqrt(sqr(x[0])+sqr(x[1])))+
		5.0*exp(cos(3.0*x[0])+sin(3.0*x[1])) );
}

static const CTestFn TestFnAckley3 = { "Ackley3", 2, -32.0, 32.0,
	-195.6290282622794, &calcAckley3 };

static double calcAckley4( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += exp(1.0/pow(sqr(x[i])+sqr(x[i+1]),5.0))+
			3*(cos(2.0*x[i])+sin(2.0*x[i+1]));
	}

	return( s );
}

static const CTestFn TestFnAckley4 = { "Ackley4", 2, -32.0, 32.0, -5.0,
	&calcAckley4 };

static double calcRosenbrock( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += 100.0*sqr(sqr(x[i])-x[i+1])+sqr(x[i]-1.0);
	}

	return( s );
}

static const CTestFn TestFnRosenbrock = { "Rosenbrock", 0, -5.0, 10.0, 0.0,
	&calcRosenbrock };

static double calcBeale( const double* const x, const int N )
{
	return( sqr(1.5-x[0]+x[0]*x[1])+sqr(2.25-x[0]+x[0]*sqr(x[1]))+
		sqr(2.625-x[0]+x[0]*pow(x[1],3.0) ));
}

static const CTestFn TestFnBeale = { "Beale", 2, -10.0, 10.0, 0.0,
	&calcBeale };

static double calcBohachevsky1( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(x[i])+2.0*sqr(x[i+1])-0.3*cos(3.0*M_PI*x[i])-
			0.4*cos(4.0*M_PI*x[i+1])+0.7;
	}

	return( s );
}

static const CTestFn TestFnBohachevsky1 = { "Bohachevsky1", 0, -15.0, 15.0,
	0.0, &calcBohachevsky1 };

static double calcBohachevsky2( const double* const x, const int N )
{
	return( sqr(x[0])+2.0*sqr(x[1])-0.3*cos(3.0*M_PI*x[0])*
		cos(4.0*M_PI*x[1])+0.3 );
}

static const CTestFn TestFnBohachevsky2 = { "Bohachevsky2", 2, -100.0, 100.0,
	0.0, &calcBohachevsky2 };

static double calcBohachevsky3( const double* const x, const int N )
{
	return( sqr(x[0])+2.0*sqr(x[1])-0.3*cos(3.0*M_PI*x[0]+4.0*M_PI*x[1])+0.3 );
}

static const CTestFn TestFnBohachevsky3 = { "Bohachevsky3", 2, -100.0, 100.0,
	0.0, &calcBohachevsky3 };

static double calcEasomN( const double* const x, const int N )
{
	const double a = 20.0;
	const double b = 0.2;
	const double c = 2.0 * M_PI;
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr(x[i]);
		s2 += cos(c*x[i]);
	}

	return( a-a/exp(b*sqrt(s1/N))+exp(1.0)-exp(s2/N));
}

static const CTestFn TestFnEasomN = { "EasomN", 0, -100.0, 100.0, 0.0,
	&calcEasomN };

static double calcEasom( const double* const x, const int N )
{
	return( -cos(x[0])*cos(x[1])*exp(-sqr(x[0]-M_PI)-sqr(x[1]-M_PI)));
}

static const CTestFn TestFnEasom = { "Easom", 2, -100.0, 100.0, -1.0,
	&calcEasom };

static double calcCrossInTray( const double* const x, const int N )
{
	return( -0.0001*pow(fabs(sin(x[0])*sin(x[1])*
		exp(fabs(100.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI)))+1.0,0.1));
}

static const CTestFn TestFnCrossInTray = { "CrossInTray", 2, -15.0, 15.0,
	-2.0626118708227, &calcCrossInTray };

static double calcRastrigin( const double* const x, const int N )
{
	double s = 10.0*N;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i])-10.0*cos(2.0*M_PI*x[i]);
	}

	return( s );
}

static const CTestFn TestFnRastrigin = { "Rastrigin", 0, -5.12, 5.12, 0.0,
	&calcRastrigin };

static double calcDropWave( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( -(cos(12.0*sqrt(s))+1.0)/(2.0+0.5*s));
}

static const CTestFn TestFnDropWave = { "DropWave", 0, -5.12, 5.12, -1.0,
	&calcDropWave };

static double calcSumSquares( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += (i+1.0)*sqr(x[i]);
	}

	return( s );
}

static const CTestFn TestFnSumSquares = { "SumSquares", 0, -5.12, 5.12, 0.0,
	&calcSumSquares };

static double calcZacharov( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr(x[i]);
		s2 += (i+1.0)*x[i];
	}

	return( s1+sqr(0.5*s2)+sqr(sqr(0.5*s2)));
}

static const CTestFn TestFnZacharov = { "Zacharov", 0, -5.0, 10.0, 0.0,
	&calcZacharov };

static double calcRotatedHyperEllipsoid( const double* const x, const int N )
{
	double s = 0.0;
	int i;
	int j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j <= i; j++ )
		{
			s += sqr(x[j]);
		}
	}

	return( s );
}

static const CTestFn TestFnRotatedHyperEllipsoid = { "RotatHyperEllips", 0,
	-65.536, 65.536, 0.0, &calcRotatedHyperEllipsoid };

static double calcGriewank( const double* const x, const int N )
{
	double s1 = 0.0;
	double m2 = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr(x[i]);
		m2 *= cos(x[i]/sqrt(i+1.0));
	}

	return( s1/4000.0-m2+1.0 );
}

static const CTestFn TestFnGriewank = { "Griewank", 0, -600.0, 600.0, 0.0,
	&calcGriewank };

static double calcGoldsteinPrice( const double* const x, const int N )
{
	return( (1+sqr(x[0]+x[1]+1)*(19-14*x[0]+3*sqr(x[0])-14*x[1]+6*x[0]*x[1]+
		3*sqr(x[1])))*(30+sqr(2*x[0]-3*x[1])*(18-32*x[0]+12*sqr(x[0])+48*x[1]-
		36*x[0]*x[1]+27*sqr(x[1]))));
}

static const CTestFn TestFnGoldsteinPrice = { "GoldsteinPrice", 2, -2.0, 2.0,
	3.0, &calcGoldsteinPrice };

static double calcSalomon( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	s = sqrt( s );

	return( 1.0-cos(2.0*M_PI*s)+0.1*s);
}

static const CTestFn TestFnSalomon = { "Salomon", 0, -100.0, 100.0, 0.0,
	&calcSalomon };

static double calcBranin01( const double* const x, const int N )
{
	return( sqr(-1.275*sqr(x[0])/(M_PI*M_PI)+5.0*x[0]/M_PI+x[1]-6.0)+
		(10.0-5.0/(4.0*M_PI))*cos(x[0])+10.0 );
}

static const CTestFn TestFnBranin01 = { "Branin01", 2, -5.0, 10.0,
	0.3978873577297, &calcBranin01 };

static double calcBranin02( const double* const x, const int N )
{
	return( sqr(-1.275*sqr(x[0])/(M_PI*M_PI)+5.0*x[0]/M_PI+x[1]-6.0)+
		(10.0-5.0/(4.0*M_PI))*cos(x[0])*cos(x[1])+
		log(sqr(x[0])+sqr(x[1])+1.0)+10.0 );
}

static const CTestFn TestFnBranin02 = { "Branin02", 2, -5.0, 15.0,
	5.5589144038938, &calcBranin02 };

static double calcTrefethen( const double* const x, const int N )
{
	return( 0.25*sqr(x[0])+0.25*sqr(x[1])+exp(sin(50.0*x[0]))-
		sin(10.0*x[0]+10.0*x[1])+sin(60.0*exp(x[1]))+sin(70.0*sin(x[0]))+
		sin(sin(80.0*x[1])) );
}

static const CTestFn TestFnTrefethen = { "Trefethen", 2, -10.0, 10.0,
	-3.3068686474752, &calcTrefethen };

static double calcWhitley( const double* const x, const int N )
{
	double s = 0.0;
	int j;
	int i;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < N; j++ )
		{
			const double v = 100.0*sqr(sqr(x[i])-x[j])+sqr(1.0-x[j]);
			s += sqr(v)/4000.0-cos(v)+1.0;
		}
	}

	return( s );
}

static const CTestFn TestFnWhitley = { "Whitley", 0, -10.24, 10.24, 0.0,
	&calcWhitley };

static double calcPrice02( const double* const x, const int N )
{
	return( 1.0+sqr(sin(x[0]))+sqr(sin(x[1]))-0.1*exp(-sqr(x[0])-sqr(x[1])) );
}

static const CTestFn TestFnPrice02 = { "Price02", 2, -10.0, 10.0, 0.9,
	&calcPrice02 };

static double calcWavy( const double* const x, const int N )
{
	const double k = 10.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += cos(k*x[i])*exp(-0.5*sqr(x[i]));
	}

	return( 1.0-s/N );
}

static const CTestFn TestFnWavy = { "Wavy", 0, -M_PI, M_PI, 0.0, &calcWavy };

static double calcShubert01( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 1; i <= 5; i++ )
	{
		s1 += i*cos((i+1)*x[0]+i);
		s2 += i*cos((i+1)*x[1]+i);
	}

	return( s1*s2 );
}

static const CTestFn TestFnShubert01 = { "Shubert01", 2, -10.0, 10.0,
	-186.7309088310240, &calcShubert01 };

static double calcShubert03( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int j;
	int i;

	for( i = 0; i < N; i++ )
	{
		s2 = 0.0;

		for( j = 1; j <= 5; j++ )
		{
			s2 += j*sin((j+1)*x[i]+j);
		}

		s1 += s2;
	}

	return( s1 );
}

static const CTestFn TestFnShubert03 = { "Shubert03", 2, -10.0, 10.0,
	-29.6759000514212, &calcShubert03 };

static double calcShubert04( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int j;
	int i;

	for( i = 0; i < N; i++ )
	{
		s2 = 0.0;

		for( j = 1; j <= 5; j++ )
		{
			s2 += j*cos((j+1)*x[i]+j);
		}

		s1 += s2;
	}

	return( s1 );
}

static const CTestFn TestFnShubert04 = { "Shubert04", 2, -10.0, 10.0,
	-25.7417709954514, &calcShubert04 };

static double calcWeierstrass( const double* const x, const int N )
{
	const double a = 0.5;
	const double b = 3.0;
	const int kmax = 20;
	double ak[ kmax + 1 ];
	double bk[ kmax + 1 ];
	double s = 0.0;
	double s2 = 0.0;
	int i;
	int k;

	for( k = 0; k <= kmax; k++ )
	{
		ak[ k ] = pow(a,k);
		bk[ k ] = pow(b,k);
		s2 += ak[k]*cos(M_PI*bk[k]);
	}

	s2 *= N;

	for( i = 0; i < N; i++ )
	{
		double s1 = 0.0;

		for( k = 0; k <= kmax; k++ )
		{
			s1 += ak[k]*cos(2.0*M_PI*bk[k]*(x[i]+0.5));
		}

		s += s1-s2;
	}

	return( s );
}

static double calcWeierstrass_p( double* const minv, double* const maxv,
	const int N )
{
	double x[ N ];
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -0.5;
		maxv[ i ] = 0.5;
		x[ i ] = 0.0;
	}

	return( calcWeierstrass( x, N ));
}

static const CTestFn TestFnWeierstrass = { "Weierstrass", 0, 0.0, 0.0, 0.0,
	&calcWeierstrass, &calcWeierstrass_p };

static double calcFreudensteinRoth( const double* const x, const int N )
{
	return( sqr(x[0]-13.0+((5.0-x[1])*x[1]-2.0)*x[1]) +
		sqr(x[0]-29.0+((x[1]+1.0)*x[1]-14.0)*x[1]) );
}

static const CTestFn TestFnFreudensteinRoth = { "FreudensteinRoth", 2,
	-10.0, 10.0, 0.0, &calcFreudensteinRoth };

static double calcTrigonometric02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		const double sp = sqr(x[i]-0.9);
		s += 8.0*sqr(sin(7.0*sp))+6.0*sqr(sin(14.0*sp))+sp;
	}

	return( 1.0+s );
}

static const CTestFn TestFnTrigonometric02 = { "Trigonometric02", 0,
	-500.0, 500.0, 1.0, &calcTrigonometric02 };

static double calcBird( const double* const x, const int N )
{
	return( sqr(x[0]-x[1])+exp(sqr(1.0-sin(x[0])))*cos(x[1])+
		exp(sqr(1.0-cos(x[1])))*sin(x[0]) );
}

static const CTestFn TestFnBird = { "Bird", 2, -2.0 * M_PI, 2.0 * M_PI,
	-106.7645367492648, &calcBird };

static double calcTreccani( const double* const x, const int N )
{
	return( sqr(sqr(x[0]))+4.0*sqr(x[0])*x[0]+4.0*sqr(x[0])+sqr(x[1]) );
}

static const CTestFn TestFnTreccani = { "Treccani", 2, -5.0, 5.0, 0.0,
	&calcTreccani };

static double calcXinSheYang02( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += fabs(x[i]);
		s2 += sin(sqr(x[i]));
	}

	return( s1*exp(-s2) );
}

static const CTestFn TestFnXinSheYang02 = { "XinSheYang02", 0,
	-2.0 * M_PI, 2.0 * M_PI, 0.0, &calcXinSheYang02 };

static double calcXinSheYang03( const double* const x, const int N )
{
	const double m = 5.0;
	const double Beta = 15.0;
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += pow(x[i]/Beta,2.0*m);
		s2 += sqr(x[i]);
		s3 *= sqr(cos(x[i]));
	}

	return( exp(-s1)-2.0*exp(-s2)*s3 );
}

static const CTestFn TestFnXinSheYang03 = { "XinSheYang03", 0,
	-2.0 * M_PI, 2.0 * M_PI, -1.0, &calcXinSheYang03 };

static double calcXinSheYang04( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( sin( x[i]));
		s2 += sqr( x[i]);
		s3 += sqr( sin( sqrt( fabs( x[i]))));
	}

	return( (s1-exp(-s2))*exp(-s3) );
}

static const CTestFn TestFnXinSheYang04 = { "XinSheYang04", 0, -10.0, 10.0,
	-1.0, &calcXinSheYang04 };

static double calcBiggsEXP2( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i <= 9; i++ )
	{
		s += sqr(exp(-i*x[0]/10.0)-5.0*exp(-i*x[1]/10.0)-exp(-i/10.0)+
			5.0*exp(-(double)i));
	}

	return( s );
}

static const CTestFn TestFnBiggsEXP2 = { "BiggsEXP2", 2, 0.0, 20.0, 0.0,
	&calcBiggsEXP2 };

static double calcSchwefel06( const double* const x, const int N )
{
	double v1 = fabs(x[0]+2.0*x[1]-7.0);
	double v2 = fabs(2.0*x[0]+x[1]-5.0);

	return( v1 > v2 ? v1 : v2 );
}

static const CTestFn TestFnSchwefel06 = { "Schwefel06", 2, -100.0, 100.0, 0.0,
	&calcSchwefel06 };

static double calcChichinadze( const double* const x, const int N )
{
	return( sqr(x[0])-12.0*x[0]+8.0*sin(5.0/2.0*M_PI*x[0])+
		10.0*cos(0.5*M_PI*x[0])+11.0-0.2*sqrt(5.0)/exp(0.5*sqr(x[1]-0.5)) );
}

static const CTestFn TestFnChichinadze = { "Chichinadze", 2, -30.0, 30.0,
	-42.9443870189910, &calcChichinadze };

static double calcEggHolder( const double* const x, const int N )
{
	return( -x[0]*sin(sqrt(fabs(x[0]-x[1]-47.0)))-
		(x[1]+47.0)*sin(sqrt(fabs(0.5*x[0]+x[1]+47.0))) );
}

static const CTestFn TestFnEggHolder = { "EggHolder", 2, -512.0, 512.0,
	-959.6406627208510, &calcEggHolder };

static double calcHolderTable1( const double* const x, const int N )
{
	return( -fabs(exp(fabs(1.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI))*cos(x[0])*
		cos(x[1])) );
}

static const CTestFn TestFnHolderTable1 = { "HolderTable1", 2, -10.0, 10.0,
	-26.9203355555189, &calcHolderTable1 };

static double calcHolderTable2( const double* const x, const int N )
{
	return( -fabs(exp(fabs(1.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI))*sin(x[0])*
		cos(x[1])) );
}

static const CTestFn TestFnHolderTable2 = { "HolderTable2", 2, -10.0, 10.0,
	-19.2085025678868, &calcHolderTable2 };

static double calcComplexHolderTable( const double* const x, const int N )
{
	double x0 = x[ 0 ];
	double x1 = x[ 1 ];
	double sign = 1;
	double j;

	for( j = -4; j < 9; j += 0.5 )
	{
		if( j < x0 && x0 < j+0.5 )
		{
			x0 += sign*0.25;
		}

		sign *= -1;
	}

	return( -(fabs(sin(x0)*cos(x1)*exp(fabs(1-sqrt(x0*x0+x1*x1)/M_PI)))-
		(x0+x1)/10-sin(x0*10)*cos(x1*10)) );
}

static const CTestFn TestFnComplexHolderTable = { "ComplexHolderTable", 2,
	-10.0, 10.0, -21.9210396635469, &calcComplexHolderTable };

static double calcPowellSum( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow(fabs(x[i]),i+2.0);
	}

	return( s );
}

static const CTestFn TestFnPowellSum = { "PowellSum", 0, -1.0, 1.0, 0.0,
	&calcPowellSum };

static double calcPrice01( const double* const x, const int N )
{
	return( sqr(fabs(x[0])-5.0)+sqr(fabs(x[1])-5.0) );
}

static const CTestFn TestFnPrice01 = { "Price01", 2, -500.0, 500.0, 0.0,
	&calcPrice01 };

static double calcPrice03( const double* const x, const int N )
{
	return( 100.0*sqr(x[1]-sqr(x[0]))+sqr(6.4*sqr(x[1]-0.5)-x[0]-0.6) );
}

static const CTestFn TestFnPrice03 = { "Price03", 2, -50.0, 50.0, 0.0,
	&calcPrice03 };

static double calcPrice04( const double* const x, const int N )
{
	return( sqr(2.0*pow(x[0],3.0)*x[1]-pow(x[1],3.0))+
		sqr(6.0*x[0]-sqr(x[1])+x[1]) );
}

static const CTestFn TestFnPrice04 = { "Price04", 2, -50.0, 50.0, 0.0,
	&calcPrice04 };

static double calcBrown( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += pow(sqr(x[i]),sqr(x[i+1])+1.0) +
			pow(sqr(x[i+1]),sqr(x[i]+1.0));
	}

	return( s );
}

static const CTestFn TestFnBrown = { "Brown", 0, -1.0, 4.0, 0.0, &calcBrown };

static double calcBrent( const double* const x, const int N )
{
	return( sqr(x[0]+10.0)+sqr(x[1]+10.0)+exp(-sqr(x[0])-sqr(x[1])) );
}

static const CTestFn TestFnBrent = { "Brent", 2, -10.0, 10.0, 0.0,
	&calcBrent };

static double calcNewFunction01( const double* const x, const int N )
{
	return( sqrt(fabs(cos(sqrt(fabs(sqr(x[0])+x[1])))))+(x[0]+x[1])/100.0 );
}

static const CTestFn TestFnNewFunction01 = { "NewFunction01", 2, -10.0, 10.0,
	-0.17894509347721144, &calcNewFunction01 };

static double calcNewFunction02( const double* const x, const int N )
{
	return( sqrt(fabs(sin(sqrt(fabs(sqr(x[0])+x[1])))))+(x[0]+x[1])/100.0 );
}

static const CTestFn TestFnNewFunction02 = { "NewFunction02", 2, -10.0, 10.0,
	-0.1971881059905, &calcNewFunction02 };

static double calcLevy05( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 1; i <= 5; i++ )
	{
		s1 += i*cos((i-1)*x[0]+i);
		s2 += i*cos((i+1)*x[1]+i);
	}

	return( s1*s2+sqr(x[0]+1.42513)+sqr(x[1]+0.80032) );
}

static const CTestFn TestFnLevy05 = { "Levy05", 2, -10.0, 10.0,
	-176.1375780016295, &calcLevy05 };

static double calcDamavandi( const double* const x, const int N )
{
	return( (1.0-pow(fabs(sin(M_PI*(x[0]-2.0))*sin(M_PI*(x[1]-2.0))/
		(M_PI*M_PI*(x[0]-2.0)*(x[1]-2.0))),5.0))*
		(2.0+sqr(x[0]-7.0)+2.0*sqr(x[1]-7.0)) );
}

static const CTestFn TestFnDamavandi = { "Damavandi", 2, 0.0, 14.0, 0.0,
	&calcDamavandi };

static double calcPowerSum( const double* const x, const int N )
{
	const double b[ 4 ] = { 8.0, 18.0, 44.0, 114.0 };
	double s1 = 0.0;
	int i;
	int k;

	for( k = 0; k < 4; k++ )
	{
		double s2 = 0.0;

		for( i = 0; i < 4; i++ )
		{
			s2 += pow(x[i],k+1.0);
		}

		s1 += sqr(s2-b[ k ]);
	}

	return( s1 );
}

static const CTestFn TestFnPowerSum = { "PowerSum", 4, 0.0, 4.0, 0.0,
	&calcPowerSum };

static double calcPowell( const double* const x, const int N )
{
	return( sqr(x[2]+10.0*x[0])+5.0*sqr(x[1]-x[3])+sqr(sqr(x[0]-2.0*x[1]))+
		10.0*sqr(sqr(x[2]-x[3])) );
}

static const CTestFn TestFnPowell = { "Powell", 4, -4.0, 5.0, 0.0,
	&calcPowell };

static double calcPowellBadlyScaled( const double* const x, const int N )
{
	return( sqr(10000.0*x[0]*x[1]-1.0)+sqr(exp(-x[0])+exp(-x[1])-1.0001) );
}

static const CTestFn TestFnPowellBadlyScaled = { "PowellBadlyScaled", 2,
	-10.0, 10.0, 0.0, &calcPowellBadlyScaled };

static double calcPaviani( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 1.0;
	int i;

	for( i = 0; i < 10; i++ )
	{
		s1 += sqr(log(10.0-x[i]))+sqr(log(x[i]-2.0));
		s2 *= x[i];
	}

	return( s1-pow(s2,0.2) );
}

static const CTestFn TestFnPaviani = { "Paviani", 10, 2.001, 9.999,
	-45.7784697074463, &calcPaviani };

static double calcDolan( const double* const x, const int N )
{
	return( (x[0]+1.7*x[1])*sin(x[0])-1.5*x[2]-
		0.1*x[3]*cos(x[4]+x[3]-x[0])+0.2*sqr(x[4])-x[1]-1.0 );
}

static const CTestFn TestFnDolan = { "Dolan", 5, -100.0, 100.0,
	-529.8714387324576, &calcDolan };

static double calcTrid6( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < 6; i++ )
	{
		s1 += sqr(x[i]-1.0);
	}

	for( i = 1; i < 6; i++ )
	{
		s2 += x[i]*x[i-1];
	}

	return( s1-s2 );
}

static const CTestFn TestFnTrid6 = { "Trid6", 6, -20.0, 20.0, -50.0,
	&calcTrid6 };

static double calcTrid10( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < 10; i++ )
	{
		s1 += sqr(x[i]-1.0);
	}

	for( i = 1; i < 10; i++ )
	{
		s2 += x[i]*x[i-1];
	}

	return( s1-s2 );
}

static const CTestFn TestFnTrid10 = { "Trid10", 10, -100.0, 100.0, -210.0,
	&calcTrid10 };

static double calcMieleCantrell( const double* const x, const int N )
{
	return( sqr(sqr(exp(-x[0])-x[1]))+100.0*pow(x[1]-x[2],6.0)+
		sqr(sqr(tan(x[2]-x[3])))+pow(x[0],8.0) );
}

static const CTestFn TestFnMieleCantrell = { "MieleCantrell", 4, -1.0, 1.0,
	0.0, &calcMieleCantrell };

static double calcColville( const double* const x, const int N )
{
	return( 100.0*sqr(x[0]-sqr(x[1]))+sqr(1.0-x[0])+90.0*sqr(x[3]-sqr(x[2]))+
		sqr(1.0-x[2])+10.1*(sqr(x[1]-1.0)+sqr(x[3]-1.0))+
		19.8*(x[1]-1.0)*(x[3]-1.0) );
}

static const CTestFn TestFnColville = { "Colville", 4, -10.0, 10.0, 0.0,
	&calcColville };

static double calcWolfe( const double* const x, const int N )
{
	return( 4.0/3.0*pow(sqr(x[0])+sqr(x[1])-x[0]*x[1],0.75)+x[2] );
}

static const CTestFn TestFnWolfe = { "Wolfe", 3, 0.0, 2.0, 0.0, &calcWolfe };

static double calcAlpine1( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(x[i]*sin(x[i])+0.1*x[i]);
	}

	return( s );
}

static const CTestFn TestFnAlpine1 = { "Alpine1", 0, -10.0, 10.0, 0.0,
	&calcAlpine1 };

static double calcBartelsConn( const double* const x, const int N )
{
	return( fabs(sqr(x[0])+sqr(x[1])+x[0]*x[1])+
		fabs(sin(x[0]))+fabs(cos(x[1])) );
}

static const CTestFn TestFnBartelsConn = { "BartelsConn", 2, -500.0, 500.0,
	1.0, &calcBartelsConn };

static double calcSixHumpCamel( const double* const x, const int N )
{
	return( (4.0-2.1*sqr(x[0])+pow(x[0],4.0)/3.0)*sqr(x[0])+x[0]*x[1]+
		(4.0*sqr(x[1])-4.0)*sqr(x[1]) );
}

static const CTestFn TestFnSixHumpCamel = { "SixHumpCamel", 2, -5.0, 5.0,
	-1.0316284534899, &calcSixHumpCamel };

static double calcChenV( const double* const x, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(sqr(x[0])+sqr(x[1])-1.0))-
		b/(sqr(b)+sqr(sqr(x[0])+sqr(x[1])-0.5))-
		b/(sqr(b)+sqr(sqr(x[0])-sqr(x[1]))) );
}

static const CTestFn TestFnChenV = { "ChenV", 2, -500.0, 500.0,
	-2000.003999984000, &calcChenV };

static double calcChenBird( const double* const x, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(x[0]-0.4*x[1]-0.1))-
		b/(sqr(b)+sqr(2.0*x[0]+x[1]-1.5)) );
}

static const CTestFn TestFnChenBird = { "ChenBird", 2, -500.0, 500.0, -2000.0,
	&calcChenBird };

static double calcChungReynolds( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( s );
}

static const CTestFn TestFnChungReynolds = { "ChungReynolds", 0,
	-100.0, 100.0, 0.0, &calcChungReynolds };

static double calcCube( const double* const x, const int N )
{
	return( 100.0*pow(x[1]-pow(x[0],3.0),2.0)+pow(1.0-x[0],2.0) );
}

static const CTestFn TestFnCube = { "Cube", 2, -10.0, 10.0, 0.0, &calcCube };

static double calcDeckkersAarts( const double* const x, const int N )
{
	return( pow(10.0,5.0)*sqr(x[0])+sqr(x[1])-sqr(sqr(x[0])+sqr(x[1]))+
		pow(10.0,-5.0)*pow(sqr(x[0])+sqr(x[1]),4.0) );
}

static const CTestFn TestFnDeckkersAarts = { "DeckkersAarts", 2, -20.0, 20.0,
	-24776.5183423176930, &calcDeckkersAarts };

static double calcEggCrate( const double* const x, const int N )
{
	return( sqr(x[0])+sqr(x[1])+25.0*(sqr(sin(x[0]))+sqr(sin(x[1]))) );
}

static const CTestFn TestFnEggCrate = { "EggCrate", 2, -5.0, 5.0, 0.0,
	&calcEggCrate };

static double calcKeane( const double* const x, const int N )
{
	return( -sqr(sin(x[0]-x[1]))*sqr(sin(x[0]+x[1]))/
		sqrt(sqr(x[0])+sqr(x[1])) );
}

static const CTestFn TestFnKeane = { "Keane", 2, 0.00001, 10.0,
	-0.6736675209904, &calcKeane };

static double calcLeon( const double* const x, const int N )
{
	return( 100.0*sqr(x[1]-pow(x[0],3.0))+sqr(1.0-x[0]) );
}

static const CTestFn TestFnLeon = { "Leon", 2, 0.0, 10.0, 0.0, &calcLeon };

static double calcQing( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sqr(x[i])-(i+1.0));
	}

	return( s );
}

static const CTestFn TestFnQing = { "Qing", 0, -500.0, 500.0, 0.0,
	&calcQing };

static double calcSchwefel220( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(x[i]);
	}

	return( s );
}

static const CTestFn TestFnSchwefel220 = { "Schwefel220", 0, -100.0, 100.0,
	0.0, &calcSchwefel220 };

static double calcSchwefel221( const double* const x, const int N )
{
	double s = fabs(x[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		if( fabs(x[i]) > s )
		{
			s = fabs(x[i]);
		}
	}

	return( s );
}

static const CTestFn TestFnSchwefel221 = { "Schwefel221", 0, -100.0, 100.0,
	0.0, &calcSchwefel221 };

static double calcSchwefel222( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 1.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s1 += fabs(x[i]);
		s2 *= fabs(x[i]);
	}

	return( s1+s2 );
}

static const CTestFn TestFnSchwefel222 = { "Schwefel222", 0, -100.0, 100.0,
	0.0, &calcSchwefel222 };

static double calcSchwefel225( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s += sqr(x[i]-1.0)+sqr(x[0]-sqr(x[i]));
	}

	return( s );
}

static const CTestFn TestFnSchwefel225 = { "Schwefel225", 0, 0.0, 10.0,
	0.0, &calcSchwefel225 };

static double calcTestTubeHolder( const double* const x, const int N )
{
	return( -4.0*fabs(exp(fabs(cos(sqr(x[0])/200.0+sqr(x[1])/200.0)))*
		sin(x[0])*cos(x[1])) );
}

static const CTestFn TestFnTestTubeHolder = { "TestTubeHolder", 2,
	-10.0, 10.0, -10.8723001056227, &calcTestTubeHolder };

static double calcWayburnSeader02( const double* const x, const int N )
{
	return( sqr(1.613-4.0*sqr(x[0]-0.3125)-4.0*sqr(x[1]-1.625))+
		sqr(x[1]-1.0) );
}

static const CTestFn TestFnWayburnSeader02 = { "WayburnSeader02", 2,
	-500.0, 500.0, 0.0, &calcWayburnSeader02 };

static double calcSawtoothxy( const double* const x, const int N )
{
	const double t = atan2(x[1],x[0]);
	const double r = sqrt(sqr(x[0])+sqr(x[1]));
	const double ht = 0.5*cos(2*t-0.5)+cos(t)+2.0;
	const double gr = (sin(r)-sin(2.0*r)/2.0+sin(3.0*r)/3.0-sin(4.0*r)/4.0+
		4.0)*(r*r/(r+1.0));

	return( gr*ht );
}

static const CTestFn TestFnSawtoothxy = { "Sawtoothxy", 2, -20.0, 20.0, 0.0,
	&calcSawtoothxy };

static double calcMishra04( const double* const x, const int N )
{
	return( sqrt(fabs(sin(sqrt(fabs(sqr(x[0])+x[1])))))+(x[0]+x[1])/100.0 );
}

static const CTestFn TestFnMishra04 = { "Mishra04", 2, -10.00, 10.0,
	-0.1994114689025, &calcMishra04 };

static double calcZeroSum( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += x[i];
	}

	return( s == 0.0 ? 0.0 : 1.0+pow(10000.0*fabs(s),0.5) );
}

static const CTestFn TestFnZeroSum = { "ZeroSum", 0, -10.0, 10.0, 0.0,
	&calcZeroSum };

static double calcSchwefel( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += x[i]*sin(sqrt(fabs(x[i])));
	}

	return( 418.9828872724339*N-s );
}

static const CTestFn TestFnSchwefel = { "Schwefel", 0, -500.0, 500.0,
	0.0, &calcSchwefel };

static double calcAdjiman( const double* const x, const int N )
{
	return( cos(x[0])*sin(x[1])-x[0]/(sqr(x[1])+1.0) );
}

static double calcAdjiman_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -1.0;
	maxv[ 0 ] = 2.0;
	minv[ 1 ] = -1.0;
	maxv[ 1 ] = 1.0;
	return( -2.0218067833598 );
}

static const CTestFn TestFnAdjiman = { "Adjiman", 2, 0.0, 0.0, 0.0,
	&calcAdjiman, &calcAdjiman_p };

static double calcAlpine2( const double* const x, const int N )
{
	double s = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s *= sqrt(x[i])*sin(x[i]);
	}

	return( -s );
}

static double calcAlpine2_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = 0.0;
		maxv[ i ] = 10.0;
	}

	return( -pow( 2.808, (double) N ));
}

static const CTestFn TestFnAlpine2 = { "Alpine2", 0, 0.0, 0.0, 0.0,
	&calcAlpine2, &calcAlpine2_p };

static double calcBukin2( const double* const x, const int N )
{
	return( 100.0*sqr(x[1]-0.01*sqr(x[0])+1.0)+0.01*sqr(x[0]+10.0) );
}

static double calcBukin2_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -15.0;
	maxv[ 0 ] = -5.0;
	minv[ 1 ] = -3.0;
	maxv[ 1 ] = 3.0;
	return( 0.0 );
}

static const CTestFn TestFnBukin2 = { "Bukin2", 2, 0.0, 0.0, 0.0, &calcBukin2,
	&calcBukin2_p };

static double calcBukin4( const double* const x, const int N )
{
	return( 100.0*sqr(x[1])+0.01*fabs(x[0]+10.0) );
}

static double calcBukin4_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -15.0;
	maxv[ 0 ] = -5.0;
	minv[ 1 ] = -3.0;
	maxv[ 1 ] = 3.0;
	return( 0.0 );
}

static const CTestFn TestFnBukin4 = { "Bukin4", 2, 0.0, 0.0, 0.0, &calcBukin4,
	&calcBukin4_p };

static double calcBukin6( const double* const x, const int N )
{
	return( 100.0*sqrt(fabs(x[1]-0.01*sqr(x[0])))+0.01*fabs(x[0]+10.0) );
}

static double calcBukin6_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -15.0;
	maxv[ 0 ] = -5.0;
	minv[ 1 ] = -3.0;
	maxv[ 1 ] = 3.0;
	return( 0.0 );
}

static const CTestFn TestFnBukin6 = { "Bukin6", 2, 0.0, 0.0, 0.0, &calcBukin6,
	&calcBukin6_p };

static double calcStyblinskiTank( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sqr(x[i]))-16.0*sqr(x[i])+5.0*x[i];
	}

	return( 0.5*s );
}

static double calcStyblinskiTank_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -10.0;
		maxv[ i ] = 10.0;
	}

	return( -39.1661657037714*N );
}

static const CTestFn TestFnStyblinskiTank = { "StyblinskiTank", 0, 0.0, 0.0,
	0.0, &calcStyblinskiTank, &calcStyblinskiTank_p };

static double calcMcCormick( const double* const x, const int N )
{
	return( sin(x[0]+x[1])+sqr(x[0]-x[1])-1.5*x[0]+2.5*x[1]+1.0 );
}

static double calcMcCormick_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -1.5;
	maxv[ 0 ] = 4.0;
	minv[ 1 ] = -3.0;
	maxv[ 1 ] = 3.0;
	return( -1.9132229549810 );
}

static const CTestFn TestFnMcCormick = { "McCormick", 2, 0.0, 0.0, 0.0,
	&calcMcCormick, &calcMcCormick_p };

static double calcHimmelblau( const double* const x, const int N )
{
	return( sqr(sqr(x[0])+x[1]-11.0)+sqr(x[0]+sqr(x[1])-7.0) );
}

static const CTestFn TestFnHimmelblau = { "Himmelblau", 2, -6.0, 6.0, 0.0,
	&calcHimmelblau };

static double calcMichalewicz( const double* const x, const int N )
{
	const double m = 10.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sin(x[i])*pow(sin((i+1)*sqr(x[i])/M_PI),2.0*m);
	}

	return( -s );
}

static const CTestFn TestFnMichalewicz = { "Michalewicz", 2, 0.0, M_PI,
	-1.8013034100986, &calcMichalewicz };

static double calcBoxBettsExpQuadSum( const double* const x, const int N )
{
	double s = 0.0;
	int j;

	for( j = 1; j <= 10; j++ )
	{
		double g = exp(-0.1*j*x[0])-exp(-0.1*j*x[1])-
			(exp(-0.1*j)-exp(-(double)j))*x[2];

		s += g*g;
	}

	return( s );
}

static double calcBoxBettsExpQuadSum_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = 0.9;
	maxv[ 0 ] = 1.2;
	minv[ 1 ] = 9.0;
	maxv[ 1 ] = 11.2;
	minv[ 2 ] = 0.9;
	maxv[ 2 ] = 1.2;
	return( 0.0 );
}

static const CTestFn TestFnBoxBettsExpQuadSum = { "BoxBettsExpQuadSum", 3,
	0.0, 0.0, 0.0, &calcBoxBettsExpQuadSum, &calcBoxBettsExpQuadSum_p };

static double calcPowellQuartic( const double* const x, const int N )
{
	return( sqr(x[0]+10.0*x[1])+5.0*sqr(x[2]-x[3])+sqr(sqr(x[1]-2.0*x[2]))+
		10.0*sqr(sqr(x[0]-x[3])) );
}

static const CTestFn TestFnPowellQuartic = { "PowellQuartic", 4, -10.0, 10.0,
	0.0, &calcPowellQuartic };

static double calcMishra05( const double* const x, const int N )
{
	const double f1 = sqr(sin(sqr(cos(x[0])+cos(x[1]))));
	const double f2 = sqr(cos(sqr(sin(x[0])+sin(x[1]))));

	return( sqr(f1+f2+x[0])+0.01*x[0]+0.1*x[1] );
}

static const CTestFn TestFnMishra05 = { "Mishra05", 2, -10.0, 10.0,
	-1.0198295199309, &calcMishra05 };

static double calcMishra06( const double* const x, const int N )
{
	const double f1 = sqr(sin(sqr(cos(x[0])+cos(x[1]))));
	const double f2 = sqr(cos(sqr(sin(x[0])+sin(x[1]))));
	const double f3 = 0.1*(sqr(x[0]-1.0)+sqr(x[1]-1.0));

	return( -log(sqr(f1-f2+x[0]))+f3 );
}

static const CTestFn TestFnMishra06 = { "Mishra06", 2, -10.0, 10.0,
	-2.2839498384748, &calcMishra06 };

static double calcMishra09( const double* const x, const int N )
{
	const double f1 = 2.0*pow(x[0],3.0)+5.0*x[0]*x[1]+4.0*x[2]-
		2.0*sqr(x[0])*x[2]-18.0;

	const double f2 = x[0]+pow(x[1],3.0)+x[0]*sqr(x[1])+x[0]*sqr(x[2])-22.0;
	const double f3 = 8.0*sqr(x[0])+2.0*x[1]*x[2]+2.0*sqr(x[1])+
		3.0*pow(x[1],3.0)-52.0;

	return( sqr(f1*sqr(f2)*f3+f1*f2*sqr(f3)+sqr(f2)+sqr(x[0]+x[1]-x[2])) );
}

static const CTestFn TestFnMishra09 = { "Mishra09", 3, -10.0, 10.0, 0.0,
	&calcMishra09 };

static double calcZirilli( const double* const x, const int N )
{
	return( 0.25*sqr(sqr(x[0]))-0.5*sqr(x[0])+0.1*x[0]+0.5*sqr(x[1]) );
}

static const CTestFn TestFnZirilli = { "Zirilli", 2, -10.0, 10.0,
	-0.3523860738000, &calcZirilli };

static double calcCamel( const double* const x, const int N )
{
	return( -(-sqr(sqr(x[0]))+4.5*sqr(x[0])+2.0)/exp(2.0*sqr(x[1])) );
}

static const CTestFn TestFnCamel = { "Camel", 2, -2.0, 2.0, -7.0625000000000,
	&calcCamel };

static double calcComplex( const double* const x, const int N )
{
	return( sqr(pow(x[0],3.0)-3.0*x[0]*sqr(x[1])-1.0)+
		sqr(3.0*x[1]*sqr(x[0])-pow(x[1],3.0)) );
}

static const CTestFn TestFnComplex = { "Complex", 2, -2.0, 2.0, 0.0,
	&calcComplex };

static double calcDavis( const double* const x, const int N )
{
	return( pow(sqr(x[0])+sqr(x[1]),0.25)*(sqr(sin(50.0*
		pow(3.0*sqr(x[0])+sqr(x[1]),0.1)))+1.0) );
}

static const CTestFn TestFnDavis = { "Davis", 2, -100.0, 100.0, 0.0,
	&calcDavis };

static double calcDownhillStep( const double* const x, const int N )
{
	return( floor(10.0*(10.0-exp(-sqr(x[0])-3.0*sqr(x[1]))))/10.0 );
}

static const CTestFn TestFnDownhillStep = { "DownhillStep", 2, -10.0, 10.0,
	9.0, &calcDownhillStep };

static double calcEngvall( const double* const x, const int N )
{
	return( sqr(sqr(x[0]))+sqr(sqr(x[1]))+2.0*sqr(x[0])*sqr(x[1])-4.0*x[0]+
		3.0 );
}

static const CTestFn TestFnEngvall = { "Engvall", 2, -2000.0, 2000.0, 0.0,
	&calcEngvall };

static double calcCrossLegTable( const double* const x, const int N )
{
	return( -1.0/pow(fabs(sin(x[0])*sin(x[1])*
		exp(fabs(100.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI)))+1.0,0.1));
}

static const CTestFn TestFnCrossLegTable = { "CrossLegTable", 2, -10.0, 10.0,
	-1.0, &calcCrossLegTable };

static double calcCrownedCross( const double* const x, const int N )
{
	return( 0.0001*pow(fabs(sin(x[0])*sin(x[1])*
		exp(fabs(100.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI)))+1.0,0.1) );
}

static const CTestFn TestFnCrownedCross = { "CrownedCross", 2, -10.0, 10.0,
	0.0001, &calcCrownedCross };

static double calcGramacyLee02( const double* const x, const int N )
{
	return( x[0]*exp(-sqr(x[0])-sqr(x[1])) );
}

static const CTestFn TestFnGramacyLee02 = { "GramacyLee02", 2, -1.5, 1.5,
	-0.4288819424804, &calcGramacyLee02 };

static double calcGramacyLee03( const double* const x, const int N )
{
	const double f1 = exp(-sqr(x[0]-1.0))+exp(-0.8*sqr(x[0]+1.0))-
		0.05*sin(8.0*x[0]+0.8);

	const double f2 = exp(-sqr(x[1]-1.0))+exp(-0.8*sqr(x[1]+1.0))-
		0.05*sin(8.0*x[1]+0.8);

	return( -f1*f2 );
}

static const CTestFn TestFnGramacyLee03 = { "GramacyLee03", 2, -1.5, 1.5,
	-1.1268717457863, &calcGramacyLee03 };

static double calcGiunta( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sin(1.0-16.0/15.0*x[i]))-1.0/50.0*sin(4.0-64.0/15.0*x[i])-
			sin(1.0-16.0/15.0*x[i]);
	}

	return( 0.6+s );
}

static const CTestFn TestFnGiunta = { "Giunta", 2, -1.0, 1.0, 0.0644704205369,
	&calcGiunta };

static double calcHosaki( const double* const x, const int N )
{
	return( (1.0-8.0*x[0]+7.0*sqr(x[0])-7.0/3.0*sqr(x[0])*x[0]+
		1.0/4.0*sqr(sqr(x[0])))*sqr(x[1])*exp(-x[1]) );
}

static const CTestFn TestFnHosaki = { "Hosaki", 2, 0.0, 10.0,
	-2.3458115761013, &calcHosaki };

static double calcKearfott( const double* const x, const int N )
{
	return( sqr(sqr(x[0])+sqr(x[1])-2.0)+sqr(sqr(x[0])-sqr(x[1])-1.0) );
}

static const CTestFn TestFnKearfott = { "Kearfott", 2, -3.0, 4.0, 0.0,
	&calcKearfott };

static double calcJennrichSampson( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= 10; i++ )
	{
		s += sqr(2.0+2*i-(exp(i*x[0])+exp(i*x[1])));
	}

	return( s );
}

static const CTestFn TestFnJennrichSampson = { "JennrichSampson", 2,
	-1.0, 1.0, 124.3621823556148, &calcJennrichSampson };

static double calcTsoulos( const double* const x, const int N )
{
	return( sqr(x[0])+sqr(x[1])-cos(18.0*x[0])-cos(18.0*x[1]) );
}

static const CTestFn TestFnTsoulos = { "Tsoulos", 2, -1.0, 1.0, -2.0,
	&calcTsoulos };

static double calcUrsemWaves( const double* const x, const int N )
{
	return( -pow(0.3*x[0],3.0)+(sqr(x[1])-4.5*sqr(x[1]))*x[0]*x[1]+
		4.7*cos(3.0*x[0]-sqr(x[1])*(2.0+x[0]))*sin(2.5*M_PI*x[0]) );
}

static double calcUrsemWaves_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -0.9;
	maxv[ 0 ] = 1.2;
	minv[ 1 ] = -1.2;
	maxv[ 1 ] = 1.2;
	return( -7.3069987313245 );
}

static const CTestFn TestFnUrsemWaves = { "UrsemWaves", 2, 0.0, 0.0, 0.0,
	&calcUrsemWaves, &calcUrsemWaves_p };

static double calcGulfResearchProblem( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= 99; i++ )
	{
		const double ti = i/100.0;
		const double ui = 25.0+pow(-50.0*log(ti),1.0/1.5);
		s += sqr(exp(-pow(ui-x[1],x[2])/x[0])-ti);
	}

	return( s );
}

static double calcGulfResearchProblem_p( double* const minv,
	double* const maxv, const int N )
{
	minv[ 0 ] = 0.1;
	maxv[ 0 ] = 100.0;
	minv[ 1 ] = 0.0;
	maxv[ 1 ] = 25.6;
	minv[ 2 ] = 0.0;
	maxv[ 2 ] = 5.0;
	return( 0.0 );
}

static const CTestFn TestFnGulfResearchProblem = { "GulfRsrchProblem", 3,
	0.0, 0.0, 0.0, &calcGulfResearchProblem, &calcGulfResearchProblem_p };

static double calcYaoLiu04( const double* const x, const int N )
{
	double m = fabs(x[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		const double v = fabs(x[i]);
		m = ( v > m ? v : m );
	}

	return( m );
}

static const CTestFn TestFnYaoLiu04 = { "YaoLiu04", 0, -10.0, 10.0, 0.0,
	&calcYaoLiu04 };

static double calcBentCigar( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( sqr(x[0])+1e6*s );
}

static const CTestFn TestFnBentCigar = { "BentCigar", 0, -100.0, 100.0, 0.0,
	&calcBentCigar };

static double calcDeflCorrSpring( const double* const x, const int N )
{
	const double a = 5.0;
	const double k = 5.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]-a);
	}

	return( 0.1*s-cos(k*sqrt(s)) );
}

static const CTestFn TestFnDeflCorrSpring = { "DeflCorrSpring", 0,
	0.0, 10.0, -1.0, &calcDeflCorrSpring };

static double calcHolzman( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += (i+1)*sqr(sqr(x[i]));
	}

	return( s );
}

static const CTestFn TestFnHolzman = { "Holzman", 0, -10.0, 10.0, 0.0,
	&calcHolzman };

static double calcHyperGrid( const double* const x, const int N )
{
	const double a = 6.0;
	const double c = 5.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow(sin(c*M_PI*x[i]),a);
	}

	return( -s/N );
}

static const CTestFn TestFnHyperGrid = { "HyperGrid", 0, 0.0, 1.0, -1.0,
	&calcHyperGrid };

static double calcQuintic( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(pow(x[i],5.0)-3.0*pow(x[i],4.0)+4.0*pow(x[i],3.0)+
			2.0*sqr(x[i])-10.0*x[i]-4.0);
	}

	return( s );
}

static const CTestFn TestFnQuintic = { "Quintic", 0, -10.0, 10.0, 0.0,
	&calcQuintic };

static double calcVincent( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sin(10.0*log(x[i]));
	}

	return( -s/N );
}

static const CTestFn TestFnVincent = { "Vincent", 0, 0.25, 10.0, -1.0,
	&calcVincent };

static double calcStep01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += floor(fabs(x[i]));
	}

	return( s );
}

static const CTestFn TestFnStep01 = { "Step01", 0, -100, 100.0, 0.0,
	&calcStep01 };

static double calcStep02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(floor(fabs(x[i]+0.5)));
	}

	return( s );
}

static const CTestFn TestFnStep02 = { "Step02", 0, -100, 100.0, 0.0,
	&calcStep02 };

static double calcStep03( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += floor(sqr(x[i]));
	}

	return( s );
}

static const CTestFn TestFnStep03 = { "Step03", 0, -100, 100.0, 0.0,
	&calcStep03 };

static double calcHelicalValley( const double* const x, const int N )
{
	double theta;

	if( x[ 0 ] >= 0.0 )
	{
		theta = (M_PI/2.0)/tan(x[1]/x[0]);
	}
	else
	{
		theta = (M_PI/2.0)/tan(x[1]/x[0]+M_PI);
	}

	return( 100.0*(sqr(x[2]-10.0*theta)+sqr(sqrt(sqr(x[0])+sqr(x[1]))-1.0))+
		sqr(x[2]) );
}

static const CTestFn TestFnHelicalValley = { "HelicalValley", 3,
	-10.0, 10.0, 0.0, &calcHelicalValley };

static double calcDixonPrice( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s += (i+1.0)*sqr(2.0*sqr(x[i])-x[i-1]);
	}

	return( sqr(x[0]-1.0)+s );
}

static const CTestFn TestFnDixonPrice = { "DixonPrice", 0, -10.0, 10.0, 0.0,
	&calcDixonPrice };

static double calcHartman3( const double* const x, const int N )
{
	static const double a[ 4 ][ 3 ] =
		{
			{ 3, 10, 30 },
			{ 0.1, 10, 35 },
			{ 3, 10, 30 },
			{ 0.1, 10, 35 }
		};

	static const double c[ 4 ] = { 1, 1.2, 3, 3.2 };
	static const double px[ 4 ][ 3 ] =
		{
			{ 0.3689, 0.1170, 0.2673 },
			{ 0.4699, 0.4387, 0.7470 },
			{ 0.1091, 0.8732, 0.5547 },
			{ 0.0381, 0.5743, 0.8828 }
		};

	double s = 0.0;
	int i;

	for( i = 0; i < 4; i++ )
	{
		double s2 = 0.0;
		int j;

		for( j = 0; j < N; j++ )
		{
			s2 += a[i][j]*sqr(x[j]-px[i][j]);
		}

		s += c[i]*exp(-s2);
	}

	return( -s );
}

static const CTestFn TestFnHartman3 = { "Hartman3", 3, 0.0, 1.0,
	-3.8627797873327, &calcHartman3 };

static double calcHartman6( const double* const x, const int N )
{
	static const double a[ 4 ][ 6 ] =
		{
			{ 10.0, 3.0, 17.0, 3.5, 1.7, 8.0 },
			{ 0.05, 10.0, 17.0, 0.1, 8.0, 14.0 },
			{ 3.0, 3.5, 1.7, 10.0, 17.0, 8.0 },
			{ 17.0, 8.0, 0.05, 10.0, 0.1 , 14.0 }
		};

	static const double c[ 4 ] = { 1, 1.2, 3, 3.2 };
	static const double px[ 4 ][ 6 ] =
		{
			{ 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886 },
			{ 0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991 },
			{ 0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650 },
			{ 0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381 }
		};

	double s = 0.0;
	int i;

	for( i = 0; i < 4; i++ )
	{
		double s2 = 0.0;
		int j;

		for( j = 0; j < N; j++ )
		{
			s2 += a[i][j]*sqr(x[j]-px[i][j]);
		}

		s += c[i]*exp(-s2);
	}

	return( -s );
}

static const CTestFn TestFnHartman6 = { "Hartman6", 6, 0.0, 1.0,
	-3.3223680114155, &calcHartman6 };

static double calcDeVilliersGlasser02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= 24; i++ )
	{
		const double ti = 0.1*(i-1);
		const double yi = 53.81*pow(1.27,ti)*tanh(3.012*ti+sin(2.13*ti))*
			cos(exp(0.507)*ti);

		s += sqr(x[0]*pow(x[1],ti)*tanh(x[2]*ti+sin(x[3]*ti))*
			cos(ti*exp(x[4]))-yi);
	}

	return( s );
}

static const CTestFn TestFnDeVilliersGlasser02 = { "DeVilliersGlasser02", 5,
	1.0, 60.0, 0.0, &calcDeVilliersGlasser02 };

static double calcBuecheRastrigin( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += cos(2.0*M_PI*x[i]);
		s2 += sqr(x[i]);
	}

	return( 10.0*(N-s1)+s2 );
}

static const CTestFn TestFnBuecheRastrigin = { "BuecheRastrigin", 0,
	-5.0, 5.0, 0.0, &calcBuecheRastrigin };

static double calcDifferentPowers( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++)
	{
		const double exponent = 2.0+(4.0*i)/(N-1);
		s += pow(fabs(x[i]),exponent);
	}

	return( sqrt(s) );
}

static const CTestFn TestFnDifferentPowers = { "DifferentPowers", 0,
	-5.0, 5.0, 0.0, &calcDifferentPowers };

static double calcDiscus( const double* const x, const int N )
{
	static const double condition = 1.0e6;
	double s = condition * sqr(x[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( s );
}

static const CTestFn TestFnDiscus = { "Discus", 0, -5.0, 5.0, 0.0,
	&calcDiscus };

static double calcEllipsoid( const double* const x, const int N )
{
	static const double condition = 1.0e6;
	double s = sqr(x[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		const double exponent = 1.0*i/(N-1);
		s += pow(condition,exponent)*sqr(x[i]);
	}

	return( s );
}

static const CTestFn TestFnEllipsoid = { "Ellipsoid", 0, -5.0, 5.0, 0.0,
	&calcEllipsoid };

static double calcGriewankRosenbrock( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		const double c1 = sqr(x[i])-x[i+1];
		const double c2 = 1.0-x[i];
		const double tmp = 100.0*c1*c1+c2*c2;
		s += tmp/4000.0-cos(tmp);
	}

	return( 10.0+10.0*s/(N-1));
}

static const CTestFn TestFnGriewankRosenbrock = { "GriewankRosenbrock", 0,
	-5.0, 5.0, 1.0, &calcGriewankRosenbrock };

static double calcSchaffer07( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		const double tmp = sqr(x[i])+sqr(x[i+1]);
		s += pow(tmp,0.25)*(1.0+pow(sin(50.0*pow(tmp,0.1)),2.0));
	}

	return( sqr(s/(N-1)) );
}

static const CTestFn TestFnSchaffer07 = { "Schaffer07", 0, -5.0, 5.0, 0.0,
	&calcSchaffer07 };

static double calcKatsuura( const double* const x, const int N )
{
	double s = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		double tmp = 0.0;
		int j;

		for( j = 1; j <= 32; j++ )
		{
			const double tmp2 = pow(2.0,(double) j);
			tmp += roundf(tmp2*x[i])/tmp2;
		}

		s *= 1.0+(i+1)*tmp;
	}

	return( s );
}

static const CTestFn TestFnKatsuura = { "Katsuura", 0, 0.0, 100.0, 1.0,
	&calcKatsuura };

static double calcRotatedEllipse01( const double* const x, const int N )
{
	return( 7.0*sqr(x[0])-6.0*sqrt(3.0)*x[0]*x[1]+13.0*sqr(x[1]) );
}

static const CTestFn TestFnRotatedEllipse01 = { "RotatedEllipse01", 2,
	-500.0, 500.0, 0.0, &calcRotatedEllipse01 };

static double calcRotatedEllipse02( const double* const x, const int N )
{
	return( sqr(x[0])-x[0]*x[1]+sqr(x[1]) );
}

static const CTestFn TestFnRotatedEllipse02 = { "RotatedEllipse02", 2,
	-500.0, 500.0, 0.0, &calcRotatedEllipse02 };

static double calcTrigonometric01( const double* const x, const int N )
{
	double s0 = 0.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s0 += cos(x[i]);
	}

	for( i = 0; i < N; i++ )
	{
		s += sqr(N-s0+(i+1)*(1.0-cos(x[i])-sin(x[i])));
	}

	return( s );
}

static const CTestFn TestFnTrigonometric01 = { "Trigonometric01", 0,
	0.0, M_PI, 0.0, &calcTrigonometric01 };

static double calcExponential( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( -exp(-0.5*s) );
}

static const CTestFn TestFnExponential = { "Exponential", 0,
	-1.0, 1.0, -1.0, &calcExponential };

static double calcUrsem01( const double* const x, const int N )
{
	return( -sin(2.0*x[0]-0.5*M_PI)-3.0*cos(x[1])-0.5*x[0] );
}

static double calcUrsem01_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -2.5;
	maxv[ 0 ] = 3.0;
	minv[ 1 ] = -2.0;
	maxv[ 1 ] = 2.0;
	return( -4.8168140637348 );
}

static const CTestFn TestFnUrsem01 = { "Ursem01", 2, 0.0, 0.0, 0.0,
	&calcUrsem01, &calcUrsem01_p };

static double calcQuadratic( const double* const x, const int N )
{
	return( -3803.84-138.08*x[0]-232.92*x[1]+128.08*sqr(x[0])+203.64*sqr(x[1])+
		182.25*x[0]*x[1] );
}

static const CTestFn TestFnQuadratic = { "Quadratic", 2,
	-10.0, 10.0, -3873.7241821862713, &calcQuadratic };

static double calcSchwefel01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]);
	}

	return( pow(s,sqrt(M_PI)));
}

static const CTestFn TestFnSchwefel01 = { "Schwefel01", 0, -100.0, 100.0, 0.0,
	&calcSchwefel01 };

static double calcSchwefel02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= N; i++ )
	{
		double s2 = 0.0;
		int j;

		for( j = 1; j <= i; j++ )
		{
			s2 += x[i-1];
		}

		s += sqr(s2);
	}

	return( s );
}

static const CTestFn TestFnSchwefel02 = { "Schwefel02", 0, -100.0, 100.0, 0.0,
	&calcSchwefel02 };

static double calcSchwefel04( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i]-1.0)+sqr(x[0]-sqr(x[i]));
	}

	return( s );
}

static const CTestFn TestFnSchwefel04 = { "Schwefel04", 0, 0.0, 10.0, 0.0,
	&calcSchwefel04 };

static double calcSchwefel236( const double* const x, const int N )
{
	return( -x[0]*x[1]*(72.0-2.0*x[0]-2.0*x[1]) );
}

static const CTestFn TestFnSchwefel236 = { "Schwefel236", 2, 0.0, 500.0,
	-3456.0, &calcSchwefel236 };

static double calcShekel_inner( const double* const x, const int N,
	const int m )
{
	static const double a[ 10 ][ 4 ] =
	{
		{ 4.0, 4.0, 4.0, 4.0 },
		{ 1.0, 1.0, 1.0, 1.0 },
		{ 8.0, 8.0, 8.0, 8.0 },
		{ 6.0, 6.0, 6.0, 6.0 },
		{ 3.0, 7.0, 3.0, 7.0 },
		{ 2.0, 9.0, 2.0, 9.0 },
		{ 5.0, 3.0, 5.0, 3.0 },
		{ 8.0, 1.0, 8.0, 1.0 },
		{ 6.0, 2.0, 6.0, 2.0 },
		{ 7.0, 3.0, 7.0, 3.0 }
	};

	static const double c[ 10 ] = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3,
		0.7, 0.5, 0.5 };

	double s = 0.0;
	int i;

	for( i = 0; i < m; i++ )
	{
		double s2 = 0.0;
		int j;

		for( j = 0; j < N; j++ )
		{
			s2 += sqr(x[j]-a[i][j]);
		}

		s += 1.0/(c[i]+s2);
	}

	return( -s );
}

static double calcShekel05( const double* const x, const int N )
{
	return( calcShekel_inner( x, N, 5 ));
}

static const CTestFn TestFnShekel05 = { "Shekel05", 4, 0.0, 10.0,
	-10.1531996790582, &calcShekel05 };

static double calcShekel07( const double* const x, const int N )
{
	return( calcShekel_inner( x, N, 7 ));
}

static const CTestFn TestFnShekel07 = { "Shekel07", 4, 0.0, 10.0,
	-10.4029153367777, &calcShekel07 };

static double calcShekel10( const double* const x, const int N )
{
	return( calcShekel_inner( x, N, 10 ));
}

static const CTestFn TestFnShekel10 = { "Shekel10", 4, 0.0, 10.0,
	-10.5320872211865, &calcShekel10 };

static double calcMishra01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += x[i];
	}

	s = N-s;

	return( pow(1.0+s,s) );
}

static const CTestFn TestFnMishra01 = { "Mishra01", 2, 0.0, 1.0, 2.0,
	&calcMishra01 };

static double calcMishra02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += (x[i]+x[i+1])/2.0;
	}

	s = N-s;

	return( pow(1.0+s,s) );
}

static const CTestFn TestFnMishra02 = { "Mishra02", 2, 0.0, 1.0, 2.0,
	&calcMishra02 };

static double calcMishra07( const double* const x, const int N )
{
	double s = 1.0;
	int nf = 1;
	int i;

	for( i = 0; i < N; i++ )
	{
		s *= x[i];
		nf *= i+1;
	}

	return( sqr(s-nf) );
}

static const CTestFn TestFnMishra07 = { "Mishra07", 0, -10.0, 10.0, 0.0,
	&calcMishra07 };

static double calcZettl( const double* const x, const int N )
{
	return( 1.0/4.0*x[0]+sqr(sqr(x[0])-2.0*x[0]+sqr(x[1])) );
}

static const CTestFn TestFnZettl = { "Zettl", 2, -1.0, 5.0,
	-0.003791237220468656, &calcZettl };

static double calcMultiModal( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += fabs(x[i]);
		s2 *= fabs(x[i]);
	}

	return( s1*s2 );
}

static const CTestFn TestFnMultiModal = { "MultiModal", 0, -10.0, 10.0, 0.0,
	&calcMultiModal };

static double calcParsopoulos( const double* const x, const int N )
{
	return( sqr(cos(x[0]))+sqr(sin(x[1])) );
}

static const CTestFn TestFnParsopoulos = { "Parsopoulos", 2, -5.0, 5.0, 0.0,
	&calcParsopoulos };

static double calcDeb01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow(sin(5.0*M_PI*x[i]),6.0);
	}

	return( -s/N );
}

static const CTestFn TestFnDeb01 = { "Deb01", 0, -1.0, 1.0, -1.0, &calcDeb01 };

static double calcDeb02( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow(sin(5.0*M_PI*(pow(x[i],3.0/4.0)-0.05)),6.0);
	}

	return( -s/N );
}

static const CTestFn TestFnDeb02 = { "Deb02", 0, 0.0, 1.0, -1.0, &calcDeb02 };

static double calcCarromTable( const double* const x, const int N )
{
	return( -1.0/30.0*exp(2.0*fabs(1.0-sqrt(sqr(x[0])+sqr(x[1]))/M_PI))*
		sqr(cos(x[0]))*sqr(cos(x[1])) );
}

static const CTestFn TestFnCarromTable = { "CarromTable", 2, -10.0, 10.0,
	-24.1568155473913, &calcCarromTable };

static double calcNewFunction03( const double* const x, const int N )
{
	return( 0.01*x[0]+0.1*x[1]+sqr(x[0]+sqr(sin(sqr(cos(x[0])+cos(x[1]))))+
		sqr(cos(sqr(sin(x[0])+sin(x[1]))))) );
}

static const CTestFn TestFnNewFunction03 = { "NewFunction03", 2, -10.0, 10.0,
	-1.0198295199309, &calcNewFunction03 };

static double calcLevy03( const double* const x, const int N )
{
	double y[ N ];
	int i;

	for( i = 0; i < N; i++ )
	{
		y[ i ] = 1.0+(x[i]-1.0)/4.0;
	}

	double s = 0.0;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(y[i]-1.0)*(1.0+10.0*sqr(sin(M_PI*y[i+1])))+sqr(y[N-1]-1.0);
	}

	return( sqr(sin(M_PI*y[0]))+s );
}

static const CTestFn TestFnLevy03 = { "Levy03", 0, -10.0, 10.0, 0.0,
	&calcLevy03 };

static double calcGear( const double* const x, const int N )
{
	return( sqr(1.0/6.931-floor(x[0])*floor(x[1])/(floor(x[2])*floor(x[3]))) );
}

static const CTestFn TestFnGear = { "Gear", 4, 12.0, 60.0,
	2.7e-12, &calcGear };

static double calcStretchedV( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		const double t = sqr(x[i+1])+sqr(x[i]);
		s += pow(t,1.0/4.0)*sqr(sin(50.0*pow(t,0.1))+1.0);
	}

	return( s );
}

static const CTestFn TestFnStretchedV = { "StretchedV", 0, -10.0, 10.0, 0.0,
	&calcStretchedV };

static double calcUrsem04( const double* const x, const int N )
{
	return( -3.0*sin(0.5*M_PI*x[0]+0.5*M_PI)*
		(2.0-sqrt(sqr(x[0])+sqr(x[1])))/4.0 );
}

static const CTestFn TestFnUrsem04 = { "Ursem04", 2, -2.0, 2.0,
	-1.5, &calcUrsem04 };

static double calcJudge( const double* const x, const int N )
{
	static const double A[ 20 ] = { 4.284,4.149,3.877,0.533,2.211,2.389,2.145,
		3.231,1.998,1.379,2.106,1.428,1.011,2.179,2.858,1.388,1.651,1.593,
		1.046,2.152 };

	static const double B[ 20 ] = { 0.286,0.973,0.384,0.276,0.973,0.543,0.957,
		0.948,0.543,0.797,0.936,0.889,0.006,0.828,0.399,0.617,0.939,0.784,
		0.072,0.889 };

	static const double C[ 20 ] = { 0.645,0.585,0.310,0.058,0.455,0.779,0.259,
		0.202,0.028,0.099,0.142,0.296,0.175,0.180,0.842,0.039,0.103,0.620,
		0.158,0.704 };

	double s = 0.0;
	int i;

	for( i = 0; i < 20; i++ )
	{
		s += sqr(x[0]+B[i]*x[1]+C[i]*sqr(x[1])-A[i]);
	}

	return( s );
}

static const CTestFn TestFnJudge = { "Judge", 2, -10.0, 10.0,
	16.0817301329604, &calcJudge };

static double calcWayburnSeader01( const double* const x, const int N )
{
	return( sqr(pow(x[0],6.0)+pow(x[1],4.0)-17.0)+sqr(2.0*x[0]+x[1]-4.0) );
}

static const CTestFn TestFnWayburnSeader01 = { "WayburnSeader01", 2,
	-5.0, 5.0, 0.0, &calcWayburnSeader01 };

static double calcWayburnSeader03( const double* const x, const int N )
{
	return( 2.0/3.0*sqr(x[0])*x[0]-8.0*sqr(x[0])+33.0*x[0]-x[0]*x[1]+5.0+
		sqr(sqr(x[0]-4.0)+sqr(x[1]-5.0)-4.0) );
}

static const CTestFn TestFnWayburnSeader03 = { "WayburnSeader03", 2,
	-500.0, 500.0, 19.1058797945680, &calcWayburnSeader03 };

static double calcVenterSobiSobieski( const double* const x, const int N )
{
	return( sqr(x[0])-100.0*sqr(cos(x[0]))-100.0*cos(sqr(x[0])/30.0)+sqr(x[1])-
		100.0*sqr(cos(x[1]))-100.0*cos(sqr(x[1])/30.0) );
}

static const CTestFn TestFnVenterSobiSobieski = { "VenterSobiSobieski", 2,
	-50.0, 50.0, -400.0, &calcVenterSobiSobieski };

static double calcElAttarVidyasDutta( const double* const x, const int N )
{
	return( sqr(sqr(x[0])+x[1]-10.0)+sqr(x[0]+sqr(x[1])-7.0)+
		sqr(sqr(x[0])+pow(x[1],3.0)-1.0) );
}

static const CTestFn TestFnElAttarVidyasDutta = { "ElAttarVidyasDutta", 2,
	-100.0, 100.0, 1.7127803548622, &calcElAttarVidyasDutta };

static double calcModifiedRosenbrock( const double* const x, const int N )
{
	return( 74.0+100.0*sqr(x[1]-sqr(x[0]))+sqr(1.0-x[0])-
		400.0*exp(-(sqr(x[0]+1.0)+sqr(x[1]+1.0))/0.1) );
}

static const CTestFn TestFnModifiedRosenbrock = { "ModifiedRosenbrock", 2,
	-2.0, 2.0, 34.0402431066405, &calcModifiedRosenbrock };

static double calcStochastic( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += rnd.getRndValue()*fabs(x[i]-1.0/(i+1));
	}

	return( s );
}

static const CTestFn TestFnStochastic = { "Stochastic", 0, -5.0, 5.0, 0.0,
	&calcStochastic };

static double calcXinSheYang01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += rnd.getRndValue()*pow(fabs(x[i]),i+1);
	}

	return( s );
}

static const CTestFn TestFnXinSheYang01 = { "XinSheYang01", 0, -5.0, 5.0, 0.0,
	&calcXinSheYang01 };

static double calcKowalik( const double* const x, const int N )
{
	static const double a[ 11 ] = { 4.0, 2.0, 1.0, 1.0/2, 1.0/4, 1.0/6,
		1.0/8, 1.0/10, 1.0/12, 1.0/14, 1.0/16 };

	static const double b[ 11 ] = { 0.1957, 0.1947, 0.1735, 0.1600, 0.0844,
		0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246 };

	double s = 0.0;
	int i;

	for( i = 0; i <= 10; i++ )
	{
		s += sqr(b[i]-(x[0]*(sqr(a[i])+a[i]*x[1]))/
			(sqr(a[i])+a[i]*x[2]+x[3]));
	}

	return( s );
}

static const CTestFn TestFnKowalik = { "Kowalik", 4, -5.0, 5.0,
	0.0003074859878, &calcKowalik };

static double calcTripod( const double* const x, const int N )
{
	const double px = ( x[0] >= 0.0 ? 1.0 : 0.0 );
	const double py = ( x[1] >= 0.0 ? 1.0 : 0.0 );

	return( py*(1.0+px)+fabs(x[0]+50.0*py*(1.0-2.0*px))+
		fabs(x[1]+50.0*(1.0-2.0*py)) );
}

static const CTestFn TestFnTripod = { "Tripod", 2, -100.0, 100.0, 0.0,
	&calcTripod };

static double calcPathological( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += 0.5+(sqr(sin(sqrt(100.0*sqr(x[i+1])+sqr(x[i]))))-0.5)/
			(1.0+0.001*sqr(sqr(x[i])-2.0*x[i]*x[i+1]+sqr(x[i+1])));
	}

	return( s );
}

static const CTestFn TestFnPathological = { "Pathological", 0, -100.0, 100.0,
	0.0, &calcPathological };

static double calcYaoLiu09( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(x[i])-10.0*cos(2.0*M_PI*x[i])+10.0;
	}

	return( s );
}

static const CTestFn TestFnYaoLiu09 = { "YaoLiu09", 0, -5.12, 5.12, 0.0,
	&calcYaoLiu09 };

static double calcUrsem03( const double* const x, const int N )
{
	return( -sin(2.2*M_PI*x[0]+0.5*M_PI)*(2.0-fabs(x[0]))/2.0*
		(3.0-fabs(x[0]))/2.0-sin(2.2*M_PI*x[1]+0.5*M_PI)*
		(2.0-fabs(x[1]))/2.0*(3.0-fabs(x[1]))/2.0 );
}

static double calcUrsem03_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -2.0;
	maxv[ 0 ] = 2.0;
	minv[ 1 ] = -1.5;
	maxv[ 1 ] = 1.5;
	return( -3.0 );
}

static const CTestFn TestFnUrsem03 = { "Ursem03", 2, 0.0, 0.0, 0.0,
	&calcUrsem03, &calcUrsem03_p };

static double calcMishra08( const double* const x, const int N )
{
	const double F1 = fabs(pow(x[0],10.0)-20.0*pow(x[0],9.0)+
		180.0*pow(x[0],8.0)-960.0*pow(x[0],7.0)+3360.0*pow(x[0],6.0)-
		8064.0*pow(x[0],5.0)+13340.0*pow(x[0],4.0)-15360.0*pow(x[0],3.0)+
		11520.0*pow(x[0],2.0)-5120.0*x[0]+2624.0);
	const double F2 = fabs(pow(x[1],4.0)+12.0*pow(x[1],3.0)+
		54.0*pow(x[1],2.0)+108.0*x[1]+81.0);

	return( 0.001*sqr(F1+F2));
}

static const CTestFn TestFnMishra08 = { "Mishra08", 2, -10.0, 10.0, 0.0,
	&calcMishra08 };

static double calcAluffiPentini( const double* const x, const int N )
{
	return( 0.25*sqr(sqr(x[0]))-0.5*sqr(x[0])+0.1*x[0]+0.5*sqr(x[1]) );
}

static const CTestFn TestFnAluffiPentini = { "AluffiPentini", 2, -10.0, 10.0,
	-0.3523860738000, &calcAluffiPentini };

static double calcBeckerLago( const double* const x, const int N )
{
	return( sqr(fabs(x[0])-5.0)+sqr(fabs(x[1])-5.0) );
}

static const CTestFn TestFnBeckerLago = { "BeckerLago", 2, -10.0, 10.0,
	0.0, &calcBeckerLago };

static double calcCosineMixture( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += cos(5.0*M_PI*x[i]);
		s2 += sqr(x[i]);
	}

	return( -0.1*s1+s2 );
}

static double calcCosineMixture_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -1.0;
		maxv[ i ] = 1.0;
	}

	return( -0.1*N );
}

static const CTestFn TestFnCosineMixture = { "CosineMixture", 0, 0.0, 0.0,
	0.0, &calcCosineMixture, &calcCosineMixture_p };

static double calcMeyerRoth( const double* const x, const int N )
{
	static const double t[ 5 ] = { 1.0, 2.0, 1.0, 2.0, 0.1 };
	static const double v[ 5 ] = { 1.0, 1.0, 2.0, 2.0, 0.0 };
	static const double y[ 5 ] = { 0.126, 0.219, 0.076, 0.126, 0.186 };
	double s = 0.0;
	int i;

	for( i = 0; i < 5; i++ )
	{
		s += sqr(x[0]*x[2]*t[i]/(1.0+x[0]*t[i]+x[1]*v[i])-y[i]);
	}

	return( s );
}

static const CTestFn TestFnMeyerRoth = { "MeyerRoth", 3, -20.0, 20.0,
	0.0000435526619, &calcMeyerRoth };

static double calcMultiGaussian( const double* const x, const int N )
{
	static const double a[ 5 ] = { 0.5, 1.2, 1.0, 1.0, 1.2 };
	static const double b[ 5 ] = { 0.0, 1.0, 0.0, -0.5, 0.0 };
	static const double c[ 5 ] = { 0.0, 0.0, -0.5, 0.0, 1.0 };
	static const double d[ 5 ] = { 0.1, 0.5, 0.5, 0.5, 0.5 };
	double s = 0.0;
	int i;

	for( i = 0; i < 5; i++ )
	{
		s += a[i]*exp(-(sqr(x[0]-b[i])+sqr(x[1]-c[i]))/sqr(d[i]));
	}

	return( -s );
}

static const CTestFn TestFnMultiGaussian = { "MultiGaussian", 2, -2.0, 2.0,
	-1.2969540459538, &calcMultiGaussian };

static double calcPeriodic( const double* const x, const int N )
{
	return( 1.0+sqr(sin(x[0]))+sqr(sin(x[1]))-0.1*exp(-sqr(x[0])-sqr(x[1])) );
}

static const CTestFn TestFnPeriodic = { "Periodic", 2, -10.0, 10.0,
	0.9, &calcPeriodic };

static double calcLevyMontalvo2( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(x[i]-1.0)*(1.0+sqr(sin(3.0*M_PI*x[i+1])));
	}

	return( 0.1*(sqr(sin(3.0*M_PI*x[0]))+s+sqr(x[N-1]-1.0)*
		(1.0+sqr(sin(2.0*M_PI*x[N-1])))) );
}

static const CTestFn TestFnLevyMontalvo2 = { "LevyMontalvo2", 0, -5.0, 5.0,
	0.0, &calcLevyMontalvo2 };

static double calcLangermann( const double* const x, const int N )
{
	static const double a[ 5 ] = {3,5,2,1,7};
	static const double b[ 5 ] = {5,2,1,4,9};
	static const double c[ 5 ] = {1,2,5,2,3};
	double s = 0.0;
	int i;

	for( i = 0; i < 5; i++ )
	{
		s += c[i]*exp(-(1.0/M_PI)*(sqr(x[0]-a[i]) +
			sqr(x[1]-b[i])))*cos(M_PI*(sqr(x[0]-a[i])+sqr(x[1]-b[i])));
	}

	return( -s );
}

static const CTestFn TestFnLangermann = { "Langermann", 2, 0.0, 10.0,
	-5.1621261599640, &calcLangermann };

static double calcLangerman5( const double* const x, const int N )
{
	const double a[ 5 ][ 10 ] = {
		{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
		{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
		{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
		{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
		{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567}
	};

	const double c[ 5 ] = { 0.806, 0.517, 0.100, 0.908, 0.965 };

	double s2 = 0.0;
	int i;

	for( i = 0; i < 5; i++ )
	{
		double s = 0.0;
		int j;

		for( j = 0; j < N; j++ )
		{
			s += sqr(x[j]-a[i][j]);
		}

		s2 += c[i]*exp((-1.0/M_PI)*s)*cos(M_PI*s);
	}

	return( -s2 );
}

static const CTestFn TestFnLangerman5 = { "Langerman5", 5, 0.0, 10.0,
	-0.9649999197933, &calcLangerman5 };

static double calcMishra10( const double* const x, const int N )
{
	return( sqr(x[0]+x[1]-x[0]*x[1]) );
}

static const CTestFn TestFnMishra10 = { "Mishra10", 2, -10.0, 10.0,
	0.0, &calcMishra10 };

static double calcMishra10b( const double* const x, const int N )
{
	return( fabs(x[0]+x[1]-x[0]*x[1]) );
}

static const CTestFn TestFnMishra10b = { "Mishra10b", 2, -10.0, 10.0,
	0.0, &calcMishra10b };

static double calcXor( const double* const x, const int N )
{
	const double F11 = x[6]/(1.0+exp(-x[0]-x[1]-x[4]));
	const double F12 = x[7]/(1.0+exp(-x[2]-x[3]-x[5]));
	const double F1 = 1.0/sqr(1.0+exp(-F11-F12-x[8]));
	const double F21 = x[6]/(1.0+exp(-x[4]));
	const double F22 = x[7]/(1.0+exp(-x[5]));
	const double F2 = 1.0/sqr(1.0+exp(-F21-F22-x[8]));
	const double F31 = x[6]/(1.0+exp(-x[0]-x[4]));
	const double F32 = x[7]/(1.0+exp(-x[2]-x[5]));
	const double F3 = sqr(1.0-1.0/(1.0+exp(-F31-F32-x[8])));
	const double F41 = x[6]/(1.0+exp(-x[1]-x[4]));
	const double F42 = x[7]/(1.0+exp(-x[3]-x[5]));
	const double F4 = sqr(1.0-1.0/(1.0+exp(-F41-F42-x[8])));

	return( F1+F2+F3+F4 );
}

static const CTestFn TestFnXor = { "Xor", 9, -1.0, 1.0,
	0.9597587570120, &calcXor };

static double calcRana( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		const double E = x[i]+1.0;
		const double x0 = x[i];

		s += E*cos(sqrt(fabs(E-x0)))*sin(sqrt(fabs(E+x0)))+
			x0*cos(sqrt(fabs(E+x0)))*sin(sqrt(fabs(E-x0)));
	}

	return( s );
}

static double calcRana_p( double* const minv, double* const maxv,
	const int N )
{
	double x[ N ];
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -500.000001;
		maxv[ i ] = 500.000001;
		x[ i ] = -500.0;
	}

	return( calcRana( x, N ));
}

static const CTestFn TestFnRana = { "Rana", 0, 0.0, 0.0, 0.0, &calcRana,
	&calcRana_p };

static double calcPenaltyU( const double xi, const double a,
	const double k, const double m )
{
	if( xi > a )
	{
		return( k * pow( xi - a, m ));
	}

	if( xi < -a )
	{
		return( k * pow( -xi - a, m ));
	}

	return( 0.0 );
}

static double calcPenalty01( const double* const x, const int N )
{
	double y[ N ];
	double penn = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		y[ i ] = 1.0+(x[i]+1.0)/4.0;
		penn += calcPenaltyU( x[ i ], 10, 100, 4 );
	}

	double s = 0.0;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(y[i]-1.0)*(1.0+10.0*sqr(sin(M_PI*y[i+1])));
	}

	return( M_PI/30.0*(10.0*sqr(sin(M_PI*y[0]))+s+sqr(y[N-1]-1.0))+penn );
}

static const CTestFn TestFnPenalty01 = { "Penalty01", 0, -50.0, 50.0,
	0.0, &calcPenalty01 };

static double calcPenalty02( const double* const x, const int N )
{
	double penn = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		penn += calcPenaltyU( x[ i ], 5, 100, 4 );
	}

	double s = 0.0;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(x[i]-1.0)*(1.0+sqr(sin(3.0*M_PI*x[i+1])));
	}

	return( 0.1*(sqr(sin(3.0*M_PI*x[0]))+s+
		sqr(x[N-1]-1.0)*(1.0+sqr(sin(2.0*M_PI*x[N-1]))))+penn );
}

static const CTestFn TestFnPenalty02 = { "Penalty02", 0, -50.0, 50.0,
	0.0, &calcPenalty02 };

static double calcPeaks( const double* const x, const int N )
{
	const double F1 = 3.0*sqr(1.0-x[0])*exp(-sqr(x[0])-sqr(x[1]+1.0));
	const double F2 = 10.0*(x[0]/5.0-pow(x[0],3.0)-pow(x[1],5.0))*
		exp(-sqr(x[0])-sqr(x[1]));

	const double F3 = 1.0/3.0*exp(-sqr(x[0]+1.0)-sqr(x[1]));

	return( F1-F2-F3 );
}

static const CTestFn TestFnPeaks = { "Peaks", 2, -4.0, 4.0,
	-6.5511333328358, &calcPeaks };

static double calcMullerBrown( const double* const x, const int N )
{
	const double A[ 4 ] = { -200, -100, -170, 15 };
	const double a[ 4 ] = { -1, -1, -6.5, 0.7 };
	const double b[ 4 ] = { 0, 0, 11, 0.6 };
	const double c[ 4 ] = { -10, -10, -6.5, 0.7 };
	const double x1[ 4 ] = { 1, 0, -0.5, -1.0 };
	const double x2[ 4 ] = { 0, 0.5, 1.5, 1.0 };
	int j;
	double s = 0.0;

	for( j = 0; j < 4; j++ )
	{
		s += A[j]*exp(a[j]*sqr(x[0]-x1[j])+b[j]*(x[0]-x1[j])*(x[1]-x2[j])+
			c[j]*sqr(x[1]-x2[j]));
	}

	return( s );
}

static double calcMullerBrown_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -1.5;
	maxv[ 0 ] = 1.0;
	minv[ 1 ] = -0.5;
	maxv[ 1 ] = 2.5;
	return( -146.6995172099541 );
}

static const CTestFn TestFnMullerBrown = { "MullerBrown", 2, 0.0, 0.0,
	0.0, &calcMullerBrown, &calcMullerBrown_p };

static double calcCorana( const double* const x, const int N )
{
	const double si = 0.2;
	const double d[ 4 ] = { 1, 1000, 10, 100 };
	double s = 0.0;
	int i;

	for( i = 0; i < 4; i++ )
	{
		const double zi = 0.2*floor(fabs(x[i]/si)+0.49999)*(x[i]<0.0?-1.0:1.0);

		if( fabs( x[i] - zi ) > 0.05 )
		{
			s += 0.15*d[i]*sqr(zi-0.05*(zi<0.0?-1.0:1.0));
		}
		else
		{
			s += d[i]*sqr(x[i]);
		}
	}

	return( s );
}

static const CTestFn TestFnCorana = { "Corana", 4, -100.0, 100.0,
	0.0, &calcCorana };

static double calcBrad( const double* const x, const int N )
{
	const double v[] = { 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37,
		0.58, 0.73, 0.96, 1.34, 2.10, 4.39 };

	double s = 0.0;
	int i;

	for( i = 0; i < 15; i++ )
	{
		const int I = i+1;
		const int J = 16-I;
		const double a = v[i]-x[0]-I;
		const double b = J*x[1]+(I<J?I:J)*x[2];
		s += sqr(a/b);
	}

	return( s );
}

static double calcBrad_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -0.25;
	maxv[ 0 ] = 0.25;
	minv[ 1 ] = 0.01;
	maxv[ 1 ] = 2.5;
	minv[ 2 ] = 0.01;
	maxv[ 2 ] = 2.5;
	return( 6.9352280697052 );
}

static const CTestFn TestFnBrad = { "Brad", 3, 0.0, 0.0,
	0.0, &calcBrad, &calcBrad_p };

static double calcLennardJones( const double* const x, const int N )
{
	int k = N / 3;
	double s = 0.0;
	int i;

	for( i = 0; i < k - 1; i++ )
	{
		int j;

		for( j = i + 1; j < k; j++ )
		{
			const double rij = sqrt(sqr(x[3*i]-x[3*j])+
				sqr(x[3*i+1]-x[3*j+1])+sqr(x[3*i+2]-x[3*j+2]));

			s += 4.0/pow(rij,12.0)-4.0/pow(rij,6.0);
		}
	}

	return( s );
}

static const CTestFn TestFnLennardJones = { "LennardJones", 6, -4.0, 4.0,
	-1.0, &calcLennardJones };

static double calcOddSquare( const double* const x, const int N )
{
	const double b[] = { 1.0, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4,
		1.0, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4 };

	int i;
	double d = sqr(x[0]-b[0]);
	double h = d;

	for( i = 1; i < N; i++ )
	{
		const double nd = sqr(x[i]-b[i]);
		d = ( nd > d ? nd : d );
		h += nd;
	}

	d *= N;

	return( -exp(-d/(2.0*M_PI))*cos(M_PI*d)*(1.0+0.02*h/(d+0.01)) );
}

static const CTestFn TestFnOddSquare = { "OddSquare", 2, -5.0 * M_PI,
	5.0 * M_PI, -1.0084672811395, &calcOddSquare };

static double calcMishra03( const double* const x, const int N )
{
	const double a = fabs(sqr(x[0])+x[1]);
	const double b = sqrt(a);
	const double c = sqrt(fabs(cos(b)));

	return( c+0.01*(x[0]+x[1]) );
}

static const CTestFn TestFnMishra03 = { "Mishra03", 2, -10.00, 10.0,
	-0.1846669920032, &calcMishra03 };

static double calcPermFunction01( const double* const x, const int N )
{
	const double Beta = 1000.0;
	double s1 = 0.0;
	int k;

	for( k = 1; k <= N; k++ )
	{
		double s2 = 0.0;
		int i;

		for( i = 1; i <= N; i++ )
		{
			s2 += ( pow( (double) i, (double) k ) + Beta ) *
				( pow( x[ i - 1 ] / i, (double) k ) - 1.0 );
		}

		s1 += s2 * s2;
	}

	return( s1 );
}

static double calcPermFunction01_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -N;
		maxv[ i ] = N + 1;
	}

	return( 0.0 );
}

static const CTestFn TestFnPermFunction01 = { "PermFunction01", 2, 0.0, 0.0,
	0.0, &calcPermFunction01, &calcPermFunction01_p };

static double calcPermFunction02( const double* const x, const int N )
{
	const double Beta = 1000.0;
	double s1 = 0.0;
	int k;

	for( k = 1; k <= N; k++ )
	{
		double s2 = 0.0;
		int i;

		for( i = 1; i <= N; i++ )
		{
			s2 += ( pow( (double) i, (double) k ) + Beta ) *
				( pow( x[ i - 1 ], (double) k ) - 1.0 / i );
		}

		s1 += s2 * s2;
	}

	return( s1 );
}

static double calcPermFunction02_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -N;
		maxv[ i ] = N + 1;
	}

	return( 0.0 );
}

static const CTestFn TestFnPermFunction02 = { "PermFunction02", 2, 0.0, 0.0,
	0.0, &calcPermFunction02, &calcPermFunction02_p };

/*static double calcWatson( const double* const x, const int N )
{
	double s1 = 0.0;
	int i;

	for( i = 0; i <= 29; i++ )
	{
		const double a = i/29.0;
		double s2 = 0.0;
		int j;

		for( j = 0; j <= 4; j++ )
		{
			s2 += (j-1)*pow(a,j)*x[j];
		}

		double s3 = 0.0;
		int k;

		for( k = 0; k <= 5; k++ )
		{
			s3 += pow(a,k)*x[k];
		}

		s1 += sqr(s2-sqr(s3)-1.0);
	}

	return( s1+sqr(x[0]) );
}

static const CTestFn TestFnWatson = { "Watson", 6, -5.00, 5.0,
	0.002288, &calcWatson };*/

static double calcDeVilliersGlasser01( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= 24; i++ )
	{
		const double ti = 0.1*(i-1);
		const double yi = 60.137*pow(1.371,ti)*sin(3.112*ti+1.761);

		s += sqr(x[0]*pow(x[1],ti)*sin(x[2]*ti+x[3])-yi);
	}

	return( s );
}

static const CTestFn TestFnDeVilliersGlasser01 = { "DeVilliersGlasser01", 4,
	1.0, 100.0, 0.0, &calcDeVilliersGlasser01 };

static double calcPinter( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		double x0;
		double x1;

		if( i == 0 )
		{
			x0 = x[ N - 1 ];
			x1 = x[ i + 1 ];
		}
		else
		if( i == N - 1 )
		{
			x0 = x[ i - 1 ];
			x1 = x[ 0 ];
		}
		else
		{
			x0 = x[ i - 1 ];
			x1 = x[ i + 1 ];
		}

		const double A = x0*sin(x[i])+sin(x1);
		const double B = sqr(x0)-2.0*x[i]+3.0*x1-cos(x[i])+1.0;

		s1 += (i+1)*sqr(x[i]);
		s2 += 20.0*(i+1)*sqr(sin(A));
		s3 += (i+1)*log(1.0+(i+1)*sqr(B))/log(10.0);
	}

	return( s1+s2+s3 );
}

static const CTestFn TestFnPinter = { "Pinter", 0, -10.0, 10.0, 0.0,
	&calcPinter };

static double calcHansen( const double* const x, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i <= 4; i++ )
	{
		s1 += (i+1)*cos(i*x[0]+i+1);
		s2 += (i+1)*cos((i+2)*x[1]+i+1);
	}

	return( s1*s2 );
}

static const CTestFn TestFnHansen = { "Hansen", 2, -10.00, 10.0,
	-176.5417931367458, &calcHansen };

static double calcSchmidtVetters( const double* const x, const int N )
{
	const double a = 1.0/(1.0+sqr(x[0]-x[1]));
	const double b = sin((M_PI*x[1]+x[2])/2.0);
	const double c = exp(sqr(((x[0]+x[1])/x[1])-2.0));

	return( a+b+c );
}

static const CTestFn TestFnSchmidtVetters = { "SchmidtVetters", 3, 0.1, 10.0,
	0.194018, &calcSchmidtVetters };

static double penaltyLE( const double x )
{
	return( x <= 0.0 ? 0.0 : 100000 + x * 9999 );
}

static double calcRosenbrockDisk( const double* const x, const int N )
{
	return( sqr(1.0-x[0])+100.0*sqr(x[1]-sqr(x[0])) +
		penaltyLE( sqr(x[0])+sqr(x[1])-2.0 ));
}

static const CTestFn TestFnRosenbrockDisk = { "RosenbrockDisk", 2, -1.5, 1.5,
	0.0, &calcRosenbrockDisk };

static double calcHougen( const double* const x, const int N )
{
	const double x1[ 13 ] =
		{ 470, 285, 470, 470, 470, 100, 100, 470, 100, 100, 100, 285, 285 };

	const double x2[ 13 ] =
		{ 300, 80, 300, 80, 80, 190, 80, 190, 300, 300, 80, 300, 190 };

	const double x3[ 13 ] =
		{ 10, 10, 120, 120, 10, 10, 65, 65, 54, 120, 120, 10, 120 };

	const double rate[ 13 ] =
		{ 8.55,3.79,4.82,0.02,2.75,14.39,2.54,4.35,13,8.5,0.05,11.32,3.13 };

	double s = 0.0;
	int i;

	for( i = 0; i < 13; i++ )
	{
		double v = (x[0]*x2[i]-x3[i]/x[4])/
			(1.0+x[1]*x1[i]+x[2]*x2[i]+x[3]*x3[i]);

		double d = rate[i]-v;
		s += d * d;
	}

	return( s );
}

static const CTestFn TestFnHougen = { "Hougen", 5, 0.01, 2.0,
	0.2989009807534, &calcHougen };

static double calcSineEnvelope( const double* const x, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		const double ss = sqr(x[i+1])+sqr(x[i]);
		s += (sqr(sin(sqrt(ss)))-0.5)/(0.001*(ss)+1.0)+0.5;
	}

	return( s );
}

static const CTestFn TestFnSineEnvelope = { "SineEnvelope", 0, -100.0, 100.0,
	0.0, &calcSineEnvelope };

static double calcZagros( const double* const x, const int N )
{
	const double a = 2.0;
	const double b = 1.0;
	const double c = 2.0;
	const double d = -1.0;
	const double m = 0.2;
	const double k = 1.0;
	return( -cos(a*x[0]+b*x[1])/(pow(sqr(a*x[0]+b*x[1]),m)+k)-
		cos(c*x[0]+d*x[1])/(pow(sqr(c*x[0]+d*x[1]),m)+k) );
}

static const CTestFn TestFnZagros = { "Zagros", 2, -10.0, 10.0, -2.0,
	&calcZagros };

// Strategy optimization corpus based on N-dimensional functions.

const CTestFn* OptCorpusND[] = { &TestFnSchwefel220, &TestFnSchwefel221,
	&TestFnSchwefel222, &TestFnQing, &TestFnSphere, &TestFnAckley,
	&TestFnAckley2, &TestFnRosenbrock, &TestFnBohachevsky1, &TestFnEasomN,
	&TestFnRastrigin, &TestFnSumSquares, &TestFnZacharov,
	&TestFnRotatedHyperEllipsoid, &TestFnWavy, &TestFnBrown, &TestFnAlpine1,
	&TestFnChungReynolds, &TestFnBentCigar, &TestFnHolzman, &TestFnHyperGrid,
	&TestFnStep01, &TestFnStep02, &TestFnStep03, &TestFnGriewank,
	&TestFnZeroSum, /*&TestFnSchwefel, */&TestFnStyblinskiTank,
	&TestFnYaoLiu04, &TestFnWeierstrass, /*&TestFnXinSheYang04, */
	&TestFnPowellSum, &TestFnAlpine2, &TestFnQuintic, &TestFnBuecheRastrigin,
	&TestFnDifferentPowers, &TestFnDiscus, &TestFnEllipsoid,
	&TestFnGriewankRosenbrock, &TestFnSchaffer07, /*&TestFnTrigonometric01,*/
	&TestFnTrigonometric02, &TestFnExponential, &TestFnSchwefel01,
	&TestFnSchwefel02, &TestFnSchwefel04, &TestFnDeb01, &TestFnLevy03,
	&TestFnYaoLiu09, &TestFnCosineMixture, &TestFnLevyMontalvo2,
	/*&TestFnVincent, &TestFnKatsuura, &TestFnDeb02,*/
	&TestFnDropWave, &TestFnSalomon, &TestFnWhitley, &TestFnXinSheYang02,
	&TestFnXinSheYang03, &TestFnDeflCorrSpring, &TestFnDixonPrice,
	&TestFnPenalty01, &TestFnPenalty02,	&TestFnPinter, &TestFnSchwefel225,
	NULL };

// Failing functions.

const CTestFn* TestCorpusFail[] = { &TestFnDamavandi, &TestFnBukin6,
	&TestFnDeVilliersGlasser02, NULL };

// Failing functions requiring more than 2000 iterations to converge.

const CTestFn* TestCorpusFailTime[] = { &TestFnTrid10, &TestFnMishra03,
	&TestFnMishra04, &TestFnDeVilliersGlasser02, &TestFnHougen, NULL };

// CPU time-consuming functions, but solving correctly:
// &TestFnDeVilliersGlasser01, &TestFnGulfResearchProblem

// Test corpus including all solvable functions.

const CTestFn* TestCorpusAll[] = { &TestFnNewFunction01, &TestFnLangerman5,
	&TestFnEggHolder, &TestFnPowerSum, &TestFnStochastic, &TestFnXor,
	&TestFnTrefethen, &TestFnBranin02, &TestFnNewFunction02, &TestFnDolan,
	&TestFnTripod, &TestFnUrsemWaves, &TestFnSchmidtVetters, &TestFnHartman6,
	&TestFnMeyerRoth, &TestFnBeale, &TestFnModifiedRosenbrock,
	&TestFnMishra06, &TestFnPeaks, &TestFnKowalik, &TestFnChenBird,
	&TestFnPrice03, &TestFnFreudensteinRoth, &TestFnShekel05, &TestFnShekel07,
	&TestFnShekel10, &TestFnGoldsteinPrice, &TestFnWhitley,
	&TestFnCrossLegTable, &TestFnRana, &TestFnHansen,
	&TestFnComplexHolderTable, &TestFnLevy05, &TestFnAlpine2,
	&TestFnSchaffer03, &TestFnSchaffer04, &TestFnXinSheYang01,
	&TestFnCrownedCross, &TestFnMullerBrown, &TestFnOddSquare,
	&TestFnPowellBadlyScaled, &TestFnHelicalValley, &TestFnAckley4,
	&TestFnHartman3, &TestFnTrid6, &TestFnZagros, &TestFnMishra09,
	&TestFnJudge, &TestFnMultiGaussian, &TestFnBranin01, &TestFnGear,
	&TestFnLangermann, &TestFnColville, &TestFnLennardJones,

	&TestFnChenV, &TestFnThreeHumpCamel, &TestFnBooth, &TestFnMatyas,
	&TestFnSphere, &TestFnLevy13, &TestFnSchaffer01, &TestFnSchaffer02,
	&TestFnSchaffer06, &TestFnAckley, &TestFnAckley2, &TestFnAckley3,
	&TestFnRosenbrock, &TestFnBohachevsky1, &TestFnBohachevsky2,
	&TestFnBohachevsky3, &TestFnEasomN, &TestFnEasom, &TestFnCrossInTray,
	&TestFnRastrigin, &TestFnDropWave, &TestFnSumSquares, &TestFnZacharov,
	&TestFnRotatedHyperEllipsoid, &TestFnGriewank, &TestFnSalomon,
	&TestFnPrice02, &TestFnWavy, &TestFnShubert01, &TestFnShubert03,
	&TestFnShubert04, &TestFnWeierstrass, &TestFnTrigonometric02, &TestFnBird,
	&TestFnTreccani, &TestFnXinSheYang02, &TestFnXinSheYang03,
	&TestFnXinSheYang04, &TestFnBiggsEXP2, &TestFnSchwefel06,
	&TestFnChichinadze, &TestFnHolderTable2, &TestFnPowellSum, &TestFnPrice01,
	&TestFnPrice04, &TestFnBrown, &TestFnBrent, &TestFnPowell, &TestFnPaviani,
	&TestFnMieleCantrell, &TestFnWolfe, &TestFnAlpine1,
	&TestFnBukin2, &TestFnBukin4, &TestFnBartelsConn, &TestFnSixHumpCamel,
	&TestFnChungReynolds, &TestFnCube, &TestFnDeckkersAarts, &TestFnEggCrate,
	&TestFnKeane, &TestFnLeon, &TestFnQing, &TestFnSchwefel220,
	&TestFnSchwefel221, &TestFnSchwefel222, &TestFnTestTubeHolder,
	&TestFnWayburnSeader02, &TestFnSawtoothxy, &TestFnSchwefel,
	&TestFnAdjiman, &TestFnStyblinskiTank, &TestFnMcCormick,
	&TestFnHimmelblau, &TestFnMichalewicz, &TestFnBoxBettsExpQuadSum,
	&TestFnPowellQuartic, &TestFnZirilli, &TestFnCamel, &TestFnComplex,
	&TestFnDavis, &TestFnDownhillStep, &TestFnEngvall, &TestFnGramacyLee02,
	&TestFnGiunta, &TestFnHosaki, &TestFnKearfott, &TestFnJennrichSampson,
	&TestFnMishra05, &TestFnTsoulos, &TestFnYaoLiu04, &TestFnBentCigar,
	&TestFnDeflCorrSpring, &TestFnHyperGrid, &TestFnQuintic, &TestFnVincent,
	&TestFnStep01, &TestFnStep02, &TestFnStep03, &TestFnDixonPrice,
	&TestFnZeroSum, &TestFnBuecheRastrigin, &TestFnDifferentPowers,
	&TestFnDiscus, &TestFnEllipsoid, &TestFnGriewankRosenbrock,
	&TestFnSchaffer07, &TestFnKatsuura, &TestFnRotatedEllipse01,
	&TestFnRotatedEllipse02, &TestFnTrigonometric01, &TestFnExponential,
	&TestFnUrsem01, &TestFnQuadratic, &TestFnSchwefel01, &TestFnSchwefel02,
	&TestFnSchwefel04, &TestFnSchwefel236, &TestFnMishra01, &TestFnMishra02,
	&TestFnMishra07, &TestFnZettl, &TestFnMultiModal, &TestFnParsopoulos,
	&TestFnDeb01, &TestFnDeb02, &TestFnCarromTable, &TestFnNewFunction03,
	&TestFnLevy03, &TestFnStretchedV, &TestFnUrsem04, &TestFnWayburnSeader01,
	&TestFnVenterSobiSobieski, &TestFnElAttarVidyasDutta, &TestFnPathological,
	&TestFnYaoLiu09, &TestFnUrsem03, &TestFnMishra08, &TestFnAluffiPentini,
	&TestFnBeckerLago, &TestFnCosineMixture, &TestFnPeriodic,
	&TestFnLevyMontalvo2, &TestFnMishra10, &TestFnMishra10b, &TestFnPenalty01,
	&TestFnPenalty02, &TestFnWayburnSeader03, &TestFnCorana, &TestFnBrad,
	&TestFnGramacyLee03, &TestFnPermFunction01, &TestFnPermFunction02,
	&TestFnPinter, &TestFnHolderTable1, &TestFnSchwefel225,
	&TestFnRosenbrockDisk, &TestFnSineEnvelope, NULL };
