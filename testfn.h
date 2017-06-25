//$ nocpp

// Test functions.
//
// Sources:
// http://infinity77.net/global_optimization/test_functions.html
// https://www.sfu.ca/~ssurjano/optimization.html
// http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page364.htm
// http://benchmarkfcns.xyz/fcns
// http://al-roomi.org/benchmarks/unconstrained/2-dimensions/120-damavandi-s-function

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

	double (*CalcFunc)( const double* const p, const int N ); ///< Calculation
		///< function.

	double (*ParamFunc)( double* const minv, double* const maxv, const int N );
};

static double calcThreeHumpCamel( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 2*x*x-1.05*sqr(sqr(x))+sqr(sqr(sqr(x)))/6+x*y+y*y );
}

static const CTestFn TestFnThreeHumpCamel = { "ThreeHumpCamel", 2, -10.0,
	10.0, 0.0, &calcThreeHumpCamel };

static double calcBooth( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x+2*y-7)+sqr(2*x+y-5));
}

static const CTestFn TestFnBooth = { "Booth", 2, -10.0, 10.0, 0.0,
	&calcBooth };

static double calcMatyas( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.26*(x*x+y*y)-0.48*x*y );
}

static const CTestFn TestFnMatyas = { "Matyas", 2, -10.0, 10.0, 0.0,
	&calcMatyas };

static double calcSphere( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr( p[ i ]);
	}

	return( s );
}

static const CTestFn TestFnSphere = { "Sphere", 0, -10.0, 10.0, 0.0,
	&calcSphere };

static double calcLevy13( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-1.0)*(sqr(sin(3.0*M_PI*y))+1.0)+
		sqr(y-1.0)*(sqr(sin(2.0*M_PI*y))+1.0)+sqr(sin(3.0*M_PI*x)));
}

static const CTestFn TestFnLevy13 = { "Levy13", 2, -10.0, 10.0, 0.0,
	&calcLevy13 };

static double calcSchaffer01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(sqr(x*x+y*y)))-0.5)/sqr(1.0+0.001*sqr(x*x+y*y)));
}

static const CTestFn TestFnSchaffer01 = { "Schaffer01", 2, -100.0, 100.0,
	0.0, &calcSchaffer01 };

static double calcSchaffer02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(x*x-y*y))-0.5)/sqr(1.0+0.001*sqr(x*x+y*y)));
}

static const CTestFn TestFnSchaffer02 = { "Schaffer02", 2, -100.0, 100.0,
	0.0, &calcSchaffer02 };

static double calcSchaffer03( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(cos(fabs(x*x-y*y))))-0.5)/
		sqr(1.0+0.001*sqr(x*x+y*y)));
}

static const CTestFn TestFnSchaffer03 = { "Schaffer03", 2, -100.0, 100.0,
	0.002456991, &calcSchaffer03 };

static double calcSchaffer04( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(cos(sin(fabs(x*x-y*y))))-0.5)/
		sqr(1.0+0.001*sqr(x*x+y*y)));
}

static const CTestFn TestFnSchaffer04 = { "Schaffer04", 2, -100.0, 100.0,
	0.292949, &calcSchaffer04 };

static double calcSchaffer06( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(sqrt(x*x+y*y)))-0.5)/
		sqr(1.0+0.001*(x*x+y*y)) );
}

static const CTestFn TestFnSchaffer06 = { "Schaffer06", 2, -100.0, 100.0, 0.0,
	&calcSchaffer06 };

static double calcAckley( const double* const p, const int N )
{
	const double a = 20.0;
	const double b = 0.2;
	const double c = 2.0 * M_PI;
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( p[ i ]);
		s2 += cos( c * p[ i ]);
	}

	return( -a*exp(-b*sqrt(s1/N))-exp(s2/N)+a+exp(1.0));
}

static const CTestFn TestFnAckley = { "Ackley", 0, -32.0, 32.0, 0.0,
	&calcAckley };

static double calcAckley2( const double* const p, const int N )
{
	return( -200.0*exp(-0.02*sqrt(sqr(p[0])+sqr(p[1]))) );
}

static const CTestFn TestFnAckley2 = { "Ackley2", 2, -32.0, 32.0, -200.0,
	&calcAckley2 };

static double calcAckley3( const double* const p, const int N )
{
	return( -200.0*exp(-0.02*sqrt(sqr(p[0])+sqr(p[1])))+
		5.0*exp(cos(3.0*p[0])+sin(3.0*p[1])) );
}

static const CTestFn TestFnAckley3 = { "Ackley3", 2, -32.0, 32.0,
	-195.629028238, &calcAckley3 };

static double calcRosenbrock( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += 100.0*sqr(sqr(p[i])-p[i+1])+sqr(p[i]-1.0);
	}

	return( s );
}

static const CTestFn TestFnRosenbrock = { "Rosenbrock", 0, -5.0, 10.0, 0.0,
	&calcRosenbrock };

static double calcBeale( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(1.5-x+x*y)+sqr(2.25-x+x*y*y)+sqr(2.625-x+x*y*y*y ));
}

static const CTestFn TestFnBeale = { "Beale", 2, -10.0, 10.0, 0.0,
	&calcBeale };

static double calcBohachevsky1( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += sqr(p[i])+2.0*sqr(p[i+1])-0.3*cos(3.0*M_PI*p[i])-
			0.4*cos(4.0*M_PI*p[i+1])+0.7;
	}

	return( s );
}

static const CTestFn TestFnBohachevsky1 = { "Bohachevsky1", 0, -15.0, 15.0,
	0.0, &calcBohachevsky1 };

static double calcBohachevsky2( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+2.0*y*y-0.3*cos(3.0*M_PI*x)*cos(4.0*M_PI*y)+0.3 );
}

static const CTestFn TestFnBohachevsky2 = { "Bohachevsky2", 2, -100.0, 100.0,
	0.0, &calcBohachevsky2 };

static double calcBohachevsky3( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+2.0*y*y-0.3*cos(3.0*M_PI*x+4.0*M_PI*y)+0.3 );
}

static const CTestFn TestFnBohachevsky3 = { "Bohachevsky3", 2, -100.0, 100.0,
	0.0, &calcBohachevsky3 };

static double calcEasomN( const double* const p, const int N )
{
	const double a = 20.0;
	const double b = 0.2;
	const double c = 2.0 * M_PI;
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( p[ i ]);
		s2 += cos( c * p[ i ]);
	}

	return( a-a/exp(b*sqrt(s1/N))+exp(1.0)-exp(s2/N));
}

static const CTestFn TestFnEasomN = { "EasomN", 0, -100.0, 100.0, 0.0,
	&calcEasomN };

static double calcEasom( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -cos(x)*cos(y)*exp(-sqr(x-M_PI)-sqr(y-M_PI)));
}

static const CTestFn TestFnEasom = { "Easom", 2, -100.0, 100.0, -1.0,
	&calcEasom };

static double calcCrossInTray( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -0.0001*pow(fabs(sin(x)*sin(y)*
		exp(fabs(100.0-sqrt(x*x+y*y)/M_PI)))+1.0,0.1));
}

static const CTestFn TestFnCrossInTray = { "CrossInTray", 2, -15.0, 15.0,
	-2.0626118708227, &calcCrossInTray };

static double calcRastrigin( const double* const p, const int N )
{
	double s = 10.0*N;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(p[i])-10.0*cos(2.0*M_PI*p[ i ]);
	}

	return( s );
}

static const CTestFn TestFnRastrigin = { "Rastrigin", 0, -5.12, 5.12, 0.0,
	&calcRastrigin };

static double calcDropWave( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr( p[ i ]);
	}

	return( -(cos(12.0*sqrt(s))+1.0)/(2.0+0.5*s));
}

static const CTestFn TestFnDropWave = { "DropWave", 0, -5.12, 5.12, -1.0,
	&calcDropWave };

static double calcSumSquares( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += (i+1.0)*sqr(p[ i ]);
	}

	return( s );
}

static const CTestFn TestFnSumSquares = { "SumSquares", 0, -5.12, 5.12, 0.0,
	&calcSumSquares };

static double calcZacharov( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( p[ i ]);
		s2 += ( i + 1.0 ) * p[ i ];
	}

	return( s1+sqr(0.5*s2)+sqr(sqr(0.5*s2)));
}

static const CTestFn TestFnZacharov = { "Zacharov", 0, -5.0, 10.0, 0.0,
	&calcZacharov };

static double calcRotatedHyperEllipsoid( const double* const p, const int N )
{
	double s = 0.0;
	int i;
	int j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j <= i; j++ )
		{
			s += sqr( p[ j ]);
		}
	}

	return( s );
}

static const CTestFn TestFnRotatedHyperEllipsoid = { "RotatHyperEllips", 0,
	-65.536, 65.536, 0.0, &calcRotatedHyperEllipsoid };

static double calcGriewank( const double* const p, const int N )
{
	double s1 = 0.0;
	double m2 = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( p[ i ]);
		m2 *= cos( p[ i ] / sqrt( i + 1.0 ));
	}

	return( s1/4000.0-m2+1.0 );
}

static const CTestFn TestFnGriewank = { "Griewank", 0, -600.0, 600.0, 0.0,
	&calcGriewank };

static double calcGoldsteinPrice( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( (1+sqr(x+y+1)*(19-14*x+3*x*x-14*y+6*x*y+3*y*y))*
		(30+sqr(2*x-3*y)*(18-32*x+12*x*x+48*y-36*x*y+27*y*y)));
}

static const CTestFn TestFnGoldsteinPrice = { "GoldsteinPrice", 2, -2.0, 2.0,
	3.0, &calcGoldsteinPrice };

static double calcSalomon( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr( p[ i ]);
	}

	s = sqrt( s );

	return( 1.0-cos(2.0*M_PI*s)+0.1*s);
}

static const CTestFn TestFnSalomon = { "Salomon", 0, -100.0, 100.0, 0.0,
	&calcSalomon };

static double calcBranin01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(-1.275*x*x/(M_PI*M_PI)+5.0*x/M_PI+y-6.0)+
		(10.0-5.0/(4.0*M_PI))*cos(x)+10.0 );
}

static const CTestFn TestFnBranin01 = { "Branin01", 2, -5.0, 10.0,
	0.39788735772973, &calcBranin01 };

static double calcBranin02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(-1.275*x*x/(M_PI*M_PI)+5.0*x/M_PI+y-6.0)+(10.0-5.0/(4.0*M_PI))*
		cos(x)*cos(y)+log(x*x+y*y+1.0)+10.0 );
}

static const CTestFn TestFnBranin02 = { "Branin02", 2, -5.0, 15.0, 5.559037,
	&calcBranin02 };

static double calcTrefethen( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.25*x*x+0.25*y*y+exp(sin(50.0*x))-sin(10.0*x+10.0*y)+
		sin(60.0*exp(y))+sin(70.0*sin(x))+sin(sin(80.0*y)) );
}

static const CTestFn TestFnTrefethen = { "Trefethen", 2, -10.0, 10.0,
	-3.3068686474, &calcTrefethen };

static double calcWhitley( const double* const p, const int N )
{
	double s = 0.0;
	int j;
	int i;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < N; j++ )
		{
			const double v = 100.0*sqr(sqr(p[i])-p[j])+sqr(1.0-p[j]);
			s += sqr(v)/4000.0-cos(v)+1.0;
		}
	}

	return( s );
}

static const CTestFn TestFnWhitley = { "Whitley", 0, -10.24, 10.24, 0.0,
	&calcWhitley };

static double calcPrice02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 1.0+sqr(sin(x))+sqr(sin(y))-0.1*exp(-x*x-y*y) );
}

static const CTestFn TestFnPrice02 = { "Price02", 2, -10.0, 10.0, 0.9,
	&calcPrice02 };

static double calcWavy( const double* const p, const int N )
{
	const double k = 10.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += cos(k*p[i])*exp(-0.5*sqr(p[i]));
	}

	return( 1.0-s/N );
}

static const CTestFn TestFnWavy = { "Wavy", 0, -M_PI, M_PI, 0.0, &calcWavy };

static double calcShubert01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 1; i <= 5; i++ )
	{
		s1 += i*cos((i+1)*x+i);
		s2 += i*cos((i+1)*y+i);
	}

	return( s1 * s2 );
}

static const CTestFn TestFnShubert01 = { "Shubert01", 2, -10.0, 10.0,
	-186.7309, &calcShubert01 };

static double calcShubert03( const double* const p, const int N )
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
			s2 += j*sin((j+1)*p[i]+j);
		}

		s1 += s2;
	}

	return( s1 );
}

static const CTestFn TestFnShubert03 = { "Shubert03", 2, -10.0, 10.0,
	-29.6733337, &calcShubert03 };

static double calcShubert04( const double* const p, const int N )
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
			s2 += j*cos((j+1)*p[i]+j);
		}

		s1 += s2;
	}

	return( s1 );
}

static const CTestFn TestFnShubert04 = { "Shubert04", 2, -10.0, 10.0,
	-25.740858, &calcShubert04 };

static double calcWeierstrass( const double* const p, const int N )
{
	const double a = 0.5;
	const double b = 3.0;
	const int kmax = 20;
	double ak[ kmax + 1 ];
	double bk[ kmax + 1 ];
	double s = 0.0;
	int i;
	int k;

	for( k = 0; k <= kmax; k++ )
	{
		ak[ k ] = pow( a, k );
		bk[ k ] = pow( b, k );
	}

	for( i = 0; i < N; i++ )
	{
		double s1 = 0.0;
		double s2 = 0.0;

		for( k = 0; k <= kmax; k++ )
		{
			s1 += ak[k]*cos(2.0*M_PI*bk[k]*(p[i]+0.5));
			s2 += ak[k]*cos(M_PI*bk[k]);
		}

		s += s1 - N*s2;
	}

	return( s );
}

static const CTestFn TestFnWeierstrass = { "Weierstrass", 0, -0.5, 0.5, 4.0,
	&calcWeierstrass };

static double calcFreudensteinRoth( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-13.0+((5.0-y)*y-2.0)*y) +
		sqr(x-29.0+((y+1.0)*y-14.0)*y) );
}

static const CTestFn TestFnFreudensteinRoth = { "FreudensteinRoth", 2,
	-10.0, 10.0, 0.0, &calcFreudensteinRoth };

static double calcTrigonometric02( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		const double sp = sqr(p[i]-0.9);
		s += 8.0*sqr(sin(7.0*sp))+6.0*sqr(sin(14.0*sp))+sp;
	}

	return( 1.0 + s );
}

static const CTestFn TestFnTrigonometric02 = { "Trigonometric02", 0,
	-500.0, 500.0, 1.0, &calcTrigonometric02 };

static double calcBird( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-y)+exp(sqr(1.0-sin(x)))*cos(y)+exp(sqr(1.0-cos(y)))*sin(x) );
}

static const CTestFn TestFnBird = { "Bird", 2, -2.0 * M_PI, 2.0 * M_PI,
	-106.7645367198034, &calcBird };

static double calcTreccani( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(sqr(x))+4.0*sqr(x)*x+4.0*sqr(x)+sqr(y) );
}

static const CTestFn TestFnTreccani = { "Treccani", 2, -5.0, 5.0, 0.0,
	&calcTreccani };

static double calcXinSheYang02( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += fabs(p[i]);
		s2 += sin(sqr(p[i]));
	}

	return( s1*exp(-s2) );
}

static const CTestFn TestFnXinSheYang02 = { "XinSheYang02", 0,
	-2.0 * M_PI, 2.0 * M_PI, 0.0, &calcXinSheYang02 };

static double calcXinSheYang03( const double* const p, const int N )
{
	const double m = 5.0;
	const double Beta = 15.0;
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += pow(p[i]/Beta,2.0*m);
		s2 += sqr(p[i]);
		s3 *= sqr(cos(p[i]));
	}

	return( exp(-s1)-2.0*exp(-s2)*s3 );
}

static const CTestFn TestFnXinSheYang03 = { "XinSheYang03", 0,
	-2.0 * M_PI, 2.0 * M_PI, -1.0, &calcXinSheYang03 };

static double calcXinSheYang04( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s1 += sqr( sin( p[i]));
		s2 += sqr( p[i]);
		s3 += sqr( sin( sqrt( fabs( p[i]))));
	}

	return( (s1-exp(-s2))*exp(-s3) );
}

static const CTestFn TestFnXinSheYang04 = { "XinSheYang04", 0, -10.0, 10.0,
	-1.0, &calcXinSheYang04 };

static double calcBiggsEXP2( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double s = 0.0;
	int i;

	for( i = 0; i <= 9; i++ )
	{
		s += sqr(exp(-i*x/10.0)-5.0*exp(-i*y/10.0)-exp(-i/10.0)+
			5.0*exp(-(double)i));
	}

	return( s );
}

static const CTestFn TestFnBiggsEXP2 = { "BiggsEXP2", 2, 0.0, 20.0, 0.0,
	&calcBiggsEXP2 };

static double calcSchwefel06( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double v1 = fabs(x+2.0*y-7.0);
	double v2 = fabs(2.0*x+y-5.0);

	return( v1 > v2 ? v1 : v2 );
}

static const CTestFn TestFnSchwefel06 = { "Schwefel06", 2, -100.0, 100.0, 0.0,
	&calcSchwefel06 };

static double calcChichinadze( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( x*x-12.0*x+8.0*sin(5.0/2.0*M_PI*x)+10.0*cos(0.5*M_PI*x)+11.0-
		0.2*sqrt(5.0)/exp(0.5*sqr(y-0.5)) );
}

static const CTestFn TestFnChichinadze = { "Chichinadze", 2, -30.0, 30.0,
	-42.944387018991, &calcChichinadze };

static double calcEggHolder( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( -x*sin(sqrt(fabs(x-y-47.0)))-
		(y+47.0)*sin(sqrt(fabs(0.5*x+y+47.0))) );
}

static const CTestFn TestFnEggHolder = { "EggHolder", 2, -512.0, 512.0,
	-959.640662711, &calcEggHolder };

static double calcHolderTable( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( -fabs(exp(fabs(1.0-sqrt(x*x+y*y)/M_PI))*sin(x)*cos(y)) );
}

static const CTestFn TestFnHolderTable = { "HolderTable", 2, -10.0, 10.0,
	-19.20850256789, &calcHolderTable };

static double calcPowellSum( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow( fabs( p[ i ]), i + 2.0 );
	}

	return( s );
}

static const CTestFn TestFnPowellSum = { "PowellSum", 0, -1.0, 1.0, 0.0,
	&calcPowellSum };

static double calcPrice01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(fabs(x)-5.0)+sqr(fabs(y)-5.0) );
}

static const CTestFn TestFnPrice01 = { "Price01", 2, -500.0, 500.0, 0.0,
	&calcPrice01 };

static double calcPrice03( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( 100.0*sqr(y-x*x)+sqr(6.4*sqr(y-0.5)-x-0.6) );
}

static const CTestFn TestFnPrice03 = { "Price03", 2, -50.0, 50.0, 0.0,
	&calcPrice03 };

static double calcPrice04( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(2.0*x*x*x*y-y*y*y)+sqr(6.0*x-y*y+y) );
}

static const CTestFn TestFnPrice04 = { "Price04", 2, -50.0, 50.0, 0.0,
	&calcPrice04 };

static double calcBrown( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N - 1; i++ )
	{
		s += pow(sqr(p[i]),sqr(p[i+1])+1.0) +
			pow(sqr(p[i+1]),sqr(p[i]+1.0));
	}

	return( s );
}

static const CTestFn TestFnBrown = { "Brown", 0, -1.0, 4.0, 0.0, &calcBrown };

static double calcBrent( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(x+10.0)+sqr(y+10.0)+exp(-x*x-y*y) );
}

static const CTestFn TestFnBrent = { "Brent", 2, -10.0, 10.0, 0.0,
	&calcBrent };

static double calcNewFunction01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqrt(fabs(cos(sqrt(fabs(x*x+y)))))+(x+y)/100.0 );
}

static const CTestFn TestFnNewFunction01 = { "NewFunction01", 2, -10.0, 10.0,
	-0.178945093, &calcNewFunction01 };

static double calcNewFunction02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqrt(fabs(sin(sqrt(fabs(x*x+y)))))+(x+y)/100.0 );
}

static const CTestFn TestFnNewFunction02 = { "NewFunction02", 2, -10.0, 10.0,
	-0.197188106, &calcNewFunction02 };

static double calcLevy05( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 1; i <= 5; i++ )
	{
		s1 += i*cos((i-1)*x+i);
		s2 += i*cos((i+1)*y+i);
	}

	return( s1*s2+sqr(x+1.42513)+sqr(y+0.80032) );
}

static const CTestFn TestFnLevy05 = { "Levy05", 2, -10.0, 10.0, -176.1375,
	&calcLevy05 };

static double calcDamavandi( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( (1.0-pow(fabs(sin(M_PI*(x-2.0))*sin(M_PI*(y-2.0))/
		(M_PI*M_PI*(x-2.0)*(y-2.0))),5.0))*(2.0+sqr(x-7.0)+2.0*sqr(y-7.0)) );
}

static const CTestFn TestFnDamavandi = { "Damavandi", 2, 0.0, 14.0, 0.0,
	&calcDamavandi };

static double calcPowerSum( const double* const p, const int N )
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
			s2 += pow( p[ i ], k + 1.0 );
		}

		s1 += sqr( s2 - b[ k ]);
	}

	return( s1 );
}

static const CTestFn TestFnPowerSum = { "PowerSum", 4, 0.0, 4.0, 0.0,
	&calcPowerSum };

static double calcPowell( const double* const p, const int N )
{
	return( sqr(p[2]+10.0*p[0])+5.0*sqr(p[1]-p[3])+sqr(sqr(p[0]-2.0*p[1]))+
		10.0*sqr(sqr(p[2]-p[3])) );
}

static const CTestFn TestFnPowell = { "Powell", 4, -4.0, 5.0, 0.0,
	&calcPowell };

static double calcPaviani( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 1.0;
	int i;

	for( i = 0; i < 10; i++ )
	{
		s1 += sqr(log(10.0-p[i]))+sqr(log(p[i]-2.0));
		s2 *= p[i];
	}

	return( s1 - pow( s2, 0.2 ));
}

static const CTestFn TestFnPaviani = { "Paviani", 10, 2.001, 9.999,
	-45.7784684040686, &calcPaviani };

static double calcDolan( const double* const p, const int N )
{
	return( (p[0]+1.7*p[1])*sin(p[0])-1.5*p[2]-
		0.1*p[3]*cos(p[4]+p[3]-p[0])+0.2*sqr(p[4])-p[1]-1.0 );
}

static const CTestFn TestFnDolan = { "Dolan", 5, -100.0, 100.0,
	-529.871438732, &calcDolan };

static double calcTrid6( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < 6; i++ )
	{
		s1 += sqr(p[i]-1.0);
	}

	for( i = 1; i < 6; i++ )
	{
		s2 += p[i]*p[i-1];
	}

	return( s1 - s2 );
}

static const CTestFn TestFnTrid6 = { "Trid6", 6, -20.0, 20.0, -50.0,
	&calcTrid6 };

static double calcTrid10( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 0.0;
	int i;

	for( i = 0; i < 10; i++ )
	{
		s1 += sqr(p[i]-1.0);
	}

	for( i = 1; i < 10; i++ )
	{
		s2 += p[i]*p[i-1];
	}

	return( s1 - s2 );
}

static const CTestFn TestFnTrid10 = { "Trid10", 10, -100.0, 100.0, -210.0,
	&calcTrid10 };

static double calcMieleCantrell( const double* const p, const int N )
{
	return( sqr(sqr(exp(-p[0])-p[1]))+100.0*pow(p[1]-p[2],6.0)+
		sqr(sqr(tan(p[2]-p[3])))+pow(p[0],8.0) );
}

static const CTestFn TestFnMieleCantrell = { "MieleCantrell", 4, -1.0, 1.0,
	0.0, &calcMieleCantrell };

static double calcColville( const double* const p, const int N )
{
	return( 100.0*sqr(p[0]-sqr(p[1]))+sqr(1.0-p[0])+90.0*sqr(p[3]-sqr(p[2]))+
		sqr(1.0-p[2])+10.1*(sqr(p[1]-1.0)+sqr(p[3]-1.0))+
		19.8*(p[1]-1.0)*(p[3]-1.0) );
}

static const CTestFn TestFnColville = { "Colville", 4, -10.0, 10.0, 0.0,
	&calcColville };

static double calcWolfe( const double* const p, const int N )
{
	return( 4.0/3.0*pow(sqr(p[0])+sqr(p[1])-p[0]*p[1],0.75)+p[2] );
}

static const CTestFn TestFnWolfe = { "Wolfe", 3, 0.0, 2.0, 0.0, &calcWolfe };

static double calcAlpine1( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(p[i]*sin(p[i]) + 0.1*p[i]);
	}

	return( s );
}

static const CTestFn TestFnAlpine1 = { "Alpine1", 0, -10.0, 10.0, 0.0,
	&calcAlpine1 };

static double calcBartelsConn( const double* const p, const int N )
{
	return( fabs(sqr(p[0])+sqr(p[1])+p[0]*p[1])+
		fabs(sin(p[0]))+fabs(cos(p[1])) );
}

static const CTestFn TestFnBartelsConn = { "BartelsConn", 2, -500.0, 500.0,
	1.0, &calcBartelsConn };

static double calcSixHumpCamel( const double* const p, const int N )
{
	return( (4.0-2.1*sqr(p[0])+pow(p[0],4.0)/3.0)*sqr(p[0])+p[0]*p[1]+
		(4.0*sqr(p[1])-4.0)*sqr(p[1]) );
}

static const CTestFn TestFnSixHumpCamel = { "SixHumpCamel", 2, -5.0, 5.0,
	-1.0316, &calcSixHumpCamel };

static double calcChenBird( const double* const p, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1])-1.0))-
		b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1])-0.5))-
		b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1]))) );
}

static const CTestFn TestFnChenBird = { "ChenBird", 2, -500.0, 500.0, -2000.0,
	&calcChenBird };

static double calcChenV( const double* const p, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(p[0]-0.4*p[1]-0.1))-
		b/(sqr(b)+sqr(2.0*p[0]+p[1]-1.5)) );
}

static const CTestFn TestFnChenV = { "ChenV", 2, -500.0, 500.0, -2000.0,
	&calcChenV };

static double calcChungReynolds( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(p[i]);
	}

	return( s );
}

static const CTestFn TestFnChungReynolds = { "ChungReynolds", 0,
	-100.0, 100.0, 0.0, &calcChungReynolds };

static double calcCube( const double* const p, const int N )
{
	return( 100.0*pow(p[1]-pow(p[0],3.0),2.0)+pow(1.0-p[0],2.0) );
}

static const CTestFn TestFnCube = { "Cube", 2, -10.0, 10.0, 0.0, &calcCube };

static double calcDeckkersAarts( const double* const p, const int N )
{
	return( pow(10.0,5.0)*sqr(p[0])+sqr(p[1])-sqr(sqr(p[0])+sqr(p[1]))+
		pow(10.0,-5.0)*pow(sqr(p[0])+sqr(p[1]),4.0) );
}

static const CTestFn TestFnDeckkersAarts = { "DeckkersAarts", 2, -20.0, 20.0,
	-24771.09375, &calcDeckkersAarts };

static double calcEggCrate( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+y*y+25.0*(sqr(sin(x))+sqr(sin(y))) );
}

static const CTestFn TestFnEggCrate = { "EggCrate", 2, -5.0, 5.0, 0.0,
	&calcEggCrate };

static double calcKeane( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -sqr(sin(x-y))*sqr(sin(x+y))/sqrt(x*x+y*y) );
}

static const CTestFn TestFnKeane = { "Keane", 2, 0.00001, 10.0,
	-0.673667521147, &calcKeane };

static double calcLeon( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 100.0*sqr(y-x*x*x)+sqr(1.0-x) );
}

static const CTestFn TestFnLeon = { "Leon", 2, 0.0, 10.0, 0.0, &calcLeon };

static double calcQing( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sqr(p[i])-(i+1.0));
	}

	return( s );
}

static const CTestFn TestFnQing = { "Qing", 0, -500.0, 500.0, 0.0,
	&calcQing };

static double calcSchwefel220( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(p[i]);
	}

	return( s );
}

static const CTestFn TestFnSchwefel220 = { "Schwefel220", 0, -100.0, 100.0,
	0.0, &calcSchwefel220 };

static double calcSchwefel221( const double* const p, const int N )
{
	double s = fabs(p[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		if( fabs(p[i]) > s )
		{
			s = fabs(p[i]);
		}
	}

	return( s );
}

static const CTestFn TestFnSchwefel221 = { "Schwefel221", 0, -100.0, 100.0,
	0.0, &calcSchwefel221 };

static double calcSchwefel222( const double* const p, const int N )
{
	double s1 = 0.0;
	double s2 = 1.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s1 += fabs(p[i]);
		s2 *= fabs(p[i]);
	}

	return( s1+s2 );
}

static const CTestFn TestFnSchwefel222 = { "Schwefel222", 0, -100.0, 100.0,
	0.0, &calcSchwefel222 };

static double calcTestTubeHolder( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -4.0*fabs(exp(fabs(cos(x*x/200.0+y*y/200.0)))*sin(x)*cos(y)) );
}

static const CTestFn TestFnTestTubeHolder = { "TestTubeHolder", 2,
	-10.0, 10.0, -10.872299901558, &calcTestTubeHolder };

static double calcWayburnSeader02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(1.613-4.0*sqr(x-0.3125)-4.0*sqr(y-1.625))+sqr(y-1.0) );
}

static const CTestFn TestFnWayburnSeader02 = { "WayburnSeader02", 2,
	-500.0, 500.0, 0.0, &calcWayburnSeader02 };

static double calcSawtoothxy( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	const double t = atan2( y, x );
	const double r = sqrt( x * x + y * y );
	const double ht = 0.5*cos(2*t-0.5)+cos(t)+2.0;
	const double gr = (sin(r)-sin(2.0*r)/2.0+sin(3.0*r)/3.0-sin(4.0*r)/4.0+
		4.0)*(r*r/(r+1.0));

	return( gr * ht );
}

static const CTestFn TestFnSawtoothxy = { "Sawtoothxy", 2, -20.0, 20.0, 0.0,
	&calcSawtoothxy };

static double calcMishra04( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqrt(fabs(sin(sqrt(fabs(x*x+y)))))+(x+y)/100.0 );
}

static const CTestFn TestFnMishra04 = { "Mishra04", 2, -10.00, 10.0,
	-0.19940697, &calcMishra04 };

static double calcZeroSum( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += p[i];
	}

	return( s == 0.0 ? 0.0 : 1.0+pow(10000.0*fabs(s),0.5) );
}

static const CTestFn TestFnZeroSum = { "ZeroSum", 0, -10.0, 10.0, 0.0,
	&calcZeroSum };

static double calcSchwefel( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += p[i]*sin(sqrt(fabs(p[i])));
	}

	return( 418.9829 * N - s );
}

static const CTestFn TestFnSchwefel = { "Schwefel", 0, -500.0, 500.0,
	0.00002546, &calcSchwefel };

static double calcAdjiman( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( cos(x)*sin(y)-x/(y*y+1.0) );
}

static double calcAdjiman_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -1.0;
	maxv[ 0 ] = 2.0;
	minv[ 1 ] = -1.0;
	maxv[ 1 ] = 1.0;
	return( -2.0218 );
}

static const CTestFn TestFnAdjiman = { "Adjiman", 2, 0.0, 0.0, 0.0,
	&calcAdjiman, &calcAdjiman_p };

static double calcAlpine2( const double* const p, const int N )
{
	double s = 1.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s *= sqrt(p[i])*sin(p[i]);
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

static double calcBukin4( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 100.0*sqr(y)+0.01*fabs(x+10.0) );
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

static double calcBukin6( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 100.0*sqrt(fabs(y-0.01*x*x))+0.01*fabs(x+10.0) );
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

static double calcStyblinskiTank( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sqr(p[i]))-16.0*sqr(p[i])+5.0*p[i];
	}

	return( 0.5*s );
}

static double calcStyblinskiTank_p( double* const minv, double* const maxv,
	const int N )
{
	int i;

	for( i = 0; i < N; i++ )
	{
		minv[ i ] = -5.0;
		maxv[ i ] = 5.0;
	}

	return( -39.16599 * N );
}

static const CTestFn TestFnStyblinskiTank = { "StyblinskiTank", 0, 0.0, 0.0,
	0.0, &calcStyblinskiTank, &calcStyblinskiTank_p };

static double calcMcCormick( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sin(x+y)+sqr(x-y)-1.5*x+2.5*y+1.0 );
}

static double calcMcCormick_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -1.5;
	maxv[ 0 ] = 4.0;
	minv[ 1 ] = -3.0;
	maxv[ 1 ] = 3.0;
	return( -1.913222887 );
}

static const CTestFn TestFnMcCormick = { "McCormick", 2, 0.0, 0.0, 0.0,
	&calcMcCormick, &calcMcCormick_p };

static double calcHimmelblau( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x*x+y-11.0)+sqr(x+y*y-7.0) );
}

static const CTestFn TestFnHimmelblau = { "Himmelblau", 2, -6.0, 6.0, 0.0,
	&calcHimmelblau };

static double calcMichalewicz( const double* const p, const int N )
{
	const double m = 10.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sin(p[i])*pow(sin((i+1)*sqr(p[i])/M_PI),2.0*m);
	}

	return( -s );
}

static const CTestFn TestFnMichalewicz = { "Michalewicz", 2, 0.0, M_PI,
	-1.8013, &calcMichalewicz };

static double calcBoxBettsExpQuadSum( const double* const p, const int N )
{
	double s = 0.0;
	int j;

	for( j = 1; j <= 10; j++ )
	{
		double g = exp(-0.1*j*p[0])-exp(-0.1*j*p[1])-
			(exp(-0.1*j)-exp(-(double)j))*p[2];

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

static double calcPowellQuartic( const double* const p, const int N )
{
	return( sqr(p[0]+10.0*p[1])+5.0*sqr(p[2]-p[3])+sqr(sqr(p[1]-2.0*p[2]))+
		10.0*sqr(sqr(p[0]-p[3])) );
}

static const CTestFn TestFnPowellQuartic = { "PowellQuartic", 4, -10.0, 10.0,
	0.0, &calcPowellQuartic };

static double calcMishra05( const double* const p, const int N )
{
	const double f1 = sqr(sin(sqr(cos(p[0])+cos(p[1]))));
	const double f2 = sqr(cos(sqr(sin(p[0])+sin(p[1]))));

	return( sqr(f1+f2+p[0])+0.01*p[0]+0.1*p[1] );
}

static const CTestFn TestFnMishra05 = { "Mishra05", 2, -10.0, 10.0,
	-1.01982951993, &calcMishra05 };

static double calcMishra06( const double* const p, const int N )
{
	const double f1 = sqr(sin(sqr(cos(p[0])+cos(p[1]))));
	const double f2 = sqr(cos(sqr(sin(p[0])+sin(p[1]))));
	const double f3 = 0.1*(sqr(p[0]-1.0)+sqr(p[1]-1.0));

	return( -log(sqr(f1-f2+p[0]))+f3 );
}

static const CTestFn TestFnMishra06 = { "Mishra06", 2, -10.0, 10.0,
	-2.28394983847, &calcMishra06 };

static double calcMishra09( const double* const p, const int N )
{
	const double f1 = 2.0*pow(p[0],3.0)+5.0*p[0]*p[1]+4.0*p[2]-
		2.0*sqr(p[0])*p[2]-18.0;

	const double f2 = p[0]+pow(p[1],3.0)+p[0]*sqr(p[1])+p[0]*sqr(p[2])-22.0;
	const double f3 = 8.0*sqr(p[0])+2.0*p[1]*p[2]+2.0*sqr(p[1])+
		3.0*pow(p[1],3.0)-52.0;

	return( sqr(f1*sqr(f2)*f3+f1*f2*sqr(f3)+sqr(f2)+sqr(p[0]+p[1]-p[2])) );
}

static const CTestFn TestFnMishra09 = { "Mishra09", 3, -10.0, 10.0, 0.0,
	&calcMishra09 };

static double calcZirilli( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.25*sqr(sqr(x))-0.5*sqr(x)+0.1*x+0.5*sqr(y) );
}

static const CTestFn TestFnZirilli = { "Zirilli", 2, -10.0, 10.0,
	-0.3523860738, &calcZirilli };

static double calcCamel( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -(-sqr(sqr(x))+4.5*sqr(x)+2.0)/exp(2.0*sqr(y)) );
}

static const CTestFn TestFnCamel = { "Camel", 2, -2.0, 2.0, -7.0625,
	&calcCamel };

static double calcComplex( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(pow(x,3.0)-3.0*x*sqr(y)-1.0)+sqr(3.0*y*sqr(x)-pow(y,3.0)) );
}

static const CTestFn TestFnComplex = { "Complex", 2, -2.0, 2.0, 0.0,
	&calcComplex };

static double calcDavis( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( pow(sqr(x)+sqr(y),0.25)*(sqr(sin(50.0*pow(3.0*sqr(x)+sqr(y),0.1)))+
		1.0) );
}

static const CTestFn TestFnDavis = { "Davis", 2, -100.0, 100.0, 0.0,
	&calcDavis };

static double calcDownhillStep( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( floor(10.0*(10.0-exp(-sqr(x)-3.0*sqr(y))))/10.0 );
}

static const CTestFn TestFnDownhillStep = { "DownhillStep", 2, -10.0, 10.0,
	9.0, &calcDownhillStep };

static double calcEngvall( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(sqr(x))+sqr(sqr(y))+2.0*sqr(x)*sqr(y)-4.0*x+3.0 );
}

static const CTestFn TestFnEngvall = { "Engvall", 2, -2000.0, 2000.0, 0.0,
	&calcEngvall };

static double calcCrossLegTable( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -1.0/pow(fabs(sin(x)*sin(y)*
		exp(fabs(100.0-sqrt(x*x+y*y)/M_PI)))+1.0,0.1));
}

static const CTestFn TestFnCrossLegTable = { "CrossLegTable", 2, -10.0, 10.0,
	-1.0, &calcCrossLegTable };

static double calcCrownedCross( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.0001*pow(fabs(sin(x)*sin(y)*
		exp(fabs(100.0-sqrt(x*x+y*y)/M_PI)))+1.0,0.1) );
}

static const CTestFn TestFnCrownedCross = { "CrownedCross", 2, -10.0, 10.0,
	0.0001, &calcCrownedCross };

static double calcGramacyLee02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*exp(-x*x-y*y) );
}

static const CTestFn TestFnGramacyLee02 = { "GramacyLee02", 2, -1.5, 1.5,
	-0.42888194248, &calcGramacyLee02 };

static double calcGiunta( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(sin(1.0-16.0/15.0*p[i]))-1.0/50.0*sin(4.0-64.0/15.0*p[i])-
			sin(1.0-16.0/15.0*p[i]);
	}

	return( 0.6+s );
}

static const CTestFn TestFnGiunta = { "Giunta", 2, -1.0, 1.0, 0.0644704205369,
	&calcGiunta };

static double calcHosaki( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( (1.0-8.0*x+7.0*sqr(x)-7.0/3.0*sqr(x)*x+1.0/4.0*sqr(sqr(x)))*
		sqr(y)*exp(-y) );
}

static const CTestFn TestFnHosaki = { "Hosaki", 2, 0.0, 10.0, -2.3458115761,
	&calcHosaki };

static double calcKearfott( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x*x+y*y-2.0)+sqr(x*x-y*y-1.0) );
}

static const CTestFn TestFnKearfott = { "Kearfott", 2, -3.0, 4.0, 0,
	&calcKearfott };

static double calcJennrichSampson( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double s = 0.0;
	int i;

	for( i = 1; i <= 10; i++ )
	{
		s += sqr(2.0+2*i-(exp(i*x)+exp(i*y)));
	}

	return( s );
}

static const CTestFn TestFnJennrichSampson = { "JennrichSampson", 2,
	-1.0, 1.0, 124.3621823556, &calcJennrichSampson };

static double calcTsoulos( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x)+sqr(y)-cos(18.0*x)-cos(18.0*y) );
}

static const CTestFn TestFnTsoulos = { "Tsoulos", 2, -1.0, 1.0, -2.0,
	&calcTsoulos };

static double calcUrsemWaves( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -pow(0.3*x,3.0)+(y*y-4.5*y*y)*x*y+4.7*cos(3.0*x-y*y*(2.0+x))*
		sin(2.5*M_PI*x) );
}

static double calcUrsemWaves_p( double* const minv, double* const maxv,
	const int N )
{
	minv[ 0 ] = -0.9;
	maxv[ 0 ] = 1.2;
	minv[ 1 ] = -1.2;
	maxv[ 1 ] = 1.2;
	return( -7.30699873132 );
}

static const CTestFn TestFnUrsemWaves = { "UrsemWaves", 2, 0.0, 0.0, 0.0,
	&calcUrsemWaves, &calcUrsemWaves_p };

static double calcGulfResearchProblem( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i <= 99; i++ )
	{
		const double ti = i/100.0;
		const double ui = 25.0+pow(-50.0*log(ti),1.0/1.5);
		s += sqr(exp(-pow(ui-p[1],p[2])/p[0])-ti);
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

static double calcYaoLiu04( const double* const p, const int N )
{
	double m = fabs(p[0]);
	int i;

	for( i = 1; i < N; i++ )
	{
		const double v = fabs(p[i]);
		m = ( v > m ? v : m );
	}

	return( m );
}

static const CTestFn TestFnYaoLiu04 = { "YaoLiu04", 0, -10.0, 10.0, 0.0,
	&calcYaoLiu04 };

static double calcBentCigar( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 1; i < N; i++ )
	{
		s += sqr( p[ i ]);
	}

	return( sqr(p[0])+1e6*s );
}

static const CTestFn TestFnBentCigar = { "BentCigar", 0, -100.0, 100.0, 0.0,
	&calcBentCigar };

static double calcDeflCorrSpring( const double* const p, const int N )
{
	const double a = 5.0;
	const double k = 5.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(p[i]-a);
	}

	return( 0.1*s-cos(k*sqrt(s)) );
}

static const CTestFn TestFnDeflCorrSpring = { "DeflCorrSpring", 0,
	0.0, 10.0, -1.0, &calcDeflCorrSpring };

static double calcHolzman( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += (i+1)*sqr(sqr(p[i]));
	}

	return( s );
}

static const CTestFn TestFnHolzman = { "Holzman", 0, -10.0, 10.0, 0.0,
	&calcHolzman };

static double calcHyperGrid( const double* const p, const int N )
{
	const double a = 6.0;
	const double c = 5.0;
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += pow(sin(c*M_PI*p[i]),a);
	}

	return( -s/N );
}

static const CTestFn TestFnHyperGrid = { "HyperGrid", 0, 0.0, 1.0, -1.0,
	&calcHyperGrid };

static double calcQuintic( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += fabs(pow(p[i],5.0)-3.0*pow(p[i],4.0)+4.0*pow(p[i],3.0)+
			2.0*sqr(p[i])-10.0*p[i]-4.0);
	}

	return( s );
}

static const CTestFn TestFnQuintic = { "Quintic", 0, -10.0, 10.0, 0.0,
	&calcQuintic };

static double calcVincent( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sin(10.0*log(p[i]));
	}

	return( -s/N );
}

static const CTestFn TestFnVincent = { "Vincent", 0, 0.25, 10.0, -1.0,
	&calcVincent };

static double calcStep01( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += floor(fabs(p[i]));
	}

	return( s );
}

static const CTestFn TestFnStep01 = { "Step01", 0, -100, 100.0, 0.0,
	&calcStep01 };

static double calcStep02( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += sqr(floor(fabs(p[i]+0.5)));
	}

	return( s );
}

static const CTestFn TestFnStep02 = { "Step02", 0, -100, 100.0, 0.0,
	&calcStep02 };

static double calcStep03( const double* const p, const int N )
{
	double s = 0.0;
	int i;

	for( i = 0; i < N; i++ )
	{
		s += floor(sqr(p[i]));
	}

	return( s );
}

static const CTestFn TestFnStep03 = { "Step03", 0, -100, 100.0, 0.0,
	&calcStep03 };

// Strategy optimization corpus based on simple N-dimensional functions.

const CTestFn* OptCorpusND[] = { &TestFnSchwefel220, &TestFnSchwefel221,
	&TestFnSchwefel222, &TestFnQing, &TestFnSphere, &TestFnAckley,
	&TestFnRosenbrock, &TestFnBohachevsky1, &TestFnEasomN, &TestFnRastrigin,
	&TestFnSumSquares, &TestFnZacharov, &TestFnRotatedHyperEllipsoid,
	&TestFnWavy, &TestFnBrown, &TestFnAlpine1, &TestFnChungReynolds,
	&TestFnBentCigar, &TestFnHolzman, &TestFnHyperGrid, &TestFnVincent,
	&TestFnStep01, &TestFnStep02, &TestFnStep03, NULL };

// Failing functions.

const CTestFn* TestCorpusFail[] = { &TestFnDamavandi, &TestFnChenBird,
	&TestFnBukin6, &TestFnCrossLegTable, &TestFnCrownedCross, NULL };

// Failing functions requiring more than 2000 iterations to converge.

const CTestFn* TestCorpusFailTime[] = {
	&TestFnTrid10, &TestFnChenV, &TestFnMishra04, &TestFnPowerSum,
	&TestFnZeroSum, NULL };

// Hard functions.

const CTestFn* TestCorpus2DHard[] = { &TestFnBeale, &TestFnGriewank,
	&TestFnGoldsteinPrice, &TestFnBranin02, &TestFnTrefethen, &TestFnWhitley,
	&TestFnPrice02, &TestFnBird, &TestFnEggHolder, &TestFnPrice03,
	&TestFnNewFunction01, &TestFnNewFunction02, &TestFnAlpine2, &TestFnDolan,
	&TestFnTestTubeHolder };

// Time-consuming function: &TestFnGulfResearchProblem

// Test corpus including all solvable functions.

const CTestFn* TestCorpusAll[] = { &TestFnThreeHumpCamel, &TestFnBooth,
	&TestFnMatyas, &TestFnSphere, &TestFnLevy13, &TestFnSchaffer01,
	&TestFnSchaffer02, &TestFnSchaffer03, &TestFnSchaffer04,
	&TestFnSchaffer06, &TestFnAckley, &TestFnAckley2, &TestFnAckley3,
	&TestFnRosenbrock, &TestFnBeale, &TestFnBohachevsky1, &TestFnBohachevsky2,
	&TestFnBohachevsky3, &TestFnEasomN, &TestFnEasom, &TestFnCrossInTray,
	&TestFnRastrigin, &TestFnDropWave, &TestFnSumSquares, &TestFnZacharov,
	&TestFnRotatedHyperEllipsoid, &TestFnGriewank, &TestFnGoldsteinPrice,
	&TestFnSalomon, &TestFnBranin01, &TestFnBranin02, &TestFnTrefethen,
	&TestFnWhitley, &TestFnPrice02, &TestFnWavy, &TestFnShubert01,
	&TestFnShubert03, &TestFnShubert04, &TestFnWeierstrass,
	&TestFnFreudensteinRoth, &TestFnTrigonometric02, &TestFnBird,
	&TestFnTreccani, &TestFnXinSheYang02, &TestFnXinSheYang03,
	&TestFnXinSheYang04, &TestFnBiggsEXP2, &TestFnSchwefel06,
	&TestFnChichinadze, &TestFnEggHolder, &TestFnHolderTable,
	&TestFnPowellSum, &TestFnPrice01, &TestFnPrice03, &TestFnPrice04,
	&TestFnBrown, &TestFnBrent, &TestFnNewFunction01, &TestFnNewFunction02,
	&TestFnLevy05, &TestFnPowell, &TestFnPaviani, &TestFnTrid6,
	&TestFnMieleCantrell, &TestFnColville, &TestFnWolfe, &TestFnAlpine1,
	&TestFnAlpine2, &TestFnBukin4, &TestFnBartelsConn, &TestFnSixHumpCamel,
	&TestFnChungReynolds, &TestFnCube, &TestFnDeckkersAarts, &TestFnEggCrate,
	&TestFnKeane, &TestFnLeon, &TestFnQing, &TestFnSchwefel220,
	&TestFnSchwefel221, &TestFnSchwefel222, &TestFnDolan,
	&TestFnTestTubeHolder, &TestFnWayburnSeader02, &TestFnSawtoothxy,
	&TestFnSchwefel, &TestFnAdjiman, &TestFnStyblinskiTank, &TestFnMcCormick,
	&TestFnHimmelblau, &TestFnMichalewicz, &TestFnBoxBettsExpQuadSum,
	&TestFnPowellQuartic, &TestFnZirilli, &TestFnCamel, &TestFnComplex,
	&TestFnDavis, &TestFnDownhillStep, &TestFnEngvall, &TestFnGramacyLee02,
	&TestFnGiunta, &TestFnHosaki, &TestFnKearfott, &TestFnJennrichSampson,
	&TestFnMishra05, &TestFnMishra06, &TestFnMishra09, &TestFnTsoulos,
	&TestFnUrsemWaves, &TestFnYaoLiu04, &TestFnBentCigar,
	&TestFnDeflCorrSpring, &TestFnHyperGrid, &TestFnQuintic, &TestFnVincent,
	&TestFnStep01, &TestFnStep02, &TestFnStep03, NULL };
