//$ nocpp

// Test functions.

#include <math.h>

#if !defined( sqr )
	#define sqr( x ) (( x ) * ( x ))
#endif // !defined( sqr )

#if !defined( M_PI )
	#define M_PI 3.14159265358979324
#endif // !defined( M_PI )

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
	double (*Calc)( const double* const p, const int N ); ///< Calculation
		///< function.
};

static double calcThreeHumpCamel( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 2*x*x-1.05*sqr(sqr(x))+sqr(sqr(sqr(x)))/6+x*y+y*y );
}

CTestFn TestFnThreeHumpCamel = { "ThreeHumpCamel", 2, -10.0, 10.0, 0.0,
	&calcThreeHumpCamel };

static double calcBooth( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x+2*y-7)+sqr(2*x+y-5));
}

CTestFn TestFnBooth = { "Booth", 2, -10.0, 10.0, 0.0, &calcBooth };

static double calcMatyas( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.26*(x*x+y*y)-0.48*x*y );
}

CTestFn TestFnMatyas = { "Matyas", 2, -10.0, 10.0, 0.0, &calcMatyas };

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

CTestFn TestFnSphere = { "Sphere", 0, -10.0, 10.0, 0.0, &calcSphere };

static double calcLevy13( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-1.0)*(sqr(sin(3.0*M_PI*y))+1.0)+
		sqr(y-1.0)*(sqr(sin(2.0*M_PI*y))+1.0)+sqr(sin(3.0*M_PI*x)));
}

CTestFn TestFnLevy13 = { "Levy13", 2, -10.0, 10.0, 0.0, &calcLevy13 };

static double calcSchaffer01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(sqr(x*x+y*y)))-0.5)/sqr(1.0+0.001*sqr(x*x+y*y)));
}

CTestFn TestFnSchaffer01 = { "Schaffer01", 2, -100.0, 100.0, 0.0,
	&calcSchaffer01 };

static double calcSchaffer02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(x*x-y*y))-0.5)/sqr(1.0+0.001*sqr(x*x+y*y)));
}

CTestFn TestFnSchaffer02 = { "Schaffer02", 2, -100.0, 100.0, 0.0,
	&calcSchaffer02 };

static double calcSchaffer03( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(sin(cos(fabs(x*x-y*y))))-0.5)/
		sqr(1.0+0.001*sqr(x*x+y*y)));
}

CTestFn TestFnSchaffer03 = { "Schaffer03", 2, -100.0, 100.0, 0.001566854526,
	&calcSchaffer03 };

static double calcSchaffer04( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.5+(sqr(cos(sin(fabs(x*x-y*y))))-0.5)/
		sqr(1.0+0.001*sqr(x*x+y*y)));
}

CTestFn TestFnSchaffer04 = { "Schaffer04", 2, -100.0, 100.0, 0.292578632,
	&calcSchaffer04 };

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

CTestFn TestFnAckley = { "Ackley", 0, -32.0, 32.0, 0.0, &calcAckley };

static double calcAckley2( const double* const p, const int N )
{
	return( -200.0*exp(-0.02*sqrt(sqr(p[0])+sqr(p[1]))) );
}

CTestFn TestFnAckley2 = { "Ackley2", 2, -32.0, 32.0, -200.0, &calcAckley2 };

static double calcAckley3( const double* const p, const int N )
{
	return( -200.0*exp(-0.02*sqrt(sqr(p[0])+sqr(p[1])))+
		5.0*exp(cos(3.0*p[0])+sin(3.0*p[1])) );
}

CTestFn TestFnAckley3 = { "Ackley3", 2, -32.0, 32.0, -195.629028238,
	&calcAckley3 };

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

CTestFn TestFnRosenbrock = { "Rosenbrock", 0, -5.0, 10.0, 0.0,
	&calcRosenbrock };

static double calcBeale( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(1.5-x+x*y)+sqr(2.25-x+x*y*y)+sqr(2.625-x+x*y*y*y ));
}

CTestFn TestFnBeale = { "Beale", 2, -10.0, 10.0, 0.0, &calcBeale };

static double calcBohachevsky( const double* const p, const int N )
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

CTestFn TestFnBohachevsky = { "Bohachevsky", 0, -15.0, 15.0, 0.0,
	&calcBohachevsky };

static double calcBohachevsky2( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+2.0*y*y-0.3*cos(3.0*M_PI*x)-0.4*cos(3.0*M_PI*y)+0.7 );
}

CTestFn TestFnBohachevsky2 = { "Bohachevsky2", 2, -100.0, 100.0, 0.0,
	&calcBohachevsky2 };

static double calcBohachevsky3( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+2.0*y*y-0.3*cos(3.0*M_PI*x)*cos(3.0*M_PI*y)+0.3 );
}

CTestFn TestFnBohachevsky3 = { "Bohachevsky3", 2, -100.0, 100.0, 0.0,
	&calcBohachevsky3 };

static double calcEasom( const double* const p, const int N )
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

CTestFn TestFnEasom = { "Easom", 0, -100.0, 100.0, 0.0, &calcEasom };

static double calcEasom2( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -cos(x)*cos(y)*exp(-sqr(x-M_PI)-sqr(y-M_PI)));
}

CTestFn TestFnEasom2 = { "Easom2", 2, -100.0, 100.0, -1.0, &calcEasom2 };

static double calcCrossInTray( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -0.0001*pow(fabs(sin(x)*sin(y)*
		exp(fabs(100.0-sqrt(x*x+y*y)/M_PI)))+1.0,0.1));
}

CTestFn TestFnCrossInTray = { "CrossInTray", 2, -15.0, 15.0, -2.0626118708227,
	&calcCrossInTray };

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

CTestFn TestFnRastrigin = { "Rastrigin", 0, -5.12, 5.12, 0.0, &calcRastrigin };

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

CTestFn TestFnDropWave = { "DropWave", 0, -5.12, 5.12, -1.0, &calcDropWave };

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

CTestFn TestFnSumSquares = { "SumSquares", 0, -5.12, 5.12, 0.0,
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

CTestFn TestFnZacharov = { "Zacharov", 0, -5.0, 10.0, 0.0, &calcZacharov };

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

CTestFn TestFnRotatedHyperEllipsoid = { "RotatHyperEllips", 0,
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

CTestFn TestFnGriewank = { "Griewank", 0, -600.0, 600.0, 0.0, &calcGriewank };

static double calcGoldsteinPrice( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( (1+sqr(x+y+1)*(19-14*x+3*x*x-14*y+6*x*y+3*y*y))*
		(30+sqr(2*x-3*y)*(18-32*x+12*x*x+48*y-36*x*y+27*y*y)));
}

CTestFn TestFnGoldsteinPrice = { "GoldsteinPrice", 2, -2.0, 2.0, 3.0,
	&calcGoldsteinPrice };

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

CTestFn TestFnSalomon = { "Salomon", 0, -100.0, 100.0, 0.0, &calcSalomon };

static double calcBranin01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(-1.275*x*x/(M_PI*M_PI)+5.0*x/M_PI+y-6.0)+
		(10.0-5.0/(4.0*M_PI))*cos(x)+10.0 );
}

CTestFn TestFnBranin01 = { "Branin01", 2, -5.0, 10.0, 0.39788735772973,
	&calcBranin01 };

static double calcBranin02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(-1.275*x*x/(M_PI*M_PI)+5.0*x/M_PI+y-6.0)+(10.0-5.0/(4.0*M_PI))*
		cos(x)*cos(y)+log(x*x+y*y+1.0)+10.0 );
}

CTestFn TestFnBranin02 = { "Branin02", 2, -5.0, 15.0, 5.559037,
	&calcBranin02 };

static double calcTrefethen( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 0.25*x*x+0.25*y*y+exp(sin(50.0*x))-sin(10.0*x+10.0*y)+
		sin(60.0*exp(y))+sin(70.0*sin(x))+sin(sin(80.0*y)) );
}

CTestFn TestFnTrefethen = { "Trefethen", 2, -10.0, 10.0, -3.3068686474,
	&calcTrefethen };

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

CTestFn TestFnWhitley = { "Whitley", 0, -10.24, 10.24, 0.0, &calcWhitley };

static double calcPrice02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 1.0+sqr(sin(x))+sqr(sin(y))-0.1*exp(-x*x-y*y) );
}

CTestFn TestFnPrice02 = { "Price02", 2, -10.0, 10.0, 0.9,
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

CTestFn TestFnWavy = { "Wavy", 0, -M_PI, M_PI, 0.0, &calcWavy };

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

CTestFn TestFnShubert01 = { "Shubert01", 2, -10.0, 10.0, -186.7309,
	&calcShubert01 };

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

CTestFn TestFnShubert03 = { "Shubert03", 2, -10.0, 10.0, -29.6733337,
	&calcShubert03 };

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

CTestFn TestFnShubert04 = { "Shubert04", 2, -10.0, 10.0, -25.740858,
	&calcShubert04 };

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

CTestFn TestFnWeierstrass = { "Weierstrass", 0, -0.5, 0.5, 4.0,
	&calcWeierstrass };

static double calcFreudensteinRoth( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-13.0+((5.0-y)*y-2.0)*y) +
		sqr(x-29.0+((y+1.0)*y-14.0)*y) );
}

CTestFn TestFnFreudensteinRoth = { "FreudensteinRoth", 2, -10.0, 10.0, 0.0,
	&calcFreudensteinRoth };

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

CTestFn TestFnTrigonometric02 = { "Trigonometric02", 0, -500.0, 500.0, 1.0,
	&calcTrigonometric02 };

static double calcBird( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(x-y)+exp(sqr(1.0-sin(x)))*cos(y)+exp(sqr(1.0-cos(y)))*sin(x) );
}

CTestFn TestFnBird = { "Bird", 2, -2.0 * M_PI, 2.0 * M_PI, -106.7645367198034,
	&calcBird };

static double calcTreccani( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(sqr(x))+4.0*sqr(x)*x+4.0*sqr(x)+sqr(y) );
}

CTestFn TestFnTreccani = { "Treccani", 2, -5.0, 5.0, 0.0, &calcTreccani };

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

CTestFn TestFnXinSheYang02 = { "XinSheYang02", 0, -2.0 * M_PI, 2.0 * M_PI,
	0.0, &calcXinSheYang02 };

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

CTestFn TestFnXinSheYang03 = { "XinSheYang03", 0, -2.0 * M_PI, 2.0 * M_PI,
	-1.0, &calcXinSheYang03 };

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

CTestFn TestFnXinSheYang04 = { "XinSheYang04", 0, -10.0, 10.0, -1.0,
	&calcXinSheYang04 };

static double calcExp2( const double* const p, const int N )
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

CTestFn TestFnExp2 = { "Exp2", 2, 0.0, 20.0, 0.0, &calcExp2 };

static double calcSchwefel06( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	double v1 = fabs(x+2.0*y-7.0);
	double v2 = fabs(2.0*x+y-5.0);

	return( v1 > v2 ? v1 : v2 );
}

CTestFn TestFnSchwefel06 = { "Schwefel06", 2, -100.0, 100.0, 0.0,
	&calcSchwefel06 };

static double calcChichinadze( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( x*x-12.0*x+8.0*sin(5.0/2.0*M_PI*x)+10.0*cos(0.5*M_PI*x)+11.0-
		0.2*sqrt(5.0)/exp(0.5*sqr(y-0.5)) );
}

CTestFn TestFnChichinadze = { "Chichinadze", 2, -30.0, 30.0, -42.944387018991,
	&calcChichinadze };

static double calcEggHolder( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( -x*sin(sqrt(fabs(x-y-47.0)))-
		(y+47.0)*sin(sqrt(fabs(0.5*x+y+47.0))) );
}

CTestFn TestFnEggHolder = { "EggHolder", 2, -512.0, 512.0, -959.640662711,
	&calcEggHolder };

static double calcHolderTable( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( -fabs(exp(fabs(1.0-sqrt(x*x+y*y)/M_PI))*sin(x)*cos(y)) );
}

CTestFn TestFnHolderTable = { "HolderTable", 2, -10.0, 10.0, -19.20850256789,
	&calcHolderTable };

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

CTestFn TestFnPowellSum = { "PowellSum", 0, -1.0, 1.0, 0.0, &calcPowellSum };

static double calcPrice01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(fabs(x)-5.0)+sqr(fabs(y)-5.0) );
}

CTestFn TestFnPrice01 = { "Price01", 2, -500.0, 500.0, 0.0,
	&calcPrice01 };

static double calcPrice03( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( 100.0*sqr(y-x*x)+sqr(6.4*sqr(y-0.5)-x-0.6) );
}

CTestFn TestFnPrice03 = { "Price03", 2, -50.0, 50.0, 0.0,
	&calcPrice03 };

static double calcPrice04( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(2.0*x*x*x*y-y*y*y)+sqr(6.0*x-y*y+y) );
}

CTestFn TestFnPrice04 = { "Price04", 2, -50.0, 50.0, 0.0,
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

CTestFn TestFnBrown = { "Brown", 0, -1.0, 4.0, 0.0, &calcBrown };

static double calcBrent( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqr(x+10.0)+sqr(y+10.0)+exp(-x*x-y*y) );
}

CTestFn TestFnBrent = { "Brent", 2, -10.0, 10.0, 0.0, &calcBrent };

static double calcNewFunction01( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqrt(fabs(cos(sqrt(fabs(x*x+y)))))+(x+y)/100.0 );
}

CTestFn TestFnNewFunction01 = { "NewFunction01", 2, -10.0, 10.0, -0.178945093,
	&calcNewFunction01 };

static double calcNewFunction02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( sqrt(fabs(sin(sqrt(fabs(x*x+y)))))+(x+y)/100.0 );
}

CTestFn TestFnNewFunction02 = { "NewFunction02", 2, -10.0, 10.0, -0.197188106,
	&calcNewFunction02 };

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

CTestFn TestFnLevy05 = { "Levy05", 2, -10.0, 10.0, -176.1375,
	&calcLevy05 };

static double calcDamavandi( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];

	return( (1.0-pow(fabs(sin(M_PI*(x-2.0))*sin(M_PI*(y-2.0))/
		(M_PI*M_PI*(x-2.0)*(y-2.0))),5.0))*(2.0+sqr(x-7.0)+2.0*sqr(y-7.0)) );
}

CTestFn TestFnDamavandi = { "Damavandi", 2, 0.0, 14.0, 0.0, &calcDamavandi };

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

CTestFn TestFnPowerSum = { "PowerSum", 4, 0.0, 4.0, 0.0, &calcPowerSum };

static double calcPowell( const double* const p, const int N )
{
	return( sqr(p[2]+10.0*p[0])+5.0*sqr(p[1]-p[3])+sqr(sqr(p[0]-2.0*p[1]))+
		10.0*sqr(sqr(p[2]-p[3])) );
}

CTestFn TestFnPowell = { "Powell", 4, -4.0, 5.0, 0.0, &calcPowell };

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

CTestFn TestFnPaviani = { "Paviani", 10, 2.001, 9.999, -45.7784684040686,
	&calcPaviani };

static double calcDolan( const double* const p, const int N )
{
	return( (p[0]+1.7*p[1])*sin(p[0])-1.5*p[2]-
		0.1*p[3]*cos(p[4]+p[3]-p[0])+0.2*sqr(p[4])-p[1]-1.0 );
}

CTestFn TestFnDolan = { "Dolan", 5, -100.0, 100.0, -529.871438732,
	&calcDolan };

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

CTestFn TestFnTrid6 = { "Trid6", 6, -20.0, 20.0, -50.0, &calcTrid6 };

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

CTestFn TestFnTrid10 = { "Trid10", 10, -100.0, 100.0, -210.0,
	&calcTrid10 };

static double calcMieleCantrell( const double* const p, const int N )
{
	return( sqr(sqr(exp(-p[0])-p[1]))+100.0*pow(p[1]-p[2],6.0)+
		sqr(sqr(tan(p[2]-p[3])))+pow(p[0],8.0) );
}

CTestFn TestFnMieleCantrell = { "MieleCantrell", 4, -1.0, 1.0, 0.0,
	&calcMieleCantrell };

static double calcColville( const double* const p, const int N )
{
	return( 100.0*sqr(p[0]-sqr(p[1]))+sqr(1.0-p[0])+90.0*sqr(p[3]-sqr(p[2]))+
		sqr(1.0-p[2])+10.1*(sqr(p[1]-1.0)+sqr(p[3]-1.0))+
		19.8*(p[1]-1.0)*(p[3]-1.0) );
}

CTestFn TestFnColville = { "Colville", 4, -10.0, 10.0, 0.0,
	&calcColville };

static double calcWolfe( const double* const p, const int N )
{
	return( 4.0/3.0*pow(sqr(p[0])+sqr(p[1])-p[0]*p[1],0.75)+p[2] );
}

CTestFn TestFnWolfe = { "Wolfe", 3, 0.0, 2.0, 0.0, &calcWolfe };

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

CTestFn TestFnAlpine1 = { "Alpine1", 0, -10.0, 10.0, 0.0, &calcAlpine1 };

static double calcBartelsConn( const double* const p, const int N )
{
	return( fabs(sqr(p[0])+sqr(p[1])+p[0]*p[1])+
		fabs(sin(p[0]))+fabs(cos(p[1])) );
}

CTestFn TestFnBartelsConn = { "BartelsConn", 2, -500.0, 500.0, 1.0,
	&calcBartelsConn };

static double calcSixHumpCamel( const double* const p, const int N )
{
	return( (4.0-2.1*sqr(p[0])+pow(p[0],4.0)/3.0)*sqr(p[0])+p[0]*p[1]+
		(4.0*sqr(p[1])-4.0)*sqr(p[1]) );
}

CTestFn TestFnSixHumpCamel = { "SixHumpCamel", 2, -5.0, 5.0, -1.0316,
	&calcSixHumpCamel };

static double calcChenBird( const double* const p, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1])-1.0))-
		b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1])-0.5))-
		b/(sqr(b)+sqr(sqr(p[0])+sqr(p[1]))) );
}

CTestFn TestFnChenBird = { "ChenBird", 2, -500.0, 500.0, -2000.0,
	&calcChenBird };

static double calcChenV( const double* const p, const int N )
{
	const double b = 0.001;
	return( -b/(sqr(b)+sqr(p[0]-0.4*p[1]-0.1))-
		b/(sqr(b)+sqr(2.0*p[0]+p[1]-1.5)) );
}

CTestFn TestFnChenV = { "ChenV", 2, -500.0, 500.0, -2000.0, &calcChenV };

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

CTestFn TestFnChungReynolds = { "ChungReynolds", 0, -100.0, 100.0, 0.0,
	&calcChungReynolds };

static double calcCube( const double* const p, const int N )
{
	return( 100.0*pow(p[1]-pow(p[0],3.0),2.0)+pow(1.0-p[0],2.0) );
}

CTestFn TestFnCube = { "Cube", 2, -10.0, 10.0, 0.0, &calcCube };

static double calcDeckkersAarts( const double* const p, const int N )
{
	return( pow(10.0,5.0)*sqr(p[0])+sqr(p[1])-sqr(sqr(p[0])+sqr(p[1]))+
		pow(10.0,-5.0)*pow(sqr(p[0])+sqr(p[1]),4.0) );
}

CTestFn TestFnDeckkersAarts = { "DeckkersAarts", 2, -20.0, 20.0, -24771.09375,
	&calcDeckkersAarts };

static double calcEggCrate( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( x*x+y*y+25.0*(sqr(sin(x))+sqr(sin(y))) );
}

CTestFn TestFnEggCrate = { "EggCrate", 2, -5.0, 5.0, 0.0, &calcEggCrate };

static double calcKeane( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -sqr(sin(x-y))*sqr(sin(x+y))/sqrt(x*x+y*y) );
}

CTestFn TestFnKeane = { "Keane", 2, 0.00001, 10.0, -0.673667521147,
	&calcKeane };

static double calcLeon( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( 100.0*sqr(y-x*x*x)+sqr(1.0-x) );
}

CTestFn TestFnLeon = { "Leon", 2, 0.0, 10.0, 0.0, &calcLeon };

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

CTestFn TestFnQing = { "Qing", 0, -500.0, 500.0, 0.0, &calcQing };

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

CTestFn TestFnSchwefel = { "Schwefel", 0, -500.0, 500.0, 0.0, &calcSchwefel };

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

CTestFn TestFnSchwefel220 = { "Schwefel220", 0, -100.0, 100.0, 0.0,
	&calcSchwefel220 };

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

CTestFn TestFnSchwefel221 = { "Schwefel221", 0, -100.0, 100.0, 0.0,
	&calcSchwefel221 };

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

CTestFn TestFnSchwefel222 = { "Schwefel222", 0, -100.0, 100.0, 0.0,
	&calcSchwefel222 };

static double calcTestTubeHolder( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( -4.0*fabs(exp(fabs(cos(x*x/200.0+y*y/200.0)))*cos(x)*cos(y)) );
}

CTestFn TestFnTestTubeHolder = { "TestTubeHolder", 2, -10.0, 10.0,
	-10.872299901558, &calcTestTubeHolder };

static double calcWayburnSeader02( const double* const p, const int N )
{
	const double x = p[ 0 ];
	const double y = p[ 1 ];
	return( sqr(1.613-4.0*sqr(x-0.3125)-4.0*sqr(y-1.625))+sqr(y-1.0) );
}

CTestFn TestFnWayburnSeader02 = { "WayburnSeader02", 2, -500.0, 500.0,
	0.0, &calcWayburnSeader02 };

// Strategy optimization corpus based on simple 2D functions.

const CTestFn* OptCorpus2D[] = { &TestFnMatyas, &TestFnThreeHumpCamel,
	&TestFnBooth, &TestFnCrossInTray, &TestFnRotatedHyperEllipsoid,
	&TestFnBohachevsky, &TestFnLevy13, &TestFnAckley,
	&TestFnSchaffer02, &TestFnRosenbrock, &TestFnEasom, &TestFnSchaffer04,
	NULL };

// Strategy optimization corpus based on simple N-dimensional functions.

const CTestFn* OptCorpusND[] = { &TestFnSchwefel220, &TestFnSchwefel221,
	&TestFnSchwefel222, &TestFnQing, &TestFnSphere, &TestFnAckley,
	&TestFnRosenbrock, &TestFnBohachevsky, &TestFnEasom, &TestFnRastrigin,
	&TestFnSumSquares, &TestFnZacharov, &TestFnRotatedHyperEllipsoid,
	&TestFnWavy, &TestFnBrown, &TestFnAlpine1, &TestFnChungReynolds,
	NULL };

// Failing functions.

const CTestFn* TestCorpusFail[] = { &TestFnDamavandi,
	&TestFnTrid10, &TestFnChenBird, &TestFnChenV, &TestFnSchaffer03,
	&TestFnSchaffer04, &TestFnSchwefel, NULL };

// Hard functions.

const CTestFn* TestCorpus2DHard[] = { &TestFnDropWave, &TestFnGriewank,
	&TestFnSalomon, &TestFnBranin02, &TestFnTrefethen, &TestFnWhitley,
	&TestFnPrice02, &TestFnEggHolder, &TestFnPrice03, &TestFnNewFunction01,
	&TestFnNewFunction02, &TestFnLevy05, &TestFnColville, NULL };

/*const CTestFn* TestCorpus2DHard[] = { &TestFnBeale, &TestFnDropWave,
	&TestFnGriewank, &TestFnGoldsteinPrice, &TestFnSalomon, &TestFnBranin02,
	&TestFnTrefethen, &TestFnWhitley, &TestFnPrice02, &TestFnFreudensteinRoth,
	&TestFnBird, &TestFnXinSheYang04, &TestFnEggHolder, &TestFnPrice03,
	&TestFnNewFunction01, &TestFnNewFunction02, &TestFnLevy05,
	&TestFnColville, NULL };*/

// Test corpus including all functions.

const CTestFn* TestCorpusAll[] = { &TestFnThreeHumpCamel, &TestFnBooth,
	&TestFnMatyas, &TestFnSphere, &TestFnLevy13, &TestFnSchaffer01,
	&TestFnSchaffer02, &TestFnAckley, &TestFnAckley2, &TestFnAckley3,
	&TestFnRosenbrock, &TestFnBeale, &TestFnBohachevsky, &TestFnBohachevsky2,
	&TestFnBohachevsky3, &TestFnEasom, &TestFnEasom2, &TestFnCrossInTray,
	&TestFnRastrigin, &TestFnDropWave, &TestFnSumSquares, &TestFnZacharov,
	&TestFnRotatedHyperEllipsoid, &TestFnGriewank, &TestFnGoldsteinPrice,
	&TestFnSalomon, &TestFnBranin01, &TestFnBranin02, &TestFnTrefethen,
	&TestFnWhitley, &TestFnPrice02, &TestFnWavy, &TestFnShubert01,
	&TestFnShubert03, &TestFnShubert04, &TestFnWeierstrass,
	&TestFnFreudensteinRoth, &TestFnTrigonometric02, &TestFnBird,
	&TestFnTreccani, &TestFnXinSheYang02, &TestFnXinSheYang03,
	&TestFnXinSheYang04, &TestFnExp2, &TestFnSchwefel06, &TestFnChichinadze,
	&TestFnEggHolder, &TestFnHolderTable, &TestFnPowellSum, &TestFnPrice01,
	&TestFnPrice03, &TestFnPrice04, &TestFnBrown, &TestFnBrent,
	&TestFnNewFunction01, &TestFnNewFunction02, &TestFnLevy05,
	&TestFnPowerSum, &TestFnPowell, &TestFnPaviani, &TestFnTrid6,
	&TestFnMieleCantrell, &TestFnColville, &TestFnWolfe, &TestFnAlpine1,
	&TestFnBartelsConn, &TestFnSixHumpCamel, &TestFnChungReynolds,
	&TestFnCube, &TestFnDeckkersAarts, &TestFnEggCrate, &TestFnKeane,
	&TestFnLeon, &TestFnQing, &TestFnSchwefel220, &TestFnSchwefel221,
	&TestFnSchwefel222, &TestFnDolan, &TestFnTestTubeHolder,
	&TestFnWayburnSeader02, NULL };

//Need to consider:Adjiman,Alpine2,Bukin N.6,Exponential??,McCormick,Quartic??,
//Styblinski-Tank
