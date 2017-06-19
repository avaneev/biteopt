#include <stdio.h>
#include "tester.h"

int main()
{
	CTester Tester;

	rnd.init( 0 );
	Tester.init( 2, TestCorpusAll, 0.000001, 500, 2000, false, true );

	Tester.run();
}
