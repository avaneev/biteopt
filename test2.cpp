#include <stdio.h>
#include "tester.h"

int main()
{
	rnd.init( 0 );

	CTester Tester;
	Tester.init( 2, TestCorpusAll, 0.000001, 500, 2000, false, true );

	Tester.run();
}
