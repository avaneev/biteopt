#include <stdio.h>
#include "tester.h"

int main()
{
	rnd.init( 0 );

	CTester Tester;
	Tester.init( 10, OptCorpusND, 0.01, 120, 150000, false, true );

	Tester.run();
}
