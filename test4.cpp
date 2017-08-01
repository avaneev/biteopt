#include <stdio.h>
#include "tester.h"

int main()
{
	CTester Tester;

	rnd.init( 0 );
	Tester.init( 10, OptCorpusND, 0.01, 15, 70000, false, true );

	Tester.run();
}
