//$ nocpp

/**
 * @file biternd.h
 *
 * @brief The inclusion file for the CBiteRnd class.
 *
 * @section license License
 * 
 * Copyright (c) 2016-2018 Aleksey Vaneev
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef BITERND_INCLUDED
#define BITERND_INCLUDED

#include <stdint.h>

/**
 * Class that implements random number generation. The default implementation
 * includes a relatively fast pseudo-random number generator (PRNG) using a
 * classic formula "seed = ( a * seed + c ) % m" (LCG). This implementation
 * uses bits 32-61 (30 bits) of the state variable which ensures at least 2^32
 * period in the lowest significant bit of the resulting pseudo-random
 * sequence. See https://en.wikipedia.org/wiki/Linear_congruential_generator
 * for more details.
 */

class CBiteRnd
{
public:
	/**
	 * Function initializes *this PRNG object.
	 *
	 * @param NewSeed New random seed value.
	 */

	void init( const int NewSeed )
	{
		seed = NewSeed;

		// Skip first values to make PRNG "settle down".

		seed = 500009 * seed + 300119;
		seed = 500009 * seed + 300119;
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	double getRndValue()
	{
		seed = 500009 * seed + 300119;

		return((( seed >> 32 ) & 0x3FFFFFFF ) * 9.31322574615478516e-10 );
	}

	/**
	 * @return Scale of "raw" random values returned by functions with the
	 * "raw" suffix.
	 */

	static int32_t getRawScale()
	{
		return( 0x40000000 );
	}

	/**
	 * @return Uniformly-distributed random number in the range [0; 1), in the
	 * "raw" scale.
	 */

	int32_t getUniformRaw()
	{
		seed = 500009 * seed + 300119;

		return( (int32_t) (( seed >> 32 ) & 0x3FFFFFFF ));
	}

	/**
	 * @return TPDF random number in the range (-1; 1), in the "raw" scale.
	 */

	int32_t getTPDFRaw()
	{
		seed = 500009 * seed + 300119;
		const int32_t v1 = (int32_t) (( seed >> 32 ) & 0x3FFFFFFF );
		seed = 500009 * seed + 300119;
		const int32_t v2 = (int32_t) (( seed >> 32 ) & 0x3FFFFFFF );

		return( v1 - v2 );
	}

private:
	uint64_t seed; ///< The current random seed value.
		///<
};

#endif // BITERND_INCLUDED
