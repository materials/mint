/* Copyright 2011-2014 Kyle Michel, Logan Ward
 *
 * Contact: Kyle Michel (kylemichel@gmail.com)
 *			Logan Ward (LoganWard2012@u.northwestern.edu)
 *
 *
 * This file is part of Mint.
 *
 * Mint is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Mint is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with Mint.  If not, see
 * <http://www.gnu.org/licenses/>.
 */



#include "num.h"
#include "multi.h"
#include "random.h"
#include "output.h"
#include "constants.h"
#include <cmath>



/* void Random::set(RandomPackage package, unsigned long int seed)
 *
 * Set the random number generator
 */

void Random::set(RandomPackage package, unsigned long int seed)
{
	
	// Set package
	if (package == RP_MERSENNE)
	{
		
		// Set method
		delete _generator;
		_generator = new Mersenne;
		_generator->seed(seed);
		return;
	}
	
	// Unknown package
	Output::newline(ERROR);
	Output::print("Attempting to set random number generator with unknown package");
	Output::quit();
}



/* int Random::integerOnNormal(int min, int max, double mean, double standardDeviation) const
 *
 * Generate a random integer on normal distribution
 */

int Random::integerOnNormal(int min, int max, double mean, double standardDeviation) const
{
    
    // Figure out range
    double range = max - mean + 0.5;
    if (mean - min > range)
        range = mean - min + 0.5;
    
    // Loop until a random number is found that is good
    double res;
    do
    {
		res = Num<double>::round(decimalOnNormal(mean, standardDeviation, range), 1);
    }
    while ((res < min) || (res > max));

    // Return result
    return (int)Num<double>::round(res, 1);
}



/* double Random::decimalOnNormal(double mean, double standardDeviation, double range) const
 *
 * Generate a random double on normal distribution
 */

double Random::decimalOnNormal(double mean, double standardDeviation, double range) const
{
	
	// Is spare is set then use it
	double res;
	if (_normalSpareSet == true)
	{
		_normalSpareSet = false;
		res = mean + _normalSpare * standardDeviation;
		if ((res >= mean - range) && (res <= mean + range))
			return res;
	}
    
    // Loop until value is in range
	double mag;
	double mult;
    double rand[2];
    do
    {
        
        // Get random numbers
		do
		{
			rand[0] = decimal(0, 1) * 2 - 1;
			rand[1] = decimal(0, 1) * 2 - 1;
			mag = rand[0]*rand[0] + rand[1]*rand[1];
		} while ((mag >= 1) || (mag == 0));

		// Save multiplier
		mult = sqrt(-2 * log(mag) / mag);
		
		// Save spare value
		_normalSpareSet = true;
		_normalSpare = mult * rand[1];
		
		// Save result
		res = mean + standardDeviation * rand[0] * mult;
    } while ((res < mean - range) || (res > mean + range));
    
    // Return result
    return res;
}



/* unsigned long int Random::readTSC(bool uniqueOnEachProcessor)
 *
 * Return the value of the time stamp counter
 */

unsigned long int Random::readTSC(bool uniqueOnEachProcessor)
{

	// Setup
    unsigned int lo;
    unsigned int hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));

	// Send value
	unsigned long int res = ((unsigned long int) hi << 32) | lo;
	
	// Broadcast value if needed
	if (!uniqueOnEachProcessor)
	{
		unsigned long int sent = res;
		Multi::broadcast(sent, 0);
		res = sent;
	}
	
	// Return result
	return res;
}

