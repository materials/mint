/* Copyright 2011-2014 Kyle Michel, Logan Ward, Christopher Wolverton
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



#ifndef RANDOM_H
#define RANDOM_H



#include "mtwist.h"
#include "randistrs.h"



// Random number generator packages
enum RandomPackage {RP_UNKNOWN, RP_MERSENNE};



// Class to deal with random number generators
class Random
{
	
	// Base class for random number generators
	class Generator
	{
	public:
		virtual ~Generator() {}
		virtual void seed(unsigned long int seed) = 0;
		virtual int integer(int min, int max) const = 0;
		virtual double decimal(double min, double max) const = 0;
	};
	
	// Mersenne generator
	class Mersenne : public Generator
	{
		mutable mt_distribution _generator;
	public:
		void seed(unsigned long int seed)				{ _generator.seed32(seed); }
		int integer(int min, int max) const				{ return _generator.iuniform(min, max+1); }
		double decimal(double min, double max) const	{ return _generator.uniform(min, max); }
	};
	
	// Variable to store generator
	Generator* _generator;
	mutable bool _normalSpareSet;
	mutable double _normalSpare;
    
public:
    
	// Constructor and destructor
	Random()	{ _generator = new Mersenne; _generator->seed(readTSC()); _normalSpareSet = false; }
	~Random()	{ delete _generator; }
	
	// Setup functions
	void set(RandomPackage package, unsigned long int seed);
	void set(RandomPackage package)		{ set(package, readTSC()); }
	void seed(unsigned long int seed)	{ _generator->seed(seed); }
    
	// Uniform distributions
	int integer(int min, int max) const				{ return _generator->integer(min, max); }
	double decimal(double min, double max) const	{ return _generator->decimal(min, max); }
    
	// Non-uniform distributions
	int integerOnNormal(int min, int max, double mean, double standardDeviation) const;
	double decimalOnNormal(double mean, double standardDeviation, double range) const;
    
	// Static member functions
	static unsigned long int readTSC(bool uniqueOnEachProcessor = false);
};



#endif
