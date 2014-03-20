/* random.h -- Random number generators
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
