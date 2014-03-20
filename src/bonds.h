/* bonds.h -- Get data about all bonds in a structure
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */


#ifndef BONDS_H
#define BONDS_H



#include "iso.h"
#include "elements.h"



// Class to store information about a single bond between elements
class Bond
{
  
	// Variables
	Element _elem1;
	Element _elem2;
	double _min;
	double _max;
    
public:

	// Functions
	void set(const Element& elem1, const Element& elem2, double min, double max)
		{ _elem1 = elem1; _elem2 = elem2; _min = min; _max = max; }
	
	// Access functions
	const Element& elem1() const	{ return _elem1; }
	const Element& elem2() const	{ return _elem2; }
	double min() const				{ return _min; }
	double max() const				{ return _max; }
};



// Class to get information about bonds
class Bonds
{   
public:
	static OList<Bond> find(const ISO& iso);
//	static OList<double>::D2 findAll(const ISO& iso);
};



#endif
