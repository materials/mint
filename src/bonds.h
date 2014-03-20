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
