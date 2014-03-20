/* elements.h -- Store data about elements
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */


#ifndef ELEMENTS_H
#define ELEMENTS_H



#include "text.h"



// Class to store information about a single element
class Element
{
    
    // Variables
	int _number;
	int _isotope;
	Word _symbol;
	Word _name;
	double _mass;
	double _radius;
    
	// Functions
	static int getNumber(const char* value, bool useNumber);
	static int getNumber(const Word& value, bool useNumber)	{ return getNumber(value.array(), useNumber); }
//	static Element getByNumber(int num);

public:
static Element getByNumber(int num);	
	// Constructors
	Element()					{ _number = _isotope = _mass = _radius = 0; }
	Element(const Element& rhs)	{ *this = rhs; }
	Element(int inNumber, int inIsotope, const char* inSymbol, const char* inName, double inMass, double radius)
		{ set(inNumber, inIsotope, inSymbol, inName, inMass, radius); }

    // Setup functions
	void clear();
	Element& operator= (const Element& rhs);
	void number(int input)			{  _number = input; }
	void isotope(int input)			{ _isotope = input; }
	void symbol(const char* input)	{  _symbol = input; }
	void symbol(const Word& input)	{  _symbol = input; }
	void name(const char* input)	{    _name = input; }
	void name(const Word& input)	{    _name = input; }
	void mass(double input)			{    _mass = input; }
	void radius(double input)		{  _radius = input; }
	void set(int inNumber, int inIsotope, const char* inSymbol, const char* inName, double inMass, double inRadius);
	void set(int inNumber, int inIsotope, const Word& inSymbol, const Word& inName, double inMass, double inRadius);

	// General functions
	bool operator== (const Element& rhs) const	{ return ((_number == rhs._number) && (_isotope == rhs._isotope)); }
	bool operator!= (const Element& rhs) const	{ return ((_number != rhs._number) || (_isotope != rhs._isotope)); }
	
	// Access functions
	int number() const			{ return _number; }
	int isotope() const			{ return _isotope; }
	const Word& symbol() const	{ return _symbol; }
	const Word& name() const	{ return _name; }
	double mass() const			{ return _mass; }
	double radius() const		{ return _radius; }
	
	// Static member functions
	static Element find(const char* value, bool useNumber = true, bool quitIfNotFound = false);
	static Element find(const Word& value, bool useNumber = true, bool quitIfNotFound = false)
		{ return find(value.array(), useNumber, quitIfNotFound); }
	static bool isElement(const char* value, bool useNumber = true) { return (getNumber(value, useNumber) != 0); }
	static bool isElement(const Word& value, bool useNumber = true) { return (getNumber(value, useNumber) != 0); }
};



/* inline void Element::clear()
 *
 * Clear data in Element object
 */

inline void Element::clear()
{
	_number = 0;
	_isotope = 0;
	_symbol.clear();
	_name.clear();
	_mass = 0;
	_radius = 0;
}



#endif
