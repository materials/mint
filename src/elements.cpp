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



#include "elements.h"
#include "output.h"



/* Element& Element::operator= (const Element& rhs)
 *
 * Copy element
 */

Element& Element::operator= (const Element& rhs)
{
	if (this != &rhs)
	{
		_number = rhs._number;
		_isotope = rhs._isotope;
		_symbol = rhs._symbol;
		_name = rhs._name;
		_mass = rhs._mass;
		_radius = rhs._radius;
	}
	return *this;
}



/* void Element::set(int inNumber, int inIsotope, const char* inSymbol, const char* inName, double inMass,
 *		double inRadius)
 *
 * Set element properties
 */

void Element::set(int inNumber, int inIsotope, const char* inSymbol, const char* inName, double inMass, \
	double inRadius)
{
	_number = inNumber;
	_isotope = inIsotope;
	_name = inName;
	_symbol = inSymbol;
	_mass = inMass;
	_radius = inRadius;
}



/* void Element::set(int inNumber, int inIsotope, const Word& inName, const Word& inSymbol, double inMass,
 *		double inRadius)
 *
 * Set element properties
 */

void Element::set(int inNumber, int inIsotope, const Word& inName, const Word& inSymbol, double inMass, \
	double inRadius)
{
	_number = inNumber;
	_isotope = inIsotope;
	_name = inName;
	_symbol = inSymbol;
	_mass = inMass;
	_radius = inRadius;
}



/* Element Element::find(const char* value, bool useNumber, bool quitIfNotFound)
 *
 * Find the element by name, symbol, or possibly number
 */

Element Element::find(const char* value, bool useNumber, bool quitIfNotFound)
{
	
	// Get the element number
	int num = getNumber(value, useNumber);
	
	// Did not find element
	if ((!num) && (quitIfNotFound))
	{
		Output::newline(ERROR);
		Output::print("\"");
		Output::print(value);
		Output::print("\" could not be identified as an element");
		Output::quit();
	}
	
	// Return element (this will be empty if element was not found)
	return getByNumber(num);
}



/* int Element::getNumber(const char* value, bool useNumber)
 *
 * Find element by name, symbol, or number
 */

int Element::getNumber(const char* value, bool useNumber)
{
	
	// Hydrogen
    if ((Text::equal(value, "H", true)) || (Text::equal(value, "hydrogen", false)) || \
        ((useNumber) && (Text::equal(value, "1"))))
		return 1;

	// Deuterium
	else if ((Text::equal(value, "D", true)) || (Text::equal(value, "deuterium", false)) || \
        	 ((useNumber) && (Text::equal(value, "1"))))
		return 1001;
	
	// Tritium
	else if ((Text::equal(value, "T", true)) || (Text::equal(value, "tritium", false)) || \
        	 ((useNumber) && (Text::equal(value, "1"))))
		return 1002;
    
    // Helium
	else if ((Text::equal(value, "He", true)) || (Text::equal(value, "helium", false)) || \
        	 ((useNumber) && (Text::equal(value, "2"))))
		return 2;
    
    // Lithium
	else if ((Text::equal(value, "Li", true)) || (Text::equal(value, "lithium", false)) || \
        	 ((useNumber) && (Text::equal(value, "3"))))
		return 3;
    
    // Beryllium
	else if ((Text::equal(value, "Be", true)) || (Text::equal(value, "beryllium", false)) || \
        	 ((useNumber) && (Text::equal(value, "4"))))
		return 4;
    
    // Boron
	else if ((Text::equal(value, "B", true)) || (Text::equal(value, "boron", false)) || \
        	 ((useNumber) && (Text::equal(value, "5"))))
		return 5;
    
    // Carbon
	else if ((Text::equal(value, "C", true)) || (Text::equal(value, "carbon", false)) || \
        	 ((useNumber) && (Text::equal(value, "6"))))
		return 6;
    
    // Nitrogen
	else if ((Text::equal(value, "N", true)) || (Text::equal(value, "nitrogen", false)) || \
        	 ((useNumber) && (Text::equal(value, "7"))))
		return 7;
    
    // Oxygen
	else if ((Text::equal(value, "O", true)) || (Text::equal(value, "oxygen", false)) || \
        	 ((useNumber) && (Text::equal(value, "8"))))
		return 8;
    
    // Fluorine
	else if ((Text::equal(value, "F", true)) || (Text::equal(value, "fluorine", false)) || \
        	 ((useNumber) && (Text::equal(value, "9"))))
		return 9;
    
    // Neon
	else if ((Text::equal(value, "Ne", true)) || (Text::equal(value, "neon", false)) || \
        	 ((useNumber) && (Text::equal(value, "10"))))
		return 10;
    
    // Sodium
	else if ((Text::equal(value, "Na", true)) || (Text::equal(value, "sodium", false)) || \
        	 ((useNumber) && (Text::equal(value, "11"))))
		return 11;

    // Magnesium
	else if ((Text::equal(value, "Mg", true)) || (Text::equal(value, "magnesium", false)) || \
        	 ((useNumber) && (Text::equal(value, "12"))))
		return 12;
    
    // Aluminum
	else if ((Text::equal(value, "Al", true)) || (Text::equal(value, "aluminum", false)) || \
        	 ((useNumber) && (Text::equal(value, "13"))))
		return 13;
    
    // Silicon
	else if ((Text::equal(value, "Si", true)) || (Text::equal(value, "silicon", false)) || \
        	 ((useNumber) && (Text::equal(value, "14"))))
		return 14;
    
    // Phosphorus
	else if ((Text::equal(value, "P", true)) || (Text::equal(value, "phosphorus", false)) || \
        	 ((useNumber) && (Text::equal(value, "15"))))
		return 15;
    
    // Sulfur
	else if ((Text::equal(value, "S", true)) || (Text::equal(value, "sulfur", false)) || \
        	 ((useNumber) && (Text::equal(value, "16"))))
		return 16;
    
    // Chlorine
	else if ((Text::equal(value, "Cl", true)) || (Text::equal(value, "chlorine", false)) || \
        	 ((useNumber) && (Text::equal(value, "17"))))
		return 17;
    
    // Argon
	else if ((Text::equal(value, "Ar", true)) || (Text::equal(value, "argon", false)) || \
        	 ((useNumber) && (Text::equal(value, "18"))))
		return 18;
    
    // Potassium
	else if ((Text::equal(value, "K", true)) || (Text::equal(value, "potassium", false)) || \
        	 ((useNumber) && (Text::equal(value, "19"))))
		return 19;
    
    // Calcium
	else if ((Text::equal(value, "Ca", true)) || (Text::equal(value, "calcium", false)) || \
        	 ((useNumber) && (Text::equal(value, "20"))))
		return 20;
    
    // Scandium
	else if ((Text::equal(value, "Sc", true)) || (Text::equal(value, "scandium", false)) || \
        	 ((useNumber) && (Text::equal(value, "21"))))
		return 21;
    
    // Titanium
	else if ((Text::equal(value, "Ti", true)) || (Text::equal(value, "titanium", false)) || \
        	 ((useNumber) && (Text::equal(value, "22"))))
		return 22;
    
    // Vanadium
	else if ((Text::equal(value, "V", true)) || (Text::equal(value, "vanadium", false)) || \
        	 ((useNumber) && (Text::equal(value, "23"))))
		return 23;
    
    // Chromium
	else if ((Text::equal(value, "Cr", true)) || (Text::equal(value, "chromium", false)) || \
        	 ((useNumber) && (Text::equal(value, "24"))))
		return 24;
    
    // Manganese
	else if ((Text::equal(value, "Mn", true)) || (Text::equal(value, "manganese", false)) || \
        	 ((useNumber) && (Text::equal(value, "25"))))
		return 25;
	
    // Iron
	else if ((Text::equal(value, "Fe", true)) || (Text::equal(value, "iron", false)) || \
        	 ((useNumber) && (Text::equal(value, "26"))))
		return 26;
    
    // Cobalt
	else if ((Text::equal(value, "Co", true)) || (Text::equal(value, "cobalt", false)) || \
        	 ((useNumber) && (Text::equal(value, "27"))))
		return 27;
    
    // Nickel
	else if ((Text::equal(value, "Ni", true)) || (Text::equal(value, "nickel", false)) || \
        	 ((useNumber) && (Text::equal(value, "28"))))
		return 28;
    
    // Copper
	else if ((Text::equal(value, "Cu", true)) || (Text::equal(value, "copper", false)) || \
        	 ((useNumber) && (Text::equal(value, "29"))))
		return 29;
    
    // Zinc
	else if ((Text::equal(value, "Zn", true)) || (Text::equal(value, "zinc", false)) || \
        	 ((useNumber) && (Text::equal(value, "30"))))
		return 30;
    
    // Gallium
	else if ((Text::equal(value, "Ga", true)) || (Text::equal(value, "gallium", false)) || \
        	 ((useNumber) && (Text::equal(value, "31"))))
		return 31;
    
    // Germanium
	else if ((Text::equal(value, "Ge", true)) || (Text::equal(value, "germanium", false)) || \
        	 ((useNumber) && (Text::equal(value, "32"))))
		return 32;
    
    // Arsenic
	else if ((Text::equal(value, "As", true)) || (Text::equal(value, "arsenic", false)) || \
        	 ((useNumber) && (Text::equal(value, "33"))))
		return 33;
    
    // Selenium
	else if ((Text::equal(value, "Se", true)) || (Text::equal(value, "selenium", false)) || \
        	 ((useNumber) && (Text::equal(value, "34"))))
		return 34;
    
    // Bromine
	else if ((Text::equal(value, "Br", true)) || (Text::equal(value, "bromine", false)) || \
        	 ((useNumber) && (Text::equal(value, "35"))))
		return 35;
    
    // Krypton
	else if ((Text::equal(value, "Kr", true)) || (Text::equal(value, "krypton", false)) || \
        	 ((useNumber) && (Text::equal(value, "36"))))
		return 36;
    
    // Rubidium
	else if ((Text::equal(value, "Rb", true)) || (Text::equal(value, "rubidium", false)) || \
        	 ((useNumber) && (Text::equal(value, "37"))))
		return 37;
    
    // Strontium
	else if ((Text::equal(value, "Sr", true)) || (Text::equal(value, "strontium", false)) || \
        	 ((useNumber) && (Text::equal(value, "38"))))
		return 38;
    
    // Yttrium
	else if ((Text::equal(value, "Y", true)) || (Text::equal(value, "yttrium", false)) || \
        	 ((useNumber) && (Text::equal(value, "39"))))
		return 39;
    
    // Zirconium
	else if ((Text::equal(value, "Zr", true)) || (Text::equal(value, "zirconium", false)) || \
        	 ((useNumber) && (Text::equal(value, "40"))))
		return 40;
    
    // Niobium
	else if ((Text::equal(value, "Nb", true)) || (Text::equal(value, "niobium", false)) || \
        	 ((useNumber) && (Text::equal(value, "41"))))
		return 41;
    
    // Molybdenum
	else if ((Text::equal(value, "Mo", true)) || (Text::equal(value, "molybdenum", false)) || \
        	 ((useNumber) && (Text::equal(value, "42"))))
		return 42;
    
    // Technetium
	else if ((Text::equal(value, "Tc", true)) || (Text::equal(value, "technetium", false)) || \
        	 ((useNumber) && (Text::equal(value, "43"))))
		return 43;
    
    // Ruthenium
	else if ((Text::equal(value, "Ru", true)) || (Text::equal(value, "ruthenium", false)) || \
        	 ((useNumber) && (Text::equal(value, "44"))))
		return 44;
    
    // Rhodium
	else if ((Text::equal(value, "Rh", true)) || (Text::equal(value, "rhodium", false)) || \
        	 ((useNumber) && (Text::equal(value, "45"))))
		return 45;
    
    // Palladium
	else if ((Text::equal(value, "Pd", true)) || (Text::equal(value, "palladium", false)) || \
        	 ((useNumber) && (Text::equal(value, "46"))))
		return 46;
    
    // Silver
	else if ((Text::equal(value, "Ag", true)) || (Text::equal(value, "silver", false)) || \
        	 ((useNumber) && (Text::equal(value, "47"))))
		return 47;
    
    // Cadmium
	else if ((Text::equal(value, "Cd", true)) || (Text::equal(value, "cadmium", false)) || \
        	 ((useNumber) && (Text::equal(value, "48"))))
		return 48;
    
    // Indium
	else if ((Text::equal(value, "In", true)) || (Text::equal(value, "indium", false)) || \
        	 ((useNumber) && (Text::equal(value, "49"))))
		return 49;
    
    // Tin
	else if ((Text::equal(value, "Sn", true)) || (Text::equal(value, "tin", false)) || \
        	 ((useNumber) && (Text::equal(value, "50"))))
		return 50;
    
    // Antimony
	else if ((Text::equal(value, "Sb", true)) || (Text::equal(value, "antimony", false)) || \
        	 ((useNumber) && (Text::equal(value, "51"))))
		return 51;
    
    // Tellurium
	else if ((Text::equal(value, "Te", true)) || (Text::equal(value, "tellurium", false)) || \
        	 ((useNumber) && (Text::equal(value, "52"))))
		return 52;
    
    // Iodine
	else if ((Text::equal(value, "I", true)) || (Text::equal(value, "iodine", false)) || \
        	 ((useNumber) && (Text::equal(value, "53"))))
		return 53;
    
    // Xenon
	else if ((Text::equal(value, "Xe", true)) || (Text::equal(value, "xenon", false)) || \
        	 ((useNumber) && (Text::equal(value, "54"))))
		return 54;
    
    // Cesium
	else if ((Text::equal(value, "Cs", true)) || (Text::equal(value, "cesium", false)) || \
        	 ((useNumber) && (Text::equal(value, "55"))))
		return 55;
    
    // Barium
	else if ((Text::equal(value, "Ba", true)) || (Text::equal(value, "barium", false)) || \
        	 ((useNumber) && (Text::equal(value, "56"))))
		return 56;
    
    // Lanthanum
	else if ((Text::equal(value, "La", true)) || (Text::equal(value, "lanthanum", false)) || \
        	 ((useNumber) && (Text::equal(value, "57"))))
		return 57;
    
    // Cerium
	else if ((Text::equal(value, "Ce", true)) || (Text::equal(value, "cerium", false)) || \
        	 ((useNumber) && (Text::equal(value, "58"))))
		return 58;
    
    // Praseodymium
	else if ((Text::equal(value, "Pr", true)) || (Text::equal(value, "praseodymium", false)) || \
        	 ((useNumber) && (Text::equal(value, "59"))))
		return 59;
    
    // Neodymium
	else if ((Text::equal(value, "Nd", true)) || (Text::equal(value, "neodymium", false)) || \
        	 ((useNumber) && (Text::equal(value, "60"))))
		return 60;
    
    // Promethium
	else if ((Text::equal(value, "Pm", true)) || (Text::equal(value, "promethium", false)) || \
        	 ((useNumber) && (Text::equal(value, "61"))))
		return 61;
    
    // Samarium
	else if ((Text::equal(value, "Sm", true)) || (Text::equal(value, "samarium", false)) || \
        	 ((useNumber) && (Text::equal(value, "62"))))
		return 62;
    
    // Europium
	else if ((Text::equal(value, "Eu", true)) || (Text::equal(value, "europium", false)) || \
        	 ((useNumber) && (Text::equal(value, "63"))))
		return 63;
    
    // Gadolinium
	else if ((Text::equal(value, "Gd", true)) || (Text::equal(value, "gadolinium", false)) || \
        	 ((useNumber) && (Text::equal(value, "64"))))
		return 64;
    
    // Terbium
	else if ((Text::equal(value, "Tb", true)) || (Text::equal(value, "terbium", false)) || \
        	 ((useNumber) && (Text::equal(value, "65"))))
		return 65;
    
    // Dysprosium
	else if ((Text::equal(value, "Dy", true)) || (Text::equal(value, "dysprosium", false)) || \
        	 ((useNumber) && (Text::equal(value, "66"))))
		return 66;
    
    // Holmium
	else if ((Text::equal(value, "Ho", true)) || (Text::equal(value, "holmium", false)) || \
        	 ((useNumber) && (Text::equal(value, "67"))))
		return 67;
    
    // Erbium
	else if ((Text::equal(value, "Er", true)) || (Text::equal(value, "erbium", false)) || \
        	 ((useNumber) && (Text::equal(value, "68"))))
		return 68;
    
    // Thulium
	else if ((Text::equal(value, "Tm", true)) || (Text::equal(value, "thulium", false)) || \
        	 ((useNumber) && (Text::equal(value, "69"))))
		return 69;
    
    // Ytterbium
	else if ((Text::equal(value, "Yb", true)) || (Text::equal(value, "ytterbium", false)) || \
        	 ((useNumber) && (Text::equal(value, "70"))))
		return 70;
    
    // Lutetium
	else if ((Text::equal(value, "Lu", true)) || (Text::equal(value, "lutetium", false)) || \
        	 ((useNumber) && (Text::equal(value, "71"))))
		return 71;
    
    // Hafnium
	else if ((Text::equal(value, "Hf", true)) || (Text::equal(value, "hafnium", false)) || \
        	 ((useNumber) && (Text::equal(value, "72"))))
		return 72;
    
    // Tantalum
	else if ((Text::equal(value, "Ta", true)) || (Text::equal(value, "tantalum", false)) || \
        	 ((useNumber) && (Text::equal(value, "73"))))
		return 73;
    
    // Tungsten
	else if ((Text::equal(value, "W", true)) || (Text::equal(value, "tungsten", false)) || \
        	 ((useNumber) && (Text::equal(value, "74"))))
		return 74;
    
    // Rhenium
	else if ((Text::equal(value, "Re", true)) || (Text::equal(value, "rhenium", false)) || \
        	 ((useNumber) && (Text::equal(value, "75"))))
		return 75;
    
    // Osmium
	else if ((Text::equal(value, "Os", true)) || (Text::equal(value, "osmium", false)) || \
        	 ((useNumber) && (Text::equal(value, "76"))))
		return 76;
    
    // Iridium
	else if ((Text::equal(value, "Ir", true)) || (Text::equal(value, "iridium", false)) || \
        	 ((useNumber) && (Text::equal(value, "77"))))
		return 77;
    
    // Platinum
	else if ((Text::equal(value, "Pt", true)) || (Text::equal(value, "platinum", false)) || \
        	 ((useNumber) && (Text::equal(value, "78"))))
		return 78;
    
    // Gold
	else if ((Text::equal(value, "Au", true)) || (Text::equal(value, "gold", false)) || \
        	 ((useNumber) && (Text::equal(value, "79"))))
		return 79;
    
    // Mercury
	else if ((Text::equal(value, "Hg", true)) || (Text::equal(value, "mercury", false)) || \
        	 ((useNumber) && (Text::equal(value, "80"))))
		return 80;
    
    // Thallium
	else if ((Text::equal(value, "Tl", true)) || (Text::equal(value, "thallium", false)) || \
        	 ((useNumber) && (Text::equal(value, "81"))))
		return 81;
    
    // Lead
	else if ((Text::equal(value, "Pb", true)) || (Text::equal(value, "lead", false)) || \
        	 ((useNumber) && (Text::equal(value, "82"))))
		return 82;
    
    // Bismuth
	else if ((Text::equal(value, "Bi", true)) || (Text::equal(value, "bismuth", false)) || \
        	 ((useNumber) && (Text::equal(value, "83"))))
		return 83;
    
    // Polonium
	else if ((Text::equal(value, "Po", true)) || (Text::equal(value, "polonium", false)) || \
        	 ((useNumber) && (Text::equal(value, "84"))))
		return 84;
    
    // Astatine
	else if ((Text::equal(value, "At", true)) || (Text::equal(value, "astatine", false)) || \
        	 ((useNumber) && (Text::equal(value, "85"))))
		return 85;
    
    // Radon
	else if ((Text::equal(value, "Rn", true)) || (Text::equal(value, "radon", false)) || \
        	 ((useNumber) && (Text::equal(value, "86"))))
		return 86;
    
    // Francium
	else if ((Text::equal(value, "Fr", true)) || (Text::equal(value, "francium", false)) || \
        	 ((useNumber) && (Text::equal(value, "87"))))
		return 87;
    
    // Radium
	else if ((Text::equal(value, "Ra", true)) || (Text::equal(value, "radium", false)) || \
        	 ((useNumber) && (Text::equal(value, "88"))))
		return 88;
    
    // Actinium
	else if ((Text::equal(value, "Ac", true)) || (Text::equal(value, "actinium", false)) || \
        	 ((useNumber) && (Text::equal(value, "89"))))
		return 89;
    
    // Thorium
	else if ((Text::equal(value, "Th", true)) || (Text::equal(value, "thorium", false)) || \
        	 ((useNumber) && (Text::equal(value, "90"))))
		return 90;
    
    // Protactinium
	else if ((Text::equal(value, "Pa", true)) || (Text::equal(value, "protactinium", false)) || \
        	 ((useNumber) && (Text::equal(value, "91"))))
		return 91;
    
    // Uranium
	else if ((Text::equal(value, "U", true)) || (Text::equal(value, "uranium", false)) || \
        	 ((useNumber) && (Text::equal(value, "92"))))
		return 92;
    
    // Neptunium
	else if ((Text::equal(value, "Np", true)) || (Text::equal(value, "neptunium", false)) || \
        	 ((useNumber) && (Text::equal(value, "93"))))
		return 93;
    
    // Plutonium
	else if ((Text::equal(value, "Pu", true)) || (Text::equal(value, "plutonium", false)) || \
        	 ((useNumber) && (Text::equal(value, "94"))))
		return 94;
    
    // Americium
	else if ((Text::equal(value, "Am", true)) || (Text::equal(value, "americium", false)) || \
        	 ((useNumber) && (Text::equal(value, "95"))))
		return 95;
    
    // Curium
	else if ((Text::equal(value, "Cm", true)) || (Text::equal(value, "curium", false)) || \
        	 ((useNumber) && (Text::equal(value, "96"))))
		return 96;
    
    // Berkelium
	else if ((Text::equal(value, "Bk", true)) || (Text::equal(value, "berkelium", false)) || \
        	 ((useNumber) && (Text::equal(value, "97"))))
		return 97;
    
    // Californium
	else if ((Text::equal(value, "Cf", true)) || (Text::equal(value, "californium", false)) || \
        	 ((useNumber) && (Text::equal(value, "98"))))
		return 98;
    
    // Einsteinium
	else if ((Text::equal(value, "Es", true)) || (Text::equal(value, "einsteinium", false)) || \
        	 ((useNumber) && (Text::equal(value, "99"))))
		return 99;
    
    // Fermium
	else if ((Text::equal(value, "Fm", true)) || (Text::equal(value, "fermium", false)) || \
        	 ((useNumber) && (Text::equal(value, "100"))))
		return 100;
    
    // Mendelevium
	else if ((Text::equal(value, "Md", true)) || (Text::equal(value, "mendelevium", false)) || \
        	 ((useNumber) && (Text::equal(value, "101"))))
		return 101;
    
    // Nobelium
	else if ((Text::equal(value, "No", true)) || (Text::equal(value, "nobelium", false)) || \
        	 ((useNumber) && (Text::equal(value, "102"))))
		return 102;
    
    // Lawrencium
	else if ((Text::equal(value, "Lr", true)) || (Text::equal(value, "lawrencium", false)) || \
        	 ((useNumber) && (Text::equal(value, "103"))))
		return 103;
    
    // Rutherfordium
	else if ((Text::equal(value, "Rf", true)) || (Text::equal(value, "rutherfordium", false)) || \
        	 ((useNumber) && (Text::equal(value, "104"))))
		return 104;
    
    // Dubnium
	else if ((Text::equal(value, "Db", true)) || (Text::equal(value, "dubnium", false)) || \
        	 ((useNumber) && (Text::equal(value, "105"))))
		return 105;
    
    // Seaborgium
	else if ((Text::equal(value, "Sg", true)) || (Text::equal(value, "seaborgium", false)) || \
        	 ((useNumber) && (Text::equal(value, "106"))))
		return 106;
    
    // Bohrium
	else if ((Text::equal(value, "Bh", true)) || (Text::equal(value, "bohrium", false)) || \
        	 ((useNumber) && (Text::equal(value, "107"))))
		return 107;
    
    // Hassium
	else if ((Text::equal(value, "Hs", true)) || (Text::equal(value, "hassium", false)) || \
        	 ((useNumber) && (Text::equal(value, "108"))))
		return 108;
    
    // Meitnerium
	else if ((Text::equal(value, "Mt", true)) || (Text::equal(value, "meitnerium", false)) || \
        	 ((useNumber) && (Text::equal(value, "109"))))
		return 109;
	
	// Unknown
	return 0;
}


/* Element Element::getByNumber(int num)
 *
 * Return the element with a given number
 */

Element Element::getByNumber(int num)
{	
	switch(num)
	{
		case 1:
			return Element(1, 0, "H", "Hydrogen", 1.00794, 0.31);
		case 1001:
			return Element(1, 1, "D", "Deuterium", 2.013553212724, 0.5);
		case 1002:
			return Element(1, 2, "T", "Tritium", 3.0160492, 1.0);
		case 2:
			return Element(2, 0, "He", "Helium", 4.002602, 0.28);
		case 3:
			return Element(3, 0, "Li", "Lithium", 6.94, 1.28);
		case 4:
			return Element(4, 0, "Be", "Beryllium", 9.012182, 0.96);
		case 5:
			return Element(5, 0, "B", "Boron", 10.81, 0.84);
		case 6:
			return Element(6, 0, "C", "Carbon", 12.011, 0.76);
		case 7:
			return Element(7, 0, "N", "Nitrogen", 14.007, 0.71);
		case 8:
			return Element(8, 0, "O", "Oxygen", 15.999, 0.66);
		case 9:
			return Element(9, 0, "F", "Fluorine", 18.9984032, 0.57);
		case 10:
			return Element(10, 0, "Ne", "Neon", 20.1797, 0.58);
		case 11:
			return Element(11, 0, "Na", "Sodium", 22.98976928, 1.67);
		case 12:
			return Element(12, 0, "Mg", "Magnesium", 24.3050, 1.42);
		case 13:
			return Element(13, 0, "Al", "Aluminum", 26.9815386, 1.21);
		case 14:
			return Element(14, 0, "Si", "Silicon", 28.085, 1.11);
		case 15:
			return Element(15, 0, "P", "Phosphorus", 30.973762, 1.07);
		case 16:
			return Element(16, 0, "S", "Sulfur", 32.06, 1.05);
		case 17:
			return Element(17, 0, "Cl", "Chlorine", 35.45, 1.02);
		case 18:
			return Element(18, 0, "Ar", "Argon", 39.948, 1.06);
		case 19:
			return Element(19, 0, "K", "Potassium", 39.0983, 2.03);
		case 20:
			return Element(20, 0, "Ca", "Calcium", 40.078, 1.76);
		case 21:
			return Element(21, 0, "Sc", "Scandium", 44.955912, 1.70);
		case 22:
			return Element(22, 0, "Ti", "Titanium", 47.867, 1.61);
		case 23:
			return Element(23, 0, "V", "Vanadium", 50.9415, 1.54);
		case 24:
			return Element(24, 0, "Cr", "Chromium", 51.9961, 1.40);
		case 25:
			return Element(25, 0, "Mn", "Manganese", 54.938045, 1.40);
		case 26:
			return Element(26, 0, "Fe", "Iron", 55.845, 1.32);
		case 27:
			return Element(27, 0, "Co", "Cobalt", 58.933195, 1.26);
		case 28:
			return Element(28, 0, "Ni", "Nickel", 58.6934, 1.24);
		case 29:
			return Element(29, 0, "Cu", "Copper", 63.546, 1.32);
		case 30:
			return Element(30, 0, "Zn", "Zinc", 65.38, 1.22);
		case 31:
			return Element(31, 0, "Ga", "Gallium", 69.723, 1.22);
		case 32:
			return Element(32, 0, "Ge", "Germanium", 72.63, 1.20);
		case 33:
			return Element(33, 0, "As", "Arsenic", 74.92160, 1.19);
		case 34:
			return Element(34, 0, "Se", "Selenium", 78.96, 1.20);
		case 35:
			return Element(35, 0, "Br", "Bromine", 79.904, 1.20);
		case 36:
			return Element(36, 0, "Kr", "Krypton", 83.798, 1.16);
		case 37:
			return Element(37, 0, "Rb", "Rubidium", 85.4678, 2.20);
		case 38:
			return Element(38, 0, "Sr", "Strontium", 87.62, 1.95);
		case 39:
			return Element(39, 0, "Y", "Yttrium", 88.90585, 1.90);
		case 40:
			return Element(40, 0, "Zr", "Zirconium", 91.224, 1.75);
		case 41:
			return Element(41, 0, "Nb", "Niobium", 92.90638, 1.65);
		case 42:
			return Element(42, 0, "Mo", "Molybdenum", 95.96, 1.55);
		case 43:
			return Element(43, 0, "Tc", "Technetium", 98.0, 1.48);
		case 44:
			return Element(44, 0, "Ru", "Ruthenium", 101.07, 1.47);
		case 45:
			return Element(45, 0, "Rh", "Rhodium", 102.90550, 1.43);
		case 46:
			return Element(46, 0, "Pd", "Palladium", 106.42, 1.40);
		case 47:
			return Element(47, 0, "Ag", "Silver", 107.8682, 1.46);
		case 48:
			return Element(48, 0, "Cd", "Cadmium", 112.411, 1.45);
		case 49:
			return Element(49, 0, "In", "Indium", 114.818, 1.43);
		case 50:
			return Element(50, 0, "Sn", "Tin", 118.710, 1.39);
		case 51:
			return Element(51, 0, "Sb", "Antimony", 121.760, 1.40);
		case 52:
			return Element(52, 0, "Te", "Tellurium", 127.60, 1.38);
		case 53:
			return Element(53, 0, "I", "Iodine", 126.90447, 1.39);
		case 54:
			return Element(54, 0, "Xe", "Xenon", 131.293, 1.41);
		case 55:
			return Element(55, 0, "Cs", "Cesium", 132.9054519, 2.44);
		case 56:
			return Element(56, 0, "Ba", "Barium", 137.327, 2.15);
		case 57:
			return Element(57, 0, "La", "Lanthanum", 138.90547, 2.08);
		case 58:
			return Element(58, 0, "Ce", "Cerium", 140.116, 2.05);
		case 59:
			return Element(59, 0, "Pr", "Praseodymium", 140.90765, 2.04);
		case 60:
			return Element(60, 0, "Nd", "Neodymium", 144.242, 2.02);
		case 61:
			return Element(61, 0, "Pm", "Promethium", 145.0, 1.99);
		case 62:
			return Element(62, 0, "Sm", "Samarium", 150.36, 1.99);
		case 63:
			return Element(63, 0, "Eu", "Europium", 151.964, 1.99);
		case 64:
			return Element(64, 0, "Gd", "Gadolinium", 157.25, 1.97);
		case 65:
			return Element(65, 0, "Tb", "Terbium", 158.92535, 1.95);
		case 66:
			return Element(66, 0, "Dy", "Dysprosium", 162.500, 1.93);
		case 67:
			return Element(67, 0, "Ho", "Holmium", 164.93032, 1.93);
		case 68:
			return Element(68, 0, "Er", "Erbium", 167.259, 1.90);
		case 69:
			return Element(69, 0, "Tm", "Thulium", 168.93421, 1.90);
		case 70:
			return Element(70, 0, "Yb", "Ytterbium", 173.054, 1.88);
		case 71:
			return Element(71, 0, "Lu", "Lutetium", 174.9668, 1.88);
		case 72:
			return Element(72, 0, "Hf", "Hafnium", 178.49, 1.75);
		case 73:
			return Element(73, 0, "Ta", "Tantalum", 180.94788, 1.71);
		case 74:
			return Element(74, 0, "W", "Tungsten", 183.84, 1.63);
		case 75:
			return Element(75, 0, "Re", "Rhenium", 186.207, 1.52);
		case 76:
			return Element(76, 0, "Os", "Osmium", 190.23, 1.44);
		case 77:
			return Element(77, 0, "Ir", "Iridium", 192.217, 1.42);
		case 78:
			return Element(78, 0, "Pt", "Platinum", 195.084, 1.37);
		case 79:
			return Element(79, 0, "Au", "Gold", 196.966569, 1.37);
		case 80:
			return Element(80, 0, "Hg", "Mercury", 200.59, 1.33);
		case 81:
			return Element(81, 0, "Tl", "Thallium", 204.38, 1.46);
		case 82:
			return Element(82, 0, "Pb", "Lead", 207.2, 1.47);
		case 83:
			return Element(83, 0, "Bi", "Bismuth", 208.98040, 1.48);
		case 84:
			return Element(84, 0, "Po", "Polonium", 209.0, 1.40);
		case 85:
			return Element(85, 0, "At", "Astatine", 210.0, 1.50);
		case 86:
			return Element(86, 0, "Rn", "Radon", 222.0, 1.50);
		case 87:
			return Element(87, 0, "Fr", "Francium", 223.0, 2.60);
		case 88:
			return Element(88, 0, "Ra", "Radium", 226.0, 2.21);
		case 89:
			return Element(89, 0, "Ac", "Actinium", 227.0, 2.15);
		case 90:
			return Element(90, 0, "Th", "Thorium", 232.03806, 2.07);
		case 91:
			return Element(91, 0, "Pa", "Protactinium", 231.03588, 2.00);
		case 92:
			return Element(92, 0, "U", "Uranium", 238.02891, 1.97);
		case 93:
			return Element(93, 0, "Np", "Neptunium", 237.0, 1.90);
		case 94:
			return Element(94, 0, "Pu", "Plutonium", 244.0, 1.87);
		case 95:
			return Element(95, 0, "Am", "Americium", 243.0, 1.81);
		case 96:
			return Element(96, 0, "Cm", "Curium", 247.0, 1.69);
		case 97:
			return Element(97, 0, "Bk", "Berkelium", 247.0, 1.5);
		case 98:
			return Element(98, 0, "Cf", "Californium", 251.0, 1.5);
		case 99:
			return Element(99, 0, "Es", "Einsteinium", 252.0, 1.5);
		case 100:
			return Element(100, 0, "Fm", "Fermium", 257.0, 1.5);
		case 101:
			return Element(101, 0, "Md", "Mendelevium", 258.0, 1.5);
		case 102:
			return Element(102, 0, "No", "Nobelium", 259.0, 1.5);
		case 103:
			return Element(103, 0, "Lr", "Lawrencium", 262.0, 1.5);
		case 104:
			return Element(104, 0, "Rf", "Rutherfordium", 267.0, 1.5);
		case 105:
			return Element(105, 0, "Db", "Dubnium", 268.0, 1.4);
		case 106:
			return Element(106, 0, "Sg", "Seaborgium", 269.0, 1.3);
		case 107:
			return Element(107, 0, "Bh", "Bohrium", 270.0, 1.3);
		case 108:
			return Element(108, 0, "Hs", "Hassium", 269.0, 1.25);
		case 109:
			return Element(109, 0, "Mt", "Meitnerium", 278.0, 1.22);
		default:
			return Element(0, 0, "", "", 0, 0);
	}
}
