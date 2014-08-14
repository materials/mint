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



#include "spaceGroup.h"
#include "language.h"
#include "output.h"



// Wyckoff letters
char Wyckoff::_letters[35][2] = {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t", \
								 "u","v","w","x","y","z","A","B","C","D","E","F","G","H","I"};



/* void Wyckoff::clear()
 *
 * Clear data in Wyckoff object
 */

void Wyckoff::clear()
{
	_rotations.clear();
	_translations.clear();
	_letter[0] = '\0';
	_name.clear();
}



/* Wyckoff& Wyckoff::operator= (const Wyckoff& rhs)
 *
 * Copy Wyckoff object
 */

Wyckoff& Wyckoff::operator= (const Wyckoff& rhs)
{
	if (this != &rhs)
	{
		_rank = rhs._rank;
		_rotations = rhs._rotations;
		_translations = rhs._translations;
		_letter[0] = rhs._letter[0];
		_letter[1] = rhs._letter[1];
		_name = rhs._name;
	}
	return *this;
}



/* void Wyckoff::set(int rotIndex, int transIndex, const OList<SymmetryOperation>& symmetry, int letterIndex)
 *
 * Set the properties of a Wyckoff position
 */

void Wyckoff::set(int rotIndex, int transIndex, const OList<SymmetryOperation>& symmetry, int letterIndex)
{

	// Clear space
	clear();
	
	// Look up rotation and translation
	Matrix3D curMat;
	Vector3D curVec;
	if (rotIndex == 0) curMat.set(1, 0, 0, 0, 1, 0, 0, 0, 1);
	else if (rotIndex == 1) curMat.set(0, 0, 0, 0, 0, 0, 0, 0, 0);
	else if (rotIndex == 2) curMat.set(0, 0, 0, 0, 1, 0, 0, 0, 0);
	else if (rotIndex == 3) curMat.set(1, 0, 0, 0, 0, 0, 0, 0, 1);
	else if (rotIndex == 4) curMat.set(0, 0, 0, 0, 0, 0, 0, 0, 1);
	else if (rotIndex == 5) curMat.set(1, 0, 0, 0, 0, 0, 0, 0, 0);
	else if (rotIndex == 6) curMat.set(0, 0, 0, 0, 1, 0, 0, 0, 1);
	else if (rotIndex == 7) curMat.set(1, 0, 0, 0, 1, 0, 0, 0, 0);
	else if (rotIndex == 8) curMat.set(1, 0, 0, 1, 0, 0, 0, 0, 0);
	else if (rotIndex == 9) curMat.set(1, 0, 0, -1, 0, 0, 0, 0, 0);
	else if (rotIndex == 10) curMat.set(1, 0, 0, 1, 0, 0, 0, 0, 1);
	else if (rotIndex == 11) curMat.set(1, 0, 0, -1, 0, 0, 0, 0, 1);
	else if (rotIndex == 12) curMat.set(1, 0, 0, 2, 0, 0, 0, 0, 0);
	else if (rotIndex == 13) curMat.set(1, 0, 0, 2, 0, 0, 0, 0, 1);
	else if (rotIndex == 14) curMat.set(1, 0, 0, 1, 0, 0, 1, 0, 0);
	else if (rotIndex == 15) curMat.set(0, 0, 0, 0, 1, 0, 0, 1, 0);
	else if (rotIndex == 16) curMat.set(0, 0, 0, 0, 1, 0, 0, -1, 0);
	else
	{
		Output::newline(ERROR);
		Output::print("Internal error: unknown Wyckoff rotation index");
		Output::quit();
	}
	if (transIndex == 0) curVec.set(0.0, 0.0, 0.0);
	else if (transIndex == 1) curVec.set(1.0/2.0, 1.0/2.0, 1.0/2.0);
	else if (transIndex == 2) curVec.set(0.0, 1.0/2.0, 1.0/2.0);
	else if (transIndex == 3) curVec.set(1.0/2.0, 0.0, 1.0/2.0);
	else if (transIndex == 4) curVec.set(1.0/2.0, 1.0/2.0, 0.0);
	else if (transIndex == 5) curVec.set(1.0/2.0, 0.0, 0.0);
	else if (transIndex == 6) curVec.set(0.0, 1.0/2.0, 0.0);
	else if (transIndex == 7) curVec.set(0.0, 0.0, 1.0/2.0);
	else if (transIndex == 8) curVec.set(0.0, 1.0/4.0, 0.0);
	else if (transIndex == 9) curVec.set(1.0/4.0, 1.0/4.0, 1.0/2.0);
	else if (transIndex == 10) curVec.set(1.0/4.0, 1.0/4.0, 0.0);
	else if (transIndex == 11) curVec.set(1.0/2.0, 0.0, 1.0/4.0);
	else if (transIndex == 12) curVec.set(0.0, 0.0, 1.0/4.0);
	else if (transIndex == 13) curVec.set(0.0, 1.0/4.0, 1.0/4.0);
	else if (transIndex == 14) curVec.set(1.0/4.0, 0.0, 1.0/4.0);
	else if (transIndex == 15) curVec.set(1.0/4.0, 1.0/4.0, 3.0/4.0);
	else if (transIndex == 16) curVec.set(1.0/4.0, 1.0/4.0, 1.0/4.0);
	else if (transIndex == 17) curVec.set(1.0/4.0, 0.0, 0.0);
	else if (transIndex == 18) curVec.set(3.0/4.0, 3.0/4.0, 3.0/4.0);
	else if (transIndex == 19) curVec.set(0.0, 1.0/2.0, 1.0/4.0);
	else if (transIndex == 20) curVec.set(1.0/2.0, 1.0/2.0, 1.0/4.0);
	else if (transIndex == 21) curVec.set(1.0/4.0, 1.0/2.0, 0.0);
	else if (transIndex == 22) curVec.set(1.0/4.0, 3.0/4.0, 0.0);
	else if (transIndex == 23) curVec.set(1.0/4.0, 0.0, 1.0/2.0);
	else if (transIndex == 24) curVec.set(5.0/8.0, 5.0/8.0, 5.0/8.0);
	else if (transIndex == 25) curVec.set(1.0/8.0, 1.0/8.0, 1.0/8.0);
	else if (transIndex == 26) curVec.set(0.0, 1.0/2.0, 3.0/4.0);
	else if (transIndex == 27) curVec.set(0.0, 1.0/4.0, 5.0/8.0);
	else if (transIndex == 28) curVec.set(0.0, 1.0/4.0, 1.0/8.0);
	else if (transIndex == 29) curVec.set(0.0, 0.0, 3.0/8.0);
	else if (transIndex == 30) curVec.set(0.0, 0.0, 3.0/4.0);
	else if (transIndex == 31) curVec.set(0.0, 0.0, 5.0/8.0);
	else if (transIndex == 32) curVec.set(1.0/4.0, 0.0, 1.0/8.0);
	else if (transIndex == 33) curVec.set(2.0/3.0, 1.0/3.0, 0.0);
	else if (transIndex == 34) curVec.set(1.0/3.0, 2.0/3.0, 0.0);
	else if (transIndex == 35) curVec.set(2.0/3.0, 1.0/3.0, 1.0/2.0);
	else if (transIndex == 36) curVec.set(1.0/3.0, 2.0/3.0, 1.0/2.0);
	else if (transIndex == 37) curVec.set(0.0, 0.0, 5.0/6.0);
	else if (transIndex == 38) curVec.set(0.0, 0.0, 1.0/3.0);
	else if (transIndex == 39) curVec.set(0.0, 0.0, 1.0/6.0);
	else if (transIndex == 40) curVec.set(0.0, 0.0, 2.0/3.0);
	else if (transIndex == 41) curVec.set(2.0/3.0, 1.0/3.0, 1.0/4.0);
	else if (transIndex == 42) curVec.set(1.0/3.0, 2.0/3.0, 1.0/4.0);
	else if (transIndex == 43) curVec.set(1.0/3.0, 2.0/3.0, 3.0/4.0);
	else if (transIndex == 44) curVec.set(1.0/8.0, 0.0, 1.0/4.0);
	else if (transIndex == 45) curVec.set(7.0/8.0, 7.0/8.0, 7.0/8.0);
	else if (transIndex == 46) curVec.set(3.0/8.0, 3.0/8.0, 3.0/8.0);
	else if (transIndex == 47) curVec.set(5.0/8.0, 0.0, 1.0/4.0);
	else if (transIndex == 48) curVec.set(7.0/8.0, 0.0, 1.0/4.0);
	else if (transIndex == 49) curVec.set(3.0/8.0, 0.0, 1.0/4.0);
	else
	{
		Output::newline(ERROR);
		Output::print("Internal error: unknown Wyckoff translation index");
		Output::quit();
	}
	
	// Save rank of position
	_rank = curMat.rank();
	
	// Loop over symmetries and add points in orbit
	int i, j, k;
	bool found;
	Matrix3D newMat;
	Vector3D newVec;
	for (i = 0; i < symmetry.length(); i++)
	{
		
		// Save rotation and loop over translations
		newMat = symmetry[i].rotation() * curMat;
		for (j = 0; j < symmetry[i].translations().length(); ++j)
		{
			
			// Save vector and check if position is known
			found = false;
			newVec = symmetry[i].rotation() * curVec + symmetry[i].translations()[j];
			for (k = 0; k < 3; ++k)
				newVec[k] = Num<double>::mod(newVec[k], 0.9999);
			for (k = 0; k < _rotations.length(); ++k)
			{
				
				// Found match
				if ((newMat == _rotations[k]) && (newVec == _translations[k]))
				{
					found = true;
					break;
				}
			}
			
			// Save point if new
			if (!found)
			{
				_rotations += newMat;
				_translations += newVec;
			}
		}
	}
	
	// Set the letter and name of the position
	_letter[0] = _letters[letterIndex][0];
	_letter[1] = '\0';
	_name = Language::numberToWord(_rotations.length()) + _letter;
}



/* void Wyckoff::print() const
 *
 * Print information about Wyckoff position
 */

void Wyckoff::print() const
{
	
	// Save current stream
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print title line
	Output::newline();
	Output::print("Orbit for Wyckoff position ");
	Output::print(_name);
	
	// Convert points to Jones Faithful representation
	int i, j;
	Words tempWords;
	Output message;
	message.addLines(_rotations.length());
	for (i = 0; i < _rotations.length(); ++i)
	{
		message.addLine();
		message.addWords(4);
		message.add("   ");
		tempWords = JonesFaithful::toString(_rotations[i], &_translations[i]);
		for (j = 0; j < 3; ++j)
			message.add(tempWords[j]);
	}
	
	// Print orbit
	Output::newline();
	Output::print(message, LEFT);
	
	// Reset method
	Output::method(origMethod);
}



/* void SpaceGroup::clear()
 *
 * Clear space group information
 */

void SpaceGroup::clear()
{
	_itcNumber.clear();
	_schoenflies.clear();
	_hermannMauguinShort.clear();
	_hermannMauguin.clear();
	_hall.clear();
	_symmetry.clear();
	_originShift = 0.0;
	_wyckoff.clear();
	_pointGroup.clear();
}



/* SpaceGroup& SpaceGroup::operator= (const SpaceGroup& rhs)
 *
 * Copy space group information
 */

SpaceGroup& SpaceGroup::operator= (const SpaceGroup& rhs)
{
	if (this != &rhs)
	{
		_itcNumber = rhs._itcNumber;
		_schoenflies = rhs._schoenflies;
		_hermannMauguinShort = rhs._hermannMauguinShort;
		_hermannMauguin = rhs._hermannMauguin;
		_hall = rhs._hall;
		_system = rhs._system;
		_centering = rhs._centering;
		_symmetry = rhs._symmetry;
		_originShift = rhs._originShift;
		_wyckoff = rhs._wyckoff;
		_pointGroup = rhs._pointGroup;
	}
	return *this;
}



/* bool SpaceGroup::set(const Word& name, bool useNumber, bool quitIfNotFound)
 *
 * Set space group by symbol
 */

bool SpaceGroup::set(const Word& name, bool useNumber, bool quitIfNotFound)
{
	
	// Clear space
	clear();
	
	// Get space group number, error if not recognized
	int number = getNumberFromName(name, useNumber);
	if (number == -1)
	{
		
		// Not quitting
		if (!quitIfNotFound)
			return false;

		// Quit
		Output::newline(ERROR);
		Output::print("Did not recognize space group");
		Output::quit();
	}
	
	// Set properties
	setByNumber(number);
	
	// Output
	Output::newline();
	Output::print("Setting space group to ");
	Output::print(name);
	Output::increase();

	// Print the result
	Output::newline();
	Output::print("ITC number: ");
	Output::print(_itcNumber);
	Output::newline();
	Output::print("Hermann-Mauguin short: ");
	Output::print(_hermannMauguinShort);
	Output::newline();
	Output::print("Hermann-Mauguin full: ");
	Output::print(_hermannMauguin);
	Output::newline();
	Output::print("Schoenflies: ");
	Output::print(_schoenflies);
	Output::newline();
	Output::print("Hall: ");
	Output::print(_hall);
	
	// Print symmetry
	int i;
	Output::newline();
	Output::print("Generated ");
	Output::print(_symmetry.length());
	Output::print(" unique symmetry operation");
	if (_symmetry.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _symmetry.length(); ++i)
	{
		Output::newline();
		Output::print(_symmetry[i].getString(), true, false);
	}
	Output::decrease();
	
	// Print centering vectors
	int j;
	Output::newline();
	Output::print("Number of centering vectors: ");
	if (_symmetry.length())
	{
		Output::print(_symmetry[0].translations().length());
		Output::increase();
		for (i = 0; i < _symmetry[0].translations().length(); ++i)
		{
			Output::newline();
			for (j = 0; j < 3; ++j)
			{
				Output::print(Language::numberToFraction(_symmetry[0].translations()[i][j]));
				Output::print(" ");
			}
		}
		Output::decrease();
	}
	else
		Output::print("0");

	// Print the Wyckoff positions
	Output::newline();
	Output::print("Generated ");
	Output::print(_wyckoff.length());
	Output::print(" Wyckoff position");
	if (_wyckoff.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _wyckoff.length(); i++)
	{

		// Print title line
		Output::newline();
		Output::print("Multiplicity of Wyckoff position ");
		Output::print(_wyckoff[i].letter());
		Output::print(": ");
		Output::print(_wyckoff[i].rotations().length());
		Output::increase();

		// Print operations
		for (j = 0; j < _wyckoff[i].rotations().length(); j++)
		{
			Output::newline();
			Output::print(_wyckoff[i].getString(j), true, false);
		}
		Output::decrease();
	}
	Output::decrease();

	// Output
	Output::decrease();

	// Return that space group was found
	return true;
}



/* int SpaceGroup::getNumberFromName(const Word& name, bool useNumber)
 *
 * Get space group number from name
 */

int SpaceGroup::getNumberFromName(const Word& name, bool useNumber)
{
	
	// Compare name to all space groups
	if (compareName(name,   "1",   "C1^1",      "P1",              "P 1",              "P 1", useNumber)) return   0;
	if (compareName(name,   "2",   "Ci^1",     "P-1",             "P -1",             "-P 1", useNumber)) return   1;
	if (compareName(name,   "3",   "C2^1",      "P2",          "P 1 2 1",             "P 2y", useNumber)) return   2;
	if (compareName(name,   "4",   "C2^2",     "P21",         "P 1 21 1",            "P 2yb", useNumber)) return   3;
	if (compareName(name,   "5",   "C2^3",      "C2",          "C 1 2 1",             "C 2y", useNumber)) return   4;
	if (compareName(name,   "6",   "Cs^1",      "Pm",          "P 1 m 1",            "P -2y", useNumber)) return   5;
	if (compareName(name,   "7",   "Cs^2",      "Pc",          "P 1 c 1",           "P -2yc", useNumber)) return   6;
	if (compareName(name,   "8",   "Cs^3",      "Cm",          "C 1 m 1",            "C -2y", useNumber)) return   7;
	if (compareName(name,   "9",   "Cs^4",      "Cc",          "C 1 c 1",           "C -2yc", useNumber)) return   8;
	if (compareName(name,  "10",  "C2h^1",    "P2/m",        "P 1 2/m 1",            "-P 2y", useNumber)) return   9;
	if (compareName(name,  "11",  "C2h^2",   "P21/m",       "P 1 21/m 1",           "-P 2yb", useNumber)) return  10;
	if (compareName(name,  "12",  "C2h^3",    "C2/m",        "C 1 2/m 1",            "-C 2y", useNumber)) return  11;
	if (compareName(name,  "13",  "C2h^4",    "P2/c",        "P 1 2/c 1",           "-P 2yc", useNumber)) return  12;
	if (compareName(name,  "14",  "C2h^5",   "P21/c",       "P 1 21/c 1",          "-P 2ybc", useNumber)) return  13;
	if (compareName(name,  "15",  "C2h^6",    "C2/c",        "C 1 2/c 1",           "-C 2yc", useNumber)) return  14;
	if (compareName(name,  "16",   "D2^1",    "P222",          "P 2 2 2",            "P 2 2", useNumber)) return  15;
	if (compareName(name,  "17",   "D2^2",   "P2221",         "P 2 2 21",           "P 2c 2", useNumber)) return  16;
	if (compareName(name,  "18",   "D2^3",  "P21212",        "P 21 21 2",          "P 2 2ab", useNumber)) return  17;
	if (compareName(name,  "19",   "D2^4", "P212121",       "P 21 21 21",        "P 2ac 2ab", useNumber)) return  18;
	if (compareName(name,  "20",   "D2^5",   "C2221",         "C 2 2 21",           "C 2c 2", useNumber)) return  19;
	if (compareName(name,  "21",   "D2^6",    "C222",          "C 2 2 2",            "C 2 2", useNumber)) return  20;
	if (compareName(name,  "22",   "D2^7",    "F222",          "F 2 2 2",            "F 2 2", useNumber)) return  21;
	if (compareName(name,  "23",   "D2^8",    "I222",          "I 2 2 2",            "I 2 2", useNumber)) return  22;
	if (compareName(name,  "24",   "D2^9", "I212121",       "I 21 21 21",          "I 2b 2c", useNumber)) return  23;
	if (compareName(name,  "25",  "C2v^1",    "Pmm2",          "P m m 2",           "P 2 -2", useNumber)) return  24;
	if (compareName(name,  "26",  "C2v^2",   "Pmc21",         "P m c 21",          "P 2c -2", useNumber)) return  25;
	if (compareName(name,  "27",  "C2v^3",    "Pcc2",          "P c c 2",          "P 2 -2c", useNumber)) return  26;
	if (compareName(name,  "28",  "C2v^4",    "Pma2",          "P m a 2",          "P 2 -2a", useNumber)) return  27;
	if (compareName(name,  "29",  "C2v^5",   "Pca21",         "P c a 21",        "P 2c -2ac", useNumber)) return  28;
	if (compareName(name,  "30",  "C2v^6",    "Pnc2",          "P n c 2",         "P 2 -2bc", useNumber)) return  29;
	if (compareName(name,  "31",  "C2v^7",   "Pmn21",         "P m n 21",         "P 2ac -2", useNumber)) return  30;
	if (compareName(name,  "32",  "C2v^8",    "Pba2",          "P b a 2",         "P 2 -2ab", useNumber)) return  31;
	if (compareName(name,  "33",  "C2v^9",   "Pna21",         "P n a 21",         "P 2c -2n", useNumber)) return  32;
	if (compareName(name,  "34", "C2v^10",    "Pnn2",          "P n n 2",          "P 2 -2n", useNumber)) return  33;
	if (compareName(name,  "35", "C2v^11",    "Cmm2",          "C m m 2",           "C 2 -2", useNumber)) return  34;
	if (compareName(name,  "36", "C2v^12",   "Cmc21",         "C m c 21",          "C 2c -2", useNumber)) return  35;
	if (compareName(name,  "37", "C2v^13",    "Ccc2",          "C c c 2",          "C 2 -2c", useNumber)) return  36;
	if (compareName(name,  "38", "C2v^14",    "Amm2",          "A m m 2",           "A 2 -2", useNumber)) return  37;
	if (compareName(name,  "39", "C2v^15",    "Abm2",          "A b m 2",          "A 2 -2c", useNumber)) return  38;
	if (compareName(name,  "40", "C2v^16",    "Ama2",          "A m a 2",          "A 2 -2a", useNumber)) return  39;
	if (compareName(name,  "41", "C2v^17",    "Aba2",          "A b a 2",         "A 2 -2ac", useNumber)) return  40;
	if (compareName(name,  "42", "C2v^18",    "Fmm2",          "F m m 2",           "F 2 -2", useNumber)) return  41;
	if (compareName(name,  "43", "C2v^19",    "Fdd2",          "F d d 2",          "F 2 -2d", useNumber)) return  42;
	if (compareName(name,  "44", "C2v^20",    "Imm2",          "I m m 2",           "I 2 -2", useNumber)) return  43;
	if (compareName(name,  "45", "C2v^21",    "Iba2",          "I b a 2",          "I 2 -2c", useNumber)) return  44;
	if (compareName(name,  "46", "C2v^22",    "Ima2",          "I m a 2",          "I 2 -2a", useNumber)) return  45;
	if (compareName(name,  "47",  "D2h^1",    "Pmmm",    "P 2/m 2/m 2/m",           "-P 2 2", useNumber)) return  46;
	if (compareName(name,  "48",  "D2h^2",    "Pnnn",    "P 2/n 2/n 2/n",        "P 2 2 -1n", useNumber)) return  47;
	if (compareName(name,  "49",  "D2h^3",    "Pccm",    "P 2/c 2/c 2/m",          "-P 2 2c", useNumber)) return  48;
	if (compareName(name,  "50",  "D2h^4",    "Pban",    "P 2/b 2/a 2/n",       "P 2 2 -1ab", useNumber)) return  49;
	if (compareName(name,  "51",  "D2h^5",    "Pmma",   "P 21/m 2/m 2/a",         "-P 2a 2a", useNumber)) return  50;
	if (compareName(name,  "52",  "D2h^6",    "Pnna",   "P 2/n 21/n 2/a",        "-P 2a 2bc", useNumber)) return  51;
	if (compareName(name,  "53",  "D2h^7",    "Pmna",   "P 2/m 2/n 21/a",         "-P 2ac 2", useNumber)) return  52;
	if (compareName(name,  "54",  "D2h^8",    "Pcca",   "P 21/c 2/c 2/a",        "-P 2a 2ac", useNumber)) return  53;
	if (compareName(name,  "55",  "D2h^9",    "Pbam",  "P 21/b 21/a 2/m",         "-P 2 2ab", useNumber)) return  54;
	if (compareName(name,  "56", "D2h^10",    "Pccn",  "P 21/c 21/c 2/n",       "-P 2ab 2ac", useNumber)) return  55;
	if (compareName(name,  "57", "D2h^11",    "Pbcm",  "P 2/b 21/c 21/m",         "-P 2c 2b", useNumber)) return  56;
	if (compareName(name,  "58", "D2h^12",    "Pnnm",  "P 21/n 21/n 2/m",          "-P 2 2n", useNumber)) return  57;
	if (compareName(name,  "59", "D2h^13",    "Pmmn",  "P 21/m 21/m 2/n",     "P 2 2ab -1ab", useNumber)) return  58;
	if (compareName(name,  "60", "D2h^14",    "Pbcn",  "P 21/b 2/c 21/n",        "-P 2n 2ab", useNumber)) return  59;
	if (compareName(name,  "61", "D2h^15",    "Pbca", "P 21/b 21/c 21/a",       "-P 2ac 2ab", useNumber)) return  60;
	if (compareName(name,  "62", "D2h^16",    "Pnma", "P 21/n 21/m 21/a",        "-P 2ac 2n", useNumber)) return  61;
	if (compareName(name,  "63", "D2h^17",    "Cmcm",   "C 2/m 2/c 21/m",          "-C 2c 2", useNumber)) return  62;
	if (compareName(name,  "64", "D2h^18",    "Cmca",   "C 2/m 2/c 21/a",         "-C 2bc 2", useNumber)) return  63;
	if (compareName(name,  "65", "D2h^19",    "Cmmm",    "C 2/m 2/m 2/m",           "-C 2 2", useNumber)) return  64;
	if (compareName(name,  "66", "D2h^20",    "Cccm",    "C 2/c 2/c 2/m",          "-C 2 2c", useNumber)) return  65;
	if (compareName(name,  "67", "D2h^21",    "Cmma",    "C 2/m 2/m 2/a",          "-C 2b 2", useNumber)) return  66;
	if (compareName(name,  "68", "D2h^22",    "Ccca",    "C 2/c 2/c 2/a",       "C 2 2 -1bc", useNumber)) return  67;
	if (compareName(name,  "69", "D2h^23",    "Fmmm",    "F 2/m 2/m 2/m",           "-F 2 2", useNumber)) return  68;
	if (compareName(name,  "70", "D2h^24",    "Fddd",    "F 2/d 2/d 2/d",        "F 2 2 -1d", useNumber)) return  69;
	if (compareName(name,  "71", "D2h^25",    "Immm",    "I 2/m 2/m 2/m",           "-I 2 2", useNumber)) return  70;
	if (compareName(name,  "72", "D2h^26",    "Ibam",    "I 2/b 2/a 2/m",          "-I 2 2c", useNumber)) return  71;
	if (compareName(name,  "73", "D2h^27",    "Ibca",    "I 2/b 2/c 2/a",         "-I 2b 2c", useNumber)) return  72;
	if (compareName(name,  "74", "D2h^28",    "Imma",    "I 2/m 2/m 2/a",          "-I 2b 2", useNumber)) return  73;
	if (compareName(name,  "75",   "C4^1",      "P4",          "P 4 1 1",              "P 4", useNumber)) return  74;
	if (compareName(name,  "76",   "C4^2",     "P41",         "P 41 1 1",             "P 4w", useNumber)) return  75;
	if (compareName(name,  "77",   "C4^3",     "P42",         "P 42 1 1",             "P 4c", useNumber)) return  76;
	if (compareName(name,  "78",   "C4^4",     "P43",         "P 43 1 1",            "P 4cw", useNumber)) return  77;
	if (compareName(name,  "79",   "C4^5",      "I4",          "I 4 1 1",              "I 4", useNumber)) return  78;
	if (compareName(name,  "80",   "C4^6",     "I41",         "I 41 1 1",            "I 4bw", useNumber)) return  79;
	if (compareName(name,  "81",   "S4^1",     "P-4",         "P -4 1 1",             "P -4", useNumber)) return  80;
	if (compareName(name,  "82",   "S4^2",     "I-4",         "I -4 1 1",             "I -4", useNumber)) return  81;
	if (compareName(name,  "83",  "C4h^1",    "P4/m",        "P 4/m 1 1",             "-P 4", useNumber)) return  82;
	if (compareName(name,  "84",  "C4h^2",   "P42/m",       "P 42/m 1 1",            "-P 4c", useNumber)) return  83;
	if (compareName(name,  "85",  "C4h^3",    "P4/n",        "P 4/n 1 1",       "P 4ab -1ab", useNumber)) return  84;
	if (compareName(name,  "86",  "C4h^4",   "P42/n",       "P 42/n 1 1",         "P 4n -1n", useNumber)) return  85;
	if (compareName(name,  "87",  "C4h^5",    "I4/m",        "I 4/m 1 1",             "-I 4", useNumber)) return  86;
	if (compareName(name,  "88",  "C4h^6",   "I41/a",       "I 41/a 1 1",       "I 4bw -1bw", useNumber)) return  87;
	if (compareName(name,  "89",   "D4^1",    "P422",          "P 4 2 2",            "P 4 2", useNumber)) return  88;
	if (compareName(name,  "90",   "D4^2",   "P4212",         "P 4 21 2",        "P 4ab 2ab", useNumber)) return  89;
	if (compareName(name,  "91",   "D4^3",   "P4122",         "P 41 2 2",          "P 4w 2c", useNumber)) return  90;
	if (compareName(name,  "92",   "D4^4",  "P41212",        "P 41 21 2",       "P 4abw 2nw", useNumber)) return  91;
	if (compareName(name,  "93",   "D4^5",   "P4222",         "P 42 2 2",           "P 4c 2", useNumber)) return  92;
	if (compareName(name,  "94",   "D4^6",  "P42212",        "P 42 21 2",          "P 4n 2n", useNumber)) return  93;
	if (compareName(name,  "95",   "D4^7",   "P4322",         "P 43 2 2",         "P 4cw 2c", useNumber)) return  94;
	if (compareName(name,  "96",   "D4^8",  "P43212",        "P 43 21 2",       "P 4nw 2abw", useNumber)) return  95;
	if (compareName(name,  "97",   "D4^9",    "I422",          "I 4 2 2",            "I 4 2", useNumber)) return  96;
	if (compareName(name,  "98",  "D4^10",   "I4122",         "I 41 2 2",        "I 4bw 2bw", useNumber)) return  97;
	if (compareName(name,  "99",  "C4v^1",    "P4mm",          "P 4 m m",           "P 4 -2", useNumber)) return  98;
	if (compareName(name, "100",  "C4v^2",    "P4bm",          "P 4 b m",         "P 4 -2ab", useNumber)) return  99;
	if (compareName(name, "101",  "C4v^3",   "P42cm",         "P 42 c m",         "P 4c -2c", useNumber)) return 100;
	if (compareName(name, "102",  "C4v^4",   "P42nm",         "P 42 n m",         "P 4n -2n", useNumber)) return 101;
	if (compareName(name, "103",  "C4v^5",    "P4cc",          "P 4 c c",          "P 4 -2c", useNumber)) return 102;
	if (compareName(name, "104",  "C4v^6",    "P4nc",          "P 4 n c",          "P 4 -2n", useNumber)) return 103;
	if (compareName(name, "105",  "C4v^7",   "P42mc",         "P 42 m c",          "P 4c -2", useNumber)) return 104;
	if (compareName(name, "106",  "C4v^8",   "P42bc",         "P 42 b c",        "P 4c -2ab", useNumber)) return 105;
	if (compareName(name, "107",  "C4v^9",    "I4mm",          "I 4 m m",           "I 4 -2", useNumber)) return 106;
	if (compareName(name, "108", "C4v^10",    "I4cm",          "I 4 c m",          "I 4 -2c", useNumber)) return 107;
	if (compareName(name, "109", "C4v^11",   "I41md",         "I 41 m d",         "I 4bw -2", useNumber)) return 108;
	if (compareName(name, "110", "C4v^12",   "I41cd",         "I 41 c d",        "I 4bw -2c", useNumber)) return 109;
	if (compareName(name, "111",  "D2d^1",   "P-42m",         "P -4 2 m",           "P -4 2", useNumber)) return 110;
	if (compareName(name, "112",  "D2d^2",   "P-42c",         "P -4 2 c",          "P -4 2c", useNumber)) return 111;
	if (compareName(name, "113",  "D2d^3",  "P-421m",        "P -4 21 m",         "P -4 2ab", useNumber)) return 112;
	if (compareName(name, "114",  "D2d^4",  "P-421c",        "P -4 21 c",          "P -4 2n", useNumber)) return 113;
	if (compareName(name, "115",  "D2d^5",   "P-4m2",         "P -4 m 2",          "P -4 -2", useNumber)) return 114;
	if (compareName(name, "116",  "D2d^6",   "P-4c2",         "P -4 c 2",         "P -4 -2c", useNumber)) return 115;
	if (compareName(name, "117",  "D2d^7",   "P-4b2",         "P -4 b 2",        "P -4 -2ab", useNumber)) return 116;
	if (compareName(name, "118",  "D2d^8",   "P-4n2",         "P -4 n 2",         "P -4 -2n", useNumber)) return 117;
	if (compareName(name, "119",  "D2d^9",   "I-4m2",         "I -4 m 2",          "I -4 -2", useNumber)) return 118;
	if (compareName(name, "120", "D2d^10",   "I-4c2",         "I -4 c 2",         "I -4 -2c", useNumber)) return 119;
	if (compareName(name, "121", "D2d^11",   "I-42m",         "I -4 2 m",           "I -4 2", useNumber)) return 120;
	if (compareName(name, "122", "D2d^12",   "I-42d",         "I -4 2 d",         "I -4 2bw", useNumber)) return 121;
	if (compareName(name, "123",  "D4h^1",  "P4/mmm",    "P 4/m 2/m 2/m",           "-P 4 2", useNumber)) return 122;
	if (compareName(name, "124",  "D4h^2",  "P4/mcc",    "P 4/m 2/c 2/c",          "-P 4 2c", useNumber)) return 123;
	if (compareName(name, "125",  "D4h^3",  "P4/nbm",    "P 4/n 2/b 2/m",       "P 4 2 -1ab", useNumber)) return 124;
	if (compareName(name, "126",  "D4h^4",  "P4/nnc",    "P 4/n 2/n 2/c",        "P 4 2 -1n", useNumber)) return 125;
	if (compareName(name, "127",  "D4h^5",  "P4/mbm",     "P 4/m 21/b m",         "-P 4 2ab", useNumber)) return 126;
	if (compareName(name, "128",  "D4h^6",  "P4/mnc",     "P 4/m 21/n c",          "-P 4 2n", useNumber)) return 127;
	if (compareName(name, "129",  "D4h^7",  "P4/nmm",     "P 4/n 21/m m",   "P 4ab 2ab -1ab", useNumber)) return 128;
	if (compareName(name, "130",  "D4h^8",  "P4/ncc",     "P 4/n 21/c c",    "P 4ab 2n -1ab", useNumber)) return 129;
	if (compareName(name, "131",  "D4h^9", "P42/mmc",   "P 42/m 2/m 2/c",          "-P 4c 2", useNumber)) return 130;
	if (compareName(name, "132", "D4h^10", "P42/mcm",   "P 42/m 2/c 2/m",         "-P 4c 2c", useNumber)) return 131;
	if (compareName(name, "133", "D4h^11", "P42/nbc",   "P 42/n 2/b 2/c",      "P 4n 2c -1n", useNumber)) return 132;
	if (compareName(name, "134", "D4h^12", "P42/nnm",   "P 42/n 2/n 2/m",       "P 4n 2 -1n", useNumber)) return 133;
	if (compareName(name, "135", "D4h^13", "P42/mbc",  "P 42/m 21/b 2/c",        "-P 4c 2ab", useNumber)) return 134;
	if (compareName(name, "136", "D4h^14", "P42/mnm",  "P 42/m 21/n 2/m",         "-P 4n 2n", useNumber)) return 135;
	if (compareName(name, "137", "D4h^15", "P42/nmc",  "P 42/n 21/m 2/c",      "P 4n 2n -1n", useNumber)) return 136;
	if (compareName(name, "138", "D4h^16", "P42/ncm",  "P 42/n 21/c 2/m",     "P 4n 2ab -1n", useNumber)) return 137;
	if (compareName(name, "139", "D4h^17",  "I4/mmm",    "I 4/m 2/m 2/m",           "-I 4 2", useNumber)) return 138;
	if (compareName(name, "140", "D4h^18",  "I4/mcm",    "I 4/m 2/c 2/m",          "-I 4 2c", useNumber)) return 139;
	if (compareName(name, "141", "D4h^19", "I41/amd",   "I 41/a 2/m 2/d",   "I 4bw 2bw -1bw", useNumber)) return 140;
	if (compareName(name, "142", "D4h^20", "I41/acd",   "I 41/a 2/c 2/d",   "I 4bw 2aw -1bw", useNumber)) return 141;
	if (compareName(name, "143",   "C3^1",      "P3",          "P 3 1 1",              "P 3", useNumber)) return 142;
	if (compareName(name, "144",   "C3^2",     "P31",         "P 31 1 1",             "P 31", useNumber)) return 143;
	if (compareName(name, "145",   "C3^3",     "P32",         "P 32 1 1",             "P 32", useNumber)) return 144;
	if (compareName(name, "146",   "C3^4",      "R3",          "R 3 1 1",              "R 3", useNumber)) return 145;
	if (compareName(name, "147",  "C3i^1",     "P-3",         "P -3 1 1",             "-P 3", useNumber)) return 146;
	if (compareName(name, "148",  "C3i^2",     "R-3",         "R -3 1 1",             "-R 3", useNumber)) return 147;
	if (compareName(name, "149",   "D3^1",    "P312",          "P 3 1 2",            "P 3 2", useNumber)) return 148;
	if (compareName(name, "150",   "D3^2",    "P321",          "P 3 2 1",          "P 3 2\"", useNumber)) return 149;
	if (compareName(name, "151",   "D3^3",   "P3112",         "P 31 1 2",  "P 31 2c (0 0 1)", useNumber)) return 150;
	if (compareName(name, "152",   "D3^4",   "P3121",         "P 31 2 1",         "P 31 2\"", useNumber)) return 151;
	if (compareName(name, "153",   "D3^5",   "P3212",         "P 32 1 2", "P 32 2c (0 0 -1)", useNumber)) return 152;
	if (compareName(name, "154",   "D3^6",   "P3221",         "P 32 2 1",         "P 32 2\"", useNumber)) return 153;
	if (compareName(name, "155",   "D3^7",     "R32",          "R 3 2 1",          "R 3 2\"", useNumber)) return 154;
	if (compareName(name, "156",  "C3v^1",    "P3m1",          "P 3 m 1",         "P 3 -2\"", useNumber)) return 155;
	if (compareName(name, "157",  "C3v^2",    "P31m",          "P 3 1 m",           "P 3 -2", useNumber)) return 156;
	if (compareName(name, "158",  "C3v^3",    "P3c1",          "P 3 c 1",        "P 3 -2\"c", useNumber)) return 157;
	if (compareName(name, "159",  "C3v^4",    "P31c",          "P 3 1 c",          "P 3 -2c", useNumber)) return 158;
	if (compareName(name, "160",  "C3v^5",     "R3m",          "R 3 m 1",         "R 3 -2\"", useNumber)) return 159;
	if (compareName(name, "161",  "C3v^6",     "R3c",          "R 3 c 1",        "R 3 -2\"c", useNumber)) return 160;
	if (compareName(name, "162",  "D3d^1",   "P-31m",       "P -3 1 2/m",           "-P 3 2", useNumber)) return 161;
	if (compareName(name, "163",  "D3d^2",   "P-31c",       "P -3 1 2/c",          "-P 3 2c", useNumber)) return 162;
	if (compareName(name, "164",  "D3d^3",   "P-3m1",       "P -3 2/m 1",         "-P 3 2\"", useNumber)) return 163;
	if (compareName(name, "165",  "D3d^4",   "P-3c1",       "P -3 2/c 1",        "-P 3 2\"c", useNumber)) return 164;
	if (compareName(name, "166",  "D3d^5",    "R-3m",       "R -3 2/m 1",         "-R 3 2\"", useNumber)) return 165;
	if (compareName(name, "167",  "D3d^6",    "R-3c",       "R -3 2/c 1",        "-R 3 2\"c", useNumber)) return 166;
	if (compareName(name, "168",   "C6^1",      "P6",          "P 6 1 1",              "P 6", useNumber)) return 167;
	if (compareName(name, "169",   "C6^2",     "P61",         "P 61 1 1",             "P 61", useNumber)) return 168;
	if (compareName(name, "170",   "C6^3",     "P65",         "P 65 1 1",             "P 65", useNumber)) return 169;
	if (compareName(name, "171",   "C6^4",     "P62",         "P 62 1 1",             "P 62", useNumber)) return 170;
	if (compareName(name, "172",   "C6^5",     "P64",         "P 64 1 1",             "P 64", useNumber)) return 171;
	if (compareName(name, "173",   "C6^6",     "P63",         "P 63 1 1",             "P 6c", useNumber)) return 172;
	if (compareName(name, "174",  "C3h^1",     "P-6",         "P -6 1 1",             "P -6", useNumber)) return 173;
	if (compareName(name, "175",  "C6h^1",    "P6/m",        "P 6/m 1 1",             "-P 6", useNumber)) return 174;
	if (compareName(name, "176",  "C6h^2",   "P63/m",       "P 63/m 1 1",            "-P 6c", useNumber)) return 175;
	if (compareName(name, "177",   "D6^1",    "P622",          "P 6 2 2",            "P 6 2", useNumber)) return 176;
	if (compareName(name, "178",   "D6^2",   "P6122",         "P 61 2 2",  "P 61 2 (0 0 -1)", useNumber)) return 177;
	if (compareName(name, "179",   "D6^3",   "P6522",         "P 65 2 2",   "P 65 2 (0 0 1)", useNumber)) return 178;
	if (compareName(name, "180",   "D6^4",   "P6222",         "P 62 2 2",  "P 62 2c (0 0 1)", useNumber)) return 179;
	if (compareName(name, "181",   "D6^5",   "P6422",         "P 64 2 2", "P 64 2c (0 0 -1)", useNumber)) return 180;
	if (compareName(name, "182",   "D6^6",   "P6322",         "P 63 2 2",          "P 6c 2c", useNumber)) return 181;
	if (compareName(name, "183",  "C6v^1",    "P6mm",          "P 6 m m",           "P 6 -2", useNumber)) return 182;
	if (compareName(name, "184",  "C6v^2",    "P6cc",          "P 6 c c",          "P 6 -2c", useNumber)) return 183;
	if (compareName(name, "185",  "C6v^3",   "P63cm",         "P 63 c m",          "P 6c -2", useNumber)) return 184;
	if (compareName(name, "186",  "C6v^4",   "P63mc",         "P 63 m c",         "P 6c -2c", useNumber)) return 185;
	if (compareName(name, "187",  "D3h^1",   "P-6m2",         "P -6 m 2",           "P -6 2", useNumber)) return 186;
	if (compareName(name, "188",  "D3h^2",   "P-6c2",         "P -6 c 2",          "P -6c 2", useNumber)) return 187;
	if (compareName(name, "189",  "D3h^3",   "P-62m",         "P -6 2 m",          "P -6 -2", useNumber)) return 188;
	if (compareName(name, "190",  "D3h^4",   "P-62c",         "P -6 2 c",        "P -6c -2c", useNumber)) return 189;
	if (compareName(name, "191",  "D6h^1",  "P6/mmm",    "P 6/m 2/m 2/m",           "-P 6 2", useNumber)) return 190;
	if (compareName(name, "192",  "D6h^2",  "P6/mcc",    "P 6/m 2/c 2/c",          "-P 6 2c", useNumber)) return 191;
	if (compareName(name, "193",  "D6h^3", "P63/mcm",   "P 63/m 2/c 2/m",          "-P 6c 2", useNumber)) return 192;
	if (compareName(name, "194",  "D6h^4", "P63/mmc",   "P 63/m 2/m 2/c",         "-P 6c 2c", useNumber)) return 193;
	if (compareName(name, "195",    "T^1",     "P23",          "P 2 3 1",          "P 2 2 3", useNumber)) return 194;
	if (compareName(name, "196",    "T^2",     "F23",          "F 2 3 1",          "F 2 2 3", useNumber)) return 195;
	if (compareName(name, "197",    "T^3",     "I23",          "I 2 3 1",          "I 2 2 3", useNumber)) return 196;
	if (compareName(name, "198",    "T^4",    "P213",         "P 21 3 1",      "P 2ac 2ab 3", useNumber)) return 197;
	if (compareName(name, "199",    "T^5",    "I213",         "I 21 3 1",        "I 2b 2c 3", useNumber)) return 198;
	if (compareName(name, "200",   "Th^1",    "Pm-3",       "P 2/m -3 1",         "-P 2 2 3", useNumber)) return 199;
	if (compareName(name, "201",   "Th^2",    "Pn-3",       "P 2/n -3 1",      "P 2 2 3 -1n", useNumber)) return 200;
	if (compareName(name, "202",   "Th^3",    "Fm-3",       "F 2/m -3 1",         "-F 2 2 3", useNumber)) return 201;
	if (compareName(name, "203",   "Th^4",    "Fd-3",       "F 2/d -3 1",      "F 2 2 3 -1d", useNumber)) return 202;
	if (compareName(name, "204",   "Th^5",    "Im-3",       "I 2/m -3 1",         "-I 2 2 3", useNumber)) return 203;
	if (compareName(name, "205",   "Th^6",    "Pa-3",      "P 21/a -3 1",     "-P 2ac 2ab 3", useNumber)) return 204;
	if (compareName(name, "206",   "Th^7",    "Ia-3",      "I 21/a -3 1",       "-I 2b 2c 3", useNumber)) return 205;
	if (compareName(name, "207",    "O^1",    "P432",          "P 4 3 2",          "P 4 2 3", useNumber)) return 206;
	if (compareName(name, "208",    "O^2",   "P4232",         "P 42 3 2",         "P 4n 2 3", useNumber)) return 207;
	if (compareName(name, "209",    "O^3",    "F432",          "F 4 3 2",          "F 4 2 3", useNumber)) return 208;
	if (compareName(name, "210",    "O^4",   "F4132",         "F 41 3 2",         "F 4d 2 3", useNumber)) return 209;
	if (compareName(name, "211",    "O^5",    "I432",          "I 4 3 2",          "I 4 2 3", useNumber)) return 210;
	if (compareName(name, "212",    "O^6",   "P4332",         "P 43 3 2",     "P 4acd 2ab 3", useNumber)) return 211;
	if (compareName(name, "213",    "O^7",   "P4132",         "P 41 3 2",      "P 4bd 2ab 3", useNumber)) return 212;
	if (compareName(name, "214",    "O^8",   "I4132",         "I 41 3 2",       "I 4bd 2c 3", useNumber)) return 213;
	if (compareName(name, "215",   "Td^1",   "P-43m",         "P -4 3 m",         "P -4 2 3", useNumber)) return 214;
	if (compareName(name, "216",   "Td^2",   "F-43m",         "F -4 3 m",         "F -4 2 3", useNumber)) return 215;
	if (compareName(name, "217",   "Td^3",   "I-43m",         "I -4 3 m",         "I -4 2 3", useNumber)) return 216;
	if (compareName(name, "218",   "Td^4",   "P-43n",         "P -4 3 n",        "P -4n 2 3", useNumber)) return 217;
	if (compareName(name, "219",   "Td^5",   "F-43c",         "F -4 3 c",        "F -4c 2 3", useNumber)) return 218;
	if (compareName(name, "220",   "Td^6",   "I-43d",         "I -4 3 d",      "I -4bd 2c 3", useNumber)) return 219;
	if (compareName(name, "221",   "Oh^1",   "Pm-3m",     "P 4/m -3 2/m",         "-P 4 2 3", useNumber)) return 220;
	if (compareName(name, "222",   "Oh^2",   "Pn-3n",     "P 4/n -3 2/n",      "P 4 2 3 -1n", useNumber)) return 221;
	if (compareName(name, "223",   "Oh^3",   "Pm-3n",    "P 42/m -3 2/n",        "-P 4n 2 3", useNumber)) return 222;
	if (compareName(name, "224",   "Oh^4",   "Pn-3m",    "P 42/n -3 2/m",     "P 4n 2 3 -1n", useNumber)) return 223;
	if (compareName(name, "225",   "Oh^5",   "Fm-3m",     "F 4/m -3 2/m",         "-F 4 2 3", useNumber)) return 224;
	if (compareName(name, "226",   "Oh^6",   "Fm-3c",     "F 4/m -3 2/c",        "-F 4c 2 3", useNumber)) return 225;
	if (compareName(name, "227",   "Oh^7",   "Fd-3m",    "F 41/d -3 2/m",     "F 4d 2 3 -1d", useNumber)) return 226;
	if (compareName(name, "228",   "Oh^8",   "Fd-3c",    "F 41/d -3 2/c",    "F 4d 2 3 -1cd", useNumber)) return 227;
	if (compareName(name, "229",   "Oh^9",   "Im-3m",     "I 4/m -3 2/m",         "-I 4 2 3", useNumber)) return 228;
	if (compareName(name, "230",  "Oh^10",   "Ia-3d",    "I 41/a -3 2/d",      "-I 4bd 2c 3", useNumber)) return 229;
	
	// Did not recognize space group name
	return -1;
}



/* bool SpaceGroup::compareName(const Word& name, const char* number, const char* schoe, const char* hmShort,
 *		const char* hmFull, const char* hall, bool useNumber)
 *
 * Compare name to space group symbols
 */

bool SpaceGroup::compareName(const Word& name, const char* number, const char* schoe, const char* hmShort, \
	const char* hmFull, const char* hall, bool useNumber)
{
	if ((name.equal(schoe, false, ' '))  || (name.equal(hmShort, false, ' ')) || \
		(name.equal(hmFull, false, ' ')) || (name.equal(hall, false, ' ')))
		return true;
	if (useNumber)
	{
		if (name.equal(number))
			return true;
	}
	return false;
}



/* void SpaceGroup::set(const ISO& iso, double tol)
 *
 * Get space group of a structure
 */

void SpaceGroup::set(const ISO& iso, double tol)
{
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Determing the space group of the current structure");
	Output::increase();
	
	// Get the point group
	_pointGroup.set(iso, tol);
	
	// Create conversion matrices
	OList<Matrix3D> conversions = Matrix3D::identity();
	if ((_pointGroup.crystalSystem() == CS_MONOCLINIC) || (_pointGroup.crystalSystem() == CS_ORTHORHOMBIC))
	{
		Matrix3D R2;
		Matrix3D R3;
		if (_pointGroup.crystalSystem() == CS_MONOCLINIC)
		{
			R2.set(0, 0, 1, 0, -1, 0, 1, 0, 0);
			R3.set(-1, 0, 1, 0, 1, 0, -1, 0, 0);
		}
		else
		{
			R2.set(0, 1, 0, 1, 0, 0, 0, 0, -1);
			R3.set(0, 0, 1, 1, 0, 0, 0, 1, 0);
		}
		conversions += R3;
		conversions += R3 * R3;
		conversions += R2;
		conversions += R2 * R3;
		conversions += R2 * R3 * R3;
	}
	
	// Loop over conversion matrices to get space group
	int i, j;
	int index;
	int number = -1;
	List<int> generators;
	Matrix3D Q;
	Matrix3D P;
	Matrix3D tempMat;
	List<int> numbers;
	OList<Vector3D> translations;
	OList<Matrix3D> goodConversions;
	for (i = 0; i < conversions.length(); ++i)
	{
		
		// Set current transformation matrices
		P = conversions[i].transpose();
		Q = P.inverse();
		
		// Get generator indices
		generators.clear();
		translations.clear();
		for (j = 0; j < _pointGroup.convSymmetry().length(); ++j)
		{
			
			// Get current generators
			tempMat  = P;
			tempMat *= _pointGroup.convSymmetry()[j].rotation();
			tempMat *= Q;
			index = getGeneratorIndex(tempMat);
			if (index != -1)
			{
				generators += index;
				translations += _pointGroup.convSymmetry()[j].translations()[0];
				translations.last() *= Q;
				ISO::moveIntoCell(translations.last());
			}
		}
		
		// Get the space group number from generators
		number = generatorsToNumber(_pointGroup.convSymmetry().length(), generators, translations, conversions[i]);
		
		// Found a match
		if (number != -1)
		{
			numbers += number;
			goodConversions += conversions[i];
		}
	}
	
	// Failed to match space group
	if (goodConversions.length() == 0)
	{
		Output::newline(ERROR);
		Output::print("Failed to find space group - try changing the tolerance");
		Output::quit();
	}
	
	// Initialize result
	number = numbers[0];
	_unitToConv = goodConversions[0] * _pointGroup.unitToConv();
	
	// Loop over conversions and save the one that gives the "best" cell
	double cur;
	double best = 0;
	Matrix3D convBasis;
	for (i = 0; i < goodConversions.length(); ++i)
	{
		convBasis = goodConversions[i] * _pointGroup.unitToConv() * iso.basis().vectors();
		cur = Num<double>::abs((Num<double>::angle(convBasis[1], convBasis[2]) + Num<double>::angle(convBasis[0], convBasis[2]) + \
			   Num<double>::angle(convBasis[0], convBasis[1]))/3 - Constants::pi/2);
		if (i == 0)
			best = cur;
		else if (cur < best)
		{
			best = cur;
			number = numbers[i];
			_unitToConv = goodConversions[i] * _pointGroup.unitToConv();
		}
	}
	
	// Set space group properties
	setByNumber(number);
	
	// Print results
	Output::newline();
	Output::print("ITC number: ");
	Output::print(_itcNumber);
	Output::newline();
	Output::print("Hermann-Mauguin short: ");
	Output::print(_hermannMauguinShort);
	Output::newline();
	Output::print("Hermann-Mauguin full: ");
	Output::print(_hermannMauguin);
	Output::newline();
	Output::print("Schoenflies: ");
	Output::print(_schoenflies);
	Output::newline();
	Output::print("Hall: ");
	Output::print(_hall);
	
	// Print conversion to conventional cell
	Output::newline();
	Output::print("Conversion from unit cell to conventional cell");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(_unitToConv(i, j));
			Output::print(" ");
		}
	}
	
	// Print origin shift
	Output::newline();
	Output::print("Origin shift in conventional cell:");
	for (i = 0; i < 3; ++i)
	{
		Output::print(" ");
		Output::print(Language::numberToFraction(_originShift[i]));
	}
	
	// Output
	Output::decrease();
}



/* int SpaceGroup::getGeneratorIndex(const Matrix3D& rotation)
 *
 * Return the generator number given a matrix
 */

int SpaceGroup::getGeneratorIndex(const Matrix3D& rotation)
{
	double tol = 1e-4;
	if ((Num<double>::abs(rotation(0, 0) - 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) - 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 0;
	if ((Num<double>::abs(rotation(0, 0) + 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) + 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 1;
	if ((Num<double>::abs(rotation(0, 0) + 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) - 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 2;
	if ((Num<double>::abs(rotation(0, 0) - 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) + 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 3;
	if ((Num<double>::abs(rotation(0, 0) + 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) + 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 4;
	if ((Num<double>::abs(rotation(0, 0) - 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) + 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 5;
	if ((Num<double>::abs(rotation(0, 0) + 1) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0)) < tol) && \
		(Num<double>::abs(rotation(1, 1) - 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 6;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 7;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) - 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) + 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 8;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1) + 1) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 9;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) + 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 10;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) - 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 11;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) + 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 12;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1) - 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 13;
	if ((Num<double>::abs(rotation(0, 0) - 1) < tol) && (Num<double>::abs(rotation(0, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) - 1) < tol))
		return 14;
	if ((Num<double>::abs(rotation(0, 0) + 1) < tol) && (Num<double>::abs(rotation(0, 1) - 1) < tol) && \
		(Num<double>::abs(rotation(0, 2)) < tol) && (Num<double>::abs(rotation(1, 0) + 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1)) < tol) && \
		(Num<double>::abs(rotation(2, 2) + 1) < tol))
		return 15;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2) - 1) < tol) && (Num<double>::abs(rotation(1, 0) - 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1) - 1) < tol) && \
		(Num<double>::abs(rotation(2, 2)) < tol))
		return 16;
	if ((Num<double>::abs(rotation(0, 0)) < tol) && (Num<double>::abs(rotation(0, 1)) < tol) && \
		(Num<double>::abs(rotation(0, 2) + 1) < tol) && (Num<double>::abs(rotation(1, 0) + 1) < tol) && \
		(Num<double>::abs(rotation(1, 1)) < tol) && (Num<double>::abs(rotation(1, 2)) < tol) && \
		(Num<double>::abs(rotation(2, 0)) < tol) && (Num<double>::abs(rotation(2, 1) + 1) < tol) && \
		(Num<double>::abs(rotation(2, 2)) < tol))
		return 17;
	return -1;
}



/* int SpaceGroup::generatorsToNumber(int numOps, const List<int>& gens, const OList<Vector3D >& trans,
 *		const Matrix3D& conv)
 *
 * Get the space group number from list of generators
 */

int SpaceGroup::generatorsToNumber(int numOps, const List<int>& gens, const OList<Vector3D >& trans, \
	const Matrix3D& conv)
{
	if (compGen(numOps, conv, gens, trans, 1, 0, 0, -1, -1, -1, -1, -1, -1, LS_TRICLINIC, LC_P)) return 0;
	if (compGen(numOps, conv, gens, trans, 2, 1, 0, -1, -1, -1, -1, -1, -1, LS_TRICLINIC, LC_P)) return 1;
	if (compGen(numOps, conv, gens, trans, 2, 2, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 2;
	if (compGen(numOps, conv, gens, trans, 2, 2, 1, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 3;
	if (compGen(numOps, conv, gens, trans, 2, 2, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C)) return 4;
	if (compGen(numOps, conv, gens, trans, 2, 3, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 5;
	if (compGen(numOps, conv, gens, trans, 2, 3, 2, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 6;
	if (compGen(numOps, conv, gens, trans, 2, 3, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C)) return 7;
	if (compGen(numOps, conv, gens, trans, 2, 3, 2, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C)) return 8;
	if (compGen(numOps, conv, gens, trans, 4, 2, 0, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 9;
	if (compGen(numOps, conv, gens, trans, 4, 2, 1, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 10;
	if (compGen(numOps, conv, gens, trans, 4, 2, 0, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_C)) return 11;
	if (compGen(numOps, conv, gens, trans, 4, 2, 2, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 12;
	if (compGen(numOps, conv, gens, trans, 4, 2, 3, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P)) return 13;
	if (compGen(numOps, conv, gens, trans, 4, 2, 2, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_C)) return 14;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 15;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 16;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 5, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 17;
	if (compGen(numOps, conv, gens, trans, 4, 4, 5, 5, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 18;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 19;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 20;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F)) return 21;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 22;
	if (compGen(numOps, conv, gens, trans, 4, 4, 1, 5, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 23;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 24;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 25;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 26;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 27;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 6, 5, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 28;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 3, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 29;
	if (compGen(numOps, conv, gens, trans, 4, 4, 5, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 30;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 31;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 6, 7, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 32;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 7, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 33;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 34;
	if (compGen(numOps, conv, gens, trans, 4, 4, 2, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 35;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 36;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A)) return 37;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A)) return 38;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A)) return 39;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 5, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A)) return 40;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F)) return 41;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 8, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F)) return 42;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 43;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 44;
	if (compGen(numOps, conv, gens, trans, 4, 4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 45;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 46;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 7, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 47;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 48;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 4, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 49;
	if (compGen(numOps, conv, gens, trans, 8, 4, 6, 5, 6, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 50;
	if (compGen(numOps, conv, gens, trans, 8, 4, 6, 5, 3, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 51;
	if (compGen(numOps, conv, gens, trans, 8, 4, 5, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 52;
	if (compGen(numOps, conv, gens, trans, 8, 4, 6, 5, 5, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 53;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 54;
	if (compGen(numOps, conv, gens, trans, 8, 4, 4, 5, 5, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 55;
	if (compGen(numOps, conv, gens, trans, 8, 4, 2, 5, 1, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 56;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 7, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 57;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 4, 1, 4, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 58;
	if (compGen(numOps, conv, gens, trans, 8, 4, 7, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 59;
	if (compGen(numOps, conv, gens, trans, 8, 4, 5, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 60;
	if (compGen(numOps, conv, gens, trans, 8, 4, 5, 5, 7, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P)) return 61;
	if (compGen(numOps, conv, gens, trans, 8, 4, 2, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 62;
	if (compGen(numOps, conv, gens, trans, 8, 4, 3, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 63;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 64;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 65;
	if (compGen(numOps, conv, gens, trans, 8, 4, 1, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 66;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 3, -1, -1, LS_ORTHORHOMBIC, LC_C)) return 67;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_F)) return 68;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 8, -1, -1, LS_ORTHORHOMBIC, LC_F)) return 69;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 70;
	if (compGen(numOps, conv, gens, trans, 8, 4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 71;
	if (compGen(numOps, conv, gens, trans, 8, 4, 1, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 72;
	if (compGen(numOps, conv, gens, trans, 8, 4, 1, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I)) return 73;
	if (compGen(numOps, conv, gens, trans, 4, 7, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 74;
	if (compGen(numOps, conv, gens, trans, 4, 7, 9, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 75;
	if (compGen(numOps, conv, gens, trans, 4, 7, 2, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 76;
	if (compGen(numOps, conv, gens, trans, 4, 7, 10, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 77;
	if (compGen(numOps, conv, gens, trans, 4, 7, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 78;
	if (compGen(numOps, conv, gens, trans, 4, 7, 11, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 79;
	if (compGen(numOps, conv, gens, trans, 4, 8, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 80;
	if (compGen(numOps, conv, gens, trans, 4, 8, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 81;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 82;
	if (compGen(numOps, conv, gens, trans, 8, 7, 2, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 83;
	if (compGen(numOps, conv, gens, trans, 8, 7, 4, 1, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 84;
	if (compGen(numOps, conv, gens, trans, 8, 7, 7, 1, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 85;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 86;
	if (compGen(numOps, conv, gens, trans, 8, 7, 11, 1, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 87;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 88;
	if (compGen(numOps, conv, gens, trans, 8, 7, 4, 5, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 89;
	if (compGen(numOps, conv, gens, trans, 8, 7, 9, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 90;
	if (compGen(numOps, conv, gens, trans, 8, 7, 12, 5, 13, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 91;
	if (compGen(numOps, conv, gens, trans, 8, 7, 2, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 92;
	if (compGen(numOps, conv, gens, trans, 8, 7, 7, 5, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 93;
	if (compGen(numOps, conv, gens, trans, 8, 7, 10, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 94;
	if (compGen(numOps, conv, gens, trans, 8, 7, 13, 5, 12, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 95;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 96;
	if (compGen(numOps, conv, gens, trans, 8, 7, 11, 5, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 97;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 98;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 99;
	if (compGen(numOps, conv, gens, trans, 8, 7, 2, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 100;
	if (compGen(numOps, conv, gens, trans, 8, 7, 7, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 101;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 102;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 103;
	if (compGen(numOps, conv, gens, trans, 8, 7, 2, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 104;
	if (compGen(numOps, conv, gens, trans, 8, 7, 2, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 105;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 106;
	if (compGen(numOps, conv, gens, trans, 8, 7, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 107;
	if (compGen(numOps, conv, gens, trans, 8, 7, 11, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 108;
	if (compGen(numOps, conv, gens, trans, 8, 7, 11, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 109;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 110;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 111;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 112;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 113;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 114;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 115;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 116;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P)) return 117;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 118;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 119;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 120;
	if (compGen(numOps, conv, gens, trans, 8, 8, 0, 5, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I)) return 121;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 122;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 123;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 0, 1, 4, -1, -1, LS_TETRAGONAL, LC_P)) return 124;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 0, 1, 7, -1, -1, LS_TETRAGONAL, LC_P)) return 125;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 4, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 126;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 7, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 127;
	if (compGen(numOps, conv, gens, trans, 16, 7, 4, 5, 4, 1, 4, -1, -1, LS_TETRAGONAL, LC_P)) return 128;
	if (compGen(numOps, conv, gens, trans, 16, 7, 4, 5, 7, 1, 4, -1, -1, LS_TETRAGONAL, LC_P)) return 129;
	if (compGen(numOps, conv, gens, trans, 16, 7, 2, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 130;
	if (compGen(numOps, conv, gens, trans, 16, 7, 2, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 131;
	if (compGen(numOps, conv, gens, trans, 16, 7, 7, 5, 2, 1, 7, -1, -1, LS_TETRAGONAL, LC_P)) return 132;
	if (compGen(numOps, conv, gens, trans, 16, 7, 7, 5, 0, 1, 7, -1, -1, LS_TETRAGONAL, LC_P)) return 133;
	if (compGen(numOps, conv, gens, trans, 16, 7, 2, 5, 4, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 134;
	if (compGen(numOps, conv, gens, trans, 16, 7, 7, 5, 7, 1, 0, -1, -1, LS_TETRAGONAL, LC_P)) return 135;
	if (compGen(numOps, conv, gens, trans, 16, 7, 7, 5, 7, 1, 7, -1, -1, LS_TETRAGONAL, LC_P)) return 136;
	if (compGen(numOps, conv, gens, trans, 16, 7, 7, 5, 4, 1, 7, -1, -1, LS_TETRAGONAL, LC_P)) return 137;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_I)) return 138;
	if (compGen(numOps, conv, gens, trans, 16, 7, 0, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_I)) return 139;
	if (compGen(numOps, conv, gens, trans, 16, 7, 11, 5, 11, 1, 11, -1, -1, LS_TETRAGONAL, LC_I)) return 140;
	if (compGen(numOps, conv, gens, trans, 16, 7, 11, 5, 14, 1, 11, -1, -1, LS_TETRAGONAL, LC_I)) return 141;
	if (compGen(numOps, conv, gens, trans, 3, 9, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 142;
	if (compGen(numOps, conv, gens, trans, 3, 9, 15, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 143;
	if (compGen(numOps, conv, gens, trans, 3, 9, 16, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 144;
	if (compGen(numOps, conv, gens, trans, 3, 9, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_R)) return 145;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 146;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R)) return 147;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 148;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 149;
	if (compGen(numOps, conv, gens, trans, 6, 9, 15, 10, 16, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 150;
	if (compGen(numOps, conv, gens, trans, 6, 9, 15, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 151;
	if (compGen(numOps, conv, gens, trans, 6, 9, 16, 10, 15, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 152;
	if (compGen(numOps, conv, gens, trans, 6, 9, 16, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 153;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R)) return 154;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 12, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 155;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 156;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 12, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 157;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 158;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 12, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R)) return 159;
	if (compGen(numOps, conv, gens, trans, 6, 9, 0, 12, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_R)) return 160;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 161;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 162;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 11, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 163;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 11, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 164;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 11, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_R)) return 165;
	if (compGen(numOps, conv, gens, trans, 12, 9, 0, 11, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_R)) return 166;
	if (compGen(numOps, conv, gens, trans, 6, 14, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 167;
	if (compGen(numOps, conv, gens, trans, 6, 14, 17, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 168;
	if (compGen(numOps, conv, gens, trans, 6, 14, 18, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 169;
	if (compGen(numOps, conv, gens, trans, 6, 14, 15, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 170;
	if (compGen(numOps, conv, gens, trans, 6, 14, 16, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 171;
	if (compGen(numOps, conv, gens, trans, 6, 14, 2, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 172;
	if (compGen(numOps, conv, gens, trans, 6, 15, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 173;
	if (compGen(numOps, conv, gens, trans, 12, 14, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 174;
	if (compGen(numOps, conv, gens, trans, 12, 14, 2, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 175;
	if (compGen(numOps, conv, gens, trans, 12, 14, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 176;
	if (compGen(numOps, conv, gens, trans, 12, 14, 17, 10, 18, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 177;
	if (compGen(numOps, conv, gens, trans, 12, 14, 18, 10, 17, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 178;
	if (compGen(numOps, conv, gens, trans, 12, 14, 15, 10, 16, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 179;
	if (compGen(numOps, conv, gens, trans, 12, 14, 16, 10, 15, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 180;
	if (compGen(numOps, conv, gens, trans, 12, 14, 2, 10, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 181;
	if (compGen(numOps, conv, gens, trans, 12, 14, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 182;
	if (compGen(numOps, conv, gens, trans, 12, 14, 0, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 183;
	if (compGen(numOps, conv, gens, trans, 12, 14, 2, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 184;
	if (compGen(numOps, conv, gens, trans, 12, 14, 2, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 185;
	if (compGen(numOps, conv, gens, trans, 12, 15, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 186;
	if (compGen(numOps, conv, gens, trans, 12, 15, 2, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 187;
	if (compGen(numOps, conv, gens, trans, 12, 15, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 188;
	if (compGen(numOps, conv, gens, trans, 12, 15, 2, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P)) return 189;
	if (compGen(numOps, conv, gens, trans, 24, 14, 0, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 190;
	if (compGen(numOps, conv, gens, trans, 24, 14, 0, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 191;
	if (compGen(numOps, conv, gens, trans, 24, 14, 2, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 192;
	if (compGen(numOps, conv, gens, trans, 24, 14, 2, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P)) return 193;
	if (compGen(numOps, conv, gens, trans, 12, 4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 194;
	if (compGen(numOps, conv, gens, trans, 12, 4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F)) return 195;
	if (compGen(numOps, conv, gens, trans, 12, 4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 196;
	if (compGen(numOps, conv, gens, trans, 12, 4, 5, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 197;
	if (compGen(numOps, conv, gens, trans, 12, 4, 1, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 198;
	if (compGen(numOps, conv, gens, trans, 24, 4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P)) return 199;
	if (compGen(numOps, conv, gens, trans, 24, 4, 0, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P)) return 200;
	if (compGen(numOps, conv, gens, trans, 24, 4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F)) return 201;
	if (compGen(numOps, conv, gens, trans, 24, 4, 0, 5, 0, 16, 0, 17, 8, LS_CUBIC, LC_F)) return 202;
	if (compGen(numOps, conv, gens, trans, 24, 4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_I)) return 203;
	if (compGen(numOps, conv, gens, trans, 24, 4, 5, 5, 4, 16, 0, 1, 0, LS_CUBIC, LC_P)) return 204;
	if (compGen(numOps, conv, gens, trans, 24, 4, 1, 5, 2, 16, 0, 1, 0, LS_CUBIC, LC_I)) return 205;
	if (compGen(numOps, conv, gens, trans, 24, 7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 206;
	if (compGen(numOps, conv, gens, trans, 24, 7, 7, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 207;
	if (compGen(numOps, conv, gens, trans, 24, 7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F)) return 208;
	if (compGen(numOps, conv, gens, trans, 24, 7, 8, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F)) return 209;
	if (compGen(numOps, conv, gens, trans, 24, 7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 210;
	if (compGen(numOps, conv, gens, trans, 24, 7, 19, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 211;
	if (compGen(numOps, conv, gens, trans, 24, 7, 20, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 212;
	if (compGen(numOps, conv, gens, trans, 24, 7, 20, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 213;
	if (compGen(numOps, conv, gens, trans, 24, 8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 214;
	if (compGen(numOps, conv, gens, trans, 24, 8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F)) return 215;
	if (compGen(numOps, conv, gens, trans, 24, 8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 216;
	if (compGen(numOps, conv, gens, trans, 24, 8, 7, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P)) return 217;
	if (compGen(numOps, conv, gens, trans, 24, 8, 2, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F)) return 218;
	if (compGen(numOps, conv, gens, trans, 24, 8, 20, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I)) return 219;
	if (compGen(numOps, conv, gens, trans, 48, 7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P)) return 220;
	if (compGen(numOps, conv, gens, trans, 48, 7, 0, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P)) return 221;
	if (compGen(numOps, conv, gens, trans, 48, 7, 7, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P)) return 222;
	if (compGen(numOps, conv, gens, trans, 48, 7, 7, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P)) return 223;
	if (compGen(numOps, conv, gens, trans, 48, 7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F)) return 224;
	if (compGen(numOps, conv, gens, trans, 48, 7, 2, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F)) return 225;
	if (compGen(numOps, conv, gens, trans, 48, 7, 8, 5, 0, 16, 0, 17, 8, LS_CUBIC, LC_F)) return 226;
	if (compGen(numOps, conv, gens, trans, 48, 7, 8, 5, 0, 16, 0, 17, 21, LS_CUBIC, LC_F)) return 227;
	if (compGen(numOps, conv, gens, trans, 48, 7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_I)) return 228;
	if (compGen(numOps, conv, gens, trans, 48, 7, 20, 5, 2, 16, 0, 1, 0, LS_CUBIC, LC_I)) return 229;
	return -1;
}



/* bool SpaceGroup::compGen(int numOps, const Matrix3D& conv, const List<int>& gens,
 *		const OList<Vector3D>& translations, int numSGOps, int gr1, int gt1, int gr2, int gt2, int gr3,
 *		int gt3, int gr4, int gt4, LatticeSystem system, LatticeCentering centering)
 *
 * Compare generators of space group and current structure
 */

bool SpaceGroup::compGen(int numOps, const Matrix3D& conv, const List<int>& gens, \
	const OList<Vector3D>& translations, int numSGOps, int gr1, int gt1, int gr2, int gt2, int gr3, \
	int gt3, int gr4, int gt4, LatticeSystem system, LatticeCentering centering)
{
	
	// Check number of symmetries, lattice system, and centering
	if (numOps != numSGOps)
		return false;
	if ((system == LS_CUBIC) && (_pointGroup.crystalSystem() != CS_CUBIC))
		return false;
	if ((system == LS_TETRAGONAL) && (_pointGroup.crystalSystem() != CS_TETRAGONAL))
		return false;
	if ((system == LS_HEXAGONAL) && (_pointGroup.crystalSystem() != CS_HEXAGONAL) && \
		(_pointGroup.crystalSystem() != CS_TRIGONAL))
		return false;
	if ((system == LS_RHOMBOHEDRAL) && (_pointGroup.crystalSystem() != CS_TRIGONAL))
		return false;
	if ((system == LS_ORTHORHOMBIC) && (_pointGroup.crystalSystem() != CS_ORTHORHOMBIC))
		return false;
	if ((system == LS_MONOCLINIC) && (_pointGroup.crystalSystem() != CS_MONOCLINIC))
		return false;
	if ((system == LS_TRICLINIC) && (_pointGroup.crystalSystem() != CS_TRICLINIC))
		return false;
	if ((centering == LC_P) && (_pointGroup.centeringType() != CT_PRIMITIVE))
		return false;
	if ((centering == LC_I) && (_pointGroup.centeringType() != CT_BODY))
		return false;
	if ((centering == LC_R) && (_pointGroup.centeringType() != CT_RHOMBOHEDRAL))
		return false;
	if ((centering == LC_F) && (_pointGroup.centeringType() != CT_ALL_FACE))
		return false;
	if ((centering == LC_A) || (centering == LC_B) || (centering == LC_C))
	{
		
		// Check that centering is one face
		if (_pointGroup.centeringType() != CT_ONE_FACE)
			return false;
		
		// Save centering vectors for current cell
		Matrix3D tempMat = conv * _pointGroup.redToConv();
		OList<Vector3D > centVecs = ISO::getLatticePoints(tempMat);
		Matrix3D Q = tempMat.inverse().transpose();
		for (int i = 0; i < centVecs.length(); ++i)
			centVecs[i] *= Q;

		// Make sure that centering of current cell is correct for current space group
		if (centering == LC_A)
		{
			if (!PointGroup::findTranslation(centVecs, 0, 0.5, 0.5))
				return false;
		}
		if (centering == LC_B)
		{
			if (!PointGroup::findTranslation(centVecs, 0.5, 0, 0.5))
				return false;
		}
		if (centering == LC_C)
		{
			if (!PointGroup::findTranslation(centVecs, 0.5, 0.5, 0))
				return false;
		}
	}
	
	// Make generator lists
	List<int> sgRotGen;
	List<int> sgTransGen;
	if (gr1 != -1) { sgRotGen += gr1; sgTransGen += gt1; }
	if (gr2 != -1) { sgRotGen += gr2; sgTransGen += gt2; }
	if (gr3 != -1) { sgRotGen += gr3; sgTransGen += gt3; }
	if (gr4 != -1) { sgRotGen += gr4; sgTransGen += gt4; }
	
	// Make sure that each operator is in space group
	int i, j;
	bool found;
	OList<Vector3D> symTrans (sgTransGen.length());
	for (i = 0; i < sgRotGen.length(); ++i)
	{
		found = false;
		for (j = 0; j < gens.length(); ++j)
		{
			if (sgRotGen[i] == gens[j])
			{
				found = true;
				symTrans[i] = translations[j];
				break;
			}
		}
		if (!found)
			return false;
	}
	
	// Look up rotations
	OList<Matrix3D > sgRot (sgRotGen.length());
	for (i = 0; i < sgRot.length(); ++i)
		sgRot[i] = generatorToMatrix(sgRotGen[i]);
	
	// Look up translations
	OList<Vector3D > sgTrans (sgTransGen.length());
	for (i = 0; i < sgTransGen.length(); ++i)
		sgTrans[i] = generatorToTranslation(sgTransGen[i]);
	
	// Transformations
	Matrix3D Q = (conv * _pointGroup.redToConv()).transpose();
	Matrix3D P = Q.inverse();
	
	// Convert everything to primitive cell
	double tol = 1e-4;
	OList<Matrix3D> primSGRot (sgRot.length());
	OList<Vector3D> primSGTrans (sgTrans.length());
	OList<Vector3D> primSymTrans (sgTrans.length());
	for (i = 0; i < sgRot.length(); ++i)
	{
		primSGRot[i] = P;
		primSGRot[i] *= sgRot[i];
		primSGRot[i] *= Q;
		primSGTrans[i] = sgTrans[i];
		primSGTrans[i] *= Q;
		primSymTrans[i] = symTrans[i];
		primSymTrans[i] *= Q;
		for (j = 0; j < 3; ++j)
		{
			primSGTrans[i][j] = Num<double>::mod(primSGTrans[i][j], 1-tol);
			primSymTrans[i][j] = Num<double>::mod(primSymTrans[i][j], 1-tol);
		}
	}
	
	// Solving for Ax = b, make A and b
	int k;
	Vector b (3 * primSGRot.length());
	Matrix A (3 * primSGTrans.length(), 3);
	for (i = 0; i < primSGRot.length(); ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			b[3*i + j] = primSymTrans[i][j] - primSGTrans[i][j];
			for (k = 0; k < 3; ++k)
				A(3*i + j, k) = primSGRot[i](j, k);
			A(3*i + j, j) -= 1;
		}
	}
	
	// Initialize identity matrices
	Matrix IP = Matrix::identity(A.numRows());
	Matrix IQ = Matrix::identity(3);
	
	// Loop until in Smith normal form
	int count;
	const int maxCount = 1000;
	Matrix M = A;
	for (count = 0; count < maxCount; ++count)
	{
		
		// Convert M to row echelon form
		M = M.rowEchelon(0, &IP, true);
		
		// M is diagonal
		if (M.isDiagonal())
			break;
		
		// Convert transpose of M to row echelon form
		M.makeTranspose();
		M = M.rowEchelon(0, &IQ, true);
		
		// Make M its transpose
		M.makeTranspose();
		
		// M is diagonal
		if (M.isDiagonal())
			break;
	}
	
	// Failed to find Smith normal form
	if (count == maxCount)
	{
		Output::newline(ERROR);
		Output::print("Failed to find Smith normal form of matrix");
		Output::quit();
	}
	
	// Make results vector (D = PAQ already is equal to M)
	Vector v = IP * b;
	for (i = 0; i < v.length(); ++i)
		v[i] = Num<double>::mod(v[i], 1-tol);
	
	// Make sure that there is a valid solution
	for (i = 0; i < 3; ++i)
	{
		if ((Num<double>::abs(M(i, i)) < tol) && (Num<double>::abs(v[i]) > tol))
			return false;
	}
	for (i = 3; i < v.length(); ++i)
	{
		if (Num<double>::abs(v[i]) > tol)
			return false;
	}
	
	// Make Q into 3x3
	Matrix3D IQ33;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			IQ33(i, j) = IQ(j, i);
	}
	
	// Save shift in primitive cell
	for (i = 0; i < 3; i++)
		_originShift[i] = (Num<double>::abs(M(i, i)) < tol) ? 0 : -v[i] / M(i, i);
	_originShift *= IQ33;
	
	// Convert shift to conventional cell
	_originShift *= P;
	ISO::moveIntoCell(_originShift);
	
	// Return that space group matched
	return true;
}



/* Matrix3D SpaceGroup::generatorToMatrix(int index)
 *
 * Return generator matrix
 */

Matrix3D SpaceGroup::generatorToMatrix(int index)
{
	if (index == 0) return Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
	if (index == 1) return Matrix3D(-1, 0, 0, 0, -1, 0, 0, 0, -1);
	if (index == 2) return Matrix3D(-1, 0, 0, 0, 1, 0, 0, 0, -1);
	if (index == 3) return Matrix3D(1, 0, 0, 0, -1, 0, 0, 0, 1);
	if (index == 4) return Matrix3D(-1, 0, 0, 0, -1, 0, 0, 0, 1);
	if (index == 5) return Matrix3D(1, 0, 0, 0, -1, 0, 0, 0, -1);
	if (index == 6) return Matrix3D(-1, 0, 0, 0, 1, 0, 0, 0, 1);
	if (index == 7) return Matrix3D(0, -1, 0, 1, 0, 0, 0, 0, 1);
	if (index == 8) return Matrix3D(0, 1, 0, -1, 0, 0, 0, 0, -1);
	if (index == 9) return Matrix3D(0, -1, 0, 1, -1, 0, 0, 0, 1);
	if (index == 10) return Matrix3D(0, -1, 0, -1, 0, 0, 0, 0, -1);
	if (index == 11) return Matrix3D(0, 1, 0, 1, 0, 0, 0, 0, -1);
	if (index == 12) return Matrix3D(0, -1, 0, -1, 0, 0, 0, 0, 1);
	if (index == 13) return Matrix3D(0, 1, 0, 1, 0, 0, 0, 0, 1);
	if (index == 14) return Matrix3D(1, -1, 0, 1, 0, 0, 0, 0, 1);
	if (index == 15) return Matrix3D(-1, 1, 0, -1, 0, 0, 0, 0, -1);
	if (index == 16) return Matrix3D(0, 0, 1, 1, 0, 0, 0, 1, 0);
	if (index == 17) return Matrix3D(0, 0, -1, -1, 0, 0, 0, -1, 0);
	Output::newline(ERROR);
	Output::print("Internal error: Unknown rotation generator index");
	Output::quit();
	return Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0);
}



/* Vector3D SpaceGroup::generatorToTranslation(int index)
 *
 * Return generator translation
 */

Vector3D SpaceGroup::generatorToTranslation(int index)
{
	if (index == 0) return Vector3D(0.0, 0.0, 0.0);
	if (index == 1) return Vector3D(0.0, 1.0/2.0, 0.0);
	if (index == 2) return Vector3D(0.0, 0.0, 1.0/2.0);
	if (index == 3) return Vector3D(0.0, 1.0/2.0, 1.0/2.0);
	if (index == 4) return Vector3D(1.0/2.0, 1.0/2.0, 0.0);
	if (index == 5) return Vector3D(1.0/2.0, 0.0, 1.0/2.0);
	if (index == 6) return Vector3D(1.0/2.0, 0.0, 0.0);
	if (index == 7) return Vector3D(1.0/2.0, 1.0/2.0, 1.0/2.0);
	if (index == 8) return Vector3D(1.0/4.0, 1.0/4.0, 1.0/4.0);
	if (index == 9) return Vector3D(0.0, 0.0, 1.0/4.0);
	if (index == 10) return Vector3D(0.0, 0.0, 3.0/4.0);
	if (index == 11) return Vector3D(0.0, 1.0/2.0, 1.0/4.0);
	if (index == 12) return Vector3D(1.0/2.0, 1.0/2.0, 1.0/4.0);
	if (index == 13) return Vector3D(1.0/2.0, 1.0/2.0, 3.0/4.0);
	if (index == 14) return Vector3D(1.0/2.0, 0.0, 1.0/4.0);
	if (index == 15) return Vector3D(0.0, 0.0, 1.0/3.0);
	if (index == 16) return Vector3D(0.0, 0.0, 2.0/3.0);
	if (index == 17) return Vector3D(0.0, 0.0, 1.0/6.0);
	if (index == 18) return Vector3D(0.0, 0.0, 5.0/6.0);
	if (index == 19) return Vector3D(3.0/4.0, 1.0/4.0, 3.0/4.0);
	if (index == 20) return Vector3D(1.0/4.0, 3.0/4.0, 1.0/4.0);
	if (index == 21) return Vector3D(1.0/4.0, 1.0/4.0, 3.0/4.0);
	Output::newline(ERROR);
	Output::print("Internal error: Unknown translation generator index");
	Output::quit();
	return Vector3D(0, 0, 0);
}



/* void SpaceGroup::setByNumber(int number)
 *
 * Set space group properties by number
 */

void SpaceGroup::setByNumber(int number)
{
	
	// Set properties of each space group
	if (number == 0)
	{
		setName("1", "C1^1", "P1", "P 1", "P 1");
		setSymmetry(0, 0, -1, -1, -1, -1, -1, -1, LS_TRICLINIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 1)
	{
		setName("2", "Ci^1", "P-1", "P -1", "-P 1");
		setSymmetry(1, 0, -1, -1, -1, -1, -1, -1, LS_TRICLINIC, LC_P);
		setWyckoff(0, 0, 1, 1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 2)
	{
		setName("3", "C2^1", "P2", "P 1 2 1", "P 2y");
		setSymmetry(2, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 2, 3, 2, 5, 2, 7, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 3)
	{
		setName("4", "C2^2", "P21", "P 1 21 1", "P 2yb");
		setSymmetry(2, 1, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 4)
	{
		setName("5", "C2^3", "C2", "C 1 2 1", "C 2y");
		setSymmetry(2, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C);
		setWyckoff(0, 0, 2, 7, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 5)
	{
		setName("6", "Cs^1", "Pm", "P 1 m 1", "P -2y");
		setSymmetry(3, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 6)
	{
		setName("7", "Cs^2", "Pc", "P 1 c 1", "P -2yc");
		setSymmetry(3, 2, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 7)
	{
		setName("8", "Cs^3", "Cm", "C 1 m 1", "C -2y");
		setSymmetry(3, 0, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C);
		setWyckoff(0, 0, 3, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 8)
	{
		setName("9", "Cs^4", "Cc", "C 1 c 1", "C -2yc");
		setSymmetry(3, 2, -1, -1, -1, -1, -1, -1, LS_MONOCLINIC, LC_C);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 9)
	{
		setName("10", "C2h^1", "P2/m", "P 1 2/m 1", "-P 2y");
		setSymmetry(2, 0, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, 2, 3, 2, 7, 2, 5, 2, 0, 1, 1, 1, 3, 1, 2, \
			1, 4, 1, 5, 1, 7, 1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 10)
	{
		setName("11", "C2h^2", "P21/m", "P 1 21/m 1", "-P 2yb");
		setSymmetry(2, 1, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 3, 8, 1, 3, 1, 7, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 11)
	{
		setName("12", "C2h^3", "C2/m", "C 1 2/m 1", "-C 2y");
		setSymmetry(2, 0, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_C);
		setWyckoff(0, 0, 3, 0, 2, 7, 2, 0, 1, 9, 1, 10, 1, 2, 1, 7, 1, 6, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 12)
	{
		setName("13", "C2h^4", "P2/c", "P 1 2/c 1", "-P 2yc");
		setSymmetry(2, 2, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 2, 11, 2, 12, 1, 5, 1, 6, 1, 4, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 13)
	{
		setName("14", "C2h^5", "P21/c", "P 1 21/c 1", "-P 2ybc");
		setSymmetry(2, 3, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_P);
		setWyckoff(0, 0, 1, 3, 1, 7, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 14)
	{
		setName("15", "C2h^6", "C2/c", "C 1 2/c 1", "-C 2yc");
		setSymmetry(2, 2, 1, 0, -1, -1, -1, -1, LS_MONOCLINIC, LC_C);
		setWyckoff(0, 0, 2, 12, 1, 9, 1, 10, 1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 15)
	{
		setName("16", "D2^1", "P222", "P 2 2 2", "P 2 2");
		setSymmetry(4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 4, 4, 6, 4, 5, 4, 0, 2, 3, 2, 5, 2, 7, 2, 0, 5, 2, \
			5, 6, 5, 7, 5, 0, 1, 1, 1, 2, 1, 3, 1, 4, 1, 7, 1, 6, \
			1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 16)
	{
		setName("17", "D2^2", "P2221", "P 2 2 21", "P 2c 2");
		setSymmetry(4, 2, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 2, 11, 2, 12, 5, 6, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 17)
	{
		setName("18", "D2^3", "P21212", "P 21 21 2", "P 2 2ab");
		setSymmetry(4, 0, 5, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 18)
	{
		setName("19", "D2^4", "P212121", "P 21 21 21", "P 2ac 2ab");
		setSymmetry(4, 5, 5, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 19)
	{
		setName("20", "D2^5", "C2221", "C 2 2 21", "C 2c 2");
		setSymmetry(4, 2, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 2, 12, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 20)
	{
		setName("21", "D2^6", "C222", "C 2 2 2", "C 2 2");
		setSymmetry(4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 4, 10, 4, 6, 4, 0, 2, 7, 2, 0, 5, 7, 5, 0, 1, 7, 1, 3, \
			1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 21)
	{
		setName("22", "D2^7", "F222", "F 2 2 2", "F 2 2");
		setSymmetry(4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F);
		setWyckoff(0, 0, 5, 13, 2, 14, 4, 10, 4, 0, 2, 0, 5, 0, 1, 15, 1, 16, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 22)
	{
		setName("23", "D2^8", "I222", "I 2 2 2", "I 2 2");
		setSymmetry(4, 0, 5, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 4, 6, 4, 0, 2, 5, 2, 0, 5, 7, 5, 0, 1, 6, 1, 7, 1, 5, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 23)
	{
		setName("24", "D2^9", "I212121", "I 21 21 21", "I 2b 2c");
		setSymmetry(4, 1, 5, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 4, 8, 2, 17, 5, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 24)
	{
		setName("25", "C2v^1", "Pmm2", "P m m 2", "P 2 -2");
		setSymmetry(4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 5, 6, 0, 3, 6, 3, 0, 4, 4, 4, 5, 4, 6, 4, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 25)
	{
		setName("26", "C2v^2", "Pmc21", "P m c 21", "P 2c -2");
		setSymmetry(4, 2, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 5, 6, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 26)
	{
		setName("27", "C2v^3", "Pcc2", "P c c 2", "P 2 -2c");
		setSymmetry(4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 4, 4, 5, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 27)
	{
		setName("28", "C2v^4", "Pma2", "P m a 2", "P 2 -2a");
		setSymmetry(4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 17, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 28)
	{
		setName("29", "C2v^5", "Pca21", "P c a 21", "P 2c -2ac");
		setSymmetry(4, 2, 6, 5, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 29)
	{
		setName("30", "C2v^6", "Pnc2", "P n c 2", "P 2 -2bc");
		setSymmetry(4, 0, 6, 3, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 30)
	{
		setName("31", "C2v^7", "Pmn21", "P m n 21", "P 2ac -2");
		setSymmetry(4, 5, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 31)
	{
		setName("32", "C2v^8", "Pba2", "P b a 2", "P 2 -2ab");
		setSymmetry(4, 0, 6, 4, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 32)
	{
		setName("33", "C2v^9", "Pna21", "P n a 21", "P 2c -2n");
		setSymmetry(4, 2, 6, 7, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 33)
	{
		setName("34", "C2v^10", "Pnn2", "P n n 2", "P 2 -2n");
		setSymmetry(4, 0, 6, 7, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 34)
	{
		setName("35", "C2v^11", "Cmm2", "C m m 2", "C 2 -2");
		setSymmetry(4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 6, 0, 3, 0, 4, 10, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 35)
	{
		setName("36", "C2v^12", "Cmc21", "C m c 21", "C 2c -2");
		setSymmetry(4, 2, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 6, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 36)
	{
		setName("37", "C2v^13", "Ccc2", "C c c 2", "C 2 -2c");
		setSymmetry(4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 4, 10, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 37)
	{
		setName("38", "C2v^14", "Amm2", "A m m 2", "A 2 -2");
		setSymmetry(4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A);
		setWyckoff(0, 0, 6, 5, 6, 0, 3, 0, 4, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 38)
	{
		setName("39", "C2v^15", "Abm2", "A b m 2", "A 2 -2c");
		setSymmetry(4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A);
		setWyckoff(0, 0, 3, 8, 4, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 39)
	{
		setName("40", "C2v^16", "Ama2", "A m a 2", "A 2 -2a");
		setSymmetry(4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A);
		setWyckoff(0, 0, 6, 17, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 40)
	{
		setName("41", "C2v^17", "Aba2", "A b a 2", "A 2 -2ac");
		setSymmetry(4, 0, 6, 5, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_A);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 41)
	{
		setName("42", "C2v^18", "Fmm2", "F m m 2", "F 2 -2");
		setSymmetry(4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F);
		setWyckoff(0, 0, 3, 0, 6, 0, 4, 10, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 42)
	{
		setName("43", "C2v^19", "Fdd2", "F d d 2", "F 2 -2d");
		setSymmetry(4, 0, 6, 8, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_F);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 43)
	{
		setName("44", "C2v^20", "Imm2", "I m m 2", "I 2 -2");
		setSymmetry(4, 0, 6, 0, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 6, 0, 3, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 44)
	{
		setName("45", "C2v^21", "Iba2", "I b a 2", "I 2 -2c");
		setSymmetry(4, 0, 6, 2, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 45)
	{
		setName("46", "C2v^22", "Ima2", "I m a 2", "I 2 -2a");
		setSymmetry(4, 0, 6, 6, -1, -1, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 6, 17, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 46)
	{
		setName("47", "D2h^1", "Pmmm", "P 2/m 2/m 2/m", "-P 2 2");
		setSymmetry(4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 3, 6, 3, 0, 6, 5, 6, 0, 4, 4, 4, 5, 4, 6, \
			4, 0, 2, 3, 2, 5, 2, 7, 2, 0, 5, 2, 5, 6, 5, 7, 5, 0, \
			1, 1, 1, 2, 1, 4, 1, 6, 1, 3, 1, 7, 1, 5, 1, 0);
	}
	else if (number == 47)
	{
		setName("48", "D2h^2", "Pnnn", "P 2/n 2/n 2/n", "P 2 2 -1n");
		setSymmetry(4, 0, 5, 0, 1, 7, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, 2, 5, 2, 0, 5, 7, 5, 0, 1, 18, 1, 16, 1, 6, \
			1, 7, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 48)
	{
		setName("49", "D2h^3", "Pccm", "P 2/c 2/c 2/m", "-P 2 2c");
		setSymmetry(4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 7, 0, 4, 5, 4, 6, 4, 4, 4, 0, 2, 11, 2, 12, 5, 19, 5, 12, \
			1, 20, 1, 19, 1, 11, 1, 12, 1, 5, 1, 6, 1, 4, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 49)
	{
		setName("50", "D2h^4", "Pban", "P 2/b 2/a 2/n", "P 2 2 -1ab");
		setSymmetry(4, 0, 5, 0, 1, 4, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, 2, 7, 2, 0, 5, 7, 5, 0, 1, 9, 1, 10, 1, 7, \
			1, 3, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 50)
	{
		setName("51", "D2h^5", "Pmma", "P 21/m 2/m 2/a", "-P 2a 2a");
		setSymmetry(4, 6, 5, 6, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 17, 3, 6, 3, 0, 2, 7, 2, 0, 4, 21, 4, 17, 1, 2, 1, 7, \
			1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 51)
	{
		setName("52", "D2h^6", "Pnna", "P 2/n 21/n 2/a", "-P 2a 2bc");
		setSymmetry(4, 6, 5, 3, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 5, 13, 4, 17, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 52)
	{
		setName("53", "D2h^7", "Pmna", "P 2/m 2/n 21/a", "-P 2ac 2");
		setSymmetry(4, 5, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 6, 0, 2, 14, 5, 6, 5, 0, 1, 6, 1, 4, 1, 5, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 53)
	{
		setName("54", "D2h^8", "Pcca", "P 21/c 2/c 2/a", "-P 2a 2ac");
		setSymmetry(4, 6, 5, 5, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 21, 4, 17, 2, 12, 1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 54)
	{
		setName("55", "D2h^9", "Pbam", "P 21/b 21/a 2/m", "-P 2 2ab");
		setSymmetry(4, 0, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 4, 6, 4, 0, 1, 2, 1, 6, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 55)
	{
		setName("56", "D2h^10", "Pccn", "P 21/c 21/c 2/n", "-P 2ab 2ac");
		setSymmetry(4, 4, 5, 5, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 4, 22, 4, 10, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 56)
	{
		setName("57", "D2h^11", "Pbcm", "P 2/b 21/c 21/m", "-P 2c 2b");
		setSymmetry(4, 2, 5, 1, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 7, 12, 5, 8, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 57)
	{
		setName("58", "D2h^12", "Pnnm", "P 21/n 21/n 2/m", "-P 2 2n");
		setSymmetry(4, 0, 5, 7, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 7, 0, 4, 6, 4, 0, 1, 2, 1, 6, 1, 7, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 58)
	{
		setName("59", "D2h^13", "Pmmn", "P 21/m 21/m 2/n", "P 2 2ab -1ab");
		setSymmetry(4, 0, 5, 4, 1, 4, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 3, 0, 6, 0, 1, 9, 1, 10, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 59)
	{
		setName("60", "D2h^14", "Pbcn", "P 21/b 2/c 21/n", "-P 2n 2ab");
		setSymmetry(4, 7, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 2, 12, 1, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 60)
	{
		setName("61", "D2h^15", "Pbca", "P 21/b 21/c 21/a", "-P 2ac 2ab");
		setSymmetry(4, 5, 5, 4, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 61)
	{
		setName("62", "D2h^16", "Pnma", "P 21/n 21/m 21/a", "-P 2ac 2n");
		setSymmetry(4, 5, 5, 7, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_P);
		setWyckoff(0, 0, 3, 8, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 62)
	{
		setName("63", "D2h^17", "Cmcm", "C 2/m 2/c 21/m", "-C 2c 2");
		setSymmetry(4, 2, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 7, 12, 6, 0, 5, 0, 1, 10, 2, 12, 1, 6, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 63)
	{
		setName("64", "D2h^18", "Cmca", "C 2/m 2/c 21/a", "-C 2bc 2");
		setSymmetry(4, 3, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 6, 0, 2, 14, 5, 0, 1, 10, 1, 5, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 64)
	{
		setName("65", "D2h^19", "Cmmm", "C 2/m 2/m 2/m", "-C 2 2");
		setSymmetry(4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 7, 7, 7, 0, 3, 0, 6, 0, 4, 10, 4, 6, 4, 0, 2, 7, 2, 0, \
			5, 7, 5, 0, 1, 9, 1, 10, 1, 7, 1, 3, 1, 5, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 65)
	{
		setName("66", "D2h^20", "Cccm", "C 2/c 2/c 2/m", "-C 2 2c");
		setSymmetry(4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 7, 0, 4, 10, 4, 6, 4, 0, 2, 12, 5, 12, 1, 22, 1, 10, 1, 6, \
			1, 0, 1, 19, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 66)
	{
		setName("67", "D2h^21", "Cmma", "C 2/m 2/m 2/a", "-C 2b 2");
		setSymmetry(4, 1, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 3, 8, 6, 0, 4, 17, 2, 23, 2, 17, 5, 7, 5, 0, 4, 8, 1, 9, \
			1, 10, 1, 7, 1, 0, 1, 23, 1, 17, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 67)
	{
		setName("68", "D2h^22", "Ccca", "C 2/c 2/c 2/a", "C 2 2 -1bc");
		setSymmetry(4, 0, 5, 0, 1, 3, -1, -1, LS_ORTHORHOMBIC, LC_C);
		setWyckoff(0, 0, 4, 10, 4, 0, 2, 0, 5, 0, 1, 13, 1, 14, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 68)
	{
		setName("69", "D2h^23", "Fmmm", "F 2/m 2/m 2/m", "-F 2 2");
		setSymmetry(4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_F);
		setWyckoff(0, 0, 7, 0, 3, 0, 6, 0, 5, 13, 2, 14, 4, 10, 4, 0, 2, 0, 5, 0, \
			1, 16, 1, 10, 1, 14, 1, 13, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 69)
	{
		setName("70", "D2h^24", "Fddd", "F 2/d 2/d 2/d", "F 2 2 -1d");
		setSymmetry(4, 0, 5, 0, 1, 8, -1, -1, LS_ORTHORHOMBIC, LC_F);
		setWyckoff(0, 0, 4, 0, 2, 0, 5, 0, 1, 24, 1, 25, 1, 7, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 70)
	{
		setName("71", "D2h^25", "Immm", "I 2/m 2/m 2/m", "-I 2 2");
		setSymmetry(4, 0, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 7, 0, 3, 0, 6, 0, 1, 16, 4, 5, 4, 0, 2, 7, 2, 0, 5, 6, \
			5, 0, 1, 3, 1, 4, 1, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 71)
	{
		setName("72", "D2h^26", "Ibam", "I 2/b 2/a 2/m", "-I 2 2c");
		setSymmetry(4, 0, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 7, 0, 4, 6, 4, 0, 2, 12, 5, 12, 1, 16, 1, 5, 1, 0, 1, 11, \
			1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 72)
	{
		setName("73", "D2h^27", "Ibca", "I 2/b 2/c 2/a", "-I 2b 2c");
		setSymmetry(4, 1, 5, 2, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 4, 8, 2, 17, 5, 12, 1, 16, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 73)
	{
		setName("74", "D2h^28", "Imma", "I 2/m 2/m 2/a", "-I 2b 2");
		setSymmetry(4, 1, 5, 0, 1, 0, -1, -1, LS_ORTHORHOMBIC, LC_I);
		setWyckoff(0, 0, 3, 8, 6, 0, 2, 14, 5, 0, 4, 8, 1, 15, 1, 16, 1, 7, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 74)
	{
		setName("75", "C4^1", "P4", "P 4 1 1", "P 4");
		setSymmetry(7, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 75)
	{
		setName("76", "C4^2", "P41", "P 41 1 1", "P 4w");
		setSymmetry(7, 9, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 76)
	{
		setName("77", "C4^3", "P42", "P 42 1 1", "P 4c");
		setSymmetry(7, 2, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 77)
	{
		setName("78", "C4^4", "P43", "P 43 1 1", "P 4cw");
		setSymmetry(7, 10, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 78)
	{
		setName("79", "C4^5", "I4", "I 4 1 1", "I 4");
		setSymmetry(7, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 79)
	{
		setName("80", "C4^6", "I41", "I 41 1 1", "I 4bw");
		setSymmetry(7, 11, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 80)
	{
		setName("81", "S4^1", "P-4", "P -4 1 1", "P -4");
		setSymmetry(8, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, 1, 1, 1, 4, 1, 7, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 81)
	{
		setName("82", "S4^2", "I-4", "I -4 1 1", "I -4");
		setSymmetry(8, 0, -1, -1, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 4, 6, 4, 0, 1, 26, 1, 19, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 82)
	{
		setName("83", "C4h^1", "P4/m", "P 4/m 1 1", "-P 4");
		setSymmetry(7, 0, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 4, 6, 4, 4, 4, 0, 1, 2, 1, 6, 1, 1, 1, 4, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 83)
	{
		setName("84", "C4h^2", "P42/m", "P 42/m 1 1", "-P 4c");
		setSymmetry(7, 2, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 4, 6, 4, 4, 4, 0, 1, 20, 1, 12, 1, 2, 1, 6, 1, 4, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 84)
	{
		setName("85", "C4h^3", "P4/n", "P 4/n 1 1", "P 4ab -1ab");
		setSymmetry(7, 4, 1, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 0, 1, 9, 1, 10, 4, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 85)
	{
		setName("86", "C4h^4", "P42/n", "P 42/n 1 1", "P 4n -1n");
		setSymmetry(7, 7, 1, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 0, 4, 6, 1, 15, 1, 16, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 86)
	{
		setName("87", "C4h^5", "I4/m", "I 4/m 1 1", "-I 4");
		setSymmetry(7, 0, 1, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 7, 0, 4, 6, 1, 16, 4, 0, 1, 19, 1, 6, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 87)
	{
		setName("88", "C4h^6", "I41/a", "I 41/a 1 1", "I 4bw -1bw");
		setSymmetry(7, 11, 1, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 4, 0, 1, 27, 1, 28, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 88)
	{
		setName("89", "D4^1", "P422", "P 4 2 2", "P 4 2");
		setSymmetry(7, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 5, 6, 5, 7, 5, 2, 5, 0, 8, 7, 8, 0, 4, 6, 4, 4, 4, 0, \
			1, 3, 1, 5, 1, 1, 1, 4, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 89)
	{
		setName("90", "D4^2", "P4212", "P 4 21 2", "P 4ab 2ab");
		setSymmetry(7, 4, 5, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 7, 8, 0, 4, 0, 4, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 90)
	{
		setName("91", "D4^3", "P4122", "P 41 2 2", "P 4w 2c");
		setSymmetry(7, 9, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 29, 2, 5, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 91)
	{
		setName("92", "D4^4", "P41212", "P 41 21 2", "P 4abw 2nw");
		setSymmetry(7, 12, 5, 13, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 92)
	{
		setName("93", "D4^5", "P4222", "P 42 2 2", "P 4c 2");
		setSymmetry(7, 2, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 30, 8, 12, 5, 6, 5, 7, 5, 2, 5, 0, 4, 6, 4, 4, 4, 0, \
			1, 20, 1, 12, 1, 2, 1, 6, 1, 4, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 93)
	{
		setName("94", "D4^6", "P42212", "P 42 21 2", "P 4n 2n");
		setSymmetry(7, 7, 5, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 7, 8, 0, 4, 6, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 94)
	{
		setName("95", "D4^7", "P4322", "P 43 2 2", "P 4cw 2c");
		setSymmetry(7, 10, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 31, 2, 5, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 95)
	{
		setName("96", "D4^8", "P43212", "P 43 21 2", "P 4nw 2abw");
		setSymmetry(7, 13, 5, 12, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 96)
	{
		setName("97", "D4^9", "I422", "I 4 2 2", "I 4 2");
		setSymmetry(7, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 8, 19, 5, 7, 5, 0, 8, 0, 4, 6, 4, 0, 1, 19, 1, 6, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 97)
	{
		setName("98", "D4^10", "I4122", "I 41 2 2", "I 4bw 2bw");
		setSymmetry(7, 11, 5, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 5, 28, 9, 0, 8, 0, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 98)
	{
		setName("99", "C4v^1", "P4mm", "P 4 m m", "P 4 -2");
		setSymmetry(7, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, 10, 0, 4, 5, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 99)
	{
		setName("100", "C4v^2", "P4bm", "P 4 b m", "P 4 -2ab");
		setSymmetry(7, 0, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 4, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 100)
	{
		setName("101", "C4v^3", "P42cm", "P 42 c m", "P 4c -2c");
		setSymmetry(7, 2, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 4, 6, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 101)
	{
		setName("102", "C4v^4", "P42nm", "P 42 n m", "P 4n -2n");
		setSymmetry(7, 7, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 102)
	{
		setName("103", "C4v^5", "P4cc", "P 4 c c", "P 4 -2c");
		setSymmetry(7, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 103)
	{
		setName("104", "C4v^6", "P4nc", "P 4 n c", "P 4 -2n");
		setSymmetry(7, 0, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 104)
	{
		setName("105", "C4v^7", "P42mc", "P 42 m c", "P 4c -2");
		setSymmetry(7, 2, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, 4, 6, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 105)
	{
		setName("106", "C4v^8", "P42bc", "P 42 b c", "P 4c -2ab");
		setSymmetry(7, 2, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 106)
	{
		setName("107", "C4v^9", "I4mm", "I 4 m m", "I 4 -2");
		setSymmetry(7, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 3, 0, 10, 0, 4, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 107)
	{
		setName("108", "C4v^10", "I4cm", "I 4 c m", "I 4 -2c");
		setSymmetry(7, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 10, 6, 4, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 108)
	{
		setName("109", "C4v^11", "I41md", "I 41 m d", "I 4bw -2");
		setSymmetry(7, 11, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 6, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 109)
	{
		setName("110", "C4v^12", "I41cd", "I 41 c d", "I 4bw -2c");
		setSymmetry(7, 11, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 110)
	{
		setName("111", "D2d^1", "P-42m", "P -4 2 m", "P -4 2");
		setSymmetry(8, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 4, 6, 5, 6, 5, 7, 5, 2, 5, 0, 4, 4, 4, 0, 1, 3, \
			1, 5, 1, 4, 1, 7, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 111)
	{
		setName("112", "D2d^2", "P-42c", "P -4 2 c", "P -4 2c");
		setSymmetry(8, 0, 5, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, 2, 12, 5, 19, 2, 11, 5, 12, 1, 4, 1, 0, \
			1, 19, 1, 20, 1, 11, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 112)
	{
		setName("113", "D2d^3", "P-421m", "P -4 21 m", "P -4 2ab");
		setSymmetry(8, 0, 5, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 4, 0, 4, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 113)
	{
		setName("114", "D2d^4", "P-421c", "P -4 21 c", "P -4 2n");
		setSymmetry(8, 0, 5, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 114)
	{
		setName("115", "D2d^5", "P-4m2", "P -4 m 2", "P -4 -2");
		setSymmetry(8, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, 8, 7, 8, 0, 4, 6, 4, 4, 4, 0, 1, 7, 1, 1, \
			1, 4, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 115)
	{
		setName("116", "D2d^6", "P-4c2", "P -4 c 2", "P -4 -2c");
		setSymmetry(8, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 4, 4, 4, 0, 8, 30, 8, 12, 1, 4, 1, 0, 1, 20, 1, 12, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 116)
	{
		setName("117", "D2d^7", "P-4b2", "P -4 b 2", "P -4 -2ab");
		setSymmetry(8, 0, 6, 4, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 2, 8, 6, 4, 6, 4, 0, 1, 2, 1, 6, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 117)
	{
		setName("118", "D2d^8", "P-4n2", "P -4 n 2", "P -4 -2n");
		setSymmetry(8, 0, 6, 7, -1, -1, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 4, 6, 8, 19, 9, 19, 4, 0, 1, 26, 1, 19, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 118)
	{
		setName("119", "D2d^9", "I-4m2", "I -4 m 2", "I -4 -2");
		setSymmetry(8, 0, 6, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 3, 0, 8, 19, 8, 0, 4, 6, 4, 0, 1, 26, 1, 19, 1, 7, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 119)
	{
		setName("120", "D2d^10", "I-4c2", "I -4 c 2", "I -4 -2c");
		setSymmetry(8, 0, 6, 2, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 8, 6, 4, 6, 4, 0, 8, 12, 1, 6, 1, 19, 1, 0, 1, 12, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 120)
	{
		setName("121", "D2d^11", "I-42m", "I -4 2 m", "I -4 2");
		setSymmetry(8, 0, 5, 0, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 10, 0, 4, 6, 5, 7, 5, 0, 4, 0, 1, 19, 1, 6, 1, 7, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 121)
	{
		setName("122", "D2d^12", "I-42d", "I -4 2 d", "I -4 2bw");
		setSymmetry(8, 0, 5, 11, -1, -1, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 5, 28, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 122)
	{
		setName("123", "D4h^1", "P4/mmm", "P 4/m 2/m 2/m", "-P 4 2");
		setSymmetry(7, 0, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 3, 6, 3, 0, 10, 0, 7, 7, 7, 0, 5, 2, 5, 6, 5, 7, 5, 0, \
			8, 7, 8, 0, 4, 6, 4, 4, 4, 0, 1, 6, 1, 2, 1, 1, 1, 4, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 123)
	{
		setName("124", "D4h^2", "P4/mcc", "P 4/m 2/c 2/c", "-P 4 2c");
		setSymmetry(7, 0, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 5, 19, 5, 12, 8, 12, 4, 6, 4, 4, 4, 0, 1, 19, 1, 6, \
			1, 4, 1, 20, 1, 0, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 124)
	{
		setName("125", "D4h^3", "P4/nbm", "P 4/n 2/b 2/m", "P 4 2 -1ab");
		setSymmetry(7, 0, 5, 0, 1, 4, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 5, 7, 5, 0, 8, 7, 8, 0, 4, 6, 4, 0, 1, 9, 1, 10, \
			1, 2, 1, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 125)
	{
		setName("126", "D4h^4", "P4/nnc", "P 4/n 2/n 2/c", "P 4 2 -1n");
		setSymmetry(7, 0, 5, 0, 1, 7, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 5, 7, 5, 0, 8, 0, 4, 5, 1, 16, 4, 0, 1, 11, 1, 5, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 126)
	{
		setName("127", "D4h^5", "P4/mbm", "P 4/m 21/b m", "-P 4 2ab");
		setSymmetry(7, 0, 5, 4, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 7, 7, 7, 0, 8, 2, 8, 6, 4, 6, 4, 0, 1, 6, 1, 2, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 127)
	{
		setName("128", "D4h^6", "P4/mnc", "P 4/m 21/n c", "-P 4 2n");
		setSymmetry(7, 0, 5, 7, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 8, 19, 4, 6, 4, 0, 1, 19, 1, 6, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 128)
	{
		setName("129", "D4h^7", "P4/nmm", "P 4/n 21/m m", "P 4ab 2ab -1ab");
		setSymmetry(7, 4, 5, 4, 1, 4, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 6, 0, 8, 7, 8, 0, 4, 0, 1, 9, 1, 10, 4, 6, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 129)
	{
		setName("130", "D4h^8", "P4/ncc", "P 4/n 21/c c", "P 4ab 2n -1ab");
		setSymmetry(7, 4, 5, 7, 1, 4, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 12, 4, 0, 1, 10, 4, 6, 1, 0, 1, 12, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 130)
	{
		setName("131", "D4h^9", "P42/mmc", "P 42/m 2/m 2/c", "-P 4c 2");
		setSymmetry(7, 2, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 6, 5, 6, 0, 8, 12, 5, 6, 5, 7, 5, 2, 5, 0, 4, 6, \
			4, 4, 4, 0, 1, 20, 1, 12, 1, 2, 1, 6, 1, 4, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 131)
	{
		setName("132", "D4h^10", "P42/mcm", "P 42/m 2/c 2/m", "-P 4c 2c");
		setSymmetry(7, 2, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 7, 0, 5, 19, 5, 12, 4, 6, 8, 7, 8, 0, 4, 4, 4, 0, \
			1, 6, 1, 19, 1, 20, 1, 4, 1, 12, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 132)
	{
		setName("133", "D4h^11", "P42/nbc", "P 42/n 2/b 2/c", "P 4n 2c -1n");
		setSymmetry(7, 7, 5, 2, 1, 7, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 8, 6, 5, 30, 5, 12, 4, 0, 4, 6, 1, 16, 1, 0, 1, 6, 1, 12, \
			1, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 133)
	{
		setName("134", "D4h^12", "P42/nnm", "P 42/n 2/n 2/m", "P 4n 2 -1n");
		setSymmetry(7, 7, 5, 0, 1, 7, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 8, 26, 8, 19, 5, 7, 5, 0, 4, 6, 4, 0, 1, 18, 1, 16, \
			1, 19, 1, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 134)
	{
		setName("135", "D4h^13", "P42/mbc", "P 42/m 21/b 2/c", "-P 4c 2ab");
		setSymmetry(7, 2, 5, 4, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 8, 19, 4, 6, 4, 0, 1, 19, 1, 6, 1, 12, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 135)
	{
		setName("136", "D4h^14", "P42/mnm", "P 42/m 21/n 2/m", "-P 4n 2n");
		setSymmetry(7, 7, 5, 7, 1, 0, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 0, 7, 0, 4, 6, 9, 0, 8, 0, 4, 0, 1, 19, 1, 6, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 136)
	{
		setName("137", "D4h^15", "P42/nmc", "P 42/n 21/m 2/c", "P 4n 2n -1n");
		setSymmetry(7, 7, 5, 7, 1, 7, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 6, 0, 8, 0, 1, 16, 4, 6, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 137)
	{
		setName("138", "D4h^16", "P42/ncm", "P 42/n 21/c 2/m", "P 4n 2ab -1n");
		setSymmetry(7, 7, 5, 4, 1, 7, -1, -1, LS_TETRAGONAL, LC_P);
		setWyckoff(0, 0, 10, 6, 8, 30, 8, 12, 4, 0, 4, 6, 1, 15, 1, 16, 1, 0, 1, 12, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 138)
	{
		setName("139", "D4h^17", "I4/mmm", "I 4/m 2/m 2/m", "-I 4 2");
		setSymmetry(7, 0, 5, 0, 1, 0, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 6, 0, 10, 0, 7, 0, 8, 19, 5, 6, 5, 0, 8, 0, 4, 6, 1, 16, \
			4, 0, 1, 19, 1, 6, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 139)
	{
		setName("140", "D4h^18", "I4/mcm", "I 4/m 2/c 2/m", "-I 4 2c");
		setSymmetry(7, 0, 5, 2, 1, 0, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 10, 6, 7, 0, 5, 12, 8, 12, 8, 6, 4, 6, 4, 0, 1, 16, 1, 6, \
			1, 0, 1, 19, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 140)
	{
		setName("141", "D4h^19", "I41/amd", "I 41/a 2/m 2/d", "I 4bw 2bw -1bw");
		setSymmetry(7, 11, 5, 11, 1, 11, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 6, 0, 8, 0, 5, 28, 4, 0, 1, 27, 1, 28, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 141)
	{
		setName("142", "D4h^20", "I41/acd", "I 41/a 2/c 2/d", "I 4bw 2aw -1bw");
		setSymmetry(7, 11, 5, 14, 1, 11, -1, -1, LS_TETRAGONAL, LC_I);
		setWyckoff(0, 0, 8, 12, 2, 32, 4, 0, 1, 28, 1, 12, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 142)
	{
		setName("143", "C3^1", "P3", "P 3 1 1", "P 3");
		setSymmetry(9, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 33, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 143)
	{
		setName("144", "C3^2", "P31", "P 31 1 1", "P 31");
		setSymmetry(9, 15, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 144)
	{
		setName("145", "C3^3", "P32", "P 32 1 1", "P 32");
		setSymmetry(9, 16, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 145)
	{
		setName("146", "C3^4", "R3", "R 3 1 1", "R 3");
		setSymmetry(9, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 146)
	{
		setName("147", "C3i^1", "P-3", "P -3 1 1", "-P 3");
		setSymmetry(9, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 1, 3, 1, 5, 4, 34, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 147)
	{
		setName("148", "C3i^2", "R-3", "R -3 1 1", "-R 3");
		setSymmetry(9, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 1, 5, 1, 3, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 148)
	{
		setName("149", "D3^1", "P312", "P 3 1 2", "P 3 2");
		setSymmetry(9, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 9, 7, 9, 0, 4, 33, 4, 34, 4, 0, 1, 35, 1, 33, 1, 36, 1, 34, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 149)
	{
		setName("150", "D3^2", "P321", "P 3 2 1", "P 3 2\"");
		setSymmetry(9, 0, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 5, 7, 5, 0, 4, 34, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 150)
	{
		setName("151", "D3^3", "P3112", "P 31 1 2", "P 31 2c (0 0 1)");
		setSymmetry(9, 15, 10, 16, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 9, 37, 9, 38, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 151)
	{
		setName("152", "D3^4", "P3121", "P 31 2 1", "P 31 2\"");
		setSymmetry(9, 15, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 5, 37, 5, 38, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 152)
	{
		setName("153", "D3^5", "P3212", "P 32 1 2", "P 32 2c (0 0 -1)");
		setSymmetry(9, 16, 10, 15, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 9, 39, 9, 40, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 153)
	{
		setName("154", "D3^6", "P3221", "P 32 2 1", "P 32 2\"");
		setSymmetry(9, 16, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 5, 39, 5, 40, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 154)
	{
		setName("155", "D3^7", "R32", "R 3 2 1", "R 3 2\"");
		setSymmetry(9, 0, 11, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 5, 7, 5, 0, 4, 0, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 155)
	{
		setName("156", "C3v^1", "P3m1", "P 3 m 1", "P 3 -2\"");
		setSymmetry(9, 0, 12, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 11, 0, 4, 33, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 156)
	{
		setName("157", "C3v^2", "P31m", "P 3 1 m", "P 3 -2");
		setSymmetry(9, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 3, 0, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 157)
	{
		setName("158", "C3v^3", "P3c1", "P 3 c 1", "P 3 -2\"c");
		setSymmetry(9, 0, 12, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 33, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 158)
	{
		setName("159", "C3v^4", "P31c", "P 3 1 c", "P 3 -2c");
		setSymmetry(9, 0, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 159)
	{
		setName("160", "C3v^5", "R3m", "R 3 m 1", "R 3 -2\"");
		setSymmetry(9, 0, 12, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 11, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 160)
	{
		setName("161", "C3v^6", "R3c", "R 3 c 1", "R 3 -2\"c");
		setSymmetry(9, 0, 12, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 161)
	{
		setName("162", "D3d^1", "P-31m", "P -3 1 2/m", "-P 3 2");
		setSymmetry(9, 0, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 3, 0, 9, 7, 9, 0, 4, 34, 1, 3, 1, 5, 4, 0, 1, 36, 1, 34, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 162)
	{
		setName("163", "D3d^2", "P-31c", "P -3 1 2/c", "-P 3 2c");
		setSymmetry(9, 0, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 9, 12, 1, 5, 4, 34, 4, 0, 1, 41, 1, 42, 1, 0, 1, 12, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 163)
	{
		setName("164", "D3d^3", "P-3m1", "P -3 2/m 1", "-P 3 2\"");
		setSymmetry(9, 0, 11, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 11, 0, 5, 7, 5, 0, 1, 3, 1, 5, 4, 34, 4, 0, 1, 7, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 164)
	{
		setName("165", "D3d^4", "P-3c1", "P -3 2/c 1", "-P 3 2\"c");
		setSymmetry(9, 0, 11, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 5, 12, 1, 5, 4, 34, 4, 0, 1, 0, 1, 12, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 165)
	{
		setName("166", "D3d^5", "R-3m", "R -3 2/m 1", "-R 3 2\"");
		setSymmetry(9, 0, 11, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 11, 0, 5, 7, 5, 0, 1, 5, 1, 3, 4, 0, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 166)
	{
		setName("167", "D3d^6", "R-3c", "R -3 2/c 1", "-R 3 2\"c");
		setSymmetry(9, 0, 11, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_R);
		setWyckoff(0, 0, 5, 12, 1, 5, 4, 0, 1, 0, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 167)
	{
		setName("168", "C6^1", "P6", "P 6 1 1", "P 6");
		setSymmetry(14, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 5, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 168)
	{
		setName("169", "C6^2", "P61", "P 61 1 1", "P 61");
		setSymmetry(14, 17, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 169)
	{
		setName("170", "C6^3", "P65", "P 65 1 1", "P 65");
		setSymmetry(14, 18, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 170)
	{
		setName("171", "C6^4", "P62", "P 62 1 1", "P 62");
		setSymmetry(14, 15, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 171)
	{
		setName("172", "C6^5", "P64", "P 64 1 1", "P 64");
		setSymmetry(14, 16, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 4, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 172)
	{
		setName("173", "C6^6", "P63", "P 63 1 1", "P 6c");
		setSymmetry(14, 2, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 173)
	{
		setName("174", "C3h^1", "P-6", "P -6 1 1", "P -6");
		setSymmetry(15, 0, -1, -1, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 4, 33, 4, 34, 4, 0, 1, 35, 1, 33, 1, 36, 1, 34, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 174)
	{
		setName("175", "C6h^1", "P6/m", "P 6/m 1 1", "-P 6");
		setSymmetry(14, 0, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 4, 5, 4, 34, 1, 3, 1, 5, 4, 0, 1, 36, 1, 34, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 175)
	{
		setName("176", "C6h^2", "P63/m", "P 63/m 1 1", "-P 6c");
		setSymmetry(14, 2, 1, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 12, 1, 5, 4, 34, 4, 0, 1, 41, 1, 42, 1, 0, 1, 12, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 176)
	{
		setName("177", "D6^1", "P622", "P 6 2 2", "P 6 2");
		setSymmetry(14, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 9, 7, 9, 0, 5, 7, 5, 0, 4, 5, 4, 34, 1, 3, 1, 5, 4, 0, \
			1, 36, 1, 34, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 177)
	{
		setName("178", "D6^2", "P6122", "P 61 2 2", "P 61 2 (0 0 -1)");
		setSymmetry(14, 17, 10, 18, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 12, 12, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 178)
	{
		setName("179", "D6^3", "P6522", "P 65 2 2", "P 65 2 (0 0 1)");
		setSymmetry(14, 18, 10, 17, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 12, 30, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 179)
	{
		setName("180", "D6^4", "P6222", "P 62 2 2", "P 62 2c (0 0 1)");
		setSymmetry(14, 15, 10, 16, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 12, 7, 12, 0, 5, 7, 5, 0, 4, 5, 4, 0, 1, 3, 1, 5, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 180)
	{
		setName("181", "D6^5", "P6422", "P 64 2 2", "P 64 2c (0 0 -1)");
		setSymmetry(14, 16, 10, 15, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 12, 7, 12, 0, 5, 7, 5, 0, 4, 5, 4, 0, 1, 3, 1, 5, 1, 7, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 181)
	{
		setName("182", "D6^6", "P6322", "P 63 2 2", "P 6c 2c");
		setSymmetry(14, 2, 10, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 12, 12, 5, 0, 4, 34, 4, 0, 1, 43, 1, 42, 1, 12, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 182)
	{
		setName("183", "C6v^1", "P6mm", "P 6 m m", "P 6 -2");
		setSymmetry(14, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 11, 0, 3, 0, 4, 5, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 183)
	{
		setName("184", "C6v^2", "P6cc", "P 6 c c", "P 6 -2c");
		setSymmetry(14, 0, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 4, 5, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 184)
	{
		setName("185", "C6v^3", "P63cm", "P 63 c m", "P 6c -2");
		setSymmetry(14, 2, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 3, 0, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 185)
	{
		setName("186", "C6v^4", "P63mc", "P 63 m c", "P 6c -2c");
		setSymmetry(14, 2, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 11, 0, 4, 34, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 186)
	{
		setName("187", "D3h^1", "P-6m2", "P -6 m 2", "P -6 2");
		setSymmetry(15, 0, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 11, 0, 7, 7, 7, 0, 9, 7, 9, 0, 4, 33, 4, 34, 4, 0, 1, 35, \
			1, 33, 1, 36, 1, 34, 1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 187)
	{
		setName("188", "D3h^2", "P-6c2", "P -6 c 2", "P -6c 2");
		setSymmetry(15, 2, 10, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 12, 9, 0, 4, 33, 4, 34, 4, 0, 1, 41, 1, 33, 1, 42, 1, 34, \
			1, 12, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 188)
	{
		setName("189", "D3h^3", "P-62m", "P -6 2 m", "P -6 -2");
		setSymmetry(15, 0, 13, 0, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 3, 0, 4, 34, 5, 7, 5, 0, 4, 0, 1, 36, 1, 34, \
			1, 7, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 189)
	{
		setName("190", "D3h^4", "P-62c", "P -6 2 c", "P -6c -2c");
		setSymmetry(15, 2, 13, 2, -1, -1, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 12, 5, 0, 4, 34, 4, 0, 1, 41, 1, 42, 1, 12, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 190)
	{
		setName("191", "D6h^1", "P6/mmm", "P 6/m 2/m 2/m", "-P 6 2");
		setSymmetry(14, 0, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 7, 7, 0, 13, 0, 3, 0, 12, 7, 12, 0, 5, 7, 5, 0, 4, 5, \
			4, 34, 1, 3, 1, 5, 4, 0, 1, 36, 1, 34, 1, 7, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 191)
	{
		setName("192", "D6h^2", "P6/mcc", "P 6/m 2/c 2/c", "-P 6 2c");
		setSymmetry(14, 0, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 7, 0, 12, 12, 5, 12, 4, 5, 4, 34, 1, 5, 1, 11, 4, 0, 1, 34, \
			1, 42, 1, 0, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 192)
	{
		setName("193", "D6h^3", "P63/mcm", "P 63/m 2/c 2/m", "-P 6c 2");
		setSymmetry(14, 2, 10, 0, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 3, 0, 7, 12, 12, 0, 4, 34, 5, 12, 1, 5, 4, 0, 1, 34, 1, 42, \
			1, 0, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 193)
	{
		setName("194", "D6h^4", "P63/mmc", "P 63/m 2/m 2/c", "-P 6c 2c");
		setSymmetry(14, 2, 10, 2, 1, 0, -1, -1, LS_HEXAGONAL, LC_P);
		setWyckoff(0, 0, 13, 0, 7, 12, 5, 0, 12, 12, 1, 5, 4, 34, 4, 0, 1, 43, 1, 42, \
			1, 12, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 194)
	{
		setName("195", "T^1", "P23", "P 2 3 1", "P 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 5, 2, 5, 6, 5, 7, 5, 0, 14, 0, 1, 5, 1, 2, 1, 1, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 195)
	{
		setName("196", "T^2", "F23", "F 2 3 1", "F 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 5, 13, 5, 0, 14, 0, 1, 18, 1, 16, 1, 1, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 196)
	{
		setName("197", "T^3", "I23", "I 2 3 1", "I 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 5, 6, 5, 0, 14, 0, 1, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 197)
	{
		setName("198", "T^4", "P213", "P 21 3 1", "P 2ac 2ab 3");
		setSymmetry(4, 5, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 14, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 198)
	{
		setName("199", "T^5", "I213", "I 21 3 1", "I 2b 2c 3");
		setSymmetry(4, 1, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 5, 12, 14, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 199)
	{
		setName("200", "Th^1", "Pm-3", "P 2/m -3 1", "-P 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 6, 5, 6, 0, 14, 0, 5, 2, 5, 6, 5, 7, 5, 0, 1, 5, 1, 2, \
			1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 200)
	{
		setName("201", "Th^2", "Pn-3", "P 2/n -3 1", "P 2 2 3 -1n");
		setSymmetry(4, 0, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 5, 6, 5, 0, 14, 0, 1, 2, 1, 18, 1, 16, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 201)
	{
		setName("202", "Th^3", "Fm-3", "F 2/m -3 1", "-F 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 6, 0, 5, 13, 14, 0, 5, 0, 1, 13, 1, 16, 1, 1, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 202)
	{
		setName("203", "Th^4", "Fd-3", "F 2/d -3 1", "F 2 2 3 -1d");
		setSymmetry(4, 0, 5, 0, 16, 0, 17, 8, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 5, 0, 14, 0, 1, 24, 1, 25, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 203)
	{
		setName("204", "Th^5", "Im-3", "I 2/m -3 1", "-I 2 2 3");
		setSymmetry(4, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 6, 0, 14, 0, 5, 7, 5, 0, 1, 16, 1, 2, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 204)
	{
		setName("205", "Th^6", "Pa-3", "P 21/a -3 1", "-P 2ac 2ab 3");
		setSymmetry(4, 5, 5, 4, 16, 0, 1, 0, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 14, 0, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 205)
	{
		setName("206", "Th^7", "Ia-3", "I 21/a -3 1", "-I 2b 2c 3");
		setSymmetry(4, 1, 5, 2, 16, 0, 1, 0, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 5, 12, 14, 0, 1, 16, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 206)
	{
		setName("207", "O^1", "P432", "P 4 3 2", "P 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 15, 5, 15, 0, 5, 6, 14, 0, 5, 2, 5, 0, 1, 5, 1, 2, 1, 1, \
			1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 207)
	{
		setName("208", "O^2", "P4232", "P 42 3 2", "P 4n 2 3");
		setSymmetry(7, 7, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 15, 23, 16, 23, 5, 6, 5, 7, 5, 0, 14, 0, 1, 21, 1, 23, 1, 2, \
			1, 18, 1, 16, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 208)
	{
		setName("209", "O^3", "F432", "F 4 3 2", "F 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 5, 13, 15, 5, 15, 0, 14, 0, 5, 0, 1, 13, 1, 16, 1, 1, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 209)
	{
		setName("210", "O^4", "F4132", "F 41 3 2", "F 4d 2 3");
		setSymmetry(7, 8, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 16, 44, 5, 0, 14, 0, 1, 24, 1, 25, 1, 1, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 210)
	{
		setName("211", "O^5", "I432", "I 4 3 2", "I 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 16, 23, 15, 0, 5, 6, 14, 0, 5, 0, 1, 21, 1, 16, 1, 2, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 211)
	{
		setName("212", "O^6", "P4332", "P 43 3 2", "P 4acd 2ab 3");
		setSymmetry(7, 19, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 16, 44, 14, 0, 1, 24, 1, 25, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 212)
	{
		setName("213", "O^7", "P4132", "P 41 3 2", "P 4bd 2ab 3");
		setSymmetry(7, 20, 5, 4, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 15, 44, 14, 0, 1, 45, 1, 46, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 213)
	{
		setName("214", "O^8", "I4132", "I 41 3 2", "I 4bd 2c 3");
		setSymmetry(7, 20, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 16, 44, 15, 44, 5, 12, 14, 0, 1, 47, 1, 44, 1, 45, 1, 25, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 214)
	{
		setName("215", "Td^1", "P-43m", "P -4 3 m", "P -4 2 3");
		setSymmetry(8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 10, 0, 5, 6, 5, 2, 5, 0, 14, 0, 1, 5, 1, 2, 1, 1, 1, 0, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 215)
	{
		setName("216", "Td^2", "F-43m", "F -4 3 m", "F -4 2 3");
		setSymmetry(8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 10, 0, 5, 13, 5, 0, 14, 0, 1, 18, 1, 16, 1, 1, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 216)
	{
		setName("217", "Td^3", "I-43m", "I -4 3 m", "I -4 2 3");
		setSymmetry(8, 0, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 10, 0, 5, 6, 5, 0, 1, 21, 14, 0, 1, 2, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 217)
	{
		setName("218", "Td^4", "P-43n", "P -4 3 n", "P -4n 2 3");
		setSymmetry(8, 7, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 5, 7, 5, 6, 5, 0, 14, 0, 1, 23, 1, 21, 1, 2, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 218)
	{
		setName("219", "Td^5", "F-43c", "F -4 3 c", "F -4c 2 3");
		setSymmetry(8, 2, 5, 0, 16, 0, -1, -1, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 5, 13, 5, 0, 14, 0, 1, 17, 1, 13, 1, 16, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 219)
	{
		setName("220", "Td^6", "I-43d", "I -4 3 d", "I -4bd 2c 3");
		setSymmetry(8, 20, 5, 2, 16, 0, -1, -1, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 5, 12, 14, 0, 1, 48, 1, 49, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 220)
	{
		setName("221", "Oh^1", "Pm-3m", "P 4/m -3 2/m", "-P 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 10, 0, 6, 5, 6, 0, 15, 5, 15, 0, 5, 6, 14, 0, 5, 2, 5, 0, \
			1, 5, 1, 2, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 221)
	{
		setName("222", "Oh^2", "Pn-3n", "P 4/n -3 2/n", "P 4 2 3 -1n");
		setSymmetry(7, 0, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 15, 0, 5, 7, 14, 0, 5, 0, 1, 23, 1, 16, 1, 2, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 222)
	{
		setName("223", "Oh^3", "Pm-3n", "P 42/m -3 2/n", "-P 4n 2 3");
		setSymmetry(7, 7, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 6, 0, 15, 23, 14, 0, 5, 6, 5, 7, 5, 0, 1, 16, 1, 21, 1, 23, \
			1, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 223)
	{
		setName("224", "Oh^4", "Pn-3m", "P 42/n -3 2/m", "P 4n 2 3 -1n");
		setSymmetry(7, 7, 5, 0, 16, 0, 17, 7, LS_CUBIC, LC_P);
		setWyckoff(0, 0, 10, 0, 15, 23, 16, 23, 5, 7, 5, 0, 1, 23, 14, 0, 1, 2, 1, 18, \
			1, 16, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 224)
	{
		setName("225", "Oh^5", "Fm-3m", "F 4/m -3 2/m", "-F 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 10, 0, 6, 0, 15, 5, 15, 0, 5, 13, 14, 0, 5, 0, 1, 13, 1, 16, \
			1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 225)
	{
		setName("226", "Oh^6", "Fm-3c", "F 4/m -3 2/c", "-F 4c 2 3");
		setSymmetry(7, 2, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 6, 0, 15, 17, 14, 0, 5, 13, 5, 0, 1, 13, 1, 17, 1, 0, 1, 16, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 226)
	{
		setName("227", "Oh^7", "Fd-3m", "F 41/d -3 2/m", "F 4d 2 3 -1d");
		setSymmetry(7, 8, 5, 0, 16, 0, 17, 8, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 16, 44, 10, 0, 5, 0, 14, 0, 1, 24, 1, 25, 1, 1, 1, 0, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 227)
	{
		setName("228", "Oh^8", "Fd-3c", "F 41/d -3 2/c", "F 4d 2 3 -1cd");
		setSymmetry(7, 8, 5, 0, 16, 0, 17, 21, LS_CUBIC, LC_F);
		setWyckoff(0, 0, 16, 44, 5, 0, 14, 0, 1, 17, 1, 46, 1, 25, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 228)
	{
		setName("229", "Oh^9", "Im-3m", "I 4/m -3 2/m", "-I 4 2 3");
		setSymmetry(7, 0, 5, 0, 16, 0, 1, 0, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 10, 0, 6, 0, 16, 23, 15, 0, 5, 7, 14, 0, 5, 0, 1, 23, 1, 16, \
			1, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else if (number == 229)
	{
		setName("230", "Oh^10", "Ia-3d", "I 41/a -3 2/d", "-I 4bd 2c 3");
		setSymmetry(7, 20, 5, 2, 16, 0, 1, 0, LS_CUBIC, LC_I);
		setWyckoff(0, 0, 16, 44, 5, 12, 14, 0, 1, 49, 1, 44, 1, 25, 1, 0, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	}
	else
	{
		Output::newline(ERROR);
		Output::print("Internal error: space group number is out of range");
		Output::quit();
	}
}



/* void SpaceGroup::setName(const char* number, const char* schoe, const char* hmShort, const char* hmLong,
 *		const char* hall)
 *
 * Set space group name
 */

void SpaceGroup::setName(const char* number, const char* schoe, const char* hmShort, const char* hmLong, \
	const char* hall)
{
	_itcNumber = number;
	_schoenflies = schoe;
	_hermannMauguinShort = hmShort;
	_hermannMauguin = hmLong;
	_hall = hall;
}



/* void SpaceGroup::setSymmetry(int gen1, int vec1, int gen2, int vec2, int gen3, int vec3, int gen4, int vec4,
 *		LatticeSystem system, LatticeCentering centering)
 *
 * Set symmetry of space group
 */

void SpaceGroup::setSymmetry(int gen1, int vec1, int gen2, int vec2, int gen3, int vec3, int gen4, int vec4, \
	LatticeSystem system, LatticeCentering centering)
{
	
	// Make list of rotations and translations
	List<int> rots;
	List<int> vecs;
	if (gen1 != -1) { rots += gen1; vecs += vec1; }
	if (gen2 != -1) { rots += gen2; vecs += vec2; }
	if (gen3 != -1) { rots += gen3; vecs += vec3; }
	if (gen4 != -1) { rots += gen4; vecs += vec4; }
	
	// Loop over generators and save
	int i;
	Matrix3D curMat;
	Vector3D curVec;
	for (i = 0; i < rots.length(); ++i)
	{
		
		// Match indices to rotation and translation
		if (rots[i] == 0) curMat.set(1, 0, 0, 0, 1, 0, 0, 0, 1);
		else if (rots[i] == 1) curMat.set(-1, 0, 0, 0, -1, 0, 0, 0, -1);
		else if (rots[i] == 2) curMat.set(-1, 0, 0, 0, 1, 0, 0, 0, -1);
		else if (rots[i] == 3) curMat.set(1, 0, 0, 0, -1, 0, 0, 0, 1);
		else if (rots[i] == 4) curMat.set(-1, 0, 0, 0, -1, 0, 0, 0, 1);
		else if (rots[i] == 5) curMat.set(1, 0, 0, 0, -1, 0, 0, 0, -1);
		else if (rots[i] == 6) curMat.set(-1, 0, 0, 0, 1, 0, 0, 0, 1);
		else if (rots[i] == 7) curMat.set(0, -1, 0, 1, 0, 0, 0, 0, 1);
		else if (rots[i] == 8) curMat.set(0, 1, 0, -1, 0, 0, 0, 0, -1);
		else if (rots[i] == 9) curMat.set(0, -1, 0, 1, -1, 0, 0, 0, 1);
		else if (rots[i] == 10) curMat.set(0, -1, 0, -1, 0, 0, 0, 0, -1);
		else if (rots[i] == 11) curMat.set(0, 1, 0, 1, 0, 0, 0, 0, -1);
		else if (rots[i] == 12) curMat.set(0, -1, 0, -1, 0, 0, 0, 0, 1);
		else if (rots[i] == 13) curMat.set(0, 1, 0, 1, 0, 0, 0, 0, 1);
		else if (rots[i] == 14) curMat.set(1, -1, 0, 1, 0, 0, 0, 0, 1);
		else if (rots[i] == 15) curMat.set(-1, 1, 0, -1, 0, 0, 0, 0, -1);
		else if (rots[i] == 16) curMat.set(0, 0, 1, 1, 0, 0, 0, 1, 0);
		else if (rots[i] == 17) curMat.set(0, 0, -1, -1, 0, 0, 0, -1, 0);
		else
		{
			Output::newline(ERROR);
			Output::print("Internal error: Unknown rotation generator index");
			Output::quit();
		}
		if (vecs[i] == 0) curVec.set(0.0, 0.0, 0.0);
		else if (vecs[i] == 1) curVec.set(0.0, 1.0/2.0, 0.0);
		else if (vecs[i] == 2) curVec.set(0.0, 0.0, 1.0/2.0);
		else if (vecs[i] == 3) curVec.set(0.0, 1.0/2.0, 1.0/2.0);
		else if (vecs[i] == 4) curVec.set(1.0/2.0, 1.0/2.0, 0.0);
		else if (vecs[i] == 5) curVec.set(1.0/2.0, 0.0, 1.0/2.0);
		else if (vecs[i] == 6) curVec.set(1.0/2.0, 0.0, 0.0);
		else if (vecs[i] == 7) curVec.set(1.0/2.0, 1.0/2.0, 1.0/2.0);
		else if (vecs[i] == 8) curVec.set(1.0/4.0, 1.0/4.0, 1.0/4.0);
		else if (vecs[i] == 9) curVec.set(0.0, 0.0, 1.0/4.0);
		else if (vecs[i] == 10) curVec.set(0.0, 0.0, 3.0/4.0);
		else if (vecs[i] == 11) curVec.set(0.0, 1.0/2.0, 1.0/4.0);
		else if (vecs[i] == 12) curVec.set(1.0/2.0, 1.0/2.0, 1.0/4.0);
		else if (vecs[i] == 13) curVec.set(1.0/2.0, 1.0/2.0, 3.0/4.0);
		else if (vecs[i] == 14) curVec.set(1.0/2.0, 0.0, 1.0/4.0);
		else if (vecs[i] == 15) curVec.set(0.0, 0.0, 1.0/3.0);
		else if (vecs[i] == 16) curVec.set(0.0, 0.0, 2.0/3.0);
		else if (vecs[i] == 17) curVec.set(0.0, 0.0, 1.0/6.0);
		else if (vecs[i] == 18) curVec.set(0.0, 0.0, 5.0/6.0);
		else if (vecs[i] == 19) curVec.set(3.0/4.0, 1.0/4.0, 3.0/4.0);
		else if (vecs[i] == 20) curVec.set(1.0/4.0, 3.0/4.0, 1.0/4.0);
		else if (vecs[i] == 21) curVec.set(1.0/4.0, 1.0/4.0, 3.0/4.0);
		else
		{
			Output::newline(ERROR);
			Output::print("Internal error: Unknown translation generator index");
			Output::quit();
		}
		
		// Add symmetry
		addSymmetry(curMat, curVec);
	}
	
	// Move identity to beginning
	Matrix3D identity = Matrix3D::identity();
	for (i = 0; i < _symmetry.length(); ++i)
	{
		if (_symmetry[i].rotation() == identity)
		{
			_symmetry.swap(0, i);
			break;
		}
	}
	
	// Save system and centering
	_system = system;
	_centering = centering;
	
	// Add centering vectors
	OList<Vector3D > centVecs;
	if (centering == LC_C)
		centVecs += Vector3D(0.5, 0.5, 0);
	else if (centering == LC_A)
		centVecs += Vector3D(0, 0.5, 0.5);
	else if (centering == LC_B)
		centVecs += Vector3D(0.5, 0, 0.5);
	else if (centering == LC_I)
		centVecs += Vector3D(0.5, 0.5, 0.5);
	else if (centering == LC_F)
	{
		centVecs += Vector3D(0, 0.5, 0.5);
		centVecs += Vector3D(0.5, 0, 0.5);
		centVecs += Vector3D(0.5, 0.5, 0);
	}
	else if (centering == LC_R)
	{
		centVecs += Vector3D(2.0/3.0, 1.0/3.0, 1.0/3.0);
		centVecs += Vector3D(1.0/3.0, 2.0/3.0, 2.0/3.0);
	}
	else if (centering != LC_P)
	{
		Output::newline(ERROR);
		Output::print("Unknown lattice centering in space group determination");
		Output::quit();
	}
	
	// Loop over symmetry operations
	int j;
	for (i = 0; i < _symmetry.length(); ++i)
	{
		
		// Loop over centering vectors and add new ones
		curVec = _symmetry[i].translations()[0];
		for (j = 0; j < centVecs.length(); ++j)
			_symmetry[i].addTranslation(curVec + centVecs[j]);
	}
}



/* void SpaceGroup::addSymmetry(const Matrix3D& rotation, const Vector3D& translation)
 *
 * Add symmetry operations to space group
 */

void SpaceGroup::addSymmetry(const Matrix3D& rotation, const Vector3D& translation)
{
	
	// Loop over operations to check if known
	int i;
	for (i = 0; i < _symmetry.length(); ++i)
	{
		if (rotation == _symmetry[i].rotation())
			return;
	}
	
	// Save operation
	_symmetry.add();
	_symmetry.last().setRotation(rotation);
	_symmetry.last().addTranslation(translation);
	
	// Generate new operations
	Vector3D newVec;
	for (i = 0; i < _symmetry.length(); ++i)
	{
		
		// Forward multiply
		newVec = rotation * _symmetry[i].translations()[0] + translation;
		ISO::moveIntoCell(newVec);
		addSymmetry(rotation * _symmetry[i].rotation(), newVec);
		
		// Reverse multiply
		newVec = _symmetry[i].rotation() * translation + _symmetry[i].translations()[0];
		ISO::moveIntoCell(newVec);
		addSymmetry(_symmetry[i].rotation() * rotation, newVec);
	}
}



/* void SpaceGroup::setWyckoff(int wr1, int wt1, int wr2, int wt2, int wr3, int wt3, int wr4, int wt4, int wr5,
 *		int wt5, int wr6, int wt6, int wr7, int wt7, int wr8, int wt8, int wr9, int wt9, int wr10, int wt10, int wr11,
 *		int wt11, int wr12, int wt12, int wr13, int wt13, int wr14, int wt14, int wr15, int wt15, int wr16, int wt16,
 *		int wr17, int wt17, int wr18, int wt18, int wr19, int wt19, int wr20, int wt20, int wr21, int wt21, int wr22,
 *		int wt22, int wr23, int wt23, int wr24, int wt24, int wr25, int wt25, int wr26, int wt26, int wr27, int wt27)
 *
 * Set the Wyckoff positions for a space group
 */

void SpaceGroup::setWyckoff(int wr1, int wt1, int wr2, int wt2, int wr3, int wt3, int wr4, int wt4, int wr5, \
	int wt5, int wr6, int wt6, int wr7, int wt7, int wr8, int wt8, int wr9, int wt9, int wr10, int wt10, int wr11, \
	int wt11, int wr12, int wt12, int wr13, int wt13, int wr14, int wt14, int wr15, int wt15, int wr16, int wt16, \
	int wr17, int wt17, int wr18, int wt18, int wr19, int wt19, int wr20, int wt20, int wr21, int wt21, int wr22, \
	int wt22, int wr23, int wt23, int wr24, int wt24, int wr25, int wt25, int wr26, int wt26, int wr27, int wt27)
{
	
	// Save the Wyckoff position indicies
	List<int> rotIndex;
	List<int> vecIndex;
	if (wr1  != -1) { rotIndex += wr1;  vecIndex += wt1;  }
	if (wr2  != -1) { rotIndex += wr2;  vecIndex += wt2;  }
	if (wr3  != -1) { rotIndex += wr3;  vecIndex += wt3;  }
	if (wr4  != -1) { rotIndex += wr4;  vecIndex += wt4;  }
	if (wr5  != -1) { rotIndex += wr5;  vecIndex += wt5;  }
	if (wr6  != -1) { rotIndex += wr6;  vecIndex += wt6;  }
	if (wr7  != -1) { rotIndex += wr7;  vecIndex += wt7;  }
	if (wr8  != -1) { rotIndex += wr8;  vecIndex += wt8;  }
	if (wr9  != -1) { rotIndex += wr9;  vecIndex += wt9;  }
	if (wr10 != -1) { rotIndex += wr10; vecIndex += wt10; }
	if (wr11 != -1) { rotIndex += wr11; vecIndex += wt11; }
	if (wr12 != -1) { rotIndex += wr12; vecIndex += wt12; }
	if (wr13 != -1) { rotIndex += wr13; vecIndex += wt13; }
	if (wr14 != -1) { rotIndex += wr14; vecIndex += wt14; }
	if (wr15 != -1) { rotIndex += wr15; vecIndex += wt15; }
	if (wr16 != -1) { rotIndex += wr16; vecIndex += wt16; }
	if (wr17 != -1) { rotIndex += wr17; vecIndex += wt17; }
	if (wr18 != -1) { rotIndex += wr18; vecIndex += wt18; }
	if (wr19 != -1) { rotIndex += wr19; vecIndex += wt19; }
	if (wr20 != -1) { rotIndex += wr20; vecIndex += wt20; }
	if (wr21 != -1) { rotIndex += wr21; vecIndex += wt21; }
	if (wr22 != -1) { rotIndex += wr22; vecIndex += wt22; }
	if (wr23 != -1) { rotIndex += wr23; vecIndex += wt23; }
	if (wr24 != -1) { rotIndex += wr24; vecIndex += wt24; }
	if (wr25 != -1) { rotIndex += wr25; vecIndex += wt25; }
	if (wr26 != -1) { rotIndex += wr26; vecIndex += wt26; }
	if (wr27 != -1) { rotIndex += wr27; vecIndex += wt27; }
	
	// Set Wyckoff positions
	_wyckoff.length(rotIndex.length());
	for (int i = 0; i < rotIndex.length(); i++)
		_wyckoff[i].set(rotIndex[i], vecIndex[i], _symmetry, rotIndex.length() - i - 1);
}



/* void SpaceGroup::print() const
 *
 * Print general information about space group
 */

void SpaceGroup::print() const
{
	
	// Save current stream
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);

	// Print title line
	Output::newline();
	Output::print("Space group information");

	// Make output object
	Output message;
	message.addLines(7);

	// ITC number
	message.addLine();
	message.add("    ITC number:");
	message.add(_itcNumber);

	// Short Hermann-Mauguin symbol
	message.addLine();
	message.add("    Hermann-Mauguin short:");
	message.add(_hermannMauguinShort);

	// Full Hermann-Mauguin symbol
	message.addLine();
	message.add("    Hermann-Mauguin full:");
	message.add(_hermannMauguin);

	// Schoenflies
	message.addLine();
	message.add("    Schoenflies:");
	message.add(_schoenflies);

	// Hall
	message.addLine();
	message.add("    Hall:");
	message.add(_hall);
	
	// Lattice system
	message.addLine();
	message.add("   Lattice system:");
	message.add(ISO::system(_system));
	
	// Lattice centering
	message.addLine();
	message.add("   Lattice centering:");
	message.add(ISO::centering(_centering));

	// Print information
	List<PrintAlign> align (2);
	align[0] = RIGHT;
	align[1] = LEFT;
	Output::newline();
	Output::print(message, align);
	
	// Reset method
	Output::method(origMethod);
}



/* void SpaceGroup::printComplete() const
 *
 * Print all information about space group
 */

void SpaceGroup::printComplete() const
{
	
	// Print general information
	print();
	
	// Output
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Create message to store symmetry operations
	int i;
	Output message;
	message.addLines(_symmetry.length());
	for (i = 0; i < _symmetry.length(); ++i)
	{
		message.addLine();
		message.add("   ");
		message.add(_symmetry[i].getString());
	}
	
	// Print rotations
	Output::newline();
	Output::print("Number of unique symmetry operations: ");
	Output::print(_symmetry.length());
	Output::newline();
	Output::print(message, LEFT);
	
	// Create message to store translations
	int j;
	message.clear();
	if (_symmetry.length())
	{
		message.addLines(_symmetry[0].translations().length());
		for (i = 0; i < _symmetry[0].translations().length(); ++i)
		{
			message.addLine();
			message.addWords(4);
			message.add("   ");
			for (j = 0; j < 3; ++j)
				message.add(Language::numberToFraction(_symmetry[0].translations()[i][j]));
		}
	}
	
	// Print translations
	Output::newline();
	Output::newline();
	Output::print("Number of centering vectors: ");
	Output::print(message.numLines());
	Output::newline();
	Output::print(message, RIGHT);
	
	// Print number of Wyckoff positions
	Output::newline();
	Output::newline();
	Output::print("Number of Wyckoff positions: ");
	Output::print(_wyckoff.length());
	
	// Loop over Wyckoff positions and add
	message.clear();
	message.addLines(_wyckoff.length());
	for (i = 0; i < _wyckoff.length(); i++)
	{
		message.addLine();
		message.addWords(5);
		message.add("   ");
		message.add(_wyckoff[i].name() + ":");
		message.add(_wyckoff[i].getString(0));
	}
	
	// Set alignment
	List<PrintAlign> align (6);
	align[0] = align[1] = align[2] = RIGHT;
	align[3] = align[4] = align[5] = LEFT;
	
	// Print Wyckoff positions
	Output::newline();
	Output::print(message, align);
	
	// Print all positions
	for (i = 0; i < _wyckoff.length(); ++i)
		_wyckoff[i].print();
	
	// Output
	Output::method(origMethod);
}



/* void SpaceGroup::printAll()
 *
 * Print all space groups
 */

void SpaceGroup::printAll()
{
	
	// Setup output
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print space groups
	Output::newline(); Output::print("ITC Number   HM Short            HM Full   Schoenflies               Hall");
	Output::newline(); Output::print("----------   --------            -------   -----------               ----");
	Output::newline(); Output::print("         1         P1                P 1          C1^1                P 1");
	Output::newline(); Output::print("         2        P-1               P -1          Ci^1               -P 1");
	Output::newline(); Output::print("         3         P2            P 1 2 1          C2^1               P 2y");
	Output::newline(); Output::print("         4        P21           P 1 21 1          C2^2              P 2yb");
	Output::newline(); Output::print("         5         C2            C 1 2 1          C2^3               C 2y");
	Output::newline(); Output::print("         6         Pm            P 1 m 1          Cs^1              P -2y");
	Output::newline(); Output::print("         7         Pc            P 1 c 1          Cs^2             P -2yc");
	Output::newline(); Output::print("         8         Cm            C 1 m 1          Cs^3              C -2y");
	Output::newline(); Output::print("         9         Cc            C 1 c 1          Cs^4             C -2yc");
	Output::newline(); Output::print("        10       P2/m          P 1 2/m 1         C2h^1              -P 2y");
	Output::newline(); Output::print("        11      P21/m         P 1 21/m 1         C2h^2             -P 2yb");
	Output::newline(); Output::print("        12       C2/m          C 1 2/m 1         C2h^3              -C 2y");
	Output::newline(); Output::print("        13       P2/c          P 1 2/c 1         C2h^4             -P 2yc");
	Output::newline(); Output::print("        14      P21/c         P 1 21/c 1         C2h^5            -P 2ybc");
	Output::newline(); Output::print("        15       C2/c          C 1 2/c 1         C2h^6             -C 2yc");
	Output::newline(); Output::print("        16       P222            P 2 2 2          D2^1              P 2 2");
	Output::newline(); Output::print("        17      P2221           P 2 2 21          D2^2             P 2c 2");
	Output::newline(); Output::print("        18     P21212          P 21 21 2          D2^3            P 2 2ab");
	Output::newline(); Output::print("        19    P212121         P 21 21 21          D2^4          P 2ac 2ab");
	Output::newline(); Output::print("        20      C2221           C 2 2 21          D2^5             C 2c 2");
	Output::newline(); Output::print("        21       C222            C 2 2 2          D2^6              C 2 2");
	Output::newline(); Output::print("        22       F222            F 2 2 2          D2^7              F 2 2");
	Output::newline(); Output::print("        23       I222            I 2 2 2          D2^8              I 2 2");
	Output::newline(); Output::print("        24    I212121         I 21 21 21          D2^9            I 2b 2c");
	Output::newline(); Output::print("        25       Pmm2            P m m 2         C2v^1             P 2 -2");
	Output::newline(); Output::print("        26      Pmc21           P m c 21         C2v^2            P 2c -2");
	Output::newline(); Output::print("        27       Pcc2            P c c 2         C2v^3            P 2 -2c");
	Output::newline(); Output::print("        28       Pma2            P m a 2         C2v^4            P 2 -2a");
	Output::newline(); Output::print("        29      Pca21           P c a 21         C2v^5          P 2c -2ac");
	Output::newline(); Output::print("        30       Pnc2            P n c 2         C2v^6           P 2 -2bc");
	Output::newline(); Output::print("        31      Pmn21           P m n 21         C2v^7           P 2ac -2");
	Output::newline(); Output::print("        32       Pba2            P b a 2         C2v^8           P 2 -2ab");
	Output::newline(); Output::print("        33      Pna21           P n a 21         C2v^9           P 2c -2n");
	Output::newline(); Output::print("        34       Pnn2            P n n 2        C2v^10            P 2 -2n");
	Output::newline(); Output::print("        35       Cmm2            C m m 2        C2v^11             C 2 -2");
	Output::newline(); Output::print("        36      Cmc21           C m c 21        C2v^12            C 2c -2");
	Output::newline(); Output::print("        37       Ccc2            C c c 2        C2v^13            C 2 -2c");
	Output::newline(); Output::print("        38       Amm2            A m m 2        C2v^14             A 2 -2");
	Output::newline(); Output::print("        39       Abm2            A b m 2        C2v^15            A 2 -2c");
	Output::newline(); Output::print("        40       Ama2            A m a 2        C2v^16            A 2 -2a");
	Output::newline(); Output::print("        41       Aba2            A b a 2        C2v^17           A 2 -2ac");
	Output::newline(); Output::print("        42       Fmm2            F m m 2        C2v^18             F 2 -2");
	Output::newline(); Output::print("        43       Fdd2            F d d 2        C2v^19            F 2 -2d");
	Output::newline(); Output::print("        44       Imm2            I m m 2        C2v^20             I 2 -2");
	Output::newline(); Output::print("        45       Iba2            I b a 2        C2v^21            I 2 -2c");
	Output::newline(); Output::print("        46       Ima2            I m a 2        C2v^22            I 2 -2a");
	Output::newline(); Output::print("        47       Pmmm      P 2/m 2/m 2/m         D2h^1             -P 2 2");
	Output::newline(); Output::print("        48       Pnnn      P 2/n 2/n 2/n         D2h^2          P 2 2 -1n");
	Output::newline(); Output::print("        49       Pccm      P 2/c 2/c 2/m         D2h^3            -P 2 2c");
	Output::newline(); Output::print("        50       Pban      P 2/b 2/a 2/n         D2h^4         P 2 2 -1ab");
	Output::newline(); Output::print("        51       Pmma     P 21/m 2/m 2/a         D2h^5           -P 2a 2a");
	Output::newline(); Output::print("        52       Pnna     P 2/n 21/n 2/a         D2h^6          -P 2a 2bc");
	Output::newline(); Output::print("        53       Pmna     P 2/m 2/n 21/a         D2h^7           -P 2ac 2");
	Output::newline(); Output::print("        54       Pcca     P 21/c 2/c 2/a         D2h^8          -P 2a 2ac");
	Output::newline(); Output::print("        55       Pbam    P 21/b 21/a 2/m         D2h^9           -P 2 2ab");
	Output::newline(); Output::print("        56       Pccn    P 21/c 21/c 2/n        D2h^10         -P 2ab 2ac");
	Output::newline(); Output::print("        57       Pbcm    P 2/b 21/c 21/m        D2h^11           -P 2c 2b");
	Output::newline(); Output::print("        58       Pnnm    P 21/n 21/n 2/m        D2h^12            -P 2 2n");
	Output::newline(); Output::print("        59       Pmmn    P 21/m 21/m 2/n        D2h^13       P 2 2ab -1ab");
	Output::newline(); Output::print("        60       Pbcn    P 21/b 2/c 21/n        D2h^14          -P 2n 2ab");
	Output::newline(); Output::print("        61       Pbca   P 21/b 21/c 21/a        D2h^15         -P 2ac 2ab");
	Output::newline(); Output::print("        62       Pnma   P 21/n 21/m 21/a        D2h^16          -P 2ac 2n");
	Output::newline(); Output::print("        63       Cmcm     C 2/m 2/c 21/m        D2h^17            -C 2c 2");
	Output::newline(); Output::print("        64       Cmca     C 2/m 2/c 21/a        D2h^18           -C 2bc 2");
	Output::newline(); Output::print("        65       Cmmm      C 2/m 2/m 2/m        D2h^19             -C 2 2");
	Output::newline(); Output::print("        66       Cccm      C 2/c 2/c 2/m        D2h^20            -C 2 2c");
	Output::newline(); Output::print("        67       Cmma      C 2/m 2/m 2/a        D2h^21            -C 2b 2");
	Output::newline(); Output::print("        68       Ccca      C 2/c 2/c 2/a        D2h^22         C 2 2 -1bc");
	Output::newline(); Output::print("        69       Fmmm      F 2/m 2/m 2/m        D2h^23             -F 2 2");
	Output::newline(); Output::print("        70       Fddd      F 2/d 2/d 2/d        D2h^24          F 2 2 -1d");
	Output::newline(); Output::print("        71       Immm      I 2/m 2/m 2/m        D2h^25             -I 2 2");
	Output::newline(); Output::print("        72       Ibam      I 2/b 2/a 2/m        D2h^26            -I 2 2c");
	Output::newline(); Output::print("        73       Ibca      I 2/b 2/c 2/a        D2h^27           -I 2b 2c");
	Output::newline(); Output::print("        74       Imma      I 2/m 2/m 2/a        D2h^28            -I 2b 2");
	Output::newline(); Output::print("        75         P4            P 4 1 1          C4^1                P 4");
	Output::newline(); Output::print("        76        P41           P 41 1 1          C4^2               P 4w");
	Output::newline(); Output::print("        77        P42           P 42 1 1          C4^3               P 4c");
	Output::newline(); Output::print("        78        P43           P 43 1 1          C4^4              P 4cw");
	Output::newline(); Output::print("        79         I4            I 4 1 1          C4^5                I 4");
	Output::newline(); Output::print("        80        I41           I 41 1 1          C4^6              I 4bw");
	Output::newline(); Output::print("        81        P-4           P -4 1 1          S4^1               P -4");
	Output::newline(); Output::print("        82        I-4           I -4 1 1          S4^2               I -4");
	Output::newline(); Output::print("        83       P4/m          P 4/m 1 1         C4h^1               -P 4");
	Output::newline(); Output::print("        84      P42/m         P 42/m 1 1         C4h^2              -P 4c");
	Output::newline(); Output::print("        85       P4/n          P 4/n 1 1         C4h^3         P 4ab -1ab");
	Output::newline(); Output::print("        86      P42/n         P 42/n 1 1         C4h^4           P 4n -1n");
	Output::newline(); Output::print("        87       I4/m          I 4/m 1 1         C4h^5               -I 4");
	Output::newline(); Output::print("        88      I41/a         I 41/a 1 1         C4h^6         I 4bw -1bw");
	Output::newline(); Output::print("        89       P422            P 4 2 2          D4^1              P 4 2");
	Output::newline(); Output::print("        90      P4212           P 4 21 2          D4^2          P 4ab 2ab");
	Output::newline(); Output::print("        91      P4122           P 41 2 2          D4^3            P 4w 2c");
	Output::newline(); Output::print("        92     P41212          P 41 21 2          D4^4         P 4abw 2nw");
	Output::newline(); Output::print("        93      P4222           P 42 2 2          D4^5             P 4c 2");
	Output::newline(); Output::print("        94     P42212          P 42 21 2          D4^6            P 4n 2n");
	Output::newline(); Output::print("        95      P4322           P 43 2 2          D4^7           P 4cw 2c");
	Output::newline(); Output::print("        96     P43212          P 43 21 2          D4^8         P 4nw 2abw");
	Output::newline(); Output::print("        97       I422            I 4 2 2          D4^9              I 4 2");
	Output::newline(); Output::print("        98      I4122           I 41 2 2         D4^10          I 4bw 2bw");
	Output::newline(); Output::print("        99       P4mm            P 4 m m         C4v^1             P 4 -2");
	Output::newline(); Output::print("       100       P4bm            P 4 b m         C4v^2           P 4 -2ab");
	Output::newline(); Output::print("       101      P42cm           P 42 c m         C4v^3           P 4c -2c");
	Output::newline(); Output::print("       102      P42nm           P 42 n m         C4v^4           P 4n -2n");
	Output::newline(); Output::print("       103       P4cc            P 4 c c         C4v^5            P 4 -2c");
	Output::newline(); Output::print("       104       P4nc            P 4 n c         C4v^6            P 4 -2n");
	Output::newline(); Output::print("       105      P42mc           P 42 m c         C4v^7            P 4c -2");
	Output::newline(); Output::print("       106      P42bc           P 42 b c         C4v^8          P 4c -2ab");
	Output::newline(); Output::print("       107       I4mm            I 4 m m         C4v^9             I 4 -2");
	Output::newline(); Output::print("       108       I4cm            I 4 c m        C4v^10            I 4 -2c");
	Output::newline(); Output::print("       109      I41md           I 41 m d        C4v^11           I 4bw -2");
	Output::newline(); Output::print("       110      I41cd           I 41 c d        C4v^12          I 4bw -2c");
	Output::newline(); Output::print("       111      P-42m           P -4 2 m         D2d^1             P -4 2");
	Output::newline(); Output::print("       112      P-42c           P -4 2 c         D2d^2            P -4 2c");
	Output::newline(); Output::print("       113     P-421m          P -4 21 m         D2d^3           P -4 2ab");
	Output::newline(); Output::print("       114     P-421c          P -4 21 c         D2d^4            P -4 2n");
	Output::newline(); Output::print("       115      P-4m2           P -4 m 2         D2d^5            P -4 -2");
	Output::newline(); Output::print("       116      P-4c2           P -4 c 2         D2d^6           P -4 -2c");
	Output::newline(); Output::print("       117      P-4b2           P -4 b 2         D2d^7          P -4 -2ab");
	Output::newline(); Output::print("       118      P-4n2           P -4 n 2         D2d^8           P -4 -2n");
	Output::newline(); Output::print("       119      I-4m2           I -4 m 2         D2d^9            I -4 -2");
	Output::newline(); Output::print("       120      I-4c2           I -4 c 2        D2d^10           I -4 -2c");
	Output::newline(); Output::print("       121      I-42m           I -4 2 m        D2d^11             I -4 2");
	Output::newline(); Output::print("       122      I-42d           I -4 2 d        D2d^12           I -4 2bw");
	Output::newline(); Output::print("       123     P4/mmm      P 4/m 2/m 2/m         D4h^1             -P 4 2");
	Output::newline(); Output::print("       124     P4/mcc      P 4/m 2/c 2/c         D4h^2            -P 4 2c");
	Output::newline(); Output::print("       125     P4/nbm      P 4/n 2/b 2/m         D4h^3         P 4 2 -1ab");
	Output::newline(); Output::print("       126     P4/nnc      P 4/n 2/n 2/c         D4h^4          P 4 2 -1n");
	Output::newline(); Output::print("       127     P4/mbm       P 4/m 21/b m         D4h^5           -P 4 2ab");
	Output::newline(); Output::print("       128     P4/mnc       P 4/m 21/n c         D4h^6            -P 4 2n");
	Output::newline(); Output::print("       129     P4/nmm       P 4/n 21/m m         D4h^7     P 4ab 2ab -1ab");
	Output::newline(); Output::print("       130     P4/ncc       P 4/n 21/c c         D4h^8      P 4ab 2n -1ab");
	Output::newline(); Output::print("       131    P42/mmc     P 42/m 2/m 2/c         D4h^9            -P 4c 2");
	Output::newline(); Output::print("       132    P42/mcm     P 42/m 2/c 2/m        D4h^10           -P 4c 2c");
	Output::newline(); Output::print("       133    P42/nbc     P 42/n 2/b 2/c        D4h^11        P 4n 2c -1n");
	Output::newline(); Output::print("       134    P42/nnm     P 42/n 2/n 2/m        D4h^12         P 4n 2 -1n");
	Output::newline(); Output::print("       135    P42/mbc    P 42/m 21/b 2/c        D4h^13          -P 4c 2ab");
	Output::newline(); Output::print("       136    P42/mnm    P 42/m 21/n 2/m        D4h^14           -P 4n 2n");
	Output::newline(); Output::print("       137    P42/nmc    P 42/n 21/m 2/c        D4h^15        P 4n 2n -1n");
	Output::newline(); Output::print("       138    P42/ncm    P 42/n 21/c 2/m        D4h^16       P 4n 2ab -1n");
	Output::newline(); Output::print("       139     I4/mmm      I 4/m 2/m 2/m        D4h^17             -I 4 2");
	Output::newline(); Output::print("       140     I4/mcm      I 4/m 2/c 2/m        D4h^18            -I 4 2c");
	Output::newline(); Output::print("       141    I41/amd     I 41/a 2/m 2/d        D4h^19     I 4bw 2bw -1bw");
	Output::newline(); Output::print("       142    I41/acd     I 41/a 2/c 2/d        D4h^20     I 4bw 2aw -1bw");
	Output::newline(); Output::print("       143         P3            P 3 1 1          C3^1                P 3");
	Output::newline(); Output::print("       144        P31           P 31 1 1          C3^2               P 31");
	Output::newline(); Output::print("       145        P32           P 32 1 1          C3^3               P 32");
	Output::newline(); Output::print("       146         R3            R 3 1 1          C3^4                R 3");
	Output::newline(); Output::print("       147        P-3           P -3 1 1         C3i^1               -P 3");
	Output::newline(); Output::print("       148        R-3           R -3 1 1         C3i^2               -R 3");
	Output::newline(); Output::print("       149       P312            P 3 1 2          D3^1              P 3 2");
	Output::newline(); Output::print("       150       P321            P 3 2 1          D3^2             P 3 2\"");
	Output::newline(); Output::print("       151      P3112           P 31 1 2          D3^3    P 31 2c (0 0 1)");
	Output::newline(); Output::print("       152      P3121           P 31 2 1          D3^4            P 31 2\"");
	Output::newline(); Output::print("       153      P3212           P 32 1 2          D3^5   P 32 2c (0 0 -1)");
	Output::newline(); Output::print("       154      P3221           P 32 2 1          D3^6            P 32 2\"");
	Output::newline(); Output::print("       155        R32            R 3 2 1          D3^7             R 3 2\"");
	Output::newline(); Output::print("       156       P3m1            P 3 m 1         C3v^1            P 3 -2\"");
	Output::newline(); Output::print("       157       P31m            P 3 1 m         C3v^2             P 3 -2");
	Output::newline(); Output::print("       158       P3c1            P 3 c 1         C3v^3           P 3 -2\"c");
	Output::newline(); Output::print("       159       P31c            P 3 1 c         C3v^4            P 3 -2c");
	Output::newline(); Output::print("       160        R3m            R 3 m 1         C3v^5            R 3 -2\"");
	Output::newline(); Output::print("       161        R3c            R 3 c 1         C3v^6           R 3 -2\"c");
	Output::newline(); Output::print("       162      P-31m         P -3 1 2/m         D3d^1             -P 3 2");
	Output::newline(); Output::print("       163      P-31c         P -3 1 2/c         D3d^2            -P 3 2c");
	Output::newline(); Output::print("       164      P-3m1         P -3 2/m 1         D3d^3            -P 3 2\"");
	Output::newline(); Output::print("       165      P-3c1         P -3 2/c 1         D3d^4           -P 3 2\"c");
	Output::newline(); Output::print("       166       R-3m         R -3 2/m 1         D3d^5            -R 3 2\"");
	Output::newline(); Output::print("       167       R-3c         R -3 2/c 1         D3d^6           -R 3 2\"c");
	Output::newline(); Output::print("       168         P6            P 6 1 1          C6^1                P 6");
	Output::newline(); Output::print("       169        P61           P 61 1 1          C6^2               P 61");
	Output::newline(); Output::print("       170        P65           P 65 1 1          C6^3               P 65");
	Output::newline(); Output::print("       171        P62           P 62 1 1          C6^4               P 62");
	Output::newline(); Output::print("       172        P64           P 64 1 1          C6^5               P 64");
	Output::newline(); Output::print("       173        P63           P 63 1 1          C6^6               P 6c");
	Output::newline(); Output::print("       174        P-6           P -6 1 1         C3h^1               P -6");
	Output::newline(); Output::print("       175       P6/m          P 6/m 1 1         C6h^1               -P 6");
	Output::newline(); Output::print("       176      P63/m         P 63/m 1 1         C6h^2              -P 6c");
	Output::newline(); Output::print("       177       P622            P 6 2 2          D6^1              P 6 2");
	Output::newline(); Output::print("       178      P6122           P 61 2 2          D6^2    P 61 2 (0 0 -1)");
	Output::newline(); Output::print("       179      P6522           P 65 2 2          D6^3     P 65 2 (0 0 1)");
	Output::newline(); Output::print("       180      P6222           P 62 2 2          D6^4    P 62 2c (0 0 1)");
	Output::newline(); Output::print("       181      P6422           P 64 2 2          D6^5   P 64 2c (0 0 -1)");
	Output::newline(); Output::print("       182      P6322           P 63 2 2          D6^6            P 6c 2c");
	Output::newline(); Output::print("       183       P6mm            P 6 m m         C6v^1             P 6 -2");
	Output::newline(); Output::print("       184       P6cc            P 6 c c         C6v^2            P 6 -2c");
	Output::newline(); Output::print("       185      P63cm           P 63 c m         C6v^3            P 6c -2");
	Output::newline(); Output::print("       186      P63mc           P 63 m c         C6v^4           P 6c -2c");
	Output::newline(); Output::print("       187      P-6m2           P -6 m 2         D3h^1             P -6 2");
	Output::newline(); Output::print("       188      P-6c2           P -6 c 2         D3h^2            P -6c 2");
	Output::newline(); Output::print("       189      P-62m           P -6 2 m         D3h^3            P -6 -2");
	Output::newline(); Output::print("       190      P-62c           P -6 2 c         D3h^4          P -6c -2c");
	Output::newline(); Output::print("       191     P6/mmm      P 6/m 2/m 2/m         D6h^1             -P 6 2");
	Output::newline(); Output::print("       192     P6/mcc      P 6/m 2/c 2/c         D6h^2            -P 6 2c");
	Output::newline(); Output::print("       193    P63/mcm     P 63/m 2/c 2/m         D6h^3            -P 6c 2");
	Output::newline(); Output::print("       194    P63/mmc     P 63/m 2/m 2/c         D6h^4           -P 6c 2c");
	Output::newline(); Output::print("       195        P23            P 2 3 1           T^1            P 2 2 3");
	Output::newline(); Output::print("       196        F23            F 2 3 1           T^2            F 2 2 3");
	Output::newline(); Output::print("       197        I23            I 2 3 1           T^3            I 2 2 3");
	Output::newline(); Output::print("       198       P213           P 21 3 1           T^4        P 2ac 2ab 3");
	Output::newline(); Output::print("       199       I213           I 21 3 1           T^5          I 2b 2c 3");
	Output::newline(); Output::print("       200       Pm-3         P 2/m -3 1          Th^1           -P 2 2 3");
	Output::newline(); Output::print("       201       Pn-3         P 2/n -3 1          Th^2        P 2 2 3 -1n");
	Output::newline(); Output::print("       202       Fm-3         F 2/m -3 1          Th^3           -F 2 2 3");
	Output::newline(); Output::print("       203       Fd-3         F 2/d -3 1          Th^4        F 2 2 3 -1d");
	Output::newline(); Output::print("       204       Im-3         I 2/m -3 1          Th^5           -I 2 2 3");
	Output::newline(); Output::print("       205       Pa-3        P 21/a -3 1          Th^6       -P 2ac 2ab 3");
	Output::newline(); Output::print("       206       Ia-3        I 21/a -3 1          Th^7         -I 2b 2c 3");
	Output::newline(); Output::print("       207       P432            P 4 3 2           O^1            P 4 2 3");
	Output::newline(); Output::print("       208      P4232           P 42 3 2           O^2           P 4n 2 3");
	Output::newline(); Output::print("       209       F432            F 4 3 2           O^3            F 4 2 3");
	Output::newline(); Output::print("       210      F4132           F 41 3 2           O^4           F 4d 2 3");
	Output::newline(); Output::print("       211       I432            I 4 3 2           O^5            I 4 2 3");
	Output::newline(); Output::print("       212      P4332           P 43 3 2           O^6       P 4acd 2ab 3");
	Output::newline(); Output::print("       213      P4132           P 41 3 2           O^7        P 4bd 2ab 3");
	Output::newline(); Output::print("       214      I4132           I 41 3 2           O^8         I 4bd 2c 3");
	Output::newline(); Output::print("       215      P-43m           P -4 3 m          Td^1           P -4 2 3");
	Output::newline(); Output::print("       216      F-43m           F -4 3 m          Td^2           F -4 2 3");
	Output::newline(); Output::print("       217      I-43m           I -4 3 m          Td^3           I -4 2 3");
	Output::newline(); Output::print("       218      P-43n           P -4 3 n          Td^4          P -4n 2 3");
	Output::newline(); Output::print("       219      F-43c           F -4 3 c          Td^5          F -4c 2 3");
	Output::newline(); Output::print("       220      I-43d           I -4 3 d          Td^6        I -4bd 2c 3");
	Output::newline(); Output::print("       221      Pm-3m       P 4/m -3 2/m          Oh^1           -P 4 2 3");
	Output::newline(); Output::print("       222      Pn-3n       P 4/n -3 2/n          Oh^2        P 4 2 3 -1n");
	Output::newline(); Output::print("       223      Pm-3n      P 42/m -3 2/n          Oh^3          -P 4n 2 3");
	Output::newline(); Output::print("       224      Pn-3m      P 42/n -3 2/m          Oh^4       P 4n 2 3 -1n");
	Output::newline(); Output::print("       225      Fm-3m       F 4/m -3 2/m          Oh^5           -F 4 2 3");
	Output::newline(); Output::print("       226      Fm-3c       F 4/m -3 2/c          Oh^6          -F 4c 2 3");
	Output::newline(); Output::print("       227      Fd-3m      F 41/d -3 2/m          Oh^7       F 4d 2 3 -1d");
	Output::newline(); Output::print("       228      Fd-3c      F 41/d -3 2/c          Oh^8      F 4d 2 3 -1cd");
	Output::newline(); Output::print("       229      Im-3m       I 4/m -3 2/m          Oh^9           -I 4 2 3");
	Output::newline(); Output::print("       230      Ia-3d      I 41/a -3 2/d         Oh^10        -I 4bd 2c 3");
	
	// Reset output
	Output::method(origMethod);
}



/* Matrix3D SpaceGroup::conventionalTransformation(const ISO& iso, double tol, Vector3D* shift)
 *
 * Save transformation to conventional cell
 */

Matrix3D SpaceGroup::conventionalTransformation(const ISO& iso, double tol, Vector3D* shift)
{
	
	// Output
	Output::newline();
	Output::print("Determining conversion for structure to conventional form");
	Output::increase();
	
	// Get the space group of the structure
	SpaceGroup spaceGroup;
	spaceGroup.set(iso, tol);
	
	// Output
	Output::decrease();
	
	// Save shift
	if (shift)
		*shift = spaceGroup.originShift();
	
	// Return conversion matrix
	return spaceGroup.unitToConv();
}



/* void SpaceGroup::makeConventional(ISO& iso, double tol)
 *
 * Convert structure to conventional cell
 */

void SpaceGroup::makeConventional(ISO& iso, double tol)
{
	
	// Output
	Output::newline();
	Output::print("Converting structure to conventional form");
	Output::increase();
	
	// Get the space group of the structure
	SpaceGroup spaceGroup;
	spaceGroup.set(iso, tol);
	
	// Apply transformation to cell
	iso.transform(spaceGroup.unitToConv(), tol);
	iso.shift(spaceGroup.originShift()*-1);
	
	// Output
	Output::decrease(); 
}
