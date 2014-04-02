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



#include "pointGroup.h"
#include "output.h"



/* void PointGroup::clear()
 *
 * Clear data in PointGroup object
 */

void PointGroup::clear()
{
	_schoenflies.clear();
	_hermannMauguinShort.clear();
	_hermannMauguin.clear();
	_convSymmetry.clear();
}



/* PointGroup& PointGroup::operator= (const PointGroup& rhs)
 *
 * Assignment operator for PointGroup object
 */

PointGroup& PointGroup::operator= (const PointGroup& rhs)
{
	if (this != &rhs)
	{
		_centric = rhs._centric;
		_number = rhs._number;
		_order = rhs._order;
		_laueGroup = rhs._laueGroup;
		_centeringType = rhs._centeringType;
		_crystalSystem = rhs._crystalSystem;
		_schoenflies = rhs._schoenflies;
		_hermannMauguinShort = rhs._hermannMauguinShort;
		_hermannMauguin = rhs._hermannMauguin;
		_redToConv = rhs._redToConv;
		_unitToConv = rhs._unitToConv;
		_convSymmetry = rhs._convSymmetry;
	}
	return *this;
}



/* bool PointGroup::set(const Word& name, bool quitIfNotFound)
 *
 * Set point group by name
 */

bool PointGroup::set(const Word& name, bool quitIfNotFound)
{
	
	// Clear space
	clear();
	
	// Get point group number, error if not recognized
	int number = getNumberFromName(name);
	if (!number)
	{
		
		// Not quitting
		if (!quitIfNotFound)
			return false;

		// Quit
		Output::newline(ERROR);
		Output::print("Did not recognize point group");
		Output::quit();
	}
	
	// Set properties
	setByNumber(number);
	
	// Output
	Output::newline();
	Output::print("Setting point group to ");
	Output::print(name);
	Output::increase();

	// Print the result
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
	Output::print("Order: ");
	Output::print(_order);
	Output::newline();
	Output::print("Centric: ");
	if (_centric)
		Output::print("True");
	else
		Output::print("False");

	// Output
	Output::decrease();

	// Return that point group was found
	return true;
}



/* int PointGroup::getNumberFromName(const Word& name)
 *
 * Get the number of the point group from name
 */

int PointGroup::getNumberFromName(const Word& name)
{
	
	// Compare name to all point groups
	if (compareName(name, "C1",  "1",     "1"))         return 1;
	if (compareName(name, "Ci",  "-1",    "-1"))        return 2;
	if (compareName(name, "C2",  "2",     "2"))         return 3;
	if (compareName(name, "Cs",  "m",     "m"))         return 4;
	if (compareName(name, "C2h", "2/m",   "2/m"))       return 5;
	if (compareName(name, "D2",  "222",   "2 2 2"))     return 6;
	if (compareName(name, "C2v", "mm2",   "m m 2"))     return 7;
	if (compareName(name, "D2h", "mmm",   "2/m2/m2/m")) return 8;
	if (compareName(name, "C4",  "4",     "4"))         return 9;
	if (compareName(name, "S4",  "-4",    "-4"))        return 10;
	if (compareName(name, "C4h", "4/m",   "4/m"))       return 11;
	if (compareName(name, "D4",  "422",   "422"))       return 12;
	if (compareName(name, "C4v", "4mm",   "4mm"))       return 13;
	if (compareName(name, "D2d", "-4m2",  "-4m2"))      return 14;
	if (compareName(name, "D4h", "4/mmm", "4/m2/m2/m")) return 15;
	if (compareName(name, "C3",  "3",     "3"))         return 16;
	if (compareName(name, "C3i", "-3",    "-3"))        return 17;
	if (compareName(name, "D3",  "32",    "32"))        return 18;
	if (compareName(name, "C3v", "3m",    "3m"))        return 19;
	if (compareName(name, "D3d", "3m",    "-32/m"))     return 20;
	if (compareName(name, "C6",  "6",     "6"))         return 21;
	if (compareName(name, "C3h", "-6",    "-6"))        return 22;
	if (compareName(name, "C6h", "6/m",   "6/m"))       return 23;
	if (compareName(name, "D6",  "622",   "622"))       return 24;
	if (compareName(name, "C6v", "6mm",   "6mm"))       return 25;
	if (compareName(name, "D3h", "-6m2",  "-6m2"))      return 26;
	if (compareName(name, "D6h", "6/mmm", "6/m2/m2/m")) return 27;
	if (compareName(name, "T",   "23",    "23"))        return 28;
	if (compareName(name, "Th",  "m-3",   "2/m-3"))     return 29;
	if (compareName(name, "O",   "432",   "432"))       return 30;
	if (compareName(name, "Td",  "-43m",  "-43m"))      return 31;
	if (compareName(name, "Oh",  "m-3m",  "4/m-32/m"))  return 32;
	
	// Did not recognize point group name
	return 0;
}



/* bool PointGroup::compareName(const Word& name, const char* schoe, const char* hmShort, const char* hmFull)
 *
 * Return whether a name is that of point group
 */

bool PointGroup::compareName(const Word& name, const char* schoe, const char* hmShort, const char* hmFull)
{
	if ((name.equal(schoe, false, ' ')) || (name.equal(hmShort, false, ' ')) || (name.equal(hmFull, false, ' ')))
		return true;
	return false;
}



/* void PointGroup::set(const ISO& iso, double tol)
 *
 * Set point group by structure symmetry
 */

void PointGroup::set(const ISO& iso, double tol)
{
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Determing the point group of the current structure");
	Output::increase();
	
	// Get the transformation to reduced primitive cell
	Matrix3D unitToRed = iso.primitiveTransformation(2*tol);
	unitToRed *= Basis::reducedTransformation(unitToRed * iso.basis().vectors());
	
	// Output
	Output::newline();
	Output::print("Converting to reduced primitive cell");
	Output::increase();
	
	// Set the reduced primitive cell
	ISO reducedCell = iso;
	reducedCell.transform(unitToRed, 2*tol);
	
	// Output
	Output::decrease();
	
	// Get the symmetry of the current cell
	Symmetry redSymmetry;
	redSymmetry.set(reducedCell, tol, true);
	
	// Get point group from orders of symmetry operations
	setByNumber(getNumberFromOperations(redSymmetry));
	
	// Print the found point group
	Output::newline();
	Output::print("Structure is in point group ");
	Output::print(_hermannMauguin);
	
	// Get the transformation to the conventional cell
	setRedToConv(reducedCell, redSymmetry);
	_unitToConv = _redToConv * unitToRed;
	
	// Convert symmetry operations to conventional cell
	setConvSymmetry(redSymmetry);

	// Print type of centering vectors
	Output::newline();
	Output::print("Conventional cell has ");
	if (_centeringType == CT_PRIMITIVE)
		Output::print("primitive");
	else if (_centeringType == CT_ONE_FACE)
		Output::print("one face");
	else if (_centeringType == CT_BODY)
		Output::print("body");
	else if (_centeringType == CT_RHOMBOHEDRAL)
		Output::print("rhombohedral");
	else if (_centeringType == CT_ALL_FACE)
		Output::print("all face");
	else
		Output::print("unknown");
	Output::print(" centering");
	
	// Print the conversion to conventional cell
	int i, j;
	Output::newline();
	Output::print("Conversion from unit cell to standard cell");
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
	
	// Output
	Output::decrease();
}



/* int PointGroup::getNumberFromOperations(const Symmetry& symmetry)
 *
 * Get point group number from symmetry operations
 */

int PointGroup::getNumberFromOperations(const Symmetry& symmetry)
{
	
	// Output
	Output::newline();
	Output::print("Determining order of rotational elements of symmetry operations");
	Output::increase();
	
	// Variable to store the number of times that each rotation order occurs
	int num_1 = 0;
	int num_2 = 0;
	int num_3 = 0;
	int num_4 = 0;
	int num_6 = 0;
	int num_m1 = 0;
	int num_m2 = 0;
	int num_m3 = 0;
	int num_m4 = 0;
	int num_m6 = 0;
	
	// Count the number of times that each rotation order occurs
	int curOrder;
	for (int i = 0; i < symmetry.operations().length(); ++i)
	{
		curOrder = Symmetry::rotationOrder(symmetry.operations()[i].rotation());
		if (curOrder == 1)
			num_1++;
		else if (curOrder == 2)
			num_2++;
		else if (curOrder == 3)
			num_3++;
		else if (curOrder == 4)
			num_4++;
		else if (curOrder == 6)
			num_6++;
		else if (curOrder == -1)
			num_m1++;
		else if (curOrder == -2)
			num_m2++;
		else if (curOrder == -3)
			num_m3++;
		else if (curOrder == -4)
			num_m4++;
		else if (curOrder == -6)
			num_m6++;
	}
	
	// Print rotation types
	Output::newline();
	Output::print("Number of rotations of order 1: ");
	Output::print(num_1);
	Output::newline();
	Output::print("Number of rotations of order 2: ");
	Output::print(num_2);
	Output::newline();
	Output::print("Number of rotations of order 3: ");
	Output::print(num_3);
	Output::newline();
	Output::print("Number of rotations of order 4: ");
	Output::print(num_4);
	Output::newline();
	Output::print("Number of rotations of order 6: ");
	Output::print(num_6);
	Output::newline();
	Output::print("Number of rotations of order -1: ");
	Output::print(num_m1);
	Output::newline();
	Output::print("Number of rotations of order -2: ");
	Output::print(num_m2);
	Output::newline();
	Output::print("Number of rotations of order -3: ");
	Output::print(num_m3);
	Output::newline();
	Output::print("Number of rotations of order -4: ");
	Output::print(num_m4);
	Output::newline();
	Output::print("Number of rotations of order -6: ");
	Output::print(num_m6);
	
	// Determine the point group number
	int res = 0;
	int numSym = symmetry.operations().length();
	     if ((num_3 + num_m3 == 16) &&  (num_m1) && (numSym == 48))                                   res = 32;
	else if ((num_3 + num_m3 ==  8) && (!num_m1) && (numSym == 24) && (num_m4 == 6))                  res = 31;
	else if ((num_3 + num_m3 ==  8) && (!num_m1) && (numSym == 24) && (num_4  == 6))                  res = 30;
	else if ((num_3 + num_m3 == 16) &&  (num_m1) && (numSym == 24))                                   res = 29;
	else if ((num_3 + num_m3 ==  8) && (!num_m1) && (numSym == 12))                                   res = 28;
	else if ((num_6 + num_m6 ==  4) &&  (num_m1) && (numSym == 24))                                   res = 27;
	else if ((num_6 + num_m6 ==  2) && (!num_m1) && (numSym == 12) && (num_m6 == 2))                  res = 26;
	else if ((num_6 + num_m6 ==  2) && (!num_m1) && (numSym == 12) && (num_6  == 2) && (num_m2 == 6)) res = 25;
	else if ((num_6 + num_m6 ==  2) && (!num_m1) && (numSym == 12) && (num_6  == 2) && (num_2  == 7)) res = 24;
	else if ((num_6 + num_m6 ==  4) &&  (num_m1) && (numSym == 12))                                   res = 23;
	else if ((num_6 + num_m6 ==  2) && (!num_m1) && (numSym ==  6) && (num_m6 == 2))                  res = 22;
	else if ((num_6 + num_m6 ==  2) && (!num_m1) && (numSym ==  6) && (num_6  == 2))                  res = 21;
	else if ((num_3 + num_m3 ==  4) &&  (num_m1) && (numSym == 12))                                   res = 20;
	else if ((num_3 + num_m3 ==  2) && (!num_m1) && (numSym ==  6) && (num_m2 == 3))                  res = 19;
	else if ((num_3 + num_m3 ==  2) && (!num_m1) && (numSym ==  6) && (num_2  == 3))                  res = 18;
	else if ((num_3 + num_m3 ==  4) &&  (num_m1) && (numSym ==  6))                                   res = 17;
	else if ((num_3 + num_m3 ==  2) && (!num_m1) && (numSym ==  3))                                   res = 16;
	else if ((num_4 + num_m4 ==  4) &&  (num_m1) && (numSym == 16))                                   res = 15;
	else if ((num_4 + num_m4 ==  2) && (!num_m1) && (numSym ==  8) && (num_m4 == 2))                  res = 14;
	else if ((num_4 + num_m4 ==  2) && (!num_m1) && (numSym ==  8) && (num_4  == 2) && (num_m2 == 4)) res = 13;
	else if ((num_4 + num_m4 ==  2) && (!num_m1) && (numSym ==  8) && (num_4  == 2) && (num_2  == 5)) res = 12;
	else if ((num_4 + num_m4 ==  4) &&  (num_m1) && (numSym ==  8))                                   res = 11;
	else if ((num_4 + num_m4 ==  2) && (!num_m1) && (numSym ==  4) && (num_m4 == 2))                  res = 10;
	else if ((num_4 + num_m4 ==  2) && (!num_m1) && (numSym ==  4) && (num_4  == 2))                  res =  9;
	else if ((num_2 + num_m2 ==  6) &&  (num_m1))                                                     res =  8;
	else if ((num_2 + num_m2 ==  3) && (!num_m1) && (num_m2 ==  2))                                   res =  7;
	else if ((num_2 + num_m2 ==  3) && (!num_m1) && (num_2  ==  3))                                   res =  6;
	else if ((num_2 + num_m2 ==  2) &&  (num_m1))                                                     res =  5;
	else if ((num_2 + num_m2 ==  1) && (!num_m1) && (num_m2 ==  1))                                   res =  4;
	else if ((num_2 + num_m2 ==  1) && (!num_m1) && (num_2  ==  1))                                   res =  3;
	else if ( (num_m1) && (numSym == 2))                                                              res =  2;
	else if ((!num_m1) && (numSym == 1))                                                              res =  1;
	
	// Did not find point group
	if (!res)
	{
		Output::newline(ERROR);
		Output::print("Could not determine point group number");
		Output::quit();
	}
	
	// Output
	Output::decrease();
	
	// Return result
	return res;
}



/* void PointGroup::setByNumber(int number)
 *
 * Set point group properties by number
 */

void PointGroup::setByNumber(int number)
{
	     if (number ==  1) setProperties( 1, false,  1,  "C1",     "1",           "1",  1,    CS_TRICLINIC);
	else if (number ==  2) setProperties( 2,  true,  2,  "Ci",    "-1",          "-1",  1,    CS_TRICLINIC);
	else if (number ==  3) setProperties( 3, false,  2,  "C2",     "2",           "2",  2,   CS_MONOCLINIC);
	else if (number ==  4) setProperties( 4, false,  2,  "Cs",     "m",           "m",  2,   CS_MONOCLINIC);
	else if (number ==  5) setProperties( 5,  true,  4, "C2h",   "2/m",         "2/m",  2,   CS_MONOCLINIC);
	else if (number ==  6) setProperties( 6, false,  4,  "D2",   "222",       "2 2 2",  3, CS_ORTHORHOMBIC);
	else if (number ==  7) setProperties( 7, false,  4, "C2v",   "mm2",       "m m 2",  3, CS_ORTHORHOMBIC);
	else if (number ==  8) setProperties( 8,  true,  8, "D2h",   "mmm", "2/m 2/m 2/m",  3, CS_ORTHORHOMBIC);
	else if (number ==  9) setProperties( 9, false,  4,  "C4",     "4",           "4",  4,   CS_TETRAGONAL);
	else if (number == 10) setProperties(10, false,  4,  "S4",    "-4",          "-4",  4,   CS_TETRAGONAL);
	else if (number == 11) setProperties(11,  true,  8, "C4h",   "4/m",         "4/m",  4,   CS_TETRAGONAL);
	else if (number == 12) setProperties(12, false,  8,  "D4",   "422",       "4 2 2",  5,   CS_TETRAGONAL);
	else if (number == 13) setProperties(13, false,  8, "C4v",   "4mm",       "4 m m",  5,   CS_TETRAGONAL);
	else if (number == 14) setProperties(14, false,  8, "D2d",  "-4m2",      "-4 m 2",  5,   CS_TETRAGONAL);
	else if (number == 15) setProperties(15,  true, 16, "D4h", "4/mmm", "4/m 2/m 2/m",  5,   CS_TETRAGONAL);
	else if (number == 16) setProperties(16, false,  3,  "C3",     "3",           "3",  6,     CS_TRIGONAL);
	else if (number == 17) setProperties(17,  true,  6, "C3i",    "-3",          "-3",  6,     CS_TRIGONAL);
	else if (number == 18) setProperties(18, false,  6,  "D3",    "32",         "3 2",  7,     CS_TRIGONAL);
	else if (number == 19) setProperties(19, false,  6, "C3v",    "3m",         "3 m",  7,     CS_TRIGONAL);
	else if (number == 20) setProperties(20,  true, 12, "D3d",   "-3m",      "-3 2/m",  7,     CS_TRIGONAL);
	else if (number == 21) setProperties(21, false,  6,  "C6",     "6",           "6",  8,    CS_HEXAGONAL);
	else if (number == 22) setProperties(22, false,  6, "C3h",    "-6",          "-6",  8,    CS_HEXAGONAL);
	else if (number == 23) setProperties(23,  true, 12, "C6h",   "6/m",         "6/m",  8,    CS_HEXAGONAL);
	else if (number == 24) setProperties(24, false, 12,  "D6",   "622",       "6 2 2",  9,    CS_HEXAGONAL);
	else if (number == 25) setProperties(25, false, 12, "C6v",   "6mm",       "6 m m",  9,    CS_HEXAGONAL);
	else if (number == 26) setProperties(26, false, 12, "D3h",  "-6m2",      "-6 m 2",  9,    CS_HEXAGONAL);
	else if (number == 27) setProperties(27,  true, 24, "D6h", "6/mmm", "6/m 2/m 2/m",  9,    CS_HEXAGONAL);
	else if (number == 28) setProperties(28, false, 12,   "T",    "23",         "2 3", 10,        CS_CUBIC);
	else if (number == 29) setProperties(29,  true, 24,  "Th",   "m-3",      "2/m -3", 10,        CS_CUBIC);
	else if (number == 30) setProperties(30, false, 24,   "O",   "432",       "4 3 2", 11,        CS_CUBIC);
	else if (number == 31) setProperties(31, false, 24,  "Td",  "-43m",      "-4 3 m", 11,        CS_CUBIC);
	else if (number == 32) setProperties(32,  true, 48,  "Oh",  "m-3m",  "4/m -3 2/m", 11,        CS_CUBIC);
	else
	{
		Output::newline(ERROR);
		Output::print("Internal error: unknown point group number");
		Output::quit();
	}
}



/* void PointGroup::setProperties(int number, bool centric, int order, const char* schoe, const char* hmShort,
 *		const char* hm, int laueGroup, CrystalSystem crystalSystem)
 *
 * Set properties of point group
 */

void PointGroup::setProperties(int number, bool centric, int order, const char* schoe, const char* hmShort, \
	const char* hm, int laueGroup, CrystalSystem crystalSystem)
{
	_number = number;
	_centric = centric;
	_order = order;
	_laueGroup = laueGroup;
	_crystalSystem = crystalSystem;
	_schoenflies = schoe;
	_hermannMauguinShort = hmShort;
	_hermannMauguin = hm;
}



/* void PointGroup::setRedToConv(const ISO& reducedCell, const Symmetry& symmetry)
 *
 * Set transformation from reduced cell to conventional cell
 */

void PointGroup::setRedToConv(const ISO& reducedCell, const Symmetry& symmetry)
{
	
	// Common variables
	int i, j;
	int index;
	int curRotOrder;
	Vector3D axis;
	Matrix3D correction;
	
	// Laue group -1 is triclinic so keep current basis
	if (_laueGroup == 1)
		_redToConv = Matrix3D::identity();
	
	// Laue group 2/m, 4/m, -3, or 6/m
	else if ((_laueGroup == 2) || (_laueGroup == 4) || (_laueGroup == 6) || (_laueGroup == 8))
	{
		
		// Set the target rotation type
		int target;
		if (_laueGroup == 2)
			target = 2;
		else if (_laueGroup == 4)
			target = 4;
		else if ((_laueGroup == 6) || (_laueGroup == 8))
			target = 3;
		
		// Loop over rotations to look for target rotation
		index = -1;
		for (i = 0; i < symmetry.operations().length(); ++i)
		{
			
			// Get current rotation
			curRotOrder = Symmetry::rotationOrder(symmetry.operations()[i].rotation());
			if (!_centric)
				curRotOrder = Num<int>::abs(curRotOrder);
			
			// Found target if also a proper rotatation
			if (curRotOrder == target)
			{
				index = i;
				break;
			}
		}
		
		// Did not find rotation
		if (index < 0)
		{
			Output::newline(ERROR);
			Output::print("Did not find rotation type fixed by Laue group");
			Output::quit();
		}
		
		// Set the third basis vector to rotation axis
		_redToConv.setRow(2, Basis::rotationAxis(symmetry.operations()[index].rotation()));
		
		// Make sure rotation is proper
		Matrix3D properRotation = symmetry.operations()[index].rotation();
		if (properRotation.determinant() < 0)
			properRotation *= -1;
		
		// Calculate S matrix
		Matrix3D S = properRotation;
		Matrix3D curMat = properRotation;
		for (i = 1; i < target; ++i)
		{
			curMat *= properRotation;
			S += curMat;
		}
		
		// Reduce matrix to row echelon form
		S = S.rowEchelon(0, 0, true);
		
		// Possible combinations of vectors
		Vector3D fill[4];
		fill[0][0] =  1;
		fill[0][1] =  0;
		fill[0][2] =  0;
		fill[1][0] =  0;
		fill[1][1] =  1;
		fill[1][2] =  0;
		fill[2][0] =  1;
		fill[2][1] =  1;
		fill[2][2] =  0;
		fill[3][0] =  1;
		fill[3][1] = -1;
		fill[3][2] =  0;
		
		// Solve for each possible vector solution
		Vector3D combos[4];
		for (i = 0; i < 4; ++i)
			combos[i] = Basis::backSolve(S, &fill[i]);
		
		// Set the matrices if 2/m
		OList<Matrix3D > matrices;
		if (_laueGroup == 2)
		{
			int matIndex = 0;
			matrices.length(5);
			for (i = 0; i < 4; ++i)
			{
				for (j = i + 1; j < 4; ++j)
				{
					matrices[matIndex].setRow(0, combos[i]);
					matrices[matIndex].setRow(1, combos[j]);
					matrices[matIndex](2, 0) = _redToConv(2, 0);
					matrices[matIndex](2, 1) = _redToConv(2, 1);
					matrices[matIndex](2, 2) = _redToConv(2, 2);
					matIndex++;
				}
			}
		}
		
		// Set matrices if 4/m, -3, or 6/m
		else
		{
			matrices.length(4);
			for (i = 0; i < 4; ++i)
			{
				matrices[i].setRow(0, combos[i]);
				matrices[i].setRow(1, properRotation * combos[i]);
				matrices[i](2, 0) = _redToConv(2, 0);
				matrices[i](2, 1) = _redToConv(2, 1);
				matrices[i](2, 2) = _redToConv(2, 2);
			}
		}
		
		// Loop over possible solutions and check which gives smallest determinant
		int minIndex;
		double curDet;
		double minDet;
		for (i = 0; i < matrices.length(); i++)
		{
			curDet = Num<double>::abs(matrices[i].determinant());
			if ((curDet < minDet) || (!i))
			{
				minDet = curDet;
				minIndex = i;
			}
		}
		
		// Save conversion
		_redToConv(0, 0) = matrices[minIndex](0, 0);
		_redToConv(0, 1) = matrices[minIndex](0, 1);
		_redToConv(0, 2) = matrices[minIndex](0, 2);
		_redToConv(1, 0) = matrices[minIndex](1, 0);
		_redToConv(1, 1) = matrices[minIndex](1, 1);
		_redToConv(1, 2) = matrices[minIndex](1, 2);
	}
	
	// Laue group mmm or m-3
	else if ((_laueGroup == 3) || (_laueGroup == 10))
	{
		
		// Loop over rotations to look for target rotation
		index = 0;
		for (i = 0; i < symmetry.operations().length(); ++i)
		{
			
			// Get current rotation order
			curRotOrder = Symmetry::rotationOrder(symmetry.operations()[i].rotation());
			if (!_centric)
				curRotOrder = Num<int>::abs(curRotOrder);
			
			// Found target if also a proper rotatation
			if (curRotOrder == 2)
			{
				
				// Save rotation axis
				_redToConv.setRow(index, Basis::rotationAxis(symmetry.operations()[i].rotation()));
				
				// Increase row
				if (++index == 3)
					break;
			}
		}
		
		// Did not find all operations
		if (index != 3)
		{
			Output::newline(ERROR);
			Output::print("Did not find all 3 2-fold rotations in mmm or m-3 Laue group");
			Output::quit();
		}
	}
	
	// Laue groups 4/mmm, -3m, 6/mmm, and m-3m
	else
	{
		
		// Set target rotation orders
		int target[2] = {0, 0};
		if (_laueGroup == 5)
		{
			target[0] = 4;
			target[1] = 2;
		}
		else if ((_laueGroup == 7) || (_laueGroup == 9))
		{
			target[0] = 3;
			target[1] = 2;
		}
		else if (_laueGroup == 11)
		{
			target[0] = 4;
			target[1] = 4;
		}
		
		// Look for targets
		List<int> indexList;
		OList<Vector3D >::D2 vectors(2);
		for (i = 0; i < symmetry.operations().length(); ++i)
		{
			
			// Get current rotation
			curRotOrder = Symmetry::rotationOrder(symmetry.operations()[i].rotation());
			if (!_centric)
				curRotOrder = Num<int>::abs(curRotOrder);
			
			// Loop over targets
			for (j = 0; j < 2; ++j)
			{
				
				// Found target if also a proper rotatation
				if (curRotOrder == target[j])
				{

					// Get the rotation axis and save
					vectors[j] += Basis::rotationAxis(symmetry.operations()[i].rotation());

					// Save index if first operation
					if (!j)
						indexList += i;
				}
			}
		}
		
		// Did not find targets
		if ((!vectors[0].length()) || (!vectors[1].length()))
		{
			Output::newline(ERROR);
			Output::print("Did not find rotation types in conventional cell determination");
			Output::quit();
		}
		
		// Loop over combinations of vectors until two non-coaxial ones are found
		bool found = false;
		double angle;
		Vector3D tempVec;
		Matrix3D properRotation;
		for (i = 0; i < vectors[0].length(); ++i)
		{
			for (j = 0; j < vectors[1].length(); ++j)
			{
				
				// Found a good set of vectors
				angle = getAngle(vectors[0][i], vectors[1][j], reducedCell.basis());
				if ((angle > 75) && (angle < 105))
				{
					
					// Save vectors
					_redToConv.setRow(0, vectors[1][j]);
					_redToConv.setRow(2, vectors[0][i]);
					
					// Make sure that rotation is proper
					properRotation = symmetry.operations()[indexList[i]].rotation();
					if (properRotation.determinant() < 0)
						properRotation *= -1;
					
					// Set second vector
					tempVec.set(_redToConv(0, 0), _redToConv(0, 1), _redToConv(0, 2));
					tempVec *= properRotation;
					_redToConv(1, 0) = tempVec[0];
					_redToConv(1, 1) = tempVec[1];
					_redToConv(1, 2) = tempVec[2];
					
					// Break since found
					found = true;
					break;
				}
			}
			
			// Break if found
			if (found)
				break;
		}
	}
	
	// Swap first two rows if determinant is negative
	if (_redToConv.determinant() < 0)
		_redToConv.swapRows(0, 1);
	
	// Check if conversion gives zero volume
	if (_redToConv.determinant() < 1e-2)
	{
		Output::newline(ERROR);
		Output::print("Found a conventional cell with volume of zero");
		Output::quit();
	}
	
	// Get conversion matrix
	Matrix3D Q = _redToConv.transpose().inverse();
	
	// Get the lattice points
	LatticePoints latticePoints = ISO::getLatticePoints(_redToConv);
	
	// Could not identify centering
	if ((latticePoints.length() < 1) || (latticePoints.length() > 4))
	{
		Output::newline(ERROR);
		Output::print("Unknown number of centering vectors (");
		Output::print(latticePoints.length());
		Output::print(")");
		Output::quit();
	}
	
	// Convert centering vectors to conventional cell
	OList<Vector3D> centVecs (latticePoints.length());
	for (i = 0; i < latticePoints.length(); ++i)
	{
		centVecs[i] = Q * latticePoints[i];
		ISO::moveIntoCell(centVecs[i]);
	}
	
	// Correction for 2/m Laue group
	if (_laueGroup == 2)
	{
		correction.fill(0);
		correction(0, 1) = correction(1, 2) = correction(2, 0) = 1;
		_redToConv *= correction.inverse();
	}
	
	// Corrections for 4/m and 4/mmm
	if ((_laueGroup == 4) || (_laueGroup == 5))
	{
		
		// Taking C cell to P
		if ((centVecs.length() == 2) && (!findTranslation(centVecs, 0.5, 0.5, 0.5)))
		{
			correction(0, 0) = correction(0, 1) = correction(1, 0) = 1;
			correction(0, 2) = correction(1, 2) = correction(2, 0) = correction(2, 1) = 0;
			correction(1, 1) = correction(2, 2) = -1;
			_redToConv *= correction.inverse();
		}
		
		// Taking F cell to I
		else if (centVecs.length() == 4)
		{
			correction(0, 0) = correction(0, 1) = correction(1, 1) = correction(2, 2) = 1;
			correction(0, 2) = correction(1, 2) = correction(2, 0) = correction(2, 1) = 0;
			correction(1, 0) = -1;
			_redToConv *= correction.inverse();
		}
	}
	
	// Corrections for -3 and -3m
	if ((_laueGroup == 6) || (_laueGroup == 7))
	{
		
		// Take reverse to obverse
		if (findTranslation(centVecs, 0.5, 0.5, 0.5))
		{
			correction(2, 2) = 1;
			correction(0, 0) = correction(1, 1) = -1;
			correction(0, 1) = correction(0, 2) = correction(1, 0) = 0;
			correction(1, 2) = correction(2, 0) = correction(2, 1) = 0;
			_redToConv *= correction.inverse();
		}
	}
	
	// Correction for -3m and 6/mmm
	if ((_laueGroup == 7) || (_laueGroup == 9))
	{
		
		// Take H cell to P
		if (findTranslation(centVecs, 1.0/3.0, 2.0/3.0, 0))
		{
			correction(0, 2) = correction(1, 2) = correction(2, 0) = correction(2, 1) = 0;
			correction(0, 1) = correction(2, 2) = 1;
			correction(0, 0) = correction(1, 0) = -1;
			correction(1, 1) = -2;
			_redToConv *= correction.inverse();
		}
	}
	
	// Swap first two rows if determinant is negative
	if (_redToConv.determinant() < 0)
		_redToConv.swapRows(0, 1);
	
	// Convert symmetry rotations to conventional cell setting
	Matrix3D P = Q.inverse();
	OList<Matrix3D> convRotations(symmetry.operations().length());
	for (i = 0; i < symmetry.operations().length(); ++i)
	{
		convRotations[i]  = P;
		convRotations[i] *= symmetry.operations()[i].rotation();
		convRotations[i] *= Q;
	}
	
	// Make sure that unique axis is correct in monoclinic cell
	if (_laueGroup == 2)
	{
		int uniqueAxis = Symmetry::getConventionalCellUniqueSymmetryAxis(convRotations, LS_MONOCLINIC);
		if (uniqueAxis == 0)
			_redToConv *= Matrix3D(0, 0, 1, 1, 0, 0, 0, 1, 0);
		else if (uniqueAxis == 2)
			_redToConv *= Matrix3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
	}
	
	// Make sure that unique axis is correct in tetragonal cell
	else if ((_laueGroup == 4) || (_laueGroup == 5))
	{
		int uniqueAxis = Symmetry::getUniqueAxis(_redToConv * reducedCell.basis().vectors(), LS_TETRAGONAL);
		if (uniqueAxis == 0)
			uniqueAxis = Symmetry::getConventionalCellUniqueSymmetryAxis(convRotations, LS_TETRAGONAL);
		if (uniqueAxis == 0)
			_redToConv *= Matrix3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
		else if (uniqueAxis == 1)
			_redToConv *= Matrix3D(0, 0, 1, 1, 0, 0, 0, 1, 0);
	}
	
	// Make sure that unique axis is correct in hexagonal cell
	else if ((_laueGroup >= 6) && (_laueGroup <= 9))
	{
		int uniqueAxis = Symmetry::getUniqueAxis(_redToConv * reducedCell.basis().vectors(), LS_HEXAGONAL);
		if (uniqueAxis == 0)
			uniqueAxis = Symmetry::getConventionalCellUniqueSymmetryAxis(convRotations, LS_HEXAGONAL);
		if (uniqueAxis == 0)
			_redToConv *= Matrix3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
		else if (uniqueAxis == 1)
			_redToConv *= Matrix3D(0, 0, 1, 1, 0, 0, 0, 1, 0);
	}
	
	// Get updated conversion matrix and lattice points
	Q = _redToConv.transpose().inverse();
	latticePoints = ISO::getLatticePoints(_redToConv);
	
	// Could not identify centering
	if ((latticePoints.length() < 1) || (latticePoints.length() > 4))
	{
		Output::newline(ERROR);
		Output::print("Unknown number of centering vectors (");
		Output::print(latticePoints.length());
		Output::print(")");
		Output::quit();
	}
	
	// Convert centering vectors to conventional cell
	centVecs.length(latticePoints.length());
	for (i = 0; i < latticePoints.length(); ++i)
	{
		centVecs[i] = Q * latticePoints[i];
		ISO::moveIntoCell(centVecs[i]);
	}
	
	// Save the centering type
	_centeringType = CT_UNKNOWN;
	if (findTranslation(centVecs, 0, 0, 0))
	{
		if (centVecs.length() == 1)
			_centeringType = CT_PRIMITIVE;
		else if (centVecs.length() == 2)
		{
			if (findTranslation(centVecs, 0.5, 0.5, 0.5))
				_centeringType = CT_BODY;
			else
				_centeringType = CT_ONE_FACE;
		}
		else if (centVecs.length() == 3)
			_centeringType = CT_RHOMBOHEDRAL;
		else if (centVecs.length() == 4)
		{
			if ((findTranslation(centVecs, 0.5, 0.5, 0.0)) && (findTranslation(centVecs, 0.5, 0.0, 0.5)) && \
				(findTranslation(centVecs, 0.0, 0.5, 0.5)))
				_centeringType = CT_ALL_FACE;
		}
	}
	
	// Could not set centering vectors
	if (_centeringType == CT_UNKNOWN)
	{
		Output::newline(ERROR);
		Output::print("Unknown centering type in point group determination");
		Output::quit();
	}
}



/* double PointGroup::getAngle(const Vector3D& vec1, const Vector3D& vec2, const Basis& basis)
 *
 * Return the angle between two vectors
 */

double PointGroup::getAngle(const Vector3D& vec1, const Vector3D& vec2, const Basis& basis)
{
	Vector3D cartVec1 = basis.getCartesian(vec1);
	Vector3D cartVec2 = basis.getCartesian(vec2);
	return Num<double>::toDegrees(cartVec1.angle(cartVec2));
}



/* bool PointGroup::findTranslation(const OList<Vector3D >& centVecs, double val1, double val2, double val3)
 *
 * Return whether a translation is in list
 */

bool PointGroup::findTranslation(const OList<Vector3D >& centVecs, double val1, double val2, double val3)
{
	for (int i = 0; i < centVecs.length(); ++i)
	{
		if ((Num<double>::modEq(centVecs[i][0], val1, 1e-4)) && (Num<double>::modEq(centVecs[i][1], val2, 1e-4)) && \
			(Num<double>::modEq(centVecs[i][2], val3, 1e-4)))
			return true;
	}
	return false;
}



/* void PointGroup::setConvSymmetry(const Symmetry& symmetry)
 *
 * Set the conventional cell symmetry
 */

void PointGroup::setConvSymmetry(const Symmetry& symmetry)
{
	
	// Set conversion matrices
	Matrix3D P = _redToConv.transpose();
	Matrix3D Q = P.inverse();
	
	// Convert symmetry operations
	_convSymmetry.length(symmetry.operations().length());
	for (int i = 0; i < symmetry.operations().length(); ++i)
	{
		_convSymmetry[i].setRotation(Q * symmetry.operations()[i].rotation() * P);
		_convSymmetry[i].addTranslation(Q * symmetry.operations()[i].translations()[0]);
	}
}



/* void PointGroup::print() const
 *
 * Print point group information
 */

void PointGroup::print() const
{
	
	// Save current stream
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);

	// Print title line
	Output::newline();
	Output::print("Point group information");

	// Make output object
	Output message;
	message.addLines(5);

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

	// Order
	message.addLine();
	message.add("    Order:");
	message.add(_order);

	// Centric
	message.addLine();
	message.add("    Centric:");
	if (_centric)
		message.add("True");
	else
		message.add("False");

	// Print information
	List<PrintAlign> align (2);
	align[0] = RIGHT;
	align[1] = LEFT;
	Output::newline();
	Output::print(message, align);
	
	// Reset method
	Output::method(origMethod);
}



/* void PointGroup::printAll()
 *
 * Print table of all point groups
 */

void PointGroup::printAll()
{
	
	// Save current stream
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print groups
	Output::newline(); Output::print("HM Short       HM Full   Schoenflies   Order   Centric");
	Output::newline(); Output::print("--------       -------   -----------   -----   -------");
	Output::newline(); Output::print("       1             1           C1        1     False");
	Output::newline(); Output::print("      -1            -1           Ci        2      True");
	Output::newline(); Output::print("       2             2           C2        2     False");
	Output::newline(); Output::print("       m             m           Cs        2     False");
	Output::newline(); Output::print("     2/m           2/m           C2h       4      True");
	Output::newline(); Output::print("     222         2 2 2           D2        4     False");
	Output::newline(); Output::print("     mm2         m m 2           C2v       4     False");
	Output::newline(); Output::print("     mmm   2/m 2/m 2/m           D2h       8      True");
	Output::newline(); Output::print("       4             4           C4        4     False");
	Output::newline(); Output::print("      -4            -4           S4        4     False");
	Output::newline(); Output::print("     4/m           4/m           C4h       8      True");
	Output::newline(); Output::print("     422         4 2 2           D4        8     False");
	Output::newline(); Output::print("     4mm         4 m m           C4v       8     False");
	Output::newline(); Output::print("    -4m2        -4 m 2           D2d       8     False");
	Output::newline(); Output::print("   4/mmm   4/m 2/m 2/m           D4h      16      True");
	Output::newline(); Output::print("       3             3           C3        3     False");
	Output::newline(); Output::print("      -3            -3           C3i       6      True");
	Output::newline(); Output::print("      32           3 2           D3        6     False");
	Output::newline(); Output::print("      3m           3 m           C3v       6     False");
	Output::newline(); Output::print("     -3m        -3 2/m           D3d      12      True");
	Output::newline(); Output::print("       6             6           C6        6     False");
	Output::newline(); Output::print("      -6            -6           C3h       6     False");
	Output::newline(); Output::print("     6/m           6/m           C6h      12      True");
	Output::newline(); Output::print("     622         6 2 2           D6       12     False");
	Output::newline(); Output::print("     6mm         6 m m           C6v      12     False");
	Output::newline(); Output::print("    -6m2        -6 m 2           D3h      12     False");
	Output::newline(); Output::print("   6/mmm   6/m 2/m 2/m           D6h      24      True");
	Output::newline(); Output::print("      23           2 3           T        12     False");
	Output::newline(); Output::print("     m-3        2/m -3           Th       24      True");
	Output::newline(); Output::print("     432         4 3 2           O        24     False");
	Output::newline(); Output::print("    -43m        -4 3 m           Td       24     False");
	Output::newline(); Output::print("    m-3m    4/m -3 2/m           Oh       48      True");
	
	// Reset method
	Output::method(origMethod);
}
