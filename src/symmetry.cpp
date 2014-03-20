/* symmetry.h -- Symmetry of a crystal
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "multi.h"
#include "symmetry.h"
#include "language.h"
#include "output.h"
#include "constants.h"
#include <cmath>



/* Words JonesFaithful::toString(const Matrix3D& rotation, const Vector3D* translation)
 *
 * Convert rotation and translation to Jones-Faithful notation
 */

Words JonesFaithful::toString(const Matrix3D& rotation, const Vector3D* translation)
{
	
	// Variable to store result
	Words res(3);
	
	// Loop over directions
	int i, j;
	double abs;
	double tol = 1e-6;
	for (i = 0; i < 3; ++i)
	{
		
		// Loop over x, y, and z
		for (j = 0; j < 3; ++j)
		{
			
			// Current entry is non-zero
			abs = Num<double>::abs(rotation(i, j));
			if (abs > tol)
			{
				
				// Add + for any positive number
				if (rotation(i, j) > 0)
					res[i] += '+';
				
				// If value is not 1 or -1, turn into fraction
				if (Num<double>::abs(abs - 1) > tol)
					res[i] += Language::numberToFraction(rotation(i, j), tol);
				
				// Found -1
				else if (rotation(i, j) < 0)
					res[i] += '-';
				
				// Add x, y, or z
				if (j == 0)
					res[i] += 'x';
				else if (j == 1)
					res[i] += 'y';
				else
					res[i] += 'z';
			}
		}
		
		// Save translation
		if (translation)
		{
			abs = Num<double>::abs((*translation)[i]);
			if ((abs > tol) && (Num<double>::abs(abs - 1) > tol))
			{
				if (((*translation)[i] > 0) && (res[i].length()))
					res[i] += '+';
				res[i] += Language::numberToFraction((*translation)[i], tol);
			}
		}
		
		// Add a zero if nothing has been saved
		if (!res[i].length())
			res[i] += '0';
	}
	
	// Return result
	return res;
}



/* void JonesFaithful::fromString(Matrix3D& rotation, Vector3D& translation, const Words& input)
 *
 * Convert Jones-Faithful notation to rotation and translation
 */

void JonesFaithful::fromString(Matrix3D& rotation, Vector3D& translation, const Words& input)
{
	
	// Clear data in result
	rotation.fill(0);
	translation.fill(0);
	
	// Check whether to use commas as deliminators
	int i, j;
	int commas = 0;
	for (i = 0; i < input.length(); ++i)
	{
		for (j = 0; j < input[i].length(); ++j)
		{
			if (input[i][j] == ',')
				++commas;
		}
	}
	
	// Found 2 commas so use them as deliminators
	bool useCommas = false;
	if (commas == 2)
		useCommas = true;
	
	// Error if there are not three words
	else if (input.length() != 3)
	{
		Output::newline(ERROR);
		Output::print("Could not interpret Jones-Faithful notation");
		Output::quit();
	}
	
	// Loop over words in input
	int cur = 0;
	double mult = 1;
	Word frac;
	for (i = 0; i < input.length(); ++i)
	{
		
		// Loop over characters in word
		for (j = 0; j < input[i].length(); ++j)
		{
			
			// Found a + or -
			if (input[i][j] == '+')
			{
				if (frac.length())
					translation[cur] = mult * Language::fractionToNumber(frac);
				frac.clear();
				mult = 1;
			}
			else if (input[i][j] == '-')
			{
				if (frac.length())
					translation[cur] = mult * Language::fractionToNumber(frac);
				frac.clear();
				mult = -1;
			}
			
			// Found x, y, or z
			else if ((input[i][j] == 'x') || (input[i][j] == 'X'))
			{
				rotation(cur, 0) = mult;
				if (frac.length() > 0)
					rotation(cur, 0) *= Language::fractionToNumber(frac);
				frac.clear();
			}
			else if ((input[i][j] == 'y') || (input[i][j] == 'Y'))
			{
				rotation(cur, 1) = mult;
				if (frac.length() > 0)
					rotation(cur, 1) *= Language::fractionToNumber(frac);
				frac.clear();
			}
			else if ((input[i][j] == 'z') || (input[i][j] == 'Z'))
			{
				rotation(cur, 2) = mult;
				if (frac.length() > 0)
					rotation(cur, 2) *= Language::fractionToNumber(frac);
				frac.clear();
			}
			
			// Found a number
			else if ((input[i][j] == '0') || (input[i][j] == '1') || (input[i][j] == '2') || (input[i][j] == '3') || \
					 (input[i][j] == '4') || (input[i][j] == '5') || (input[i][j] == '6') || (input[i][j] == '7') || \
					 (input[i][j] == '8') || (input[i][j] == '9') || (input[i][j] == '.') || (input[i][j] == '/'))
				frac += input[i][j];
			
			// Found a comma
			else if ((input[i][j] == ',') && (useCommas))
			{
				if (frac.length())
					translation[cur] = mult * Language::fractionToNumber(frac);
				frac.clear();
				mult = 1;
				cur++;
			}
		}
		
		// Using words as deliminators
		if ((!useCommas) || ((useCommas) && (i == input.length() - 1)))
		{
			if (frac.length())
				translation[cur] = mult * Language::fractionToNumber(frac);
			frac.clear();
			mult = 1;
			cur++;
		}
	}
}



/* SymmetryOperation& SymmetryOperation::operator= (const SymmetryOperation& rhs)
 *
 * Assignment operator for SymmetryOperation object
 */

SymmetryOperation& SymmetryOperation::operator= (const SymmetryOperation& rhs)
{
	if (this != &rhs)
	{
		_rotation = rhs._rotation;
		_translations = rhs._translations;
	}
	return *this;
}



/* void SymmetryOperation::addTranslation(const Vector3D& translation)
 *
 * Add translation to list under a single rotation
 */

void SymmetryOperation::addTranslation(const Vector3D& translation)
{
	
	// Add the translation
	int i;
	double tol = 1e-6;
	_translations += translation;
	for (i = 0; i < 3; ++i)
		_translations.last()[i] = Num<double>::mod(_translations.last()[i], 1-tol);
	
	// Loop over translations and sort
	int j;
	for (i = _translations.length() - 2, j = _translations.length() - 1; i >= 0; --i, --j)
	{
		
		// Lower than first value
		if (_translations[j][0] < _translations[i][0] - tol)
			_translations.swap(j, i);
		
		// First values are equal
		else if (_translations[j][0] < _translations[i][0] + tol)
		{
			
			// Lower than second value
			if (_translations[j][1] < _translations[i][1] - tol)
				_translations.swap(j, i);
			
			// Second values are equal
			else if (_translations[j][1] < _translations[i][1] + tol)
			{
				
				// Lower than third value
				if (_translations[j][2] < _translations[i][2] - tol)
					_translations.swap(j, i);
				
				// Third values are equal
				else if (_translations[j][2] < _translations[i][2] + tol)
				{
					_translations.remove(j);
					return;
				}
			}
		}
	}
}



/* void Symmetry::clear()
 *
 * Clear data in symmetry object
 */

void Symmetry::clear()
{
	_operations.clear();
	_orbits.clear();
	_orbitNumbers.clear();
}



/* Symmetry& Symmetry::operator= (const Symmetry& rhs)
 *
 * Copy Symmetry object
 */

Symmetry& Symmetry::operator= (const Symmetry& rhs)
{
	if (this != &rhs)
	{
		_operations = rhs._operations;
		_orbits = rhs._orbits;
		_orbitNumbers = rhs._orbitNumbers;
		_metricMatrixConstraint = rhs._metricMatrixConstraint;
	}
	return *this;
}



/* void Symmetry::setToP1(const ISO& iso)
 *
 * Set symmetry to that of P1 space group
 */

void Symmetry::setToP1(const ISO& iso)
{
	
	// Clear space
	clear();
	
	// Save identity and zero translation
	Matrix3D identity;
	identity.makeIdentity();
	Vector3D zeroTrans = 0.0;
	
	// Add identity operation
	_operations += SymmetryOperation(identity, zeroTrans);
	
	// Add orbits for each atom
	int i, j;
	_orbits.length(iso.numAtoms());
	_orbitNumbers.length(iso.numAtoms());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			_orbits[iso.atoms()[i][j].atomNumber()].set(&iso.atoms()[i][j], identity, zeroTrans);
			_orbitNumbers[iso.atoms()[i][j].atomNumber()] = iso.atoms()[i][j].atomNumber();
		}
			
	}
	
	// Set metric matrix constraint matrix
	_metricMatrixConstraint.fill(0);
}



/* void Symmetry::set(const ISO& iso, double tol, bool isReducedPrim)
 *
 * Set the symmetry of a structure
 */

void Symmetry::set(const ISO& iso, double tol, bool isReducedPrim)
{
	
	// Clear current data
	clear();
	
	// Output
	Output::newline();
	Output::print("Determining the cell symmetry using a tolerance of ");
	Output::printSci(tol);
	Output::print(" Ang");
	Output::increase();
	
	// Get conversion to reduced cell
	Transformation unitToRed = unitToReduced(iso, 2*tol, isReducedPrim);
	
	// Get reduced primitive cell
	ISO reduced(iso);
	reduced.transform(unitToRed, 2*tol, false);
	
	// Output
	Output::newline();
	Output::print("Searching for operations");
	Output::increase();
	
	// Get the possible rotations
	Linked<Rotation> candidates;
	Basis::getPossibleRotations(candidates, reduced.basis().vectors(), tol);
	
	// Get the symmetry operations
	Linked<Rotation> redRotations;
	Linked<Translation> redTranslations;
	testOperations(reduced, tol, candidates, redRotations, redTranslations);
	
	// Output
	Output::decrease();
	Output::newline();
	Output::print("Converting primitive cell operations to unit cell");
	Output::increase();
	
	// Save unit cell operations
	setUnitOperations(redRotations, redTranslations, unitToRed);
	setMetricMatrixConstraint();
	
	// Output
	Output::decrease();
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Grouping atoms by equivalence");
	Output::increase();
	
	// Set the Wyckoff positions
	setOrbits(iso, tol);
	
	// Output
	Output::decrease();
}



/* Symmetry::Transformation Symmetry::unitToReduced(const ISO& iso, double tol, bool isReducedPrim)
 *
 * Get transformations between cell types
 */

Symmetry::Transformation Symmetry::unitToReduced(const ISO& iso, double tol, bool isReducedPrim)
{
	
	// Already in reduced primitive cell
	if (isReducedPrim)
		return Matrix3D::identity();
	
	// Output
	Output::newline();
	Output::print("Searching for transformations between cells");
	Output::increase();
	
	// Get transformation to primitive cell
	Transformation unitToPrim = iso.primitiveTransformation(tol, false);
	
	// Print transformation
	int i, j;
	Output::newline();
	Output::print("Transformation from unit cell to primitive cell");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(unitToPrim(i, j));
			Output::print(" ");
		}
	}
	
	// Get transformation to reduced primitive cell
	Transformation primToRed = Basis::reducedTransformation(unitToPrim * iso.basis().vectors());
	
	// Print transformation
	Output::newline();
	Output::print("Transformation from primitive cell to reduced primitive cell");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(primToRed(i, j));
			Output::print(" ");
		}
	}
	
	// Output
	Output::decrease();
	
	// Return result
	return primToRed * unitToPrim;
}



/* void Symmetry::testOperations(const ISO& iso, double tol, Linked<Rotation>& candidates, \
 *		Linked<Rotation>& redRotations, Linked<Translation>& redTranslations)
 *
 * Get the symmetry operations that exist for the current structure
 */

void Symmetry::testOperations(const ISO& iso, double tol, Linked<Rotation>& candidates, \
	Linked<Rotation>& redRotations, Linked<Translation>& redTranslations)
{
	
	// Get the distance of all atoms from the origin
	int i, j;
	Vector3D origin(0.0);
	List<double>::D2 distances(iso.atoms().length());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		distances[i].length(iso.atoms()[i].length());
		for (j = 0; j < iso.atoms()[i].length(); ++j)
			distances[i][j] = iso.basis().distance(iso.atoms()[i][j].fractional(), FRACTIONAL, origin, FRACTIONAL);
	}
	
	// Variables to store rotated and translated positions
	Atoms rotAtoms(iso.atoms());
	Atoms transAtoms(iso.atoms());
	
	// Save that identity is an allowed operation
	Rotation identity;
	identity.makeIdentity();
	redRotations.add(identity);
	redTranslations.add(Translation(0.0));
	
	// Check if inversion exists
	Translation curTranslation(0.0);
	Rotation inversion = identity * -1;
	if (checkOperation(inversion, iso, tol, curTranslation, rotAtoms, transAtoms, &distances))
	{
		redRotations.add(inversion);
		redTranslations.add(curTranslation);
	}
	
	// Loop until there are no operations left
	Linked<Rotation> notAllowed;
	while (candidates.length())
	{
		
		// Current operation is in space group
		if (checkOperation(*(candidates.begin()), iso, tol, curTranslation, rotAtoms, transAtoms, &distances))
			addOperation(candidates, *(candidates.begin()), curTranslation, redRotations, redTranslations, notAllowed);
		
		// Current operation is not in space group
		else
			removeOperation(candidates, *(candidates.begin()), redRotations, notAllowed);
	}
}



/* bool Symmetry::checkOperation(const Rotation& rotation, const ISO& iso, double tol, Translation& translation,
 *		Atoms& rotAtoms, Atoms& transAtoms, List<double>::D2* distances)
 *
 * Check whether a symmetry operation exists in group
 */

bool Symmetry::checkOperation(const Rotation& rotation, const ISO& iso, double tol, Translation& translation, \
	Atoms& rotAtoms, Atoms& transAtoms, List<double>::D2* distances)
{
	
	// Figure out which element occurs the least
	int i;
	int minElem = 0;
	for (i = 1; i < iso.atoms().length(); ++i)
	{
		if (iso.atoms()[i].length() < iso.atoms()[minElem].length())
			minElem = i;
	}
	
	// Transform all atoms under symmetry operation
	int j;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
			rotAtoms[i][j].fractional(rotation * iso.atoms()[i][j].fractional());
	}
	
	// Number of loops to make
	int numLoops = (int) Num<double>::ceil((double) iso.atoms()[minElem].length() / Multi::worldSize());
	
	// Loop over possible translations
	int k;
	bool areEqual;
	bool areEqualRec;
	int index;
	Translation translationRec;
	for (i = 0; i < numLoops; ++i)
	{
		
		// Figure out current atom to use in translation
		index = i * Multi::worldSize() + Multi::rank();
		
		// Check if sites map on current processor
		areEqual = false;
		if (index < rotAtoms[minElem].length())
		{
			
			// Get current translation
			translation = iso.atoms()[minElem][0].fractional() - rotAtoms[minElem][index].fractional();
			ISO::moveIntoCell(translation);
			
			// Translate all atoms by current vector
			for (j = 0; j < transAtoms.length(); ++j)
			{
				for (k = 0; k < transAtoms[j].length(); ++k)
					transAtoms[j][k].fractional(rotAtoms[j][k].fractional() + translation);
			}

			// Check if sites map
			areEqual = ISO::areSitesEqual(iso.basis(), iso.atoms(), transAtoms, tol, &translation, distances);
		}
		
		// Loop over processors
		for (j = 0; j < Multi::worldSize(); ++j)
		{
			
			// Broadcast result
			areEqualRec = areEqual;
			Multi::broadcast(areEqualRec, j);
			
			// Site do not map
			if (!areEqualRec)
				continue;
			
			// Broadcast vector
			translationRec = translation;
			Multi::broadcast(translationRec, j);
			
			// Save translation and return that generator is a valid operation
			translation = translationRec;
			return true;
		}
	}
	
	// Return that operation is not valid
	return false;
}



/* void Symmetry::addOperation(Linked<Rotation>& candidates, const Rotation& rotation, const Translation& translation,
 *		Linked<Rotation>& redRotations, Linked<Translation>& redTranslations, Linked<Rotation>& notAllowed)
 *
 * Add operation to allowed symmetries
 */

void Symmetry::addOperation(Linked<Rotation>& candidates, const Rotation& rotation, const Translation& translation, \
	Linked<Rotation>& redRotations, Linked<Translation>& redTranslations, Linked<Rotation>& notAllowed)
{
	
	// Check for operation in saved list
	Linked<Rotation>::iterator itRot;
	for (itRot = redRotations.begin(); itRot != redRotations.end(); ++itRot)
	{
		if (*itRot == rotation)
			return;
	}
	
	// Check for operation in candidates list
	for (itRot = candidates.begin(); itRot != candidates.end(); ++itRot)
	{
		if (*itRot == rotation)
			break;
	}
	
	// Did not find operation
	if (itRot == candidates.end())
		return;
	
	// Remove operation
	Rotation rotToAdd = rotation;
	candidates.remove(itRot);
	
	// Save operation
	redRotations.add(rotToAdd);
	redTranslations.add(translation);
	
	// Generate operations that cannot be in group based on operations that already are not allowed
	for (itRot = notAllowed.begin(); itRot != notAllowed.end(); ++itRot)
	{
		removeOperation(candidates, (*itRot) * rotToAdd, redRotations, notAllowed);
		removeOperation(candidates, rotToAdd * (*itRot), redRotations, notAllowed);
	}
	
	// Generate new operations
	Rotation newRotation;
	Translation newTranslation;
	Linked<Translation>::iterator itTrans = redTranslations.begin();
	for (itRot = redRotations.begin(); itRot != redRotations.end(); ++itRot, ++itTrans)
	{
		
		// Try forward multiplication
		newRotation = (*itRot) * rotToAdd;
		newTranslation = (*itRot) * translation + *itTrans;
		ISO::moveIntoCell(newTranslation);
		addOperation(candidates, newRotation, newTranslation, redRotations, redTranslations, notAllowed);
		
		// Try reverse multiplication
		newRotation = rotToAdd * (*itRot);
		newTranslation = rotToAdd * (*itTrans) + translation;
		ISO::moveIntoCell(newTranslation);
		addOperation(candidates, newRotation, newTranslation, redRotations, redTranslations, notAllowed);
	}
}



/* void Symmetry::removeOperation(Linked<Rotation>& candidates, const Rotation& rotation,
 *		Linked<Rotation>& redRotations, Linked<Rotation>& notAllowed)
 *
 * Remove operation from allowed symmetries
 */

void Symmetry::removeOperation(Linked<Rotation>& candidates, const Rotation& rotation, \
	Linked<Rotation>& redRotations, Linked<Rotation>& notAllowed)
{
	
	// Check if operation is in candidates list
	Linked<Rotation>::iterator it;
	for (it = candidates.begin(); it != candidates.end(); ++it)
	{
		if (*it == rotation)
			break;
	}
	
	// Operation has already been taken care of
	if (it == candidates.end())
		return;
	
	// Remove operation
	Rotation rotToRemove = rotation;
	candidates.remove(it);
	notAllowed.add(rotToRemove);
	
	// Remove inverse operation
	removeOperation(candidates, rotToRemove.inverse(), redRotations, notAllowed);
	
	// Loop over operations to generate those that cannot exist based on current operation
	for (it = redRotations.begin(); it != redRotations.end(); ++it)
	{
		removeOperation(candidates, *it * rotToRemove, redRotations, notAllowed);
		removeOperation(candidates, rotToRemove * *it, redRotations, notAllowed);
	}
}



/* void Symmetry::setOperations(const Linked<Rotation>& redRotations, Linked<Translation>& redTranslations,
 *		const Transformation& unitToReduced)
 *
 * Set symmetry operations in unit cell
 */

void Symmetry::setUnitOperations(const Linked<Rotation>& redRotations, Linked<Translation>& redTranslations, \
	const Transformation& unitToReduced)
{
	
	// Get conversion matrices
	Transformation Q = unitToReduced.transpose();
	Transformation P = Q.inverse();
	
	// Get the centering vectors
	int i;
	LatticePoints centeringVectors = ISO::getLatticePoints(unitToReduced.inverse());
	for (i = 0; i < centeringVectors.length(); ++i)
	{
		centeringVectors[i] *= Q;
		ISO::moveIntoCell(centeringVectors[i]);
	}
	
	// Clear current operations
	_operations.clear();
	
	// Loop over operations
	int j;
	Vector3D curTrans;
	Matrix3D curMatrix;
	Linked<Rotation>::iterator itRot = redRotations.begin();
	Linked<Translation>::iterator itTrans = redTranslations.begin();
	for (i = 0; itRot != redRotations.end(); ++itRot, ++itTrans, ++i)
	{
		
		// Save current rotation
		curMatrix = Q * *itRot * P;
		
		// Check that all components are integer
		if (curMatrix.isInteger(1e-2) == false)
			continue;
		
		// Save rotation
		_operations.add();
		_operations.last().setRotation(curMatrix);
		
		// Add translations
		curTrans = Q * *itTrans;
		for (j = 0; j < centeringVectors.length(); ++j)
			_operations.last().addTranslation(curTrans + centeringVectors[j]);
	}
	
	// Print warning if the number of rotations is less than reduced cell
	if (redRotations.length() > _operations.length())
	{
		Output::newline(WARNING);
		Output::print("Ignoring ");
		Output::print(redRotations.length() - _operations.length());
		Output::print(" non-integer rotation");
		if (redRotations.length() - _operations.length() != 1)
			Output::print("s");
		Output::print(" in the current unit cell");
	}
	
	// Print rotations
	Output::newline();
	Output::print("Found ");
	Output::print(_operations.length());
	Output::print(" unique symmetry operation");
	if (_operations.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _operations.length(); ++i)
	{
		Output::newline();
		Output::print(_operations[i].getString(), true, false);
	}
	Output::decrease();
	
	// Print centering vectors
	Output::newline();
	Output::print("Found ");
	Output::print(centeringVectors.length());
	Output::print(" centering vector");
	if (centeringVectors.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < centeringVectors.length(); ++i)
	{
		Output::newline();
		for (j = 0; j < 3; ++j)
		{
			Output::print(Language::numberToFraction(centeringVectors[i][j]));
			Output::print(" ");
		}
	}
	Output::decrease();
}



/* void Symmetry::setMetricMatrixConstraint()
 *
 * Set the constraint matrix for metric matrix
 */

void Symmetry::setMetricMatrixConstraint()
{
	
	// Make list of operations
	int i;
	Linked<const Rotation*> rotList;
	for (i = 0; i < _operations.length(); ++i)
		rotList += &(_operations[i].rotation());
	
	// Loop over symmetry operations to build full lattice symmetry constraint matrix
	int j, p, q;
	int curRow = 0;
	int curCol;
	Matrix fullMat(6 * rotList.length(), 6);
	for (Linked<const Rotation*>::iterator it = rotList.begin(); it != rotList.end(); ++it)
	{
		
		// Set components of metric matrix
		const Rotation& rot = **it;
		for (i = 0; i < 3; ++i)
		{
			for (j = i; j < 3; ++j)
			{
				curCol = 0;
				for (p = 0; p < 3; ++p)
				{
					for (q = p; q < 3; ++q)
					{
						fullMat(curRow, curCol) = (1 - Num<int>::delta(p, q)/2.0) * \
							(rot(p, i)*rot(q, j) + rot(q, i)*rot(p, j) - Num<int>::delta(p, i)*Num<int>::delta(q, j) -\
							Num<int>::delta(q, i)*Num<int>::delta(p, j));
						++curCol;
					}
				}
				++curRow;
			}
		}
	}
	
	// Save first rows of row echelon form of matrix
	if (fullMat.numRows() > 5)
	{
		fullMat = fullMat.rowEchelon();
		for (i = 0; i < 6; ++i)
		{
			for (j = 0; j < 6; ++j)
				_metricMatrixConstraint(i, j) = fullMat(i, j);
		}
	}
	else
		_metricMatrixConstraint.fill(0);
}



/* void Symmetry::setOrbits(const ISO& iso, double tol)
 *
 * Set the wyckoff positions and their orbits
 */

void Symmetry::setOrbits(const ISO& iso, double tol)
{
	
	// Set size of orbits list
	_orbitNumbers.length(iso.numAtoms());
	
	// Loop over elements
	int i, j, k;
	int count;
	double curDistance;
	Atom curAtom;
	Vector3D origin(0.0);
	Vector3D curPos;
	Vector3D curSpecialTranslation;
	Matrix3D curSpecialRotation;
	OList<Matrix3D> curPointRotations;
	OList<Vector3D> curPointTranslations;
	Linked<int> curGens;
	Linked<int>::iterator itGens;
	Linked<int> curTrans;
	Linked<int>::iterator itTrans;
	Linked<Atom*> atoms;
	Linked<Atom*>::iterator itAtom;
	Linked<Atom*> curOrbit;
	Linked<Atom*>::iterator itOrbit;
	Linked<double> distances;
	Linked<double>::iterator itDist;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Create list of atoms of current element
		atoms.clear();
		distances.clear();
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			atoms.add(&(iso.atoms()[i][j]));
			distances.add(iso.basis().distance(iso.atoms()[i][j].fractional(), FRACTIONAL, origin, FRACTIONAL));
		}
		
		// Loop until no atoms are left
		while (atoms.length())
		{
					
			// Set current atom
			curAtom  = **(atoms.begin());
			
			// Loop over symmetry operations
			count = 0;
			curSpecialRotation = 0.0;
			curSpecialTranslation = 0.0;
			curGens.clear();
			curTrans.clear();
			curOrbit.clear();
			curPointRotations.length(0);
			curPointTranslations.length(0);
			for (j = 0; j < _operations.length(); ++j)
			{
				for (k = 0; k < _operations[j].translations().length(); ++k)
				{
					
					// Generate new position
					curPos = _operations[j].rotation() * (*atoms.begin())->fractional();
					curAtom.fractional(curPos + _operations[j].translations()[k]);
					curDistance = iso.basis().distance(curAtom.fractional(), FRACTIONAL, origin, FRACTIONAL);
					
					// Loop over atoms of current element
					for (itAtom = atoms.begin(), itDist = distances.begin(); itAtom != atoms.end(); ++itAtom, ++itDist)
					{
						
						// Check if atoms are the same
						if (Num<double>::abs(curDistance - *itDist) > tol)
							continue;
						if (!curAtom.equal(**itAtom, tol))
							continue;
						
						// Returned same atom
						if (itAtom == atoms.begin())
						{
							++count;
							curSpecialRotation += _operations[j].rotation();
							curSpecialTranslation += (*atoms.begin())->fractional() - curPos;
							curPointRotations += _operations[j].rotation();
							curPointTranslations += (*atoms.begin())->fractional() - curPos;
						}
						
						// Found a new atom
						else
						{
							curGens.add(j);
							curTrans.add(k);
							curOrbit.add(*itAtom);
							atoms.remove(itAtom);
							distances.remove(itDist);
						}
						
						// Break since atom was found
						break;
					}
				}
			}
			
			// Add new orbit
			_orbits.add();
			curSpecialRotation /= count;
			curSpecialTranslation /= count;
			_orbits.last().set(*(atoms.begin()), curSpecialRotation, curSpecialTranslation, curPointRotations, \
				curPointTranslations);
			_orbitNumbers[(*atoms.begin())->atomNumber()] = _orbits.length() - 1;
			
			// Remove current atom
			atoms.remove(atoms.begin());
			distances.remove(distances.begin());
			
			// Save atoms in orbit
			itGens = curGens.begin();
			itTrans = curTrans.begin();
			for (itOrbit = curOrbit.begin(); itOrbit != curOrbit.end(); ++itOrbit, ++itGens, ++itTrans)
			{
				_orbits.last().add(*itOrbit, _operations[*itGens].rotation(), \
					_operations[*itGens].translations()[*itTrans]);
				_orbitNumbers[(*itOrbit)->atomNumber()] = _orbits.length() - 1;
			}
		}
	}
	
	// Figure out which operation is the identity
	Matrix3D identity = Matrix3D::identity();
	int Inum = 0;
	for (i = 0; i < _operations.length(); ++i)
	{
		if (_operations[i].rotation() == identity)
		{
			Inum = i;
			break;
		}
	}
	
	// For every atom, generate its full site symmetry set
	int m;
	Vector3D rotPos;
	for (i = 0; i < _orbits.length(); ++i)
	{
		for (j = 1; j < _orbits[i].atoms().length(); ++j)
		{
			
			// Loop over symmetry operations and save those that map atom onto itself
			for (k = 0; k < _operations.length(); ++k)
			{
				
				// Skip if identity
				if (k == Inum)
					continue;
				
				// Loop over translations
				rotPos  = _orbits[i].atoms()[j]->fractional();
				rotPos *= _operations[k].rotation();
				for (m = 0; m < _operations[k].translations().length(); ++m)
				{
					
					// Save current position
					curPos  = rotPos;
					curPos += _operations[k].translations()[m];
					ISO::moveIntoCell(curPos);
					
					// Skip if sites are not the same
					if (iso.basis().distance(curPos, FRACTIONAL, _orbits[i].atoms()[j]->fractional(), FRACTIONAL) > tol)
						continue;
					
					// Save and finish
					_orbits[i].specialPositions()[j].addOperation(_operations[k].rotation(), \
						_operations[k].translations()[m]);
					break;
				}
			}
		}
	}
	
	// Output
	Output::newline();
	Output::print("Found ");
	Output::print(_orbits.length());
	Output::print(" unique position");
	if (_orbits.length() != 1)
		Output::print("s");
	Output::increase();
	
	// Print all positions
	for (i = 0; i < _orbits.length(); ++i)
	{
		Output::newline();
		Output::print(_orbits[i].atoms()[0]->element().symbol());
		Output::print(" at ");
		Output::print(_orbits[i].atoms()[0]->fractional(), 8);
		Output::print(" (");
		Output::print(_orbits[i].specialPositions()[0].getString(), true, false);
		Output::print(") with multiplicity of ");
		Output::print(_orbits[i].atoms().length());
	}
	
	// Output
	Output::decrease();
}



/* void Symmetry::print(bool useJonesFaithful) const
 *
 * Print symmetry data
 */

void Symmetry::print(bool useJonesFaithful) const
{
	
	// Set output method
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Create message to store symmetry operations
	int i, j, k;
	Output message;
	List<PrintAlign> align;
	if (useJonesFaithful)
	{
		align.length(4);
		align.fill(LEFT);
		message.addLines(_operations.length());
		for (i = 0; i < _operations.length(); ++i)
		{
			message.addLine();
			message.add("   ");
			message.add(_operations[i].getString());
		}
	}
	else
	{
		align.length(5);
		align.fill(RIGHT);
		message.addLines(5 * _operations.length() - 1);
		for (i = 0; i < _operations.length(); ++i)
		{
			if (i > 0)
				message.addLine();
			for (j = 0; j < 3; ++j)
			{
				message.addLine();
				message.addWords(5);
				message.add("   ");
				for (k = 0; k < 3; ++k)
					message.add(_operations[i].rotation()(j, k));
				message.add(_operations[i].translations()[0][j]);
			}
			message.addLine();
			message.addWords(5);
			message.add("   ");
			for (k = 0; k < 3; ++k)
				message.add(0.0);
			message.add(1.0);
		}
	}
	
	// Print rotations
	Output::newline();
	Output::print("Number of unique symmetry operations: ");
	Output::print(_operations.length());
	Output::newline();
	Output::print(message, align);
	
	// Create message to store translations
	message.clear();
	if (_operations.length())
	{
		message.addLines(_operations[0].translations().length());
		for (i = 0; i < _operations[0].translations().length(); ++i)
		{
			message.addLine();
			message.addWords(4);
			message.add("   ");
			if (useJonesFaithful)
			{
				for (j = 0; j < 3; ++j)
					message.add(Language::numberToFraction(_operations[0].translations()[i][j]));
			}
			else
			{
				for (j = 0; j < 3; ++j)
					message.add(_operations[0].translations()[i][j]);
			}
		}
	}
	
	// Print translations
	Output::newline();
	Output::newline();
	Output::print("Number of centering vectors: ");
	Output::print(message.numLines());
	Output::newline();
	Output::print(message, RIGHT);
	
	// Reset output
	Output::method(origMethod);
}



/* void Symmetry::printSites() const
 *
 * Print unique sites in the structure
 */

void Symmetry::printSites() const
{
	
	// Set output method
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Make list of sites
	int i, j, k;
	OList<Words> wyckoff(_orbits.length());
	OList<Output> sites(_orbits.length());
	for (i = 0; i < _orbits.length(); ++i)
	{
		wyckoff[i] = _orbits[i].specialPositions()[0].getString();
		sites[i].addLines(_orbits[i].atoms().length());
		for (j = 0; j < _orbits[i].atoms().length(); ++j)
		{
			sites[i].addLine();
			sites[i].add("       ");
			for (k = 0; k < 3; ++k)
				sites[i].add(_orbits[i].atoms()[j]->fractional()[k], 8);
		}
	}
	
	// Print sites
	Output::newline();
	Output::print("Number of unique atoms: ");
	Output::print(_orbits.length());
	for (i = 0; i < wyckoff.length(); ++i)
	{
		Output::newline();
		Output::print("    ");
		Output::print(_orbits[i].atoms()[0]->element().symbol());
		Output::print(" at (");
		Output::print(wyckoff[i], true, false);
		Output::print(") with multiplicity ");
		Output::print(_orbits[i].atoms().length());
		Output::newline();
		Output::print(sites[i], RIGHT);
	}
	
	// Reset output
	Output::method(origMethod);
}



/* int Symmetry::rotationOrder(const Matrix3D& rotation)
 *
 * Return the rotation order for a matrix
 */

int Symmetry::rotationOrder(const Matrix3D& rotation)
{
	
	// Get determinant and trace of rotation
	double det = rotation.determinant();
	double trc = rotation.trace();
	
	// Get rotation order
	double tol = 1e-4;
	if (Num<double>::abs(det - 1) < tol)
	{
		if (Num<double>::abs(trc + 1) < tol)
			return 2;
		else if (Num<double>::abs(trc) < tol)
			return 3;
		else if (Num<double>::abs(trc - 1) < tol)
			return 4;
		else if (Num<double>::abs(trc - 2) < tol)
			return 6;
		else if (Num<double>::abs(trc - 3) < tol)
			return 1;
	}
	else if (Num<double>::abs(det + 1) < tol)
	{
		if (Num<double>::abs(trc + 3) < tol)
			return -1;
		else if (Num<double>::abs(trc + 2) < tol)
			return -6;
		else if (Num<double>::abs(trc + 1) < tol)
			return -4;
		else if (Num<double>::abs(trc) < tol)
			return -3;
		else if (Num<double>::abs(trc - 1) < tol)
			return -2;
	}
	
	// Could not get rotation order
	Output::newline(ERROR);
	Output::print("Internal error: could not identify rotation order");
	Output::quit();
	return 0;
}



/* Vector3D Symmetry::intrinsicTranslation(const Matrix3D& rotation, 
 *		const Vector3D& translation)
 *
 * Return the intrinsic part of a translation in a symmetry operation
 */

Vector3D Symmetry::intrinsicTranslation(const Matrix3D& rotation, const Vector3D& translation)
{
	
	// Get the order of the rotation
	int order = rotationOrder(rotation);
	
	// Get number of times to multiply matrix to get back to identity
	int numMult = Num<int>::abs(order);
	if ((order == -1) || (order == -3))
		numMult = -2*order;
	
	// Get translation
	int i;
	Vector3D newTrans = translation;
	for (i = 1; i < numMult; ++i)
		newTrans = rotation * newTrans + translation;
	
	// Return result
	for (i = 0; i < 3; ++i)
		newTrans[i] /= numMult;
	return newTrans;
}



/* int Symmetry::getUniqueAxis(const Matrix3D& vectors, LatticeSystem system)
 *
 * Get the unique axis for a cell
 */

int Symmetry::getUniqueAxis(const Matrix3D& vectors, LatticeSystem system)
{
	
	// Get lattice metrics
	double a = sqrt(vectors(0, 0)*vectors(0, 0) + vectors(0, 1)*vectors(0, 1) + vectors(0, 2)*vectors(0, 2));
    double b = sqrt(vectors(1, 0)*vectors(1, 0) + vectors(1, 1)*vectors(1, 1) + vectors(1, 2)*vectors(1, 2));
    double c = sqrt(vectors(2, 0)*vectors(2, 0) + vectors(2, 1)*vectors(2, 1) + vectors(2, 2)*vectors(2, 2));
	Vector3D vec0(vectors(0, 0), vectors(0, 1), vectors(0, 2));
	Vector3D vec1(vectors(1, 0), vectors(1, 1), vectors(1, 2));
	Vector3D vec2(vectors(2, 0), vectors(2, 1), vectors(2, 2));
	double alpha = vec1.angle(vec2);
    double beta  = vec0.angle(vec2);
    double gamma = vec0.angle(vec1);
	
	// Get the unique axis if tetragonal or hexagonal
    if ((system == LS_TETRAGONAL) || (system == LS_HEXAGONAL))
    {
        if ((Num<double>::abs(a - c) < Num<double>::abs(a - b)) && \
            (Num<double>::abs(a - c) < Num<double>::abs(b - c)))
			return 1;
        else if ((Num<double>::abs(a - b) < Num<double>::abs(a - c)) && \
            	 (Num<double>::abs(a - b) < Num<double>::abs(b - c)))
			return 2;
    }

    // Get the unique angle if monoclinic
    else if (system == LS_MONOCLINIC)
    {
		double ninety = Constants::pi / 2;
        if ((Num<double>::abs(beta - ninety) > Num<double>::abs(alpha - ninety)) && \
            (Num<double>::abs(beta - ninety) > Num<double>::abs(gamma - ninety)))
			return 1;
        else if ((Num<double>::abs(gamma - ninety) > Num<double>::abs(alpha - ninety)) && \
                 (Num<double>::abs(gamma - ninety) > Num<double>::abs(beta  - ninety)))
			return 2;
    }

	// Return 0 otherwise
	return 0;
}



/* int Symmetry::getConventionalCellUniqueSymmetryAxis(const OList<Matrix3D>& rotations, LatticeSystem system)
 *
 * Get the unique axis for a cell based on the symmetry operations
 */

int Symmetry::getConventionalCellUniqueSymmetryAxis(const OList<Matrix3D>& rotations, LatticeSystem system)
{
	
	// Get the unique axis if tetragonal, hexagonal, or monoclinic
    if ((system == LS_TETRAGONAL) || (system == LS_HEXAGONAL) || (system == LS_MONOCLINIC))
    {
		int i, j, k;
		int index[2];
		double tol = 1e-4;
		for (i = 0; i <= 2; ++i)
		{
			index[0] = Num<int>::next(i, 2);
			index[1] = Num<int>::next(index[0], 2);
			for (j = 0; j < rotations.length(); ++j)
			{
				for (k = 0; k < 3; ++k)
				{
					if (k == i)
					{
						if (Num<double>::abs(rotations[j](k, index[0])) > tol)
							break;
						if (Num<double>::abs(rotations[j](k, index[1])) > tol)
							break;
					}
					else
					{
						if (Num<double>::abs(rotations[j](k, i)) > tol)
							break;
					}
				}
				if (k < 3)
					break;
			}
			if (j >= rotations.length())
				return i;
		}
	}

	// Return 0 otherwise
	return 0;
}



/* void Symmetry::confineBasis(Vector3D& lengths, Vector3D& angles, LatticeSystem latticeSystem)
 *
 * Adjust cell so that it satisfies lattice symmetry
 */

void Symmetry::confineBasis(Vector3D& lengths, Vector3D& angles, LatticeSystem latticeSystem)
{
	
	// Confine to cubic system
	if (latticeSystem == LS_CUBIC)
	{
		lengths[0] = lengths[1] = lengths[2] = (lengths[0] + lengths[1] + lengths[2]) / 2;
		angles[0] = angles[1] = angles[2] = Constants::pi / 2;
	}
	
	// Confine to hexagonal system
	else if (latticeSystem == LS_HEXAGONAL)
	{
		lengths[0] = lengths[1] = (lengths[0] + lengths[1]) / 2;
		angles[0] = angles[1] = Constants::pi / 2;
		angles[2] = 2 * Constants::pi / 3;
	}
	
	// Confine to tetragonal system
	else if (latticeSystem == LS_TETRAGONAL)
	{
		lengths[0] = lengths[1] = (lengths[0] + lengths[1]) / 2;
		angles[0] = angles[1] = angles[2] = Constants::pi / 2;
	}
	
	// Confine to orthorhombic system
	else if (latticeSystem == LS_ORTHORHOMBIC)
		angles[0] = angles[1] = angles[2] = Constants::pi / 2;
	
	// Confine to monoclinic system
	else if (latticeSystem == LS_MONOCLINIC)
		angles[0] = angles[2] = Constants::pi / 2;
}



/* void Symmetry::getFullMap(List<Atom*>& map, const ISO& iso, const Matrix3D& rotation,
 *		const Vector3D& translation) const
 *
 * Get mapping of all atoms under a set symmetry operation
 */

void Symmetry::getFullMap(List<Atom*>& map, const ISO& iso, const Matrix3D& rotation, \
	const Vector3D& translation) const
{
	
	// Make sure there is enough room in map
	map.length(iso.numAtoms());
	
	// Loop over orbits
	int i, j;
	double tol = 0.01;
	double curDis;
	double nearDis;
	double origDis;
	Vector3D newPos;
	Vector3D origin(0.0);
	Linked<Atom*> atoms;
	Linked<double> distances;
	Linked<Atom*>::iterator itAtom;
	Linked<Atom*>::iterator nearAtom;
	Linked<double>::iterator itDist;
	Linked<double>::iterator nearDist;
	for (i = 0; i < _orbits.length(); ++i)
	{
		
		// Get distance of each atom in orbit from the origin
		for (j = 0; j < _orbits[i].atoms().length(); ++j)
		{
			atoms += _orbits[i].atoms()[j];
			distances += iso.basis().distance(origin, FRACTIONAL, _orbits[i].atoms()[j]->fractional(), FRACTIONAL);
		}
		
		// Loop over atoms in orbit and generate new positions
		for (j = 0; j < _orbits[i].atoms().length(); ++j)
		{
			
			// Get new position
			newPos  = rotation * _orbits[i].atoms()[j]->fractional();
			newPos += translation;
			ISO::moveIntoCell(newPos);
			
			// Get current distance from origin
			origDis = iso.basis().distance(origin, FRACTIONAL, newPos, FRACTIONAL);
			
			// Initialize nearest atom as first
			nearAtom = atoms.begin();
			nearDist = distances.begin();
			nearDis  = iso.basis().distance(newPos, FRACTIONAL, (*nearAtom)->fractional(), FRACTIONAL);
			
			// Loop over remaining positions to get nearest atom
			for (itAtom = atoms.begin() + 1, itDist = distances.begin() + 1; itAtom != atoms.end(); ++itAtom, ++itDist)
			{
				
				// Current distances from origin do not match
				if (Num<double>::abs(origDis - *itDist) > nearDis)
					continue;
				
				// Current distance between atoms do not match
				curDis = iso.basis().distance(newPos, FRACTIONAL, (*itAtom)->fractional(), FRACTIONAL);
				if (curDis > nearDis)
					continue;
				
				// Save new nearest atom
				nearDis = curDis;
				nearAtom = itAtom;
				nearDist = itDist;
				
				// Break if an exact match
				if (curDis < tol)
					break;
			}
			
			// Save match and remove atom from list
			map[_orbits[i].atoms()[j]->atomNumber()] = *nearAtom;
			atoms.remove(nearAtom);
			distances.remove(nearDist);
		}
	}
}



/* void Symmetry::refine(ISO& iso, double clusterTol) const
 *
 * Refine structure
 */

void Symmetry::refine(ISO& iso, double clusterTol) const
{
	
	// Output
	Output::newline();
	Output::print("Refining structure");
	Output::increase();
	
	// Refine basis and positions
	refineBasis(iso);
	refineAtoms(iso, clusterTol);
	
	// Output
	Output::decrease();
}



/* void Symmetry::refineBasis(ISO& iso) const
 *
 * Refine basis
 */

void Symmetry::refineBasis(ISO& iso) const
{
	
	// Output
	Output::newline();
	Output::print("Refining basis");
	Output::increase();
	
	// Refine vectors
	Matrix3D vectors = iso.basis().vectors();
	refineBasis(vectors);
	
	// Save refined basis
	iso.basis(vectors, true);
	
	// Output
	Output::decrease();
}



/* void Symmetry::refineBasis(Matrix3D& vectors) const
 *
 * Refine basis vectors
 */

void Symmetry::refineBasis(Matrix3D& vectors) const
{
	
	// Set initial metric matrix values
	Vector gOrig(6);
	double vec0[3] = {vectors(0, 0), vectors(0, 1), vectors(0, 2)};
	double vec1[3] = {vectors(1, 0), vectors(1, 1), vectors(1, 2)};
	double vec2[3] = {vectors(2, 0), vectors(2, 1), vectors(2, 2)};
	gOrig[0] = Num<double>::dot(3, vec0, vec0);
	gOrig[1] = Num<double>::dot(3, vec0, vec1);
	gOrig[2] = Num<double>::dot(3, vec0, vec2);
	gOrig[3] = Num<double>::dot(3, vec1, vec1);
	gOrig[4] = Num<double>::dot(3, vec1, vec2);
	gOrig[5] = Num<double>::dot(3, vec2, vec2);
	
	// Build refined metric matrix values
	int i, j, k;
	double tol = 1e-4;
	Vector gRef = gOrig;
	for (i = 5; i >= 0; --i)
	{
		
		// Find first non-zero value on row
		for (j = 0; j < 6; ++j)
		{
			if (Num<double>::abs(_metricMatrixConstraint(i, j)) >= tol)
			{
				
				// Set metric value
				gRef[j] = 0;
				for (k = j + 1; k < 6; ++k)
					gRef[j] -= _metricMatrixConstraint(i, k) * gRef[k];
				gRef[j] /= _metricMatrixConstraint(i, j);
				
				// Multiply values by -1 if needed
				if (((j == 0) || (j == 3) || (j == 5)) && (gRef[j] < 0))
				{
					if (gRef[j] < 0)
					{
						for (k = j; k < 6; ++k)
						{
							if (Num<double>::abs(_metricMatrixConstraint(i, k)) >= tol)
								gRef[k] *= -1;
						}
					}
				}
				
				// Break since set
				break;
			}
		}
	}
	
	// Get new lengths and angles
	Vector3D lengths(sqrt(gRef[0]), sqrt(gRef[3]), sqrt(gRef[5]));
	Vector3D angles(acos(gRef[4] / (lengths[1]*lengths[2])), acos(gRef[2] / (lengths[0]*lengths[2])), \
		acos(gRef[1] / (lengths[0]*lengths[1])));
	
	// Save refined basis
	vectors = Basis::vectors(lengths, angles);
}



/* void Symmetry::refineAtoms(ISO& iso, double clusterTol) const
 *
 * Refine positions
 */

void Symmetry::refineAtoms(ISO& iso, double clusterTol) const
{
	
	// Output
	Output::newline();
	Output::print("Refining atomic positions");
	Output::increase();
	
	// Refine the unique positions
	OList<Atom> uniqueAtoms;
	refineUniquePositions(iso, uniqueAtoms, clusterTol);
	
	// Cluster atom positions
	clusterAtoms(iso, uniqueAtoms, clusterTol);
	
	// Output
	Output::decrease();
}



/* void Symmetry::refineUniquePositions(ISO& iso, OList<Atom>& uniqueAtoms, double clusterTol) const
 *
 * Refine unique positions in the structure
 */

void Symmetry::refineUniquePositions(ISO& iso, OList<Atom>& uniqueAtoms, double clusterTol) const
{
	
	// Clear space
	uniqueAtoms.clear();
	
	// Loop over unique atoms in the structure
	int i, j, k;
	Vector3D position;
	Vector3D curPosition;
	Vector3D nearImage;
	for (i = 0; i < _orbits.length(); ++i)
	{
		
		// Average positions in orbit
		position = _orbits[i].atoms()[0]->fractional();
		for (j = 1; j < _orbits[i].atoms().length(); ++j)
		{
			
			// Generate current position from symmetry
			curPosition = _orbits[i].generators()[j].rotation().inverse() * \
				(_orbits[i].atoms()[j]->fractional() - _orbits[i].generators()[j].translations()[0]);
			ISO::moveIntoCell(curPosition);
			
			// Get nearest image of position
			iso.basis().distance(position, FRACTIONAL, curPosition, FRACTIONAL, &nearImage);
			
			// Average position
			for (k = 0; k < 3; ++k)
				position[k] = (j*position[k] + curPosition[k] + nearImage[k]) / (j + 1);
			ISO::moveIntoCell(position);
		}
		
		// Project atom onto special position
		position -= _orbits[i].specialPositions()[0].translation();
		position *= _orbits[i].specialPositions()[0].rotation();
		position += _orbits[i].specialPositions()[0].translation();
		ISO::moveIntoCell(position);
		
		// Save atom
		uniqueAtoms += *(_orbits[i].atoms()[0]);
		uniqueAtoms.last().fractional(position);
	}
}



/* void Symmetry::clusterAtoms(ISO& iso, const OList<Atom>& uniqueAtoms, double clusterTol) const
 *
 * Cluster atom positions
 */

void Symmetry::clusterAtoms(ISO& iso, const OList<Atom>& uniqueAtoms, double clusterTol) const
{
	
	// Generate full list of sites for each unique atom
	int i, j, k, m;
	bool found;
	Vector3D rotPos;
	Vector3D position;
	Vector3D nearCell;
	Linked<int> timesFound;
	Linked<int>::iterator itFound;
	Linked<Vector3D >::iterator itPos;
	OList<Linked<Vector3D > > positions(uniqueAtoms.length());
	for (i = 0; i < uniqueAtoms.length(); ++i)
	{
		
		// Loop over symmetry operations
		timesFound.clear();
		for (j = 0; j < _operations.length(); ++j)
		{
			
			// Loop over translations
			rotPos = _operations[j].rotation() * uniqueAtoms[i].fractional();
			for (k = 0; k < _operations[j].translations().length(); ++k)
			{
				
				// Generate current position
				position = rotPos + _operations[j].translations()[k];
				ISO::moveIntoCell(position);
				
				// Check whether position is already known
				found = false;
				itFound = timesFound.begin();
				itPos = positions[i].begin();
				for (; itPos != positions[i].end(); ++itPos, ++itFound)
				{
					if (iso.basis().distance(*itPos, FRACTIONAL, position, FRACTIONAL, &nearCell) < clusterTol)
					{
						for (m = 0; m < 3; ++m)
							(*itPos)[m] = (*itFound * (*itPos)[m] + position[m] + nearCell[m]) / (*itFound + 1);
						ISO::moveIntoCell(*itPos);
						++(*itFound);
						found = true;
						break;
					}
				}
				
				// Found a new position
				if (!found)
				{
					positions[i] += position;
					timesFound += 1;
				}
			}
		}
	}
	
	// Add atoms
	iso.clearAtoms();
	Atom* atom;
	for (i = 0; i < uniqueAtoms.length(); ++i)
	{
		for (itPos = positions[i].begin(); itPos != positions[i].end(); ++itPos)
		{
			
			// Set atom
			atom = iso.addAtom(uniqueAtoms[i]);
			atom->fractional(*itPos);
			
			// Output
			Output::newline();
			Output::print("Atom ");
			Output::print(atom->atomNumber() + 1);
			Output::print(" (");
			Output::print(atom->element().symbol());
			Output::print(") at ");
			Output::print(*itPos, 8);
		}
	}
}



/* void Symmetry::addAtom(ISO& iso, const Atom& atom, const OList<SymmetryOperation>& operations, double clusterTol)
 *
 * Add atom to structure and use symmetry to generate equivalent positions
 */

void Symmetry::addAtom(ISO& iso, const Atom& atom, const OList<SymmetryOperation>& operations, double clusterTol)
{
	
	// Loop over symmetry operations
	int i, j;
	bool found;
	Linked<int> timesFound;
	Vector3D rotPos;
	Vector3D position;
	Vector3D nearCell;
	Linked<Vector3D > positions;
	Linked<int>::iterator itFound;
	Linked<Vector3D>::iterator itPos;
	for (i = 0; i < operations.length(); ++i)
	{
		rotPos = operations[i].rotation() * atom.fractional();
		for (j = 0; j < operations[i].translations().length(); ++j)
		{
			
			// Get current position
			position = rotPos + operations[i].translations()[j];
			ISO::moveIntoCell(position);
			
			// Check if position is already saved
			found = false;
			for (itFound = timesFound.begin(), itPos = positions.begin(); itPos != positions.end(); ++itFound, ++itPos)
			{
				
				// Positions are the same
				if (iso.basis().distance(*itPos, FRACTIONAL, position, FRACTIONAL, &nearCell) < clusterTol)
				{
					*itPos = (*itPos * *itFound + position + nearCell) / (*itFound + 1);
					ISO::moveIntoCell(*itPos);
					++(*itFound);
					found = true;
					break;
				}
			}
			
			// Found a new position
			if (!found)
			{
				positions += position;
				timesFound += 1;
			}
		}
	}
	
	// Loop over atoms to add
	Atom* newAtom;
	for (itPos = positions.begin(); itPos != positions.end(); ++itPos)
	{
		
		// Check if position is already in the structure
		found = false;
		for (i = 0; i < iso.atoms().length(); ++i)
		{
			if (iso.atoms()[i][0].element() != atom.element())
				continue;
			for (j = 0; j < iso.atoms()[i].length(); ++j)
			{
				
				// Positions are the same
				if (iso.basis().distance(iso.atoms()[i][j].fractional(), FRACTIONAL, *itPos, FRACTIONAL) < clusterTol)
				{
					found = true;
					break;
				}
			}
			break;
		}
		if (found)
			continue;
		
		// Add atom
		newAtom = iso.addAtom(atom);
		newAtom->fractional(*itPos);
		
		// Output
		Output::newline();
		Output::print("Adding atom ");
		Output::print(newAtom->atomNumber() + 1);
		Output::print(" at ");
		Output::print(*itPos, 8);
	}
}



/* Matrix3D Symmetry::idealTransformation(const ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol)
 *
 * Return the transformation to the most ideal supercell
 */

Matrix3D Symmetry::idealTransformation(const ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol)
{
	
	// Figure out whether to run to number of atoms or minimum distance
	bool runForMinDis = (minDis < 0) ? false : true;
	
	// Output
	Output::newline();
	Output::print("Searching for the most ideal cell ");
	if (runForMinDis)
	{
		Output::print("with minimum image distance of ");
		Output::print(minDis);
		Output::print(" ang");
	}
	else
	{
		Output::print("containing at ");
		if (numAtomsAsMin)
			Output::print("least ");
		else
			Output::print("most ");
		Output::print(numAtoms);
		Output::print(" atom");
		if (numAtoms != 1)
			Output::print("s");
	}
	Output::increase();
	
	// Output
	Output::newline();
	Output::print("Converting to reduced primitive cell and calculating symmetry");
	Output::increase();
	
	// Get the transformation to the reduced primitive cell
	Matrix3D unitToPrim = iso.primitiveTransformation(tol, true);
	Matrix3D primToRed = Basis::reducedTransformation(unitToPrim * iso.basis().vectors());
	
	// Get the primitive reduced cell
	ISO primRedISO = iso;
	primRedISO.transform(primToRed * unitToPrim, tol);
	
	// Get the symmetry of the primitive reduced cell
	Symmetry primRedSymm;
	primRedSymm.set(primRedISO, tol, true);
	
	// Output
	Output::decrease();
	
	// Get the target number of cells
	double atomRatio = numAtomsAsMin == true ? Num<double>::ceil(numAtoms / primRedISO.numAtoms()) : \
		Num<double>::floor(numAtoms / primRedISO.numAtoms());
	double numCells = Num<double>::round(atomRatio, 1);
	numCells = numCells == 0 ? 1 : numCells;
	if ((runForMinDis) || ((runForMinDis == false) && (numAtomsAsMin == true)))
		numCells = 10000;
	
	// Initialize the transformation matrix
	Matrix3D resTrans = Matrix3D::identity();
	Matrix3D startTrans = resTrans;
	
	// Loop until the target number of cells is reached
	int i;
	bool first;
	int curNum;
	int bestNum;
	double range;
	double curLen;
	double minLen = 0;
	double resLen = 0;
	Vector3D vec;
	Vector3D curVec;
	Matrix3D P;
	Matrix3D Q;
	Matrix3D modMat;
	Matrix3D curTrans(0.0);
	Matrix3D bestTrans;
	Matrix3D transformedCell;
	for (double curDet = 2; curDet <= numCells; ++curDet)
	{
			
		// Loop over range from -1 to 1, and -2 to 2 if needed
		for (range = 1; range <= 2; ++range)
		{
			
			// Loop over all permutations of vectors to add
			first = true;
			for (modMat(0, 0) = -range; modMat(0, 0) <= range+1e-6; modMat(0, 0) += 1)
			{
				for (modMat(0, 1) = -range; modMat(0, 1) <= range+1e-6; modMat(0, 1) += 1)
				{
					for (modMat(0, 2) = -range; modMat(0, 2) <= range+1e-6; modMat(0, 2) += 1)
					{
						for (modMat(1, 0) = -range; modMat(1, 0) <= range+1e-6; modMat(1, 0) += 1)
						{
							for (modMat(1, 1) = -range; modMat(1, 1) <= range+1e-6; modMat(1, 1) += 1)
							{
								for (modMat(1, 2) = -range; modMat(1, 2) <= range+1e-6; modMat(1, 2) += 1)
								{
									for (modMat(2, 0) = -range; modMat(2, 0) <= range+1e-6; modMat(2, 0) += 1)
									{
										for (modMat(2, 1) = -range; modMat(2, 1) <= range+1e-6; modMat(2, 1) += 1)
										{
											for (modMat(2, 2) = -range; modMat(2, 2) <= range+1e-6; modMat(2, 2) += 1)
											{
								
												// Save transformation
												curTrans = startTrans;
												curTrans += modMat;

												// Skip if determinant is not correct
												if (Num<double>::neq(curTrans.determinant(), curDet, 1e-4))
													continue;
								
												// Get the reduced cell
												curTrans *= Basis::reducedTransformation(curTrans * \
													primRedISO.basis().vectors());
												transformedCell = curTrans * primRedISO.basis().vectors();
								
												// Get the length of the A vector (will be the shortest since reduced)
												vec[0] = transformedCell(0, 0);
												vec[1] = transformedCell(0, 1);
												vec[2] = transformedCell(0, 2);
												curLen = vec * vec;
							
												// Check whether operation preserves all symmetries as integer
												P = curTrans.transpose();
												Q = P.inverse();
												curNum = 0;
												for (i = primRedSymm.operations().length() - 1; i >= 0; --i)
												{
													if ((Q * primRedSymm.operations()[i].rotation() * \
														 P).isInteger(1e-2) == true)
														++curNum;
												}
											
												// Save if the number of symmetries improved
												if ((first) || (curNum > bestNum))
												{
													minLen = curLen;
													bestNum = curNum;
													bestTrans = curTrans;
												}
											
												// If number of symmetries was the same then check if length was better
												if ((curNum == bestNum) && (curLen > minLen))
												{
													minLen = curLen;
													bestTrans = curTrans;
												}
											
												// Save that a loop was made
												first = false;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			
			// Break condition to improve efficiency
			if (bestNum >= primRedSymm.operations().length()/2)
				break;
		}
		
		// Save best transformation for current determinant
		startTrans = bestTrans;
		
		// Only process if a valid cell was found
		if (bestNum == primRedSymm.operations().length())
		{
			
			// Not running for distance
			if (runForMinDis == false)
			{
				
				// Running for a maximum number of atoms - save only if minimum distance improved
				if (numAtomsAsMin == false)
				{
					if (sqrt(minLen) > resLen + 1e-6)
					{
						resLen = sqrt(minLen);
						resTrans = bestTrans;
					}
				}
				
				// Running for a minimum number of atoms - save no matter what
				// Break if the number of atoms was passed
				else
				{
					resTrans = bestTrans;
					if (curDet * primRedISO.numAtoms() >= numAtoms)
						break;
				}
					
			}
			
			// Save an break if running for minimum distance and it has been passed
			else if (sqrt(minLen) > minDis)
			{
				resTrans = bestTrans;
				break;
			}
		}
	}
	
	// Save transformation
	Matrix3D transformation = unitToPrim;
	transformation *= primToRed;
	transformation *= resTrans;
	
	// Output
	int j;
	Output::newline();
    Output::print("Relative to the unit cell, the most ideal basis vectors are:");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(transformation(i, j), 8);
			Output::print(" ");
		}
	}
	
	// Output
	Output::decrease();

	// Return the transformation matrix
	return transformation;
}



/* Matrix3D Symmetry::makeIdeal(ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol)
 *
 * Convert to ideal cell
 */

Matrix3D Symmetry::makeIdeal(ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol)
{
	
	// Output
	Output::newline();
	Output::print("Converting cell to ideal form");
	Output::increase();
	
	// Get the transformation to the ideal cell
	Matrix3D transformation = idealTransformation(iso, numAtoms, numAtomsAsMin, minDis, tol);
	
	// Set the new structure
	iso.transform(transformation, tol, true);
	
	// Output
	Output::decrease();
	
	// Return the transformation
	return transformation;
}
