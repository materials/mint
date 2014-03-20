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



#include "phonons.h"
#include "language.h"
#include "constants.h"
#include "output.h"
#include <cmath>
#include <cstdlib>



/* Word Phonons::generateForceConstants(const ISO& iso, const Symmetry& symmetry, const Potential& potential, 
 *		const Word& fileAppend)
 *
 * Calculate force constant matrix
 */

Word Phonons::generateForceConstants(const ISO& iso, const Symmetry& symmetry, const Potential& potential, \
	const Word& fileAppend)
{
	
	// Output
	Output::newline();
	Output::print("Generating force constants by displacing ");
	Output::print(symmetry.orbits().length());
	Output::print(" unique atom");
	if (symmetry.orbits().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Make space
	int i;
	clear();
	_forceConstants.size(3 * iso.numAtoms());
	_massFactors.size(iso.numAtoms());
	_vectors.length(iso.numAtoms());
	for (i = 0; i < _vectors.length(); ++i)
		_vectors[i].length(iso.numAtoms());
	
	// Variable to store pairs of atom types
	List<int>::D3 types;
	types.length(iso.numAtoms());
	for (i = 0; i < types.length(); ++i)
		types[i].length(iso.numAtoms());
	
	// Variable to store atom maps for each symmetry operation
	OList<List<Atom*> >::D2 atomMap;
	atomMap.length(symmetry.operations().length());
	for (i = 0; i < symmetry.operations().length(); ++i)
		atomMap[i].length(symmetry.operations()[i].translations().length());
	
	// Loop over unique atoms in the structure and build up force constant matrix
	int j, k, m, n;
	int atom1;
	int atom2;
	Atom* atom;
	List<Atom*>* curAtomMap;
	Vector curForceConstants;
	Matrix3D curMat;
	Matrix3D cartRot;
	Matrix3D cartRotTrans;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		
		// Output
		Output::newline();
		Output::print("Displacing atom ");
		Output::print(symmetry.orbits()[i].atoms()[0]->atomNumber() + 1);
		Output::print(" (");
		Output::print(symmetry.orbits()[i].atoms()[0]->element().symbol());
		Output::print(")");
		Output::increase();
	
		// Get the force constants
		atom = symmetry.orbits()[i].atoms()[0];
		for (j = 0; j < 3; ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Displacing along ");
			if (j == 0)
				Output::print("x");
			else if (j == 1)
				Output::print("y");
			else
				Output::print("z");
			Output::print(" direction");
			Output::increase();
			
			// Get force constant
			getForceConstants(curForceConstants, iso, symmetry, potential, atom, j);
			for (k = 0; k < curForceConstants.length(); ++k)
				_forceConstants(3*atom->atomNumber()+j, k) = curForceConstants[k];
			
			// Output
			Output::decrease();
		}
		
		// Output
		Output::newline();
		Output::print("Generating equivalent components of force constant matrix");
		Output::increase();
		
		// Loop over atoms that are equivalent to current
		for (j = 1; j < symmetry.orbits()[i].atoms().length(); ++j)
		{
			
			// Get current symmetry operation relative to cartesian basis
			cartRot  = iso.basis().inverseTranspose();
			cartRot *= symmetry.orbits()[i].generators()[j].rotation();
			cartRot *= iso.basis().vectorsTranspose();
			cartRotTrans = cartRot.transpose();
			
			// Get atom map for current symmetry operation
			curAtomMap = 0;
			for (k = 0; k < atomMap.length(); ++k)
			{
				if (symmetry.operations()[k].rotation() == symmetry.orbits()[i].generators()[j].rotation())
				{
					for (m = 0; m < atomMap[k].length(); ++m)
					{
						if (symmetry.operations()[k].translations()[m] == \
							symmetry.orbits()[i].generators()[j].translations()[0])
						{
							curAtomMap = &atomMap[k][m];
							break;
						}
					}
					break;
				}
			}

			// Did not find match
			if (!curAtomMap)
			{
				Output::newline(ERROR);
				Output::print("Internal error: could not match generator to symmetry operation");
				Output::quit();
			}
			
			// Get atom map if it has not been set
			if (curAtomMap->length() == 0)
				symmetry.getFullMap(*curAtomMap, iso, symmetry.orbits()[i].generators()[j].rotation(), \
					symmetry.orbits()[i].generators()[j].translations()[0]);
			
			// Loop over atoms
			atom1 = symmetry.orbits()[i].atoms()[j]->atomNumber();
			for (k = 0; k < curAtomMap->length(); ++k)
			{
				
				// Get current 3x3 for pair of atoms
				atom2 = (*curAtomMap)[k]->atomNumber();
				for (m = 0; m < 3; ++m)
				{
					for (n = 0; n < 3; ++n)
						curMat(m, n) = _forceConstants(3*atom->atomNumber()+m, 3*k+n);
				}
				
				// Convert under symmetry operation
				curMat  = curMat * cartRot;
				curMat *= cartRotTrans;
				
				// Save new 3x3 for equivalent pair
				for (m = 0; m < 3; ++m)
				{
					for (n = 0; n < 3; ++n)
						_forceConstants(3*atom1+m, 3*atom2+n) = curMat(m, n);
				}
			}
		}
		
		// Output
		Output::decrease();
		
		// Output
		Output::decrease();
	}
	
	// Make sure that force constant matrix is symmetric
	for (i = 0; i < _forceConstants.numRows(); ++i)
	{
		for (j = 0; j < i; ++j)
			_forceConstants(i, j) = _forceConstants(j, i) = (_forceConstants(i, j) + _forceConstants(j, i))/2;
	}
	
	for (i = 0; i < _forceConstants.numRows(); ++i)
		_forceConstants(i, i) += 1;
	
	// Set translation constraint
	double total;
	double absTotal;
	for (i = 0; i < _forceConstants.numRows(); ++i)
	{
		total = 0;
		absTotal = 0;
		for (j = 0; j < _forceConstants.numCols(); ++j)
		{
			total += _forceConstants(i, j);
			absTotal += Num<double>::abs(_forceConstants(i, j));
		}
		for (j = 0; j < _forceConstants.numCols(); ++j)
			_forceConstants(i, j) -= total * Num<double>::abs(_forceConstants(i, j)) / absTotal;
	}
	
	// Build up mass factor, vector, and type variables
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Loop over pairs of atoms and save values
			atom1 = iso.atoms()[i][j].atomNumber();
			for (k = 0; k < iso.atoms().length(); ++k)
			{
				for (m = 0; m < iso.atoms()[k].length(); ++m)
				{
					atom2 = iso.atoms()[k][m].atomNumber();
					types[atom1][atom2].length(2);
					types[atom1][atom2][0] = i;
					types[atom1][atom2][1] = k;
					_massFactors(atom1, atom2) = sqrt(iso.atoms()[i][j].element().mass() * \
						iso.atoms()[k][m].element().mass());
					iso.basis().distance(iso.atoms()[i][j].fractional(), FRACTIONAL, iso.atoms()[k][m].fractional(), \
						FRACTIONAL, &_vectors[atom1][atom2]);
					_vectors[atom1][atom2] += iso.atoms()[k][m].fractional();
					_vectors[atom1][atom2] -= iso.atoms()[i][j].fractional();
				}
			}
		}
	}
	
	// Print force constant file if needed
	Word fcFile;
	if (_writeFCFile)
	{
		
		// Save current output settings and open new stream
		fcFile = "fc_";
		fcFile += fileAppend;
		fcFile += ".dat";
		int origStreamID = Output::streamID();
		int newStreamID  = Output::addStream(fcFile);
		Output::setStream(newStreamID);

		// Print masses
		Output::newline();
		for (i = 0; i < iso.atoms().length(); ++i)
		{
			Output::print(iso.atoms()[i][0].element().mass());
			if (i != iso.atoms().length() - 1)
				Output::print(" ");
		}

		// Loop over pairs of atoms
		for (i = 0; i < types.length(); ++i)
		{
			for (j = 0; j < types[i].length(); ++j)
			{

				// Print types
				Output::newline();
				Output::print(types[i][j][0]);
				Output::print(" ");
				Output::print(types[i][j][1]);
				Output::print(" ");

				// Print vector components
				Output::print(_vectors[i][j], 8, false);

				// Print force constants
				for (k = 0; k < 3; ++k)
				{
					Output::newline();
					for (m = 0; m < 3; ++m)
					{
						Output::printSci(_forceConstants(3*i+k, 3*j+m), 8);
						if (m != 2)
							Output::print(" ");
					}
				}
			}
		}

		// Reset stream
		Output::setStream(origStreamID);
		Output::removeStream(newStreamID);
	}

	// Save that set
	_isSet = true;
	
	// Output
	Output::decrease();
	
	// Return name of force constant file
	return fcFile;
}



/* void Phonons::getForceConstants(Vector& constants, const ISO& iso, const Symmetry& symmetry,
 *		const Potential& potential, Atom* atom, int direction)
 *
 * Get the force constants
 */

void Phonons::getForceConstants(Vector& constants, const ISO& iso, const Symmetry& symmetry, \
	const Potential& potential, Atom* atom, int direction)
{
	
	// Make sure there is enough space for data
	constants.length(3*iso.numAtoms());
	
	// Save the original position
	Vector3D origPos = atom->cartesian();
	
	// Variables to store forces
	List<double> displacements;
	OList<OList<Vector3D > > forces;
	
	// Loop over displacements
	int i;
	Vector3D curPos;
//	Symmetry newSymmetry;
	for (double disp = -0.1; disp < 0.11; disp += 0.05)
	{
		
		// Output
		Output::newline();
		Output::print("Calculating forces at displacement of ");
		Output::print(disp);
		Output::print(" ang");
		Output::increase();
		
		// Set current position
		curPos = origPos;
		curPos[direction] += disp;
		atom->cartesian(curPos);
//		if (disp == -0.1)
//			newSymmetry.set(iso, 1e-5);
		
		// Calculate forces on atom
		forces.add();
		displacements += disp;
		potential.single(iso, /*newSymmetry, */0, &forces.last(), true, false);
		
		// Convert forces to cartesian frame
		for (i = 0; i < forces.last().length(); ++i)
			iso.basis().toCartesian(forces.last()[i]);
		
		// Output
		Output::decrease();
	}
	
	// Initialize variable to store fitting data
	List<double>::D2 data(displacements.length());
	for (i = 0; i < displacements.length(); ++i)
	{
		data[i].length(2);
		data[i][0] = displacements[i];
	}
	
	// Calculate force constants
	int j, k;
	Vector polyCoeffs;
	for (i = 0; i < iso.numAtoms(); ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			for (k = 0; k < forces.length(); ++k)
				data[k][1] = forces[k][i][j];
			polyCoeffs = Fit::polynomial(data, 0, 4);
			constants[3*i+j] = -polyCoeffs[1];
		}
	}
	
	// Reset original position
	atom->cartesian(origPos);
}



/* CVector Phonons::frequencies(const Vector3D& qFrac, CMatrix* modes) const
 *
 * Get the squared frequencies for a set reciprocal lattice vector
 */

CVector Phonons::frequencies(const Vector3D& qFrac, CMatrix* modes) const
{
	
	// Output
	Output::newline();
	Output::print("Diagonalizing dynamical matrix at q = ");
	Output::print(qFrac, 8, true);
	Output::increase();
	
	// Get the q vector
	Vector3D q = qFrac * (2 * Constants::pi);
	
	// Build up the Hessian matrix, loop over pairs of atoms
	int i, j, k, m;
	double arg;
	Complex scale;
	CMatrix H(_forceConstants.numRows(), _forceConstants.numCols());
	for (i = 0; i < _massFactors.numRows(); ++i)
	{
		for (j = 0; j < _massFactors.numCols(); ++j)
		{
			
			// Save scaling for current pair
			arg = q * _vectors[i][j];
			scale.real = cos(arg)/_massFactors(i, j);
			scale.imag = sin(arg)/_massFactors(i, j);
			
			// Loop over directions and save matrix values
			for (k = 0; k < 3; ++k)
			{
				for (m = 0; m < 3; ++m)
					H(3*i+k, 3*j+m) = _forceConstants(3*i+k, 3*j+m) * scale;
			}
		}
	}
	
	// Get the squared frequencies
	CMatrix allModes;
	CVector squaredFreqs = H.eigenvalues(&allModes, true);
	for (i = 0; i < squaredFreqs.length(); ++i)
	{
		
		// Found a frequency with imaginary component
		if (Num<double>::abs(squaredFreqs[i].imag) > 1e-6)
		{
			Output::newline(ERROR);
			Output::print("Found a squared frequency with non-zero imaginary component (");
			Output::print(squaredFreqs[i].imag);
			Output::print(")");
			Output::quit();
		}
	}
	
	// Sort the squared frequencies
	sortModes(squaredFreqs, allModes, 0, squaredFreqs.length() - 1);
	moveAcousticToStart(squaredFreqs, allModes);
	
	// Get the frequencies
	int numImaginary = 0;
	CVector freqs(squaredFreqs.length());
	for (i = 0; i < freqs.length(); ++i)
	{
		
		// Save frequency if real
		if (squaredFreqs[i].real >= 0)
		{
			freqs[i].real = sqrt(squaredFreqs[i].real) / (2 * Constants::pi);
			freqs[i].imag = 0;
		}
		
		// Save imaginary frequency and index if needed
		else
		{
			freqs[i].imag = sqrt(-squaredFreqs[i].real) / (2 * Constants::pi);
			freqs[i].real = 0;
			if (i > 2)
				++numImaginary;
		}
	}
	
	// Print warnings for imaginary frequencies
	if (numImaginary)
	{
		Output::newline(WARNING);
		Output::print("Found ");
		Output::print(numImaginary);
		Output::print(" imaginary frequenc");
		if (numImaginary == 1)
			Output::print("y");
		else
			Output::print("ies");
	}
	
	// Output
	Output::decrease();
	
	// Return the frequencies
	if (modes)
		*modes = allModes;
	return freqs;
}



/* void Phonons::sortModes(CVector& freqs, CMatrix& modes, int left, int right)
 *
 * Sort frequencies
 */

void Phonons::sortModes(CVector& freqs, CMatrix& modes, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	double pivot = freqs[pivotIndex].real;
	
	// Move pivot to end
	Num<Complex>::swap(freqs[pivotIndex], freqs[right]);
	modes.swapColumns(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (freqs[i].real < pivot)
		{
			Num<Complex>::swap(freqs[i], freqs[newPivotIndex]);
			modes.swapColumns(i, newPivotIndex);
			++newPivotIndex;
		}
	}
	
	// Move pivot to final position
	Num<Complex>::swap(freqs[newPivotIndex], freqs[right]);
	modes.swapColumns(newPivotIndex, right);
	
	// Recursive calls to next sorts
	sortModes(freqs, modes, left, newPivotIndex - 1);
	sortModes(freqs, modes, newPivotIndex + 1, right);
}



/* void Phonons::moveAcousticToStart(CVector& freqs, CMatrix& modes)
 *
 * Make sure that acoustic modes are at beginning
 */

void Phonons::moveAcousticToStart(CVector& freqs, CMatrix& modes)
{
	
	// Loop over three modes to find which are acoustic
	int i;
	List<int> indices;
	for (i = 0; i < freqs.length(); ++i)
	{
		if (isAcoustic(modes, i))
			indices += i;
	}
	
	// Did not find the correct number of acoustic modes
	if (indices.length() != 3)
	{
		Output::newline(WARNING);
		Output::print("Found ");
		Output::print(indices.length());
		Output::print(" acoustic mode");
		if (indices.length() != 1)
			Output::print("s");
		Output::print(" when expecting exactly 3");
	}
	
	// Set modes in order of lowest to highest frequency
	sortModes(freqs, indices, 0, indices.length() - 1);
	
	// Move acoustic modes to beginning
	int j;
	for (i = 0; i < Num<int>::min(3, indices.length()); ++i)
	{
		for (j = indices[i]; j > i; --j)
		{
			Num<Complex>::swap(freqs[j], freqs[j-1]);
			modes.swapColumns(j, j-1);
		}
	}
}



/* bool Phonons::isAcoustic(CMatrix& modes, int index)
 *
 * Return whether mode is acoustic
 */

bool Phonons::isAcoustic(CMatrix& modes, int index)
{
	
	// Get the max value in the mode
	int i;
	double curMag;
	double tol = Num<double>::abs(modes(0, index).real);
	for (i = 1; i < modes.numRows(); ++i)
	{
		curMag = Num<double>::abs(modes(i, index).real);
		if (curMag > tol)
			tol = curMag;
	}
	
	// Set tolerance
	tol /= modes.numRows();
	
	// Loop over directions
	int j;
	double compSign;
	for (i = 0; i < 3; ++i)
	{
		
		// Loop until first non-zero entry
		for (j = i; j < modes.numRows(); j += 3)
		{
			if (Num<double>::abs(modes(j, index).real) > tol)
			{
				compSign = Num<double>::sign(modes(j, index).real, tol);
				break;
			}
		}
		
		// Loop over remaining entries and make sure signs are correct
		for (j += 3; j < modes.numRows(); j += 3)
		{
			if (Num<double>::abs(modes(j, index).real) < tol)
				continue;
			if (Num<double>::sign(modes(j, index).real, tol) != compSign)
				return false;
		}
	}
	
	// Return that an acoustic mode if at this point
	return true;
}



/* void Phonons::sortModes(CVector& freqs, List<int>& indices, int left, int right)
 *
 * Sort acoustic phonon modes
 */

void Phonons::sortModes(CVector& freqs, List<int>& indices, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	double pivot = freqs[indices[pivotIndex]].real;
	
	// Move pivot to end
	Num<int>::swap(indices[pivotIndex], indices[right]);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (freqs[indices[i]].real < pivot)
		{
			Num<int>::swap(indices[i], indices[newPivotIndex]);
			++newPivotIndex;
		}
	}
	
	// Move pivot to final position
	Num<int>::swap(indices[newPivotIndex], indices[right]);
	
	// Recursive calls to next sorts
	sortModes(freqs, indices, left, newPivotIndex - 1);
	sortModes(freqs, indices, newPivotIndex + 1, right);
}



/* void Phonons::set(const Text& content)
 * 
 * Set force constants from file
 */

void Phonons::set(const Text& content)
{
	
	// Output
	Output::newline();
	Output::print("Reading force constant data from file");
	Output::increase();
	
	// Make sure that file is correct length
	if ((content.length() - 1) % 4 != 0)
	{
		Output::newline(ERROR);
		Output::print("Force constant file length is not correct");
		Output::quit();
	}
	
	// Clear space
	clear();
	
	// Get number of atoms
	int size = (int)sqrt((double)((content.length() - 1) / 4));
	
	// Allocate space
	int i;
	_forceConstants.size(3*size);
	_massFactors.size(size);
	_vectors.length(size);
	for (i = 0; i < size; ++i)
		_vectors[i].length(size);
	
	// Get masses
	List<double> masses(content[0].length());
	for (i = 0; i < content[0].length(); ++i)
		masses[i] = atof(content[0][i].array());
	
	// Loop over data lines
	int j, k;
	int type1;
	int type2;
	int curRow = 0;
	int curCol = 0;
	for (i = 1; i < content.length(); i += 4)
	{
		
		// Line is too short
		if (content[i].length() < 5)
		{
			Output::newline(ERROR);
			Output::print("Not enough data on force constant file line");
			Output::quit();
		}
		
		// Get masses
		type1 = atoi(content[i][0].array());
		type2 = atoi(content[i][1].array());
		if ((type1 >= masses.length()) || (type2 >= masses.length()))
		{
			Output::newline(ERROR);
			Output::print("Requesting mass of atom that has not been defined");
			Output::quit();
		}
		_massFactors(curRow, curCol) = sqrt(masses[type1]*masses[type2]);
		
		// Set vector
		_vectors[curRow][curCol].set(atof(content[i][2].array()), atof(content[i][3].array()), \
			atof(content[i][4].array()));
		
		// Get force constant values
		for (j = 0; j < 3; ++j)
		{
			for (k = 0; k < 3; ++k)
				_forceConstants(3*curRow+j, 3*curCol+k) = atof(content[i+j+1][k].array());
		}
		
		// Set new row and column
		if (++curCol == size)
		{
			curCol = 0;
			++curRow;
		}
	}
	
	// Output
	Output::decrease();
	
	// Save that set
	_isSet = true;
}



/* bool Phonons::isForceConstantFile(const Text& content)
 * 
 * Return whether file contains force constant information
 */

bool Phonons::isForceConstantFile(const Text& content)
{
	
	// Return if empty
	if (content.length() == 0)
		return false;
	if (content[0].length() == 0)
		return false;
	
	// Make sure that file is correct length
	if ((content.length() - 1) % 4 != 0)
		return false;
	
	// Make sure that all values on first line are numbers
	int i;
	for (i = 0; i < content[0].length(); ++i)
	{
		if (!Language::isNumber(content[0][i]))
			return false;
	}
	
	// Loop over remaining lines
	for (i = 1; i < content.length(); i += 4)
	{
		if ((content[i].length() != 5) || (content[i+1].length() != 3) || (content[i+2].length() != 3) || \
			(content[i+3].length() != 3))
			return false;
	}
	
	// Return that file is correct format
	return true;
}
