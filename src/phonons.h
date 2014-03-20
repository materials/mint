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



#ifndef PHONONS_H
#define PHONONS_H



#include "iso.h"
#include "symmetry.h"
#include "potential.h"
#include "text.h"
#include "fileSystem.h"
#include "num.h"
#include "list.h"



// Class for calculating phonons
class Phonons
{
	
	// Variables
	bool _isSet;
	bool _writeFCFile;
	Matrix _forceConstants;
	Matrix _massFactors;
	OList<Vector3D>::D2 _vectors;
	
	// Functions
	static void getForceConstants(Vector& constants, const ISO& iso, const Symmetry& symmetry, \
		const Potential& potential, Atom* atom, int direction);
	static void sortModes(CVector& freqs, CMatrix& modes, int left, int right);
	static void moveAcousticToStart(CVector& freqs, CMatrix& modes);
	static bool isAcoustic(CMatrix& modes, int index);
	static void sortModes(CVector& freqs, List<int>& indices, int left, int right);
	
public:
	
	// Constructor
	Phonons() { _isSet = false; _writeFCFile = true; }
	
	// Assignment
	Phonons& operator= (const Phonons& rhs);
	
	// Setup functions
	void clear();
	Word generateForceConstants(const ISO& iso, const Potential& potential, const Word& fileAppend);
	Word generateForceConstants(const ISO& iso, const Symmetry& symmetry, const Potential& potential, \
		const Word& fileAppend);
	void set(const Text& content);
	
	// Set settings
	void writeForceConstantsFile(bool input)	{ _writeFCFile = input; }
	
	// Get modes at given reciprocal lattice vector
	CVector frequencies(const Vector3D& qFrac, CMatrix* modes = 0) const;
	
	// Access functions
	bool isSet() const						{ return _isSet; }
	const Matrix& forceConstants() const	{ return _forceConstants; }
	
	// Other functions
	static bool isForceConstantFile(const Text& content);
	static bool isForceConstantFile(const Word& file)		{ return isForceConstantFile(Read::text(file)); }
};



/* inline void Phonons::clear()
 *
 * Clear data in phonons object
 */

inline void Phonons::clear()
{
	_isSet = false;
	_forceConstants.clear();
	_massFactors.clear();
	_vectors.clear();
}



/* inline Phonons& Phonons::operator= (const Phonons& rhs)
 *
 * Assignment operator for phonons
 */

inline Phonons& Phonons::operator= (const Phonons& rhs)
{
	if (this != &rhs)
	{
		clear();
		_isSet = rhs._isSet;
		_forceConstants = rhs._forceConstants;
		_massFactors = rhs._massFactors;
		_vectors = rhs._vectors;
	}
	return *this;
}



/* inline Word Phonons::generateForceConstants(const ISO& iso, const Potential& potential, const Word& fileAppend)
 *
 * Evaluate force constants without symmetry
 */

inline Word Phonons::generateForceConstants(const ISO& iso, const Potential& potential, const Word& fileAppend)
{
	Symmetry symmetry;
	symmetry.setToP1(iso);
	return generateForceConstants(iso, symmetry, potential, fileAppend);
}



#endif
