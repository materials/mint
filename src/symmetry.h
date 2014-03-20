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


#ifndef SYMMETRY_H
#define SYMMETRY_H



#include "num.h"
#include "iso.h"
#include "text.h"
#include "list.h"



// Namespace to convert between Jones-Faithful notation and rotation matrix + translation
namespace JonesFaithful
{
	Words toString(const Matrix3D& rotation, const Vector3D* translation = 0);
	void fromString(Matrix3D& rotation, Vector3D& translation, const Words& input);
}



// Class to store information about a single symmetry operation
class SymmetryOperation
{
    
    // Variables
    Matrix3D _rotation;
    OList<Vector3D> _translations;
    
public:
	
	// Constructors
	SymmetryOperation()								 {}
	SymmetryOperation(const SymmetryOperation& copy) { _rotation = copy._rotation; _translations = copy._translations; }
	SymmetryOperation(const Matrix3D& rotation, const Vector3D& translation)
		{ _rotation = rotation; _translations = translation; }

    // Setup functions
	void clear()											{ _translations.clear(); }
    void clearTranslations()								{ _translations.clear(); }
	void setRotation(const Matrix3D& inRotation)			{ _rotation = inRotation; }
	void addTranslation(const Vector3D& translation);
	SymmetryOperation& operator= (const SymmetryOperation& rhs);
	void setFromString(const Words& input)
		{ Vector3D temp; JonesFaithful::fromString(_rotation, temp, input); addTranslation(temp); }
	bool operator== (const SymmetryOperation& rhs) const
		{ return ((_rotation == rhs._rotation) && (_translations == rhs._translations)); }
    
    // Access functions
	const Matrix3D& rotation() const			{ return _rotation; }
	const OList<Vector3D>& translations() const	{ return _translations; }
    
    // Functions
	Words getString(int transNum = 0) const	{ return JonesFaithful::toString(_rotation, &_translations[transNum]); }
	Words getRotationString() const			{ return JonesFaithful::toString(_rotation); }
};



// Class to store information about a special position
class SpecialPosition
{
	
	// Variables
	Matrix3D _rotation;
	Vector3D _translation;
	OList<Matrix3D> _rotations;
	OList<Vector3D> _translations;
	
public:
	
	// Constructors
	SpecialPosition()								{}
	SpecialPosition(const SpecialPosition& copy)	{ *this = copy; }
	SpecialPosition(const Matrix3D& rotation, const Vector3D& translation);
	SpecialPosition(const Matrix3D& rotation, const Vector3D& translation, const OList<Matrix3D>& rotations, \
		const OList<Vector3D>& translations);

	// Setup functions
	void set(const Matrix3D& rotation, const Vector3D& translation, const OList<Matrix3D>& rotations, \
		const OList<Vector3D>& translations);
	void addOperation(const Matrix3D& rotation, const Vector3D& translation);
	SpecialPosition& operator= (const SpecialPosition& rhs);
	
	// Access functions
	const Matrix3D& rotation() const			{ return _rotation; }
	const Vector3D& translation() const			{ return _translation; }
	const OList<Matrix3D>& rotations() const	{ return _rotations; }
	const OList<Vector3D>& translations() const	{ return _translations; }
	Words getString() const						{ return JonesFaithful::toString(_rotation, &_translation); }
};



// Class to store Wyckoff orbit
class Orbit
{
	
	// Variables
	int _rank;
	List<Atom*> _atoms;
	OList<SymmetryOperation> _generators;
	OList<SpecialPosition> _specialPositions;
	
public:
	
	// Constructors
	Orbit() {}
	Orbit(const Orbit& rhs)
		{ _atoms = rhs._atoms; _generators = rhs._generators; _specialPositions = rhs._specialPositions; }
	
	// Setup functions
	void set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation);
	void set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation, \
		const OList<Matrix3D>& pointRotations, const OList<Vector3D>& pointTranslations);
	void add(Atom* atom, const Matrix3D& genRotation, const Vector3D& genTranslation);
	void set(int index, Atom* atom)	{ _atoms[index] = atom; }
	Orbit& operator= (const Orbit& rhs);
	
	// Access functions
	int rank() const										{ return _rank; }
	const List<Atom*>& atoms() const						{ return _atoms; }
	const OList<SymmetryOperation>& generators() const		{ return _generators; }
	const OList<SpecialPosition>& specialPositions() const	{ return _specialPositions; }
	bool anyAtomsFixed() const;
};



// Class to store information about symmetry of a structure
class Symmetry
{
	
	// Type definitions
	typedef Matrix3D Transformation;
	typedef Matrix3D Rotation;
	typedef Vector3D Translation;
	
	// Variables
	OList<SymmetryOperation> _operations;
	OList<Orbit> _orbits;
	List<int> _orbitNumbers;
	Matrix _metricMatrixConstraint;
	
	// Functions
	void setUnitOperations(const Linked<Rotation>& redRotations, Linked<Translation>& redTranslations, \
		const Transformation& unitToReduced);
	void setMetricMatrixConstraint();
	void setOrbits(const ISO& iso, double tol);

	// Initialization helper functions
	static Transformation unitToReduced(const ISO& iso, double tol, bool isReducedPrim);
	
	// Helper functions for determining symmetry operations
	static void testOperations(const ISO& iso, double tol, Linked<Rotation>& candidates, \
		Linked<Rotation>& redRotations, Linked<Translation>& redTranslations);
	static bool checkOperation(const Rotation& rotation, const ISO& iso, double tol, Translation& translation, \
		Atoms& rotAtoms, Atoms& transAtoms, List<double>::D2* distances);
	static void addOperation(Linked<Rotation>& candidates, const Rotation& rotation, const Translation& translation, \
		Linked<Rotation>& redRotations, Linked<Translation>& redTranslations, Linked<Rotation>& notAllowed);
	static void removeOperation(Linked<Rotation>& candidates, const Rotation& rotation, \
		Linked<Rotation>& redRotations, Linked<Rotation>& notAllowed);
	
	// Helper functions for refinement
	void refineUniquePositions(ISO& iso, OList<Atom>& uniqueAtoms, double clusterTol) const;
	void clusterAtoms(ISO& iso, const OList<Atom>& uniqueAtoms, double clusterTol) const;
	
	// Helper functions for ideal cell conversions
	static Matrix3D idealTransformation(const ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol);
	static Matrix3D makeIdeal(ISO& iso, int numAtoms, bool numAtomsAsMin, double minDis, double tol);

public:
	
	// Constructors
	Symmetry()								{ _metricMatrixConstraint.size(6); }
	Symmetry(const Symmetry& copy)			{ _metricMatrixConstraint.size(6); *this = copy; }
	Symmetry(const ISO& iso, double tol)	{ _metricMatrixConstraint.size(6); set(iso, tol); }
	
	// Clear data
	void clear();
	void setToP1(const ISO& iso);
	
	// Setup functions
	Symmetry& operator= (const Symmetry& rhs);
	void set(const ISO& iso, double tol, bool isReducedPrim = false);
	
	// Manual setup functions
	void operations(const OList<SymmetryOperation>& input)	{ _operations = input; setMetricMatrixConstraint(); }
	void addOrbit(const Orbit& input)						{ _orbits += input; }
	
	// Access functions
	const OList<SymmetryOperation>& operations() const	{ return _operations; }
	const OList<Orbit>& orbits() const					{ return _orbits; }
	Orbit& orbit(int index)								{ return _orbits[index]; }
	const List<int>& orbitNumbers() const				{ return _orbitNumbers; }
	void print(bool useJonesFaithful = true) const;
	void printSites() const;
	
	// Other functions
	void getFullMap(List<Atom*>& map, const ISO& iso, const Matrix3D& rotation, const Vector3D& translation) const;
	
	// Refine functions
	void refine(ISO& iso, double clusterTol) const;
	void refineBasis(ISO& iso) const;
	void refineBasis(Matrix3D& vectors) const;
	void refineAtoms(ISO& iso, double clusterTol) const;
	
	// Conversions to ideal cell
	static Matrix3D idealTransformation(const ISO& iso, int numAtoms, bool numAtomsAsMin, double tol = 1e-4);
	static Matrix3D idealTransformation(const ISO& iso, double minDis, double tol = 1e-4);
	static Matrix3D makeIdeal(ISO& iso, int numAtoms, bool numAtomsAsMin, double tol = 1e-4);
	static Matrix3D makeIdeal(ISO& iso, double minDis, double tol = 1e-4);
	
	// Static member functions
	static int rotationOrder(const Matrix3D& rotation);
	static Vector3D intrinsicTranslation(const Matrix3D& rotation, const Vector3D& translation);
	static int getUniqueAxis(const Matrix3D& vectors, LatticeSystem system);
	static int getConventionalCellUniqueSymmetryAxis(const OList<Matrix3D>& rotations, LatticeSystem system);
	static void confineBasis(Vector3D& lengths, Vector3D& angles, LatticeSystem latticeSystem);
	static void addAtom(ISO& iso, const Atom& atom, const OList<SymmetryOperation>& operations, double clusterTol);
};



// =====================================================================================================================
// SpecialPosition
// =====================================================================================================================

/* inline SpecialPosition::SpecialPosition(const Matrix3D& rotation, const Vector3D& translation)
 *
 * Constructor for SpecialPosition
 */

inline SpecialPosition::SpecialPosition(const Matrix3D& rotation, const Vector3D& translation)
{
	set(rotation, translation, Matrix3D::identity(), Vector3D(0.0));
}



/* inline SpecialPosition::SpecialPosition(const Matrix3D& rotation, const Vector3D& translation,
 *		const OList<Matrix3D>& rotations, const OList<Vector3D>& translations)
 *
 * Constructor for SpecialPosition
 */

inline SpecialPosition::SpecialPosition(const Matrix3D& rotation, const Vector3D& translation, \
	const OList<Matrix3D>& rotations, const OList<Vector3D>& translations)
{
	set(rotation, translation, rotations, translations);
}



/* inline void SpecialPosition::set(const Matrix3D& rotation, const Vector3D& translation,
 *		const OList<Matrix3D>& rotations, const OList<Vector3D>& translations)
 *
 * Set values of special position
 */

inline void SpecialPosition::set(const Matrix3D& rotation, const Vector3D& translation, \
	const OList<Matrix3D>& rotations, const OList<Vector3D>& translations)
{
	
	// Save values
	_rotation = rotation;
	_translation = translation;
	_rotations = rotations;
	_translations = translations;
	
	// Make sure translation are within cell
	ISO::moveIntoCell(_translation);
	for (int i = 0; i < _translations.length(); ++i)
		ISO::moveIntoCell(_translations[i]);
}



/* inline void SpecialPosition::addOperation(const Matrix3D& rotation, const Vector3D& translation)
 *
 * Add operation to list for special position
 */

inline void SpecialPosition::addOperation(const Matrix3D& rotation, const Vector3D& translation)
{
	_rotations += rotation;
	_translations += translation;
}



/* inline SpecialPosition& SpecialPosition::operator= (const SpecialPosition& rhs)
 *
 * Assignment operator for SpecialPosition object
 */

inline SpecialPosition& SpecialPosition::operator= (const SpecialPosition& rhs)
{
	if (this != &rhs)
		set(rhs._rotation, rhs._translation, rhs._rotations, rhs._translations);
	return *this;
}



// =====================================================================================================================
// Orbit
// =====================================================================================================================

/* inline void Orbit::set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation)
 *
 * Set first position in orbit
 */

inline void Orbit::set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation)
{
	OList<Matrix3D> mats (Matrix3D::identity());
	OList<Vector3D> vecs (Vector3D(0.0));
	set(atom, specialRotation, specialTranslation, mats, vecs);
}



/* inline void Orbit::set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation,
 *		const OList<Matrix3D>& pointRotations, const OList<Vector3D>& pointTranslations)
 *
 * Set first position in orbit
 */

inline void Orbit::set(Atom* atom, const Matrix3D& specialRotation, const Vector3D& specialTranslation, \
	const OList<Matrix3D>& pointRotations, const OList<Vector3D>& pointTranslations)
{
	_atoms = atom;
	_generators = SymmetryOperation(Matrix3D::identity(), Vector3D(0.0));
	_specialPositions = SpecialPosition(specialRotation, specialTranslation, pointRotations, pointTranslations);
	_rank = _specialPositions.last().rotation().rank();
}



/* inline void Orbit::add(Atom* atom, const Matrix3D& genRotation, const Vector3D& genTranslation)
 *
 * Add position to orbit
 */

inline void Orbit::add(Atom* atom, const Matrix3D& genRotation, const Vector3D& genTranslation)
{
	_atoms += atom;
	_generators += SymmetryOperation(genRotation, genTranslation);
	_specialPositions += SpecialPosition(genRotation * _specialPositions[0].rotation(), \
		genRotation * _specialPositions[0].translation() + genTranslation);
}



/* inline Orbit& Orbit::operator= (const Orbit& rhs)
 *
 * Assignment operator for Orbit object
 */

inline Orbit& Orbit::operator= (const Orbit& rhs)
{
	if (this != &rhs)
	{
		_rank = rhs._rank;
		_atoms = rhs._atoms;
		_generators = rhs._generators;
		_specialPositions = rhs._specialPositions;
	}
	return *this;
}



/* inline bool Orbit::anyAtomsFixed() const
 *
 * Return whether any atoms in the orbit have a fixed position
 */

inline bool Orbit::anyAtomsFixed() const
{
	for (int i = 0; i < _atoms.length(); ++i)
	{
		if ((_atoms[i]->fixed()[0]) || (_atoms[i]->fixed()[1]) || (_atoms[i]->fixed()[2]))
			return true;
	}
	return false;
}



// =====================================================================================================================
// Symmetry
// =====================================================================================================================

/* inline Matrix3D Symmetry::idealTransformation(int numAtoms, bool numAtomsAsMin, double tol) const
 *
 * Return transformation to ideal cell
 */

inline Matrix3D Symmetry::idealTransformation(const ISO& iso, int numAtoms, bool numAtomsAsMin, double tol)
{
	return idealTransformation(iso, numAtoms, numAtomsAsMin, -1, tol);
}

inline Matrix3D Symmetry::idealTransformation(const ISO& iso, double minDis, double tol)
{
	return idealTransformation(iso, 0, false, minDis, tol);
}



/* inline Matrix3D Symmetry::makeIdeal(int numAtoms, bool numAtomsAsMin, double tol)
 *
 * Convert cell to ideal form
 */

inline Matrix3D Symmetry::makeIdeal(ISO& iso, int numAtoms, bool numAtomsAsMin, double tol)
{
	return makeIdeal(iso, numAtoms, numAtomsAsMin, -1, tol);
}

inline Matrix3D Symmetry::makeIdeal(ISO& iso, double minDis, double tol)
{
	return makeIdeal(iso, 0, false, minDis, tol);
}



#endif
