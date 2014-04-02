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



#ifndef SPACEGROUP_H
#define SPACEGROUP_H



#include "num.h"
#include "iso.h"
#include "symmetry.h"
#include "pointGroup.h"
#include "text.h"
#include "list.h"



// Class for a Wyckoff position and its orbit
class Wyckoff
{
	
	// Position variables
	int _rank;
	OList<Matrix3D > _rotations;
	OList<Vector3D > _translations;
	
	// Identification variables
	char _letter[2];
	Word _name;
	
	// Static member variables
	static char _letters[35][2];
	
public:
	
	// Constructor
	Wyckoff()						{ _letter[0] = '\0'; };
	Wyckoff(const Wyckoff& copy)	{ *this = copy; }
	
	// General functions
	void clear();
	Wyckoff& operator= (const Wyckoff& rhs);
	
	// Setup functions
	void set(int rotIndex, int transIndex, const OList<SymmetryOperation>& symmetry, int letterIndex);
	
	// Print functions
	void print() const;
	
	// Access functions
	int rank() const								{ return _rank; }
	int multiplicity() const						{ return _rotations.length(); }
	const char* letter() const						{ return _letter; }
	const Word& name() const						{ return _name; }
	const OList<Matrix3D >& rotations() const		{ return _rotations; }
	const OList<Vector3D >& translations() const	{ return _translations; }
	Words getString(int index) const
		{ return JonesFaithful::toString(_rotations[index], &_translations[index]); }
};



// Class for a single space group
class SpaceGroup
{
	
	// Identification variables
	Word _itcNumber;
	Word _schoenflies;
	Word _hermannMauguinShort;
	Word _hermannMauguin;
	Word _hall;
	
	// Symmetry variables
	LatticeSystem _system;
	LatticeCentering _centering;
	OList<SymmetryOperation> _symmetry;
	Matrix3D _unitToConv;
	Vector3D _originShift;
	
	// Wyckoff data variables
	OList<Wyckoff> _wyckoff;
	
	// Point group variables
	PointGroup _pointGroup;
	
	// General functions
	void setByNumber(int number);
	void setName(const char* number, const char* schoe, const char* hmShort, const char* hmLong, const char* hall);
	void setSymmetry(int gen1, int vec1, int gen2, int vec2, int gen3, int vec3, int gen4, int vec4, \
		LatticeSystem system, LatticeCentering centering);
	void addSymmetry(const Matrix3D& rotation, const Vector3D& translation);
	void setWyckoff(int wr1, int wt1, int wr2, int wt2, int wr3, int wt3, int wr4, int wt4, int wr5, int wt5, \
		int wr6, int wt6, int wr7, int wt7, int wr8, int wt8, int wr9, int wt9, int wr10, int wt10, int wr11, \
		int wt11, int wr12, int wt12, int wr13, int wt13, int wr14, int wt14, int wr15, int wt15, int wr16, \
		int wt16, int wr17, int wt17, int wr18, int wt18, int wr19, int wt19, int wr20, int wt20, int wr21, \
		int wt21, int wr22, int wt22, int wr23, int wt23, int wr24, int wt24, int wr25, int wt25, int wr26, \
		int wt26, int wr27, int wt27);
	
	// Helper functions to get space group from name
	static int getNumberFromName(const Word& name, bool useNumber);
	static bool compareName(const Word& name, const char* number, const char* schoe, const char* hmShort, \
		const char* hmFull, const char* hall, bool useNumber);
	
	// Helper functions to get space group from structure
	int getGeneratorIndex(const Matrix3D& rotation);
	int generatorsToNumber(int numOps, const List<int>& gens, const OList<Vector3D >& trans, \
		const Matrix3D& conv);
	bool compGen(int numOps, const Matrix3D& conv, const List<int>& gens, \
		const OList<Vector3D >& translations, int numSGOps, int gr1, int gt1, int gr2, int gt2, int gr3, \
		int gt3, int gr4, int gt4, LatticeSystem system, LatticeCentering centering);
	Matrix3D generatorToMatrix(int index);
	Vector3D generatorToTranslation(int index);
	
public:
	
	// Constructor
	SpaceGroup()							{}
	SpaceGroup(const SpaceGroup& copy)		{ *this = copy; }
	SpaceGroup(const Word& name, bool useNumber = true, bool quitIfNotFound = true)
		{ set(name, useNumber, quitIfNotFound); }
	SpaceGroup(const ISO& iso, double tol)	{ set(iso, tol); }
	
	// General functions
	void clear();
	SpaceGroup& operator= (const SpaceGroup& rhs);
	
	// Setup functions
	bool set(const Word& name, bool useNumber = true, bool quitIfNotFound = true);
	void set(const ISO& iso, double tol);
	
	// Static member functions
	static bool isSpaceGroup(const Word& name, bool useNumber)	{ return (getNumberFromName(name, useNumber) != -1); }
	static Matrix3D conventionalTransformation(const ISO& iso, double tol, Vector3D* shift = 0);
	static void makeConventional(ISO& iso, double tol);
	
	// Print functions
	void print() const;
	void printComplete() const;
	static void printAll();
	
	// Access functions
	const Word& itcNumber() const						{ return _itcNumber; }
	const Word& schoenflies() const						{ return _schoenflies; }
	const Word& hermannMauguinShort() const				{ return _hermannMauguinShort; }
	const Word& hermannMauguin() const					{ return _hermannMauguin; }
	const Word& hall() const							{ return _hall; }
	LatticeSystem system() const						{ return _system; }
	LatticeCentering centering() const					{ return _centering; }
	const OList<SymmetryOperation>& symmetry() const	{ return _symmetry; }
	const Matrix3D& unitToConv() const					{ return _unitToConv; }
	const Vector3D& originShift() const					{ return _originShift; }
	const OList<Wyckoff>& wyckoff() const				{ return _wyckoff; }
	const PointGroup& pointGroup() const				{ return _pointGroup; }
};



#endif
