/* iso.h -- Internal structure object
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef ISO_H
#define ISO_H



#include "num.h"
#include "elements.h"
#include "text.h"
#include "list.h"



// Define fractional and cartesian coordinates
enum CoordinateType {FRACTIONAL, CARTESIAN};

// Crystal systems
enum CrystalSystem {CS_UNKNOWN, CS_TRICLINIC, CS_MONOCLINIC, CS_ORTHORHOMBIC, CS_TETRAGONAL, CS_TRIGONAL, \
	CS_HEXAGONAL, CS_CUBIC};

// Lattice systems
enum LatticeSystem {LS_UNKNOWN, LS_TRICLINIC, LS_MONOCLINIC, LS_ORTHORHOMBIC, LS_HEXAGONAL, LS_TETRAGONAL, \
	LS_RHOMBOHEDRAL, LS_CUBIC};

// Centering types
enum CenteringType {CT_UNKNOWN, CT_PRIMITIVE, CT_ONE_FACE, CT_BODY, CT_RHOMBOHEDRAL, CT_ALL_FACE};

// Lattice centerings
enum LatticeCentering {LC_UNKNOWN, LC_P, LC_H, LC_C, LC_A, LC_B, LC_I, LC_R, LC_F};



// Class containing information about basis vectors
class Basis
{
	
	// Variable to store lattice system if generating randomly
	mutable LatticeSystem _latticeSystem;
	
	// Variables
	double _volume;
	Matrix3D _vectors;
	Matrix3D _inverse;
	Matrix3D _metric;
	Matrix3D _vectorsTranspose;
	Matrix3D _inverseTranspose;
	Vector3D _lengths;
	Vector3D _angles;
	mutable bool _lengthFixed[3];
	mutable bool _angleFixed[3];
	
	// Reduced cell variables
	Matrix3D _reduced;
	Matrix3D _reducedTranspose;
	Matrix3D _reducedInverse;
	Matrix3D _reducedMetric;
	Matrix3D _unitToReduced;
	Matrix3D _reducedPointToUnit;
	Matrix3D _unitPointToReduced;
	
	// Helper variables
	mutable bool _first;
	mutable bool _second;
	mutable int _i, _j, _k;
	mutable int _numCellsToSearch[3];
	mutable double _curDis;
	mutable double _minDis;
	mutable double _secondDis;
	mutable double _minCell[3];
	mutable double _secondCell[3];
	mutable double _cellsToSearch[2][3];
	mutable Vector3D _frac1;
	mutable Vector3D _frac2;
	mutable Vector3D _redDif;
	mutable Vector3D _curVec;
	
	// Helper functions
	void init();
	void setVectors();
	void setMetrics();
	void finishSetup(bool showOutput);
	
	// Reduction functions
	static void NiggliParams(const Matrix3D& vectors, double& A, double& B, double& C, double& ksi, \
		double& eta, double& zeta);
	static void updateTransformation(Matrix3D& transformation, int m00, int m01, int m02, int m10, \
		int m11, int m12, int m20, int m21, int m22);
	
	// Helper function for lattice symmetry
	static bool addPossibleRotation(Linked<Matrix3D>& rotations, const Matrix3D& rotationToAdd);
	
public:
	
	// Constructors
	Basis();
	Basis(const Basis& copy)										{ init(); *this = copy; }
	Basis(const Matrix3D& vectors, bool showOutput = true)			{ init(); set(vectors, showOutput); }
	Basis(const Vector3D& lengths, const Vector3D& angles, bool showOutput = true)
		{ init(); set(lengths, angles, showOutput); }
	
	// Setup
	Basis& operator= (const Basis& rhs);
	void set(const Matrix3D& vectors, bool showOutput = true);
	void set(const Vector3D& lengths, const Vector3D& angles, bool showOutput = true);
	void latticeSystem(LatticeSystem input) const	{ _latticeSystem = input; }
	void lengthFixed(int index, bool value) const	{ _lengthFixed[index] = value; }
	void angleFixed(int index, bool value) const	{ _angleFixed[index] = value; }
	void lengthFixed(const bool* input) const
		{ _lengthFixed[0] = input[0]; _lengthFixed[1] = input[1]; _lengthFixed[2] = input[2]; }
	void angleFixed(const bool* input) const
		{ _angleFixed[0] = input[0]; _angleFixed[1] = input[1]; _angleFixed[2] = input[2]; }
	
	// Conversion functions
	void toFractional(Vector3D& input) const 					{ input *= _inverseTranspose; }
	void toCartesian(Vector3D& input) const  					{ input *= _vectorsTranspose; }
	Vector3D getFractional(const Vector3D& input) const			{ return _inverseTranspose * input; }
	Vector3D getCartesian(const Vector3D& input) const			{ return _vectorsTranspose * input; }
	
	// Distance functions
    double distance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
		CoordinateType type2, Vector3D* cell = 0) const;
	double secondDistance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
		CoordinateType type2, Vector3D* cell = 0) const;
    double absoluteDistance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
		CoordinateType type2) const;
	
	// Access functions
	double volume() const							{ return _volume; }
	const Matrix3D& vectors() const					{ return _vectors; }
	const Matrix3D& inverse() const					{ return _inverse; }
	const Matrix3D& vectorsTranspose() const		{ return _vectorsTranspose; }
	const Matrix3D& inverseTranspose() const		{ return _inverseTranspose; }
	const Vector3D& lengths() const					{ return _lengths; }
	const Vector3D& angles() const					{ return _angles; }
	LatticeSystem latticeSystem() const				{ return _latticeSystem; }
	const bool* lengthFixed() const					{ return _lengthFixed; }
	const bool* angleFixed() const					{ return _angleFixed; }
	
	// Access functions for reduced cell
	const Matrix3D& reduced() const					{ return _reduced; }
	const Matrix3D& reducedTranspose() const		{ return _reducedTranspose; }
	const Matrix3D& reducedInverse() const			{ return _reducedInverse; }
	const Matrix3D& unitToReduced() const			{ return _unitToReduced; }
	const Matrix3D& reducedPointToUnit() const		{ return _reducedPointToUnit; }
	const Matrix3D& unitPointToReduced() const		{ return _unitPointToReduced; }
	
	// Static member functions
	static double getAngle(const Vector3D& vec1, const Vector3D& vec2);
	static Vector3D rotationAxis(const Matrix3D& matrix);
	static Vector3D backSolve(const Matrix3D& rotation, const Vector3D* fill = 0);
	static Matrix3D vectors(const Vector3D& lengths, const Vector3D& angles);
	static Vector3D lengths(const Matrix3D& vectors);
	static Vector3D angles(const Matrix3D& vectors);
	static Matrix3D reducedTransformation(const Matrix3D& vectors);
	static Matrix3D reduced(const Matrix3D& vectors, Matrix3D* transformation = 0);
	static void getPossibleRotations(Linked<Matrix3D>& rotations, const Matrix3D& vectors, double tol);
};



// Class containing information about a single atom
class Atom
{
	
	// General variables
	mutable bool _assigned;
	bool _fixed[3];
	bool _interstitial;
	double _occupancy;
	Vector3D _cartesian;
	Vector3D _fractional;
	Vector3D _magneticMoment;
	Element _element;
	Words _tags;
	
	// Variables connected to ISO object
	int _atomNumber;
	const Basis* _basis;
	
	// Functions
	void basis(const Basis* input)		{ _basis = input; }
	void element(const Element& input)	{ _element = input; }
	void atomNumber(int input)			{ _atomNumber = input; }
	
public:
	
	// Constructor
	Atom();
	Atom(const Atom& copy)				{ *this = copy; }
	
	// Setup functions
	void clear();
	Atom& operator= (const Atom& rhs);
	void assigned(bool value) const		{ _assigned = value; }
	
	// Set whether atom is fixed
	void fixed(bool input)				{ _fixed[0] = input;    _fixed[1] = input;    _fixed[2] = input; }
	void fixed(const bool* input)		{ _fixed[0] = input[0]; _fixed[1] = input[1]; _fixed[2] = input[2]; }
	void fixed(int index, bool input)	{ _fixed[index] = input; }
	
	// Set whether atom is an interstitial
	void isInterstitial(bool input)		{ _interstitial = input; }
	
	// Set the occupancy
	void occupancy(double input)		{ _occupancy = input; }
	
	// Set the position
	void fractional(const double* input, bool moveIntoCell = true);
	void fractional(const Vector3D& input, bool moveIntoCell = true);
	void fractional(double x, double y, double z, bool moveIntoCell = true);
	void cartesian(const double* input, bool moveIntoCell = true);
	void cartesian(const Vector3D& input, bool moveIntoCell = true);
	void cartesian(double x, double y, double z, bool moveIntoCell = true);
	
	// Set the magnetic moment
	void magneticMoment(const double* input)			{ _magneticMoment = input; }
	void magneticMoment(const Vector3D& input)			{ _magneticMoment = input; }
	void magneticMoment(double input)					{ _magneticMoment[2] = input; }
	void magneticMoment(int index, double input)		{ _magneticMoment[index] = input; }
	
	// Add tags
	void addTag(const char* word)						{ _tags += word; }
	void addTag(const Word& word)						{ _tags += word; }
	
	// Functions
	bool equal(const Atom& rhs, double tol, double* distance = 0, Vector3D* cell = 0) const;
	
	// Access functions
	bool assigned() const							{ return _assigned; }
	const bool* fixed() const						{ return _fixed; }
	bool anyFixed() const							{ return ((_fixed[0]) || (_fixed[1]) || (_fixed[2])); }
	bool isInterstitial() const						{ return _interstitial; }
	double occupancy() const						{ return _occupancy; }
	const Vector3D& cartesian() const				{ return _cartesian; }
	const Vector3D& fractional() const				{ return _fractional; }
	const Vector3D& magneticMoment() const			{ return _magneticMoment; }
	const Element& element() const					{ return _element; }
	const Words& tags() const						{ return _tags; }
	int atomNumber() const							{ return _atomNumber; }
	
	// Friends
	friend class ISO;
	friend class CIF;
	friend class MintStructure;
};



// Types used in ISO class
typedef OList<Atom>::D2 Atoms;
typedef OList<Vector3D > LatticePoints;



// Class containing information about a structure
class ISO
{
	
	// Comment
	Words _comment;
	
	// Basis
	Basis _basis;
	
	// Atoms
	int _numAtoms;
	Atoms _atoms;
	
	// Space group when generating structure randomly
	Word _spaceGroup;
	
	// Functions
	Atoms expand(const Matrix3D& transformation);
	Atoms contract(const Matrix3D& transformation, double tol);
	
	// Helper functions
	static void fillRotations(Linked<Matrix3D>& rotations, Linked<Matrix3D>& origRotations);
	static void addRotation(Linked<Matrix3D>& rotations, const Matrix3D& curRotation);
	
public:
	
	// Constructors
	ISO()					{ _numAtoms = 0; }
	ISO(const ISO& copy)	{ *this = copy; }
	
	// General functions
	void clear();
	ISO& operator= (const ISO& rhs);
	
	// Set the comment
	void comment(const Word& input)		{ _comment += input; }
	void comment(const Words& input)	{ _comment = input; }
	
	// Conver to different cell types
	Matrix3D makeReduced(bool showOutput = true);
	Matrix3D primitiveTransformation(double tol, bool showOutput = true) const;
	Matrix3D makePrimitive(double tol, bool showOutput = true);
	
	// Change the basis
	void basis(const Basis& basis, bool showOutput = true);
	void basis(const Matrix3D& vectors, bool showOutput = true);
	void basis(const Vector3D& lengths, const Vector3D& angles, bool showOutput = true);
	void transform(const Matrix3D& transformation, double tol, bool showOutput = true);
	void rotatePositions(const Matrix3D& rotation, bool showOutput = true);
	void shift(const Vector3D& vector, bool showOutput = true);
	
	// Add, remove, or change atoms
	Atom* addAtom(const Element& element);
	Atom* addAtom(const Word& element)	{ return addAtom(Element::find(element, true)); }
	Atom* addAtom(const Atom& atom);
	void removeAtom(int index, bool showOutput = true);
	void setElement(int index, const Word& element, bool showOutput = true)
		{ setElement(index, Element::find(element, true), showOutput); }
	void setElement(int index, const Element& element, bool showOutput = true);
	void moveAtomToFirst(const Atom* atom);
	void moveAtomToLast(const Atom* atom);
	void setAtomInElement(int atomNumber, int index);
	void clearAtoms()	{ _atoms.clear(); _numAtoms = 0; }
	void orderAtomNumbers();
	
	// Get information about the structure
	List<Atom*> coordination(const Atom* atom, double fracTol = 0.1) const;
	OList<Atom>::D2 shells(const Atom* atom, double maxDistance, double tol = 1e-4) const;
	
	// Compare to another structure
	bool equivalent(const ISO& compISO, double tol, bool matchVolume = false) const;
	
	// Set the space group
	void spaceGroup(const Word& input)	{ _spaceGroup = input; }
	
	// Access functions
	const Words& comment() const	{ return _comment; }
	const Basis& basis() const		{ return _basis; }
	const Atoms& atoms() const		{ return _atoms; }
	const Word& spaceGroup() const	{ return _spaceGroup; }
	int numAtoms() const			{ return _numAtoms; }
	bool anyFixed() const;
	bool anyUnset() const;
	bool anyPartiallyOccupied() const;
	Atom* atom(int index) const;
	Atom* atom(const Element& element, int index) const;
	Atom* atom(const Vector3D& fracPos, double tol) const;
	List<Atom*> atoms(const Element& element) const;
	
	// Static member functions
	static void moveIntoCell(double& coordinate)
		{ coordinate -= Num<double>::floor(coordinate); }
	static void moveIntoCell(Vector3D& coordinates)
		{ for (int i = 0; i < 3; ++i) moveIntoCell(coordinates[i]); }
	static bool areSitesEqual(const Basis& basis, const Atoms& origAtoms, const Atoms& newAtoms, double tol, \
		Vector3D* vector = 0, List<double>::D2* origDistances = 0);
	static LatticePoints getLatticePoints(const Matrix3D& transformationToNewCell);
	static Word system(LatticeSystem input);
	static Word centering(LatticeCentering input);
	static LatticeSystem system(const Word& input);
};



// Class to store image iterator
class ImageIterator
{
	
	// Variables
	double _xRange[2];
	double _yRange[2];
	double _zRange[2];
	double _maxDistanceSquared;
	Basis _redBasis;
	Basis _unitBasis;
	
	// Variables to store current state
	bool _finished;
	double _x;
	double _y;
	double _z;
	double _distance;
	double _redMults[6];
	Vector3D _curVector;
	Vector3D _nearVector;
	
	// Functions
	double setRedMults();
	double curDistanceSquared();
	
public:
	
	// Constructor
	ImageIterator()											{}
	ImageIterator(const Basis& basis, double maxDistance)	{ setCell(basis, maxDistance); }
	
	// Functions
	void setCell(const Basis& basis, double maxDistance);
	void reset(const Vector3D& site1, const Vector3D& site2);
	double operator++ ();
	double operator++ (int)	{ return ++(*this); }
	
	// Access functions
	bool finished() const		{ return _finished; }
	double distance() const		{ return _distance; }
	Vector3D fracVector() const	{ return _unitBasis.reducedPointToUnit()*_curVector; }
	Vector3D cartVector() const	{ return _redBasis.getCartesian(_curVector); }
	Vector3D cellVector() const	{ return _unitBasis.reducedPointToUnit()*Vector3D(_x,_y,_z); }
};



// =====================================================================================================================
// Basis
// =====================================================================================================================

/* inline void Basis::init()
 *
 * Common setup
 */

inline void Basis::init()
{
	for (int i = 0; i < 3; ++i)
	{
		_lengthFixed[i] = false;
		_angleFixed[i] = false;
	}
}



/* inline Basis::Basis()
 *
 * Constructor for Basis
 */

inline Basis::Basis()
{
	init();
	_latticeSystem = LS_TRICLINIC;
	_volume = 0;
	_vectors = 0.0;
	_inverse = 0.0;
	_metric = 0.0;
	_vectorsTranspose = 0.0;
	_inverseTranspose = 0.0;
	_lengths = 0.0;
	_angles = 0.0;
	_reduced = 0.0;
	_reducedTranspose = 0.0;
	_reducedInverse = 0.0;
	_reducedMetric = 0.0;
	_unitToReduced = 0.0;
	_reducedPointToUnit = 0.0;
	_unitPointToReduced = 0.0;
}



/* inline double Basis::absoluteDistance(const Vector3D& pos1, CoordinateType type1,
 *		const Vector3D& pos2, CoordinateType type2) const
 *
 * Return the absolute distance between two positions
 */

inline double Basis::absoluteDistance(const Vector3D& pos1, CoordinateType type1, \
	const Vector3D& pos2, CoordinateType type2) const
{

	// Get fractional coordinates
	_frac1  = (type1 == FRACTIONAL) ? pos1 : getFractional(pos1);
	_curVec = (type2 == FRACTIONAL) ? pos2 : getFractional(pos2);
	
	// Return distance
	_curVec -= _frac1;
	return sqrt(_curVec * (_metric * _curVec));
}



// =====================================================================================================================
// Atom
// =====================================================================================================================

/* inline void Atom::fractional(const double* input, bool moveIntoCell)
 *
 * Set position by fractional coordinate
 */

inline void Atom::fractional(const double* input, bool moveIntoCell)
{
	_assigned = true;
	_fractional = input;
	if (moveIntoCell)
		ISO::moveIntoCell(_fractional);
	if (_basis)
		_cartesian = _basis->getCartesian(_fractional);
}



/* inline void Atom::fractional(const Vector3D& input, bool moveIntoCell)
 *
 * Set position by fractional coordinate
 */

inline void Atom::fractional(const Vector3D& input, bool moveIntoCell)
{
	_assigned = true;
	_fractional = input;
	if (moveIntoCell)
		ISO::moveIntoCell(_fractional);
	if (_basis)
		_cartesian = _basis->getCartesian(_fractional);
}



/* inline void Atom::fractional(double x, double y, double z, bool moveIntoCell)
 *
 * Set position by fractional coordinate
 */

inline void Atom::fractional(double x, double y, double z, bool moveIntoCell)
{
	fractional(Vector3D(x, y, z), moveIntoCell);
}



/* inline void Atom::cartesian(const double* input, bool moveIntoCell)
 *
 * Set position by cartesian coordinate
 */

inline void Atom::cartesian(const double* input, bool moveIntoCell)
{
	_assigned = true;
	_cartesian = input;
	if (_basis)
	{
		_fractional = _basis->getFractional(_cartesian);
		if (moveIntoCell)
			ISO::moveIntoCell(_fractional);
		_cartesian = _basis->getCartesian(_fractional);
	}
}



/* inline void Atom::cartesian(const Vector3D& input, bool moveIntoCell)
 *
 * Set position by cartesian coordinate
 */

inline void Atom::cartesian(const Vector3D& input, bool moveIntoCell)
{
	_assigned = true;
	_cartesian = input;
	if (_basis)
	{
		_fractional = _basis->getFractional(input);
		if (moveIntoCell)
			ISO::moveIntoCell(_fractional);
		_cartesian = _basis->getCartesian(_fractional);
	}
}



/* inline void Atom::cartesian(double x, double y, double z, bool moveIntoCell)
 *
 * Set position by cartesian coordinate
 */

inline void Atom::cartesian(double x, double y, double z, bool moveIntoCell)
{
	cartesian(Vector3D(x, y, z), moveIntoCell);
}



// =====================================================================================================================
// ISO
// =====================================================================================================================

/* inline void ISO::moveAtomToFirst(const Atom* atom)
 *
 * Move atom to beginning of element list
 */

inline void ISO::moveAtomToFirst(const Atom* atom)
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		if (_atoms[i][0].element() == atom->element())
		{
			for (j = 0; j < _atoms[i].length(); ++j)
			{
				if (&_atoms[i][j] == atom)
				{
					for (; j >= 1; --j)
						_atoms[i].swap(j, j-1);
					return;
				}
			}
		}
	}
}



/* inline void ISO::moveAtomToLast(const Atom* atom)
 *
 * Move atom to end of element list
 */

inline void ISO::moveAtomToLast(const Atom* atom)
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		if (_atoms[i][0].element() == atom->element())
		{
			for (j = 0; j < _atoms[i].length(); ++j)
			{
				if (&_atoms[i][j] == atom)
				{
					for (; j < _atoms[i].length() - 1; ++j)
						_atoms[i].swap(j, j+1);
					return;
				}
			}
		}
	}
}



/* inline void ISO::setAtomInElement(int atomNumber, int index)
 *
 * Move atom to set position in element list
 */

inline void ISO::setAtomInElement(int atomNumber, int index)
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].atomNumber() == atomNumber)
			{
				int step = Num<int>::sign(index - j);
				for (; j != index; j += step)
					_atoms[i].swap(j, j+step);
				return;
			}
		}
	}
}



/* inline void ISO::orderAtomNumbers()
 * 
 * Save atom numbers in current order
 */

inline void ISO::orderAtomNumbers()
{
	int i, j;
	int number = 0;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
			_atoms[i][j].atomNumber(number++);
	}
}



/* inline Atom* ISO::atom(int index) const
 *
 * Return an atom by number
 */

inline Atom* ISO::atom(int index) const
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].atomNumber() == index)
				return &_atoms[i][j];
		}
	}
	return 0;
}



/* inline Atom* ISO::atom(const Element& element, int index) const
 *
 * Return an atom by element and number
 */

inline Atom* ISO::atom(const Element& element, int index) const
{
	for (int i = 0; i < _atoms.length(); ++i)
	{
		if (_atoms[i][0].element() == element)
		{
			if (index < _atoms[i].length())
				return &_atoms[i][index];
			else
				return 0;
		}
	}
	return 0;
}



/* inline Atom* ISO::atom(const Vector3D& fracPos, double tol) const
 *
 * Return an atom by position
 */

inline Atom* ISO::atom(const Vector3D& fracPos, double tol) const
{
	
	// Return if there are no atoms
	if (_atoms.length() == 0)
		return 0;
	
	// Initialize first atom as nearest
	double minIndex[2] = {0, 0};
	double minDis = _basis.distance(fracPos, FRACTIONAL, _atoms[0][0].fractional(), FRACTIONAL);
	
	// Loop over atoms
	int i, j;
	double curDis;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = (i == 0) ? 1 : 0; j < _atoms[i].length(); ++j)
		{
			curDis = _basis.distance(fracPos, FRACTIONAL, _atoms[i][j].fractional(), FRACTIONAL);
			if (curDis < minDis)
			{
				minDis = curDis;
				minIndex[0] = i;
				minIndex[1] = j;
			}
		}
	}
	
	// Nearest atom is outside tolerance
	if (minDis > tol)
		return 0;
	
	// Return nearest atom
	return &_atoms[minIndex[0]][minIndex[1]];
}



/* inline List<Atom*> ISO::atoms(const Element& element) const
 *
 * Return all atoms of a set element
 */

inline List<Atom*> ISO::atoms(const Element& element) const
{
	int i, j;
	List<Atom*> res;
	for (i = 0; i < _atoms.length(); ++i)
	{
		if (_atoms[i][0].element() == element)
		{
			res.length(_atoms[i].length());
			for (j = 0; j < _atoms[i].length(); ++j)
				res[j] = &_atoms[i][j];
			break;
		}
	}
	return res;
}



/* inline bool ISO::anyFixed() const
 *
 * Check if any coordinates are fixed
 */

inline bool ISO::anyFixed() const
{
	int i, j, k;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			for (k = 0; k < 3; ++k)
			{
				if (_atoms[i][j].fixed()[k])
					return true;
			}
		}
	}
	return false;
}



/* inline bool ISO::anyUnset() const
 *
 * Return whether any of the coordinates have not been set
 */

inline bool ISO::anyUnset() const
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (!_atoms[i][j].assigned())
				return true;
		}
	}
	return false;
}



/* inline bool ISO::anyPartiallyOccupied() const
 *
 * Return whether any of the atoms in the structure are not fully occupied
 */

inline bool ISO::anyPartiallyOccupied() const
{
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].occupancy() < 1)
				return true;
		}
	}
	return false;
}



// =====================================================================================================================
// Image iterator
// =====================================================================================================================

/* inline void ImageIterator::sesetCelltup(const Basis& basis, double maxDistance)
 *
 * Setup image iterator
 */

inline void ImageIterator::setCell(const Basis& basis, double maxDistance)
{
	
	// Save basis
	_unitBasis = basis;
	_redBasis.set(basis.reduced(), false);
	
	// Save max distance
	_maxDistanceSquared = maxDistance*maxDistance;
	
	// Set multipliers for going from reduced cell vector to cartesian length squared
	_redMults[0] =	_redBasis.vectors()(0, 0)*_redBasis.vectors()(0, 0) + \
					_redBasis.vectors()(0, 1)*_redBasis.vectors()(0, 1) + \
				 	_redBasis.vectors()(0, 2)*_redBasis.vectors()(0, 2);
	_redMults[1] =	_redBasis.vectors()(1, 0)*_redBasis.vectors()(1, 0) + \
					_redBasis.vectors()(1, 1)*_redBasis.vectors()(1, 1) + \
					_redBasis.vectors()(1, 2)*_redBasis.vectors()(1, 2);
	_redMults[2] =	_redBasis.vectors()(2, 0)*_redBasis.vectors()(2, 0) + \
					_redBasis.vectors()(2, 1)*_redBasis.vectors()(2, 1) + \
					_redBasis.vectors()(2, 2)*_redBasis.vectors()(2, 2);
	_redMults[3] = (_redBasis.vectors()(0, 0)*_redBasis.vectors()(1, 0) + \
					_redBasis.vectors()(0, 1)*_redBasis.vectors()(1, 1) + \
					_redBasis.vectors()(0, 2)*_redBasis.vectors()(1, 2)) * 2;
	_redMults[4] = (_redBasis.vectors()(0, 0)*_redBasis.vectors()(2, 0) + \
					_redBasis.vectors()(0, 1)*_redBasis.vectors()(2, 1) + \
					_redBasis.vectors()(0, 2)*_redBasis.vectors()(2, 2)) * 2;
	_redMults[5] = (_redBasis.vectors()(1, 0)*_redBasis.vectors()(2, 0) + \
					_redBasis.vectors()(1, 1)*_redBasis.vectors()(2, 1) + \
					_redBasis.vectors()(1, 2)*_redBasis.vectors()(2, 2)) * 2;
}



/* inline void ImageIterator::initialize(const Vector3D& site1, const Vector3D& site2)
 *
 * Initialize image iterator
 */

inline void ImageIterator::reset(const Vector3D& site1, const Vector3D& site2)
{
	
	// Get points in reduced cell
	Vector3D redSite1 = _unitBasis.unitPointToReduced() * site1;
	Vector3D redSite2 = _unitBasis.unitPointToReduced() * site2;
	ISO::moveIntoCell(redSite1);
	ISO::moveIntoCell(redSite2);
	
	// Save the nearest point
	_redBasis.distance(redSite1, FRACTIONAL, redSite2, FRACTIONAL, &_nearVector);
	_nearVector += redSite2;
	_nearVector -= redSite1;
	
	// Set range for cell iterations
	_xRange[0] = Num<double>::floor(-sqrt(_maxDistanceSquared / _redMults[0]) - _nearVector[0]);
	_xRange[1] = Num<double>::ceil(sqrt(_maxDistanceSquared / _redMults[0]) - _nearVector[0]);
	_yRange[0] = Num<double>::floor(-sqrt(_maxDistanceSquared / _redMults[1]) - _nearVector[1]);
	_yRange[1] = Num<double>::ceil(sqrt(_maxDistanceSquared / _redMults[1]) - _nearVector[1]);
	_zRange[0] = Num<double>::floor(-sqrt(_maxDistanceSquared / _redMults[2]) - _nearVector[2]);
	_zRange[1] = Num<double>::ceil(sqrt(_maxDistanceSquared / _redMults[2]) - _nearVector[2]);
	
	// Set current cell
	_x = _xRange[0];
	_y = _yRange[0];
	_z = _zRange[0] - 1;
	
	// Save that a new simulation
	_finished = false;
}



/* inline double ImageIterator::operator++ ()
 *
 * Get next image distance
 */

inline double ImageIterator::operator++ ()
{
	
	// Loop until a point is found that is in range or all points have been explored
	while (!_finished)
	{
		
		// Go to next cell vector
		if (++_z > _zRange[1])
		{
			_z = _zRange[0];
			if (++_y > _yRange[1])
			{
				_y = _yRange[0];
				if (++_x > _xRange[1])
					_finished = true;
			}
		}
		
		// Break if finished
		if (_finished)
		{
			_distance = -1;
			break;
		}
		
		// Get the current distance
		_curVector = _nearVector;
		_curVector[0] += _x;
		_curVector[1] += _y;
		_curVector[2] += _z;
		_distance = curDistanceSquared();
		
		// Return distance if in range
		if (_distance <= _maxDistanceSquared)
		{
			_distance = sqrt(_distance);
			break;
		}
	}
	
	// Return distance
	return _distance;
}



/* inline double ImageIterator::curDistanceSquared()
 *
 * Return the length of the current vector squared
 */

inline double ImageIterator::curDistanceSquared()
{
	double res = _curVector[0]*_curVector[0]*_redMults[0];
	res += _curVector[1]*_curVector[1]*_redMults[1];
	res += _curVector[2]*_curVector[2]*_redMults[2];
	res += _curVector[0]*_curVector[1]*_redMults[3];
	res += _curVector[0]*_curVector[2]*_redMults[4];
	res += _curVector[1]*_curVector[2]*_redMults[5];
	return res;
}



#endif
