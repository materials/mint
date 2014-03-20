/* potential.h -- Handle energy and forces
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef POTENTIAL_H
#define POTENTIAL_H



#include "num.h"
#include "iso.h"
#include "elements.h"
#include "symmetry.h"
#include "fileSystem.h"
#include "list.h"
#include "text.h"
#include "output.h"
#include <cmath>



// Types of energy functions
enum PotentialType {PT_UNKNOWN, PT_VASP, PT_QE, PT_EWALD, PT_LENNARDJONES, PT_BUCKINGHAM, PT_POWER, PT_EXPONENTIAL, \
	PT_COVALENT};



// Declare that Potential object will be defined later
class Potential;



// Class to store a single potential object (internal potential object)
class IPO
{
	
protected:
	
	// Functions
	void readError(const OList<Word>& line);
	
public:
	
	// Destructor
	virtual ~IPO() {}
	
	// Setup from file
	virtual void add(const Text& input, PotentialType type) = 0;
	
	// Static evaluation
	virtual void single(const ISO& iso, double* energy = 0, OList<Vector3D >* forces = 0, bool restart = false, \
		bool reduce = true) const = 0;
	virtual void single(const ISO& iso, const Symmetry& symmetry, double* energy = 0, OList<Vector3D >* forces = 0, \
		bool restart = false, bool reduce = true) const = 0;
	
	// Relaxation
	virtual void relax(ISO& iso, double* energy = 0, OList<Vector3D >* forces = 0, bool restart = false, \
		bool reduce = true) const = 0;
	virtual void relax(ISO& iso, const Symmetry& symmetry, double* energy = 0, OList<Vector3D >* forces = 0, \
		bool restart = false, bool reduce = true) const = 0;
	
	// NEB
	virtual void neb(OList<ISO>& isos, double* tsEnergy = 0, ISO* tsISO = 0) const = 0;
	
	// Other functions
	virtual bool usesSymmetry() const = 0;
	virtual bool supportsNEB() const = 0;
};



// Class to store a reference potential
class Reference
{

	// Variables
	ISO _iso;
	Element _element;
	double _energyPerAtom;
	
public:
	
	// Setup functions
	void set(double value)					{ _energyPerAtom = value; }
	void set(const Element& element, bool setISO);
	void set(const Potential& potential);
	
	// Access functions
	const Element& element() const			{ return _element; }
	double energyPerAtom() const			{ return _energyPerAtom; }
};



// Class to store a potential
class Potential
{
	
	// Variables
	bool _useReferences;
	double _unphysicalCutoff;
	IPO* _ipo;
	mutable OList<Reference> _references;
	
	// Functions
	bool finish(const ISO& iso, double* energy, OList<Vector3D >* forces) const;
	double getReferenceEnergy(const ISO& iso) const;
	void setReferenceData(const Text& input);
	void errorIfNotSet() const;
	
	// Static functions
	static void initialize(const ISO& iso, double* energy, OList<Vector3D >* forces);
	static OList<Text> parseInput(const Text& input);
	static PotentialType potentialType(const Word& word);
	static Word potentialType(PotentialType type);
	static void error(PotentialType setType, PotentialType otherType);
	
public:
	
	// Destructor
	Potential()		{ _ipo = 0; _useReferences = false; _unphysicalCutoff = -5; }
	~Potential()	{ clear(); }
	
	// Setup functions
	void clear();
	
	// Setup from file
	void set(const Text& text);
	void set(const Word& file)	{ set(Read::text(file)); }
	
	// Other functions
	static bool isFormat(const Text& text)	{ return (parseInput(text).length() != 0); }
	static bool isFormat(const Word& file)	{ return (parseInput(Read::text(file)).length() != 0); }
	
	// Static evaluation
	bool single(const ISO& iso, double* energy = 0, OList<Vector3D >* forces = 0, bool restart = false, \
		bool reduce = true) const;
	bool single(const ISO& iso, const Symmetry& symmetry, double* energy = 0, OList<Vector3D >* forces = 0, \
		bool restart = false, bool reduce = true) const;
	
	// Relaxation
	bool relax(ISO& iso, double* energy = 0, OList<Vector3D >* forces = 0, bool restart = false, \
		bool reduce = true) const;
	bool relax(ISO& iso, const Symmetry& symmetry, double* energy = 0, OList<Vector3D >* forces = 0, \
		bool restart = false, bool reduce = true) const;
	
	// Nudged elastic band calculation
	void neb(OList<ISO>& isos, double* tsEnergy = 0, ISO* tsISO = 0) const;
	
	// Access functions
	bool isSet() const			{ return (_ipo != 0); }
	bool usesSymmetry() const	{ return isSet() ? _ipo->usesSymmetry() : false; }
	bool supportsNEB()  const	{ return isSet() ? _ipo->supportsNEB()  : false; }
	bool useReferences() const	{ return _useReferences; }
	
	// Friends
	friend class Reference;
};



// =====================================================================================================================
// Potential
// =====================================================================================================================

/* inline void Potential::clear()
 *
 * Destructor for Potential object
 */

inline void Potential::clear()
{
	_references.clear();
	_useReferences = false;
	_unphysicalCutoff = -5;
	if (_ipo)
		delete _ipo;
	_ipo = 0;
}



/* inline bool Potential::single(const ISO& iso, double* energy, OList<Vector3D >* forces, bool restart,
 *		bool reduce) const
 *
 * Evaluate the potential
 */

inline bool Potential::single(const ISO& iso, double* energy, OList<Vector3D >* forces, bool restart, \
	bool reduce) const
{
	errorIfNotSet();
	initialize(iso, energy, forces);
	_ipo->single(iso, energy, forces, restart, reduce);
	return finish(iso, energy, forces);
}

inline bool Potential::single(const ISO& iso, const Symmetry& symmetry, double* energy, OList<Vector3D >* forces, \
	bool restart, bool reduce) const
{
	errorIfNotSet();
	initialize(iso, energy, forces);
	_ipo->single(iso, symmetry, energy, forces, restart, reduce);
	return finish(iso, energy, forces);
}



/* inline bool Potential::relax(ISO& iso, double* energy, OList<Vector3D >* forces, bool restart, bool reduce) const
 *
 * Relax structure under potential
 */

inline bool Potential::relax(ISO& iso, double* energy, OList<Vector3D >* forces, bool restart, bool reduce) const
{
	errorIfNotSet();
	initialize(iso, energy, forces);
	_ipo->relax(iso, energy, forces, restart, reduce);
	return finish(iso, energy, forces);
}

inline bool Potential::relax(ISO& iso, const Symmetry& symmetry, double* energy, OList<Vector3D >* forces, \
	bool restart, bool reduce) const
{
	errorIfNotSet();
	initialize(iso, energy, forces);
	_ipo->relax(iso, symmetry, energy, forces, restart, reduce);
	return finish(iso, energy, forces);
}



/* inline bool Potential::finish(const ISO& iso, double* energy, OList<Vector3D >* forces) const
 *
 * Finish any calculations that are needed
 */

inline bool Potential::finish(const ISO& iso, double* energy, OList<Vector3D >* forces) const
{
	
	// Return there was no problem if energy was not set or not using references
	if ((!energy) || (!_useReferences))
		return true;
	
	// Subtract the reference energy
	*energy -= getReferenceEnergy(iso);
	
	// Warning if energy is unphysical
	if (*energy < _unphysicalCutoff)
	{
		Output::newline(WARNING);
		Output::print("Found an unphysical formation energy of ");
		Output::print(*energy/iso.numAtoms());
		Output::print(" eV/atom");
		return false;
	}
	
	// Return that run was okay
	return true;
}



/* inline void Potential::initialize(const ISO& iso, double* energy, OList<Vector3D >* forces)
 *
 * Initialize energy and forces for new evaluation
 */

inline void Potential::initialize(const ISO& iso, double* energy, OList<Vector3D >* forces)
{
	if (energy)
		*energy = 0;
	if (forces)
	{
		forces->length(iso.numAtoms());
		forces->fill(0.0);
	}
}



/* inline void Potential::error(PotentialType setType, PotentialType otherType)
 *
 * Print error message for mixing potential types
 */

inline void Potential::error(PotentialType setType, PotentialType otherType)
{
	Output::newline(ERROR);
	Output::print("Cannot mix ");
	Output::print(potentialType(setType));
	Output::print(" potential with ");
	Output::print(potentialType(otherType));
	Output::print(" potential");
	Output::quit();
}



/* inline void Potential::errorIfNotSet() const
 *
 * Exit if potential has not been set
 */

inline void Potential::errorIfNotSet() const
{
	if (_ipo == 0)
	{
		Output::newline(ERROR);
		Output::print("Cannot evaluate potential when it has not been set");
		Output::quit();
	}
}



#endif
