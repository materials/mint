/* locPotential.h -- Handle local potentials
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef LOCPOTENTIAL_H
#define LOCPOTENTIAL_H



#include "potential.h"
#include "iso.h"
#include "symmetry.h"
#include "text.h"
#include "num.h"
#include "list.h"



// Single local potential
class SingleLocalPotential
{
protected:
	
	// Functions
	void readError(const OList<Word>& line);
	
public:
	
	// Virtual functions
	virtual ~SingleLocalPotential() {}
	virtual void set(const Text& input) = 0;
	virtual void evaluate(const ISO& iso, double* energy = 0, OList<Vector3D >* forces = 0) const = 0;
	virtual void evaluate(const ISO& iso, const Symmetry& symmetry, double* energy = 0, \
		OList<Vector3D >* forces = 0) const = 0;
};



// Local potential
class LocalPotential : public IPO
{
	
	// Variables
	List<SingleLocalPotential*> _potentials;
	
	// Functions
	void initialize(const ISO& iso, double* energy, OList<Vector3D >* forces) const;
	
public:
	
	// Destructor
	~LocalPotential();
	
	// Setup from file
	void add(const Text& input, PotentialType type);
	
	// Static evaluation
	void single(const ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void single(const ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const;
	
	// Relax
	void relax(ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void relax(ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const;
	
	// NEB
	void neb(OList<ISO>& isos, double* tsEnergy = 0, ISO* tsISO = 0) const;
	
	// Other functions
	bool usesSymmetry() const	{ return true;  }
	bool supportsNEB()  const	{ return false; }
};



/* inline LocalPotential::~LocalPotential()
 *
 * Destructor for local potential object
 */

inline LocalPotential::~LocalPotential()
{
	for (int i = 0; i < _potentials.length(); ++i)
		delete _potentials[i];
}



/* inline void LocalPotential::initialize(const ISO& iso, double* energy, OList<Vector3D >* forces) const
 *
 * Initialize variables for evaluation
 */

inline void LocalPotential::initialize(const ISO& iso, double* energy, OList<Vector3D >* forces) const
{
	if (energy)
		*energy = 0;
	if (forces)
	{
		forces->length(iso.numAtoms());
		forces->fill(0.0);
	}
}



/* inline void LocalPotential::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces,
 *		bool restart, bool reduce) const
 *
 * Calculate energy and forces of structure
 */

inline void LocalPotential::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, \
	bool restart, bool reduce) const
{
	initialize(iso, totalEnergy, totalForces);
	for (int i = 0; i < _potentials.length(); ++i)
		_potentials[i]->evaluate(iso, totalEnergy, totalForces);
}



/* inline void LocalPotential::single(const ISO& iso, const Symmetry& symmetry, double* totalEnergy,
 *		OList<Vector3D >* totalForces, bool restart, bool reduce) const
 *
 * Calculate energy and forces of structure
 */

inline void LocalPotential::single(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, \
	OList<Vector3D >* totalForces, bool restart, bool reduce) const
{
	initialize(iso, totalEnergy, totalForces);
	for (int i = 0; i < _potentials.length(); ++i)
		_potentials[i]->evaluate(iso, symmetry, totalEnergy, totalForces);
}



/* inline void LocalPotential::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
 *
 * NEB calculation
 */

inline void LocalPotential::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
{
	if (tsEnergy)
		*tsEnergy = 0;
	if (tsISO)
		tsISO->clear();
}



#endif
