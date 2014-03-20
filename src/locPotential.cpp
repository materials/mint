/* locPotential.cpp -- Handle local potentials
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */

#include "locPotential.h"
#include "pairPotential.h"
#include "ewald.h"
#include "relax.h"
#include "output.h"



/* void SingleLocalPotential::readError(const OList<Word>& line)
 *
 * Error in reading file
 */

void SingleLocalPotential::readError(const OList<Word>& line)
{
	Output::newline(ERROR);
	Output::print("Did not recognize potential setting on line \"");
	for (int i = 0; i < line.length(); ++i)
	{
		Output::print(line[i]);
		if (i != line.length() - 1)
			Output::print(" ");
	}
	Output::print("\"");
	Output::quit();
}



/* void LocalPotential::add(const Text& input, PotentialType type)
 *
 * Add potential to local potential list
 */

void LocalPotential::add(const Text& input, PotentialType type)
{
	
	// Ewald potential
	if (type == PT_EWALD)
		_potentials += new Ewald;
		
	// Lennard-Jones potential
	else if (type == PT_LENNARDJONES)
		_potentials += new LennardJones;
	
	// Buckingham potential
	else if (type == PT_BUCKINGHAM)
		_potentials += new Buckingham;
	
	// Power potential
	else if (type == PT_POWER)
		_potentials += new Power;
		
	// Exponential potential
	else if (type == PT_EXPONENTIAL)
		_potentials += new Exponential;
		
	// Covalent potential
	else if (type == PT_COVALENT)
		_potentials += new Covalent;
		
	// Unknown potential
	else
	{
		Output::newline(ERROR);
		Output::print("Internal error in setting local potential");
		Output::quit();
	}
	
	// Add data for potential
	_potentials.last()->set(input);
}



/* void LocalPotential::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart,
 *		bool reduce) const
 *
 * Relax structure
 */

void LocalPotential::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, \
	bool reduce) const
{
	initialize(iso, totalEnergy, totalForces);
	Relax relax;
	relax.structure(iso, *this);
	if ((totalEnergy) || (totalForces))
		single(iso, totalEnergy, totalForces);
}



/* void LocalPotential::relax(ISO& iso, const Symmetry& symmetry, double* totalEnergy, OList<Vector3D >* totalForces,
 *		bool restart, bool reduce) const
 *
 * Relax structure
 */

void LocalPotential::relax(ISO& iso, const Symmetry& symmetry, double* totalEnergy, OList<Vector3D >* totalForces, \
	bool restart, bool reduce) const
{
	initialize(iso, totalEnergy, totalForces);
	Relax relax;
	relax.structure(iso, *this, symmetry);
	if ((totalEnergy) || (totalForces))
		single(iso, symmetry, totalEnergy, totalForces);
}
