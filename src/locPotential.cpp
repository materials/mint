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
