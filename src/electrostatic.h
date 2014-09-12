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

#ifndef POTENTIALPEGS_H
#define	POTENTIALPEGS_H

#include "locPotential.h"
#include "ewald.h"
#include "text.h"
#include "pairPotential.h"
#include <vector>

/**
 * Defines the electrostatic potential similar to that used in PEGS
 * 
 * E(R) = [columbic repulsion] + [hard-sphere potential]
 * 
 * PEGS (as implemented by Eric Majzoub) uses soft-sphere potential that is only 
 * taken into account when ions overlap. The discontinuity of this function causes
 * problems during relaxation. To overcome this problem, this implementation uses
 * the following function to represent hard-sphere interaction
 * 
 * E(R) = R &lt; R_HS ? a * (exp((R - R_HS) ^ 2) - 1) : 0
 * 
 * This function has the advantage being having continuous second and third derivatives
 *  with respect to R, which allows it to be used with most molecular-dynamics codes.
 *  (This is kind of how LAMMPS does it: <a href="http://lammps.sandia.gov/threads/msg08499.html">
 *  Thank you, Steve Plimpton</a>)
 * 
 * The parameter a is chosen to be on the order of 1000, which keeps the minimum distance for a
 *  pair of oppositely charged atoms to within 0.1 A of the hard-sphere radius for 
 *  reasonable charges.
 * 
 * Input format: &lt;element #1 Symbol, charge, radius&gt; &lt; &lt;that data for each other element&gt;
 * 
 * Options: See options for Ewald (in set)
 * 
 */
class Electrostatic : public Ewald {
	
	// Hard-sphere radius of each element (same order as inherited _elements array)
	vector<double> _radii;
	
	// Hard sphere potentials
	vector<HardSphere> _potentials;

public:
		
	void set(const Text& input);

	virtual void evaluate(const ISO& iso, double* energy, OList<Vector3D>* forces) const;
	virtual void evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, OList<Vector3D>* totalForces) const;

};

#endif	/* POTENTIALPEGS_H */

