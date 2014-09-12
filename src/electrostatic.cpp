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

#include <vector>
#include <stdlib.h>

#include "electrostatic.h"
#include "text.h"
#include "pairPotential.h"

/**
 * Read in options from file. See class documentation for those options
 * @param input Text input from file
 */
void Electrostatic::set(const Text& input) {
	
	// Read in elements
	for (int pos=1; pos < input[0].length(); pos += 3) {
		try {
			_elements += Element::find(input[0][pos], false);
			_charges += atof(input[0][pos+1].array());
			_radii.push_back(atof(input[0][pos+2].array()));
		} catch(...) {
			readError(input[0]);
		}
	}
	
	// Read in other options
	setEwaldOptions(input, true);
	
	// Define the hard-sphere potentials
	_potentials.clear();
	for (int elem1=0; elem1 < _elements.length(); elem1++) {
		for (int elem2=elem1; elem2 < _elements.length(); elem2++) {
			HardSphere hs;
			hs.setForce(1000.0);
			hs.setRadius(_radii[elem1] + _radii[elem2]);
			hs.setElementOne(_elements[elem1]);
			hs.setElementTwo(_elements[elem2]);
			_potentials.push_back(hs);
		}
	}
}

void Electrostatic::evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D>* totalForces) const {
	Ewald::evaluate(iso, totalEnergy, totalForces);
	
	for (int i=0; i<_potentials.size(); i++) {
		_potentials[i].evaluate(iso, totalEnergy, totalForces);
	}
}

void Electrostatic::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, OList<Vector3D>* totalForces) const {
	Ewald::evaluate(iso, symmetry, totalEnergy, totalForces);
	
	for (int i=0; i<_potentials.size(); i++) {
		_potentials[i].evaluate(iso, symmetry, totalEnergy, totalForces);
	}
}
