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



#ifndef EWALD_H
#define EWALD_H



#include "locPotential.h"
#include "iso.h"
#include "elements.h"
#include "symmetry.h"
#include "text.h"
#include "constants.h"
#include <cmath>



// Ewald object
class Ewald : public SingleLocalPotential
{
	
	// Elements and charges
	List<double> _charges;
	OList<Element> _elements;
	
	// General variables
	double _perm;
	double _accuracy;
	
	// Helper variables
	mutable double _alpha;
	mutable ImageIterator _realIterator;
	mutable List<double> _recipFactors;
	mutable Linked<Vector3D > _recipVectors;
	
	// Functions
	void initialize(const ISO& iso, int numUniqueAtoms) const;
	double realEnergy(const ISO& iso, Atom* atom, bool skipLowerAtoms) const;
	double recipEnergy(const ISO& iso) const;
	double selfEnergy(const ISO& iso) const;
	double chargedEnergy(const ISO& iso) const;
	double realEnergy(double distance) const	{ return erfc(_alpha * distance) / distance; }
	
	// Helper functions
	double getCharge(const Element& element) const;
	
public:
	
	// Constructor
	Ewald()	{ _perm = Constants::eps0; _accuracy = 1e-8; }
	
	// Setup by file input
	void set(const Text& input);
	
	// Evaluation functions
	void evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces) const;
	void evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, \
		OList<Vector3D >* totalForces) const;
};



#endif
