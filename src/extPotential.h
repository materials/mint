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



#ifndef EXTPOTENTIAL_H
#define EXTPOTENTIAL_H



#include "vasp.h"
#include "espresso.h"
#include "potential.h"
#include "iso.h"
#include "elements.h"
#include "symmetry.h"
#include "text.h"
#include "num.h"
#include "list.h"



// Vasp potential
class VaspPot : public IPO
{
	
	// Variables
	bool _static;
	bool _saveFiles;
	Word _executable;
	OList<Element> _elements;
	OList<Word> _potentials;
	Vasp::Accuracy _accuracy;
	OList<Word>::D2 _tags;
	
public:
	
	// Constructor
	VaspPot()	{ _saveFiles = false; _accuracy = Vasp::VA_NORMAL; _static = true; }
	
	// Setup from file
	void add(const Text& input, PotentialType type);
	
	// Static evaluation
	void single(const ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void single(const ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const { single(iso, totalEnergy, totalForces, restart, reduce); }
	
	// Relax
	void relax(ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void relax(ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const { relax(iso, totalEnergy, totalForces, restart, reduce); }
	
	// NEB
	void neb(OList<ISO>& isos, double* tsEnergy = 0, ISO* tsISO = 0) const;
	
	// Other functions
	bool usesSymmetry() const	{ return false; }
	bool supportsNEB()  const	{ return true;  }
};



// Quantum espresso potential
class QEPot : public IPO
{
	
	// Variables
	bool _static;
	bool _saveFiles;
	Word _executable;
	Word _potDirectory;
	OList<Element> _elements;
	OList<Word> _potentials;
	Espresso::Accuracy _accuracy;
	
public:
	
	// Constructor
	QEPot()	{ _saveFiles = false; _accuracy = Espresso::EA_NORMAL; _static = true; }
	
	// Setup from file
	void add(const Text& input, PotentialType type);
	
	// Static evaluation
	void single(const ISO& iso, double* energy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void single(const ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const { single(iso, totalEnergy, totalForces, restart, reduce); }
	
	// Relax
	void relax(ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, bool restart = false, \
		bool reduce = true) const;
	void relax(ISO& iso, const Symmetry& symmetry, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0, \
		bool restart = false, bool reduce = true) const { relax(iso, totalEnergy, totalForces, restart, reduce); }
	
	// NEB
	void neb(OList<ISO>& isos, double* tsEnergy = 0, ISO* tsISO = 0) const;
	
	// Other functions
	bool usesSymmetry() const	{ return false; }
	bool supportsNEB()  const	{ return false; }
};



// Class to convert structure to reduced primitive cell before external potential call
class ReduceISO
{
	
	// Variables
	Matrix3D _cellConversion;
	Matrix3D _positionConversion;
	List<Atom*>::D2 _atomMap;
	OList<Vector3D >::D2 _translations;
	
public:
	
	// Functions
	double reduce(ISO& primISO, const ISO& unitISO, double tol, bool reduce);
	void expand(ISO& unitISO, const ISO& primISO);
	void expandForces(OList<Vector3D >& totalForces, OList<Vector3D > primForces, \
		const ISO& unitISO, const ISO& primISO);
};



/* inline void QEPot::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
 *
 * NEB calculation
 */

inline void QEPot::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
{
	if (tsEnergy)
		*tsEnergy = 0;
	if (tsISO)
		tsISO->clear();
}



#endif
