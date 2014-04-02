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



#ifndef RELAX_H
#define RELAX_H



#include "iso.h"
#include "symmetry.h"
#include "locPotential.h"
#include "text.h"
#include "num.h"
#include "list.h"



// Relaxation methods
enum RelaxMethod {RM_UNKNOWN, RM_STEEPEST_DESCENT, RM_CONJUGATE_GRADIENT};



// Class to perform structure relaxation
class Relax
{
	
	// Variable to store relaxation method
	mutable void (Relax::*_getLineDirection)(OList<Vector3D >&) const;
	
	// Settings
	int _maxIterations;
	double _forceTol;
	double _lineSearchScale;
	double _maxStep;
	RelaxMethod _relaxMethod;
	
	// Storage variables
	mutable OList<Vector3D > _forces;
	mutable OList<Vector3D > _prevForces;
	mutable OList<Vector3D > _forcesNext;
	mutable OList<Vector3D >::D2 _origPositions;
	mutable OList<Vector3D > _dirNorm;
	
	// Functions
	void structure(ISO& iso, const LocalPotential& potential, const Symmetry* symmetry) const;
	
	// Minimization methods
	void SD(OList<Vector3D >& direction) const;
	void CG(OList<Vector3D >& direction) const;
	
	// Helper functions
	double lineSearch(OList<Vector3D>& direction, ISO& iso, const LocalPotential& potential, \
		const Symmetry* symmetry) const;
	void evaluateForces(OList<Vector3D>& forces, ISO& iso, const LocalPotential& potential, \
		const Symmetry* symmetry, bool print) const;
	bool areForcesConverged() const;

public:
	
	// Constructor
	Relax();
	
	// Setup functions
	void maxIterations(int input)		{ _maxIterations = input; }
	void forceTolerance(double input)	{ _forceTol = input; }
	void lineSearchScale(double input)	{ _lineSearchScale = input; }
	void maxStep(double input)			{ _maxStep = input; }

	// Functions
	void structure(ISO& iso, const LocalPotential& potential) const
		{ structure(iso, potential, 0); }
	void structure(ISO& iso, const LocalPotential& potential, const Symmetry& symmetry) const
		{ structure(iso, potential, &symmetry); }
	
	// Helper functions
	static RelaxMethod relaxMethod(const Word& method);
	static Word relaxMethod(RelaxMethod method);
};



/* inline Relax::Relax()
 *
 * Constructor for Relax object
 */

inline Relax::Relax()
{
	_maxIterations = 100;
	_forceTol = 1e-9;
	_lineSearchScale = 1;
	_maxStep = 0.01;
	_relaxMethod = RM_CONJUGATE_GRADIENT;
}



/* inline RelaxMethod Relax::relaxMethod(const Word& method)
 *
 * Convert word to relaxation method
 */

inline RelaxMethod Relax::relaxMethod(const Word& method)
{
	if (method.equal("steepest", false, 5))
		return RM_STEEPEST_DESCENT;
	if (method.equal("conjugate", false, 4))
		return RM_CONJUGATE_GRADIENT;
	return RM_UNKNOWN;
}



/* inline Word Relax::relaxMethod(RelaxMethod method)
 *
 * Convert relaxation method to word
 */

inline Word Relax::relaxMethod(RelaxMethod method)
{
	switch (method)
	{
		case RM_STEEPEST_DESCENT:
			return Word("Steepest descent");
		case RM_CONJUGATE_GRADIENT:
			return Word("Conjugate gradient");
		default:
			return Word("Unknown");
	}
}



#endif
