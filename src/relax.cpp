/* relax.cpp -- Relax structure
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "relax.h"
#include "output.h"
#include <cmath>



/* void Relax::structure(ISO& iso, const LocalPotential& potential, const Symmetry* symmetry) const
 *
 * Run relaxation
 */

void Relax::structure(ISO& iso, const LocalPotential& potential, const Symmetry* symmetry) const
{
	
	// Steepest descent run
	if (_relaxMethod == RM_STEEPEST_DESCENT)
		_getLineDirection = &Relax::SD;
	
	// Conjugate gradient run
	else if (_relaxMethod == RM_CONJUGATE_GRADIENT)
		_getLineDirection = &Relax::CG;
	
	// Unknown method
	else
	{
		Output::newline(ERROR);
		Output::print("Unknown relaxation method");
		Output::quit();
	}
	
	// Output
	Output::newline();
	Output::print("Relaxing internal coordinates using ");
	Output::print(relaxMethod(_relaxMethod).tolower());
	Output::print(" algorithm");
	Output::increase();
	
	// Loop until max loops is reached or converged
	int i, j, k;
	int loopNum;
	double stepScale;
	Vector3D newPos;
	OList<Vector3D > direction;
	for (loopNum = 0; loopNum < _maxIterations; ++loopNum)
	{
		
		// Output
		Output::newline();
		if (loopNum == 0)
			Output::print("Initial structure");
		else
		{
			Output::print("Step ");
			Output::print(loopNum);
		}
		Output::increase();

		// Evaluate the forces
		_prevForces = _forces;
		evaluateForces(_forces, iso, potential, symmetry, true);

		// Output
		Output::decrease();

		// Break if converged
		if (areForcesConverged())
			break;
		
		// Break if at max iterations
		if (loopNum == _maxIterations - 1)
			break;
		
		// Get the line search direction
		(this->*_getLineDirection)(direction);
		
		// Get the step size to make
		stepScale = lineSearch(direction, iso, potential, symmetry);
		
		// Set the new positions
		for (i = 0; i < iso.atoms().length(); ++i)
		{
			for (j = 0; j < iso.atoms()[i].length(); ++j)
			{
				newPos = _origPositions[i][j];
				for (k = 0; k < 3; ++k)
					newPos[k] += stepScale * direction[iso.atoms()[i][j].atomNumber()][k];
				iso.atoms()[i][j].cartesian(newPos);
			}
		}
	}
	
	// Did not converged
	if (loopNum >= _maxIterations - 1)
	{
		Output::newline(WARNING);
		Output::print("Failed to reach convergence criterion");
	}
	
	// Output
	Output::decrease();
}



/* void Relax::SD(OList<Vector3D >& direction) const
 *
 * Steepest descent minimization
 */

void Relax::SD(OList<Vector3D >& direction) const
{
	direction = _forces;
}



/* void Relax::CG(OList<Vector3D >& direction) const
 *
 * Conjugate gradient minimization
 */

void Relax::CG(OList<Vector3D >& direction) const
{
	
	// Get the magnitude of the previous forces
	int i;
	double prevForceMag = 0;
	for (i = 0; i < _prevForces.length(); ++i)
		prevForceMag += _prevForces[i] * _prevForces[i];
	
	// Set the scaling factor
	double scale = 0;
	if (prevForceMag != 0)
	{
		double forceMag = 0;
		for (i = 0; i < _forces.length(); ++i)
			forceMag += _forces[i] * _forces[i];
		scale = forceMag / prevForceMag;
	}
	
	// Set the direction
	if ((scale == 0) || (scale > 1))
		direction = _forces;
	else
	{
		int j;
		for (i = 0; i < _forces.length(); ++i)
		{
			for (j = 0; j < 3; ++j)
				direction[i][j] = _forces[i][j] + scale * direction[i][j];
		}
	}
}



/* double Relax::lineSearch(OList<Vector3D >& direction, ISO& iso, const LocalPotential& potential,
 *		Symmetry* symmetry) const
 *
 * Search for minimum along set direction
 */

double Relax::lineSearch(OList<Vector3D >& direction, ISO& iso, const LocalPotential& potential, \
	const Symmetry* symmetry) const
{
	
	// Save the current positions
	int i, j;
	_origPositions.length(iso.atoms().length());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		_origPositions[i].length(iso.atoms()[i].length());
		for (j = 0; j < iso.atoms()[i].length(); ++j)
			_origPositions[i][j] = iso.atoms()[i][j].cartesian();
	}
	
	// Get the max step size
	double curStep;
	double curMaxStep = 0;
	for (i = 0; i < direction.length(); ++i)
	{
		curStep = direction[i]*direction[i];
		if (curStep > curMaxStep)
			curMaxStep = curStep;
	}
	
	// Set the scaling
	curMaxStep = sqrt(curMaxStep);
	double scale = (_lineSearchScale*curMaxStep > _maxStep) ? _maxStep / curMaxStep : _lineSearchScale;	
	
	// Set the next positions
	int k;
	Vector3D newPos;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			for (k = 0; k < 3; ++k)
				newPos[k] = _origPositions[i][j][k] + scale*direction[iso.atoms()[i][j].atomNumber()][k];
			iso.atoms()[i][j].cartesian(newPos);
		}
	}
	
	// Evaluate the current forces
	evaluateForces(_forcesNext, iso, potential, symmetry, false);
	
	// Normalize direction vector
	double normalizer = 0;
	_dirNorm = direction;
	for (i = 0; i < _dirNorm.length(); ++i)
		normalizer += _dirNorm[i] * _dirNorm[i];
	normalizer = sqrt(normalizer);
	for (i = 0; i < _dirNorm.length(); ++i)
	{
		for (j = 0; j < 3; ++j)
			_dirNorm[i][j] /= normalizer;
	}
	
	// Project original and next forces onto direction vector
	double forceProjOrig = 0;
	double forceProjNext = 0;
	for (i = 0; i < _dirNorm.length(); ++i)
	{
		forceProjOrig += _forces[i] * _dirNorm[i];
		forceProjNext += _forcesNext[i] * _dirNorm[i];
	}
	
	// Get the optimal position along the line (0 = orig, 1 = next)
	double opt = forceProjOrig / (forceProjOrig - forceProjNext);

	// Make sure result is in the range of -1 to 1
	if (Num<double>::abs(opt) > 1)
		opt *= 1 / Num<double>::abs(opt);
	
	// Return the scaling factor
	return opt * scale;
}



/* void Relax::evaluateForces(OList<Vector3D >& forces, ISO& iso, const LocalPotential& potential,
 *		const Symmetry* symmetry, bool print) const
 * 
 * Evaluate the forces of the current structure
 */

void Relax::evaluateForces(OList<Vector3D >& forces, ISO& iso, const LocalPotential& potential, \
	const Symmetry* symmetry, bool print) const
{
	
	// Not using symmetry
	if (!symmetry)
		potential.single(iso, 0, &forces);
	
	// Using symmetry
	else
		potential.single(iso, *symmetry, 0, &forces);
	
	// Save cartesian forces
	int i;
	for (i = 0; i < forces.length(); ++i)
		iso.basis().toCartesian(forces[i]);
	
	// Print forces
	if (print)
	{
		int j;
		for (i = 0; i < forces.length(); ++i)
		{
			Output::newline();
			Output::print("Atom ");
			Output::print(i+1);
			Output::print(": ");
			for (j = 0; j < 3; ++j)
			{
				Output::printSci(forces[i][j], 8);
				Output::print(" ");
			}
		}
	}
}



/* bool Relax::areForcesConverged() const
 *
 * Check if all forces are below convergence tolerance
 */

bool Relax::areForcesConverged() const
{
	int i, j;
	for (i = 0; i < _forces.length(); ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if ((_forces[i][j] < -_forceTol) || (_forces[i][j] > _forceTol))
				return false;
		}
	}
	Output::newline();
	Output::print("Reached convergence criterion");
	return true;
}
