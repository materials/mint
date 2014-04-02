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



#include "extPotential.h"
#include "language.h"
#include "output.h"



/* void VaspPot::add(const Text& input, PotentialType type)
 *
 * Set vasp potential
 */

void VaspPot::add(const Text& input, PotentialType type)
{
	
	// No data
	if (!input.length())
		return;
		
	// Output
	Output::newline();
	Output::print("External vasp potential");
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 2)
		readError(input[0]);
	
	// Save executable
	_executable = input[0][1];
	
	// Loop over lines
	int i, j;
	for (i = 1; i < input.length(); ++i)
	{
		
		// Skip if too short
		if (input[i].length() < 2)
			continue;
		
		// Skip if a comment
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found element
		if (Element::isElement(input[i][0], true))
		{
			_elements += Element::find(input[i][0], true);
			_potentials += input[i][1];
		}
		
		// Found accuracy
		else if (input[i][0].equal("accuracy", false, 3))
		{
			Vasp::Accuracy tempAccuracy = Vasp::accuracy(input[i][1]);
			if (tempAccuracy != Vasp::VA_UNKNOWN)
				_accuracy = tempAccuracy;
			else
				readError(input[i]);
		}
		
		// Found save files
		else if (input[i][0].equal("save", false, 4))
		{
			if ((input[i][1][0] == 't') || (input[i][1][0] == 'T'))
				_saveFiles = true;
			else if ((input[i][1][0] == 'f') || (input[i][1][0] == 'F'))
				_saveFiles = false;
			else
				readError(input[i]);
		}
		
		// Found static
		else if (input[i][0].equal("static", false, 4))
		{
			if ((input[i][1][0] == 't') || (input[i][1][0] == 'T'))
				_static = true;
			else if ((input[i][1][0] == 'f') || (input[i][1][0] == 'F'))
				_static = false;
			else
				readError(input[i]);
		}
		
		// Found vasp setting
		else if (Vasp::Settings::isSetting(input[i][0]))
		{
			_tags.add();
			_tags.last().length(2);
			_tags.last()[0] = input[i][0];
			for (j = 1; j < input[i].length(); ++j)
			{
				if (Language::isComment(input[i][j]))
					break;
				if (j != 1)
					_tags.last()[1] += ' ';
				_tags.last()[1] += input[i][j];
			}
		}
		
		// Did not recognize setting
		else
			readError(input[i]);
	}
	
	// Print settings
	Output::newline();
	Output::print("Executable: ");
	Output::print(_executable);
	Output::newline();
	Output::print("Accuracy: ");
	Output::print(Vasp::accuracy(_accuracy));
	for (i = 0; i < _elements.length(); ++i)
	{
		Output::newline();
		Output::print(_elements[i].symbol());
		Output::print(" potential: ");
		Output::print(_potentials[i]);
	}
	
	// Output
	Output::decrease();
}



/* void VaspPot::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, 
 *		bool reduce) const
 *
 * Run static Vasp calculation
 */

void VaspPot::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, \
	bool reduce) const
{
	
	// Output
	Output::newline();
	Output::print("Performing a static VASP calculation");
	Output::increase();
	
	// Create reduced cell object
	ISO primISO;
	ReduceISO reduceISO;
	double energyScale = reduceISO.reduce(primISO, iso, 1e-6, reduce);
	
	// Create vasp job
	Vasp::Job job;
	job.saveFiles(_saveFiles);
	job.calculation(Vasp::VC_SINGLE);
	job.accuracy(_accuracy);
	job.executable(_executable);
	job.elements(_elements);
	job.potentials(_potentials);
	job.setTags(_tags);
	
	// Submit the job
	double res = 1e10;
	if (job.submit(primISO, restart, false))
	{
		
		// Save the energy
		res = energyScale * job.energy().last();
		
		// Update forces
		if (totalForces)
			reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
	}
	
	// Save energy
	if (totalEnergy)
		*totalEnergy += res;
	
	// Output
	Output::decrease();
}



/* void VaspPot::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, bool reduce) const
 *
 * Run Vasp relaxation
 */

void VaspPot::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, bool reduce) const
{
	
	// Output
	Output::newline();
	Output::print("Relaxing structure with VASP");
	Output::increase();
	
	// Create reduced cell object
	ISO primISO;
	ReduceISO reduceISO;
	double energyScale = reduceISO.reduce(primISO, iso, 1e-6, reduce);
	
	// Create vasp job
	Vasp::Job job;
	job.saveFiles(_saveFiles);
	job.calculation(Vasp::VC_SINGLE);
	job.accuracy(_accuracy);
	job.executable(_executable);
	job.elements(_elements);
	job.potentials(_potentials);
	job.setTags(_tags);
	
	// Run a static calculation to get forces
	OList<Vector3D> curForces(iso.numAtoms());
	curForces.fill(Vector3D(0.0));
	if (job.submit(primISO, restart, false))
		reduceISO.expandForces(curForces, job.forces().last(), iso, primISO);
	
	// Get the max force in the structure
	int i;
	double curForce;
	double maxForce = 0;
	for (i = 0; i < curForces.length(); ++i)
	{
		curForce = iso.basis().getCartesian(curForces[i]).magnitude();
		if (curForce > maxForce)
			maxForce = curForce;
	}
	
	// Set the relaxaton type
	Vasp::Relaxation relaxation = Vasp::VR_UNKNOWN;
/*	if (maxForce > 1)
	{
		relaxation = Vasp::VR_DAMPED_MD;
		Output::newline();
		Output::print("The forces are very large - using damped MD relaxation");
	}
	else*/ if (maxForce > 0.1)
	{
		relaxation = Vasp::VR_CONJUGATE_GRADIENT;
		Output::newline();
		Output::print("The forces are not small - using conjugate gradient relaxation");
	}
	else
	{
		relaxation = Vasp::VR_QUASI_NEWTON;
		Output::newline();
		Output::print("The forces are small - using quasi-newton relaxation");
	}
	job.relaxation(relaxation);
	
	// Check if basis is fixed
	job.calculation(Vasp::VC_RELAX_ALL);
	for (i = 0; i < 3; ++i)
	{
		if ((primISO.basis().lengthFixed()[i]) || (primISO.basis().angleFixed()[i]))
		{
			job.calculation(Vasp::VC_RELAX_POSITIONS);
			break;
		}
	}
	
	// Submit the relaxation
	double res = 1e10;
	if (job.submit(primISO, true, false))
	{
		
		// Save the new structure
		job.updateStructure(primISO);
		reduceISO.expand(iso, primISO);
		
		// Save properties if not doing a static calculation
		if (!_static)
		{
			res = energyScale * job.energy().last();
			if (totalForces)
				reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
		}
		
		// Perform static calculation
		else
		{
			job.calculation(Vasp::VC_SINGLE);
			if (job.submit(primISO, true, false))
			{
				res = energyScale * job.energy().last();
				if (totalForces)
					reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
			}
		}
	}
	
	// Save energy
	if (totalEnergy)
		*totalEnergy += res;
	
	// Output
	Output::decrease();
}



/* void VaspPot::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
 *
 * Nudged elastic band calculation
 */

void VaspPot::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
{
	
	// Output
	Output::newline();
	Output::print("Running nudged elastic band (NEB) calculation with VASP");
	Output::increase();
	
	// Create vasp job
	Vasp::Job job;
	job.saveFiles(_saveFiles);
	job.calculation(Vasp::VC_NEB);
	job.accuracy(_accuracy);
	job.executable(_executable);
	job.elements(_elements);
	job.potentials(_potentials);
	job.setTags(_tags);
	
	// Submit the calculation
	double res = 1e10;
	if (job.submitNEB(isos, false, false))
	{
		
		// Save the new structures
		job.updateStructures(isos);
		
		// Save the energy
		res = job.energy().last();
	}
	
	// Return if job did not finish
	if (!job.completed())
	{
		if (tsEnergy)
			*tsEnergy += res;
		return;
	}
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Running climbing image nudged elastic band (CINEB) calculation with VASP");
	Output::increase();
	
	// Modify vasp job to perform CINEB calculation
	job.calculation(Vasp::VC_CINEB);
	
	// Submit the job
	res = 1e10;
	if (job.submitNEB(isos, true, false))
	{
		
		// Save the new structure
		job.updateStructures(isos);
		if (tsISO)
			job.updateStructure(*tsISO);
		
		// Save the energy
		res = job.energy().last();
	}
	
	// Save energy
	if (tsEnergy)
		*tsEnergy += res;
	
	// Output
	Output::decrease();
}



/* void QEPot::add(const Text& input, PotentialType type)
 *
 * Set quantum espresso potential
 */

void QEPot::add(const Text& input, PotentialType type)
{
	
	// No data
	if (!input.length())
		return;
		
	// Output
	Output::newline();
	Output::print("External Quantum Espresso potential");
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 3)
		readError(input[0]);
	
	// Save potential directory
	_executable = input[0][1];
	_potDirectory = input[0][2];
	
	// Loop over lines
	int i;
	for (i = 1; i < input.length(); ++i)
	{
		
		// Skip if too short
		if (input[i].length() < 2)
			continue;
		
		// Skip if a comment
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found element
		if (Element::isElement(input[i][0], true))
		{
			_elements += Element::find(input[i][0], true);
			_potentials += input[i][1];
		}
		
		// Found accuracy
		else if (input[i][0].equal("accuracy", false, 3))
		{
			Espresso::Accuracy tempAccuracy = Espresso::accuracy(input[i][1]);
			if (tempAccuracy != Espresso::EA_UNKNOWN)
				_accuracy = tempAccuracy;
			else
				readError(input[i]);
		}
		
		// Found save files
		else if (input[i][0].equal("save", false, 4))
		{
			if ((input[i][1][0] == 't') || (input[i][1][0] == 'T'))
				_saveFiles = true;
			else if ((input[i][1][0] == 'f') || (input[i][1][0] == 'F'))
				_saveFiles = false;
			else
				readError(input[i]);
		}
		
		// Found static
		else if (input[i][0].equal("static", false, 4))
		{
			if ((input[i][1][0] == 't') || (input[i][1][0] == 'T'))
				_static = true;
			else if ((input[i][1][0] == 'f') || (input[i][1][0] == 'F'))
				_static = false;
			else
				readError(input[i]);
		}
		
		// Did not recognize setting
		else
			readError(input[i]);
	}
	
	// Print settings
	Output::newline();
	Output::print("Executable: ");
	Output::print(_executable);
	Output::newline();
	Output::print("Accuracy: ");
	Output::print(Espresso::accuracy(_accuracy));
	Output::newline();
	Output::print("Potential directory: ");
	Output::print(_potDirectory);
	for (i = 0; i < _elements.length(); ++i)
	{
		Output::newline();
		Output::print(_elements[i].symbol());
		Output::print(" potential: ");
		Output::print(_potentials[i]);
	}
	
	// Output
	Output::decrease();
}



/* void QEPot::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart,
 *		bool reduce) const
 *
 * Run static Quantum Espresso calculation
 */

void QEPot::single(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, \
	bool reduce) const
{
	
	// Output
	Output::newline();
	Output::print("Performing a static Quantum Espresso calculation");
	Output::increase();
	
	// Create reduced cell object
	ISO primISO;
	ReduceISO reduceISO;
	double energyScale = reduceISO.reduce(primISO, iso, Num<double>::distanceFromAU(1e-6), reduce);
	
	// Create quantum espresso job
	Espresso::Job job;
	job.saveFiles(_saveFiles);
	job.executable(_executable);
	job.potDirectory(_potDirectory);
	job.calculation(Espresso::EC_SINGLE);
	job.accuracy(_accuracy);
	job.elements(_elements);
	job.potentials(_potentials);
	
	// Submit the job
	double res = 1e10;
	if (job.submit(primISO, false))
	{
		
		// Save the energy
		res = energyScale * job.energy().last();
		
		// Update forces
		if (totalForces)
			reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
	}
	
	// Save energy
	if (totalEnergy)
		*totalEnergy += res;
	
	// Output
	Output::decrease();
}



/* void QEPot::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, bool reduce) const
 *
 * Run Quantum Espresso relaxation
 */

void QEPot::relax(ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces, bool restart, bool reduce) const
{
	
	// Output
	Output::newline();
	Output::print("Relaxing structure with Quantum Espresso");
	Output::increase();
	
	// Create reduced cell object
	ISO primISO;
	ReduceISO reduceISO;
	double energyScale = reduceISO.reduce(primISO, iso, Num<double>::distanceFromAU(1e-6), reduce);
	
	// Check if basis is fixed
	Espresso::Calculation calculation = Espresso::EC_RELAX_ALL;
	for (int i = 0; i < 3; ++i)
	{
		if ((primISO.basis().lengthFixed()[i]) || (primISO.basis().angleFixed()[i]))
		{
			calculation = Espresso::EC_RELAX_POSITIONS;
			break;
		}
	}
	
	// Create quantum espresso job
	Espresso::Job job;
	job.saveFiles(_saveFiles);
	job.executable(_executable);
	job.potDirectory(_potDirectory);
	job.calculation(calculation);
	job.accuracy(_accuracy);
	job.elements(_elements);
	job.potentials(_potentials);
	
	// Submit the job
	double res = 1e10;
	if (job.submit(primISO, false))
	{
		
		// Save the new structure
		job.updateStructure(primISO);
		reduceISO.expand(iso, primISO);
		
		// Save properties if not doing a static calculation
		if (!_static)
		{
			res = energyScale * job.energy().last();
			if (totalForces)
				reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
		}
		
		// Perform static calculation
		else
		{
			job.calculation(Espresso::EC_SINGLE);
			if (job.submit(primISO, false))
			{
				res = energyScale * job.energy().last();
				if (totalForces)
					reduceISO.expandForces(*totalForces, job.forces().last(), iso, primISO);
			}
		}
	}
	
	// Save energy
	if (totalEnergy)
		*totalEnergy += res;
	
	// Output
	Output::decrease();
}



/* double ReduceISO::reduce(ISO& primISO, const ISO& unitISO, double tol, bool reduce)
 *
 * Convert to reduced primitive cell
 */

double ReduceISO::reduce(ISO& primISO, const ISO& unitISO, double tol, bool reduce)
{
	
	// Clear primitive cell structure
	primISO.clear();
	
	// Allocate space
	int i;
	_atomMap.length(unitISO.atoms().length());
	_translations.length(unitISO.atoms().length());
	for (i = 0; i < unitISO.atoms().length(); ++i)
	{
		_atomMap[i].length(unitISO.atoms()[i].length());
		_translations[i].length(unitISO.atoms()[i].length());
	}
	
	// Get transformation to reduced primitive cell
	Matrix3D unitToPrimReduced = (reduce == true) ? unitISO.primitiveTransformation(tol, false) : Matrix3D::identity();
	if (reduce == true)
		unitToPrimReduced *= Basis::reducedTransformation(unitToPrimReduced * unitISO.basis().vectors());
	Matrix3D unitPosToPrimReduced = unitToPrimReduced.inverse().transpose();
	
	// Save reduced basis
	primISO.basis(unitToPrimReduced * unitISO.basis().vectors(), false);
	primISO.basis().lengthFixed(unitISO.basis().lengthFixed());
	primISO.basis().angleFixed(unitISO.basis().angleFixed());
	
	// Save conversions that will be used to convert back to unit cell
	_cellConversion = unitToPrimReduced.inverse();
	_positionConversion = unitPosToPrimReduced.inverse();
	
	// Loop over atoms in the unit cell
	int j, k, m;
	bool isNew;
	double curDis;
	List<int> timesFound;
	List<double> primDis;
	Vector3D newPos;
	Vector3D nearCell;
	Vector3D newMagMom;
	Vector3D origin(0.0);
	for (i = 0; i < unitISO.atoms().length(); ++i)
	{
		primDis.length(0);
		timesFound.length(0);
		for (j = 0; j < unitISO.atoms()[i].length(); ++j)
		{
			
			// Get new position
			newPos = unitPosToPrimReduced * unitISO.atoms()[i][j].fractional();
			ISO::moveIntoCell(newPos);
			curDis = primISO.basis().distance(newPos, FRACTIONAL, origin, FRACTIONAL);
			
			// Loop over atoms of current element that have already been found
			isNew = true;
			if (j != 0)
			{
				for (k = 0; k < primISO.atoms()[i].length(); ++k)
				{
					
					// Check if positions are the same
					if (Num<double>::abs(curDis - primDis[k]) > tol)
						continue;
					if (primISO.basis().distance(primISO.atoms()[i][k].fractional(), FRACTIONAL, \
						newPos, FRACTIONAL, &nearCell) > tol)
						continue;
					
					// Average values
					for (m = 0; m < 3; ++m)
					{
						newPos[m] = (newPos[m] + nearCell[m] + timesFound[k]*primISO.atoms()[i][k].fractional()[m]) /\
							(timesFound[k] + 1);
						newMagMom[m] = (unitISO.atoms()[i][j].magneticMoment()[m] + \
							timesFound[k] * primISO.atoms()[i][k].magneticMoment()[m]) / (timesFound[k] + 1);
					}
					
					// Save pointer to atom
					_atomMap[i][j] = &primISO.atoms()[i][k];

					// Update atom properties
					primISO.atoms()[i][k].fractional(newPos);
					primISO.atoms()[i][k].magneticMoment(newMagMom);
					timesFound[k]++;
					primDis[k] = primISO.basis().distance(primISO.atoms()[i][k].fractional(), FRACTIONAL, \
						origin, FRACTIONAL);
					for (m = 0; m < 3; ++m)
					{
						if (unitISO.atoms()[i][j].fixed()[m])
							primISO.atoms()[i][k].fixed(m, true);
					}

					// Break since atom was found
					isNew = false;
					break;
				}
			}
			
			// Save if a new atom
			if (isNew)
			{
				_atomMap[i][j] = primISO.addAtom(unitISO.atoms()[i][j]);
				_atomMap[i][j]->fractional(newPos);
				timesFound += 1;
				primDis += primISO.basis().distance(_atomMap[i][j]->fractional(), FRACTIONAL, origin, FRACTIONAL);
			}
		}
	}
	
	// Get the translations that will convert primitive cell atoms back to unit cell
	for (i = 0; i < unitISO.atoms().length(); ++i)
	{
		for (j = 0; j < unitISO.atoms()[i].length(); ++j)
		{
			newPos = _positionConversion * _atomMap[i][j]->fractional();
			for (k = 0; k < 3; ++k)
				_translations[i][j][k] = unitISO.atoms()[i][j].fractional()[k] - newPos[k];
		}
	}
	
	// Return the number of primitive cells in the unit cell
	return _cellConversion.determinant();
}



/* void ReduceISO::expand(ISO& unitISO, const ISO& primISO)
 *
 * Convert from reduced primitive cell
 */

void ReduceISO::expand(ISO& unitISO, const ISO& primISO)
{
	
	// Save new basis
	unitISO.basis(_cellConversion * primISO.basis().vectors(), false);
	
	// Loop over atoms in the unit cell
	int i, j;
	Vector3D newPos;
	for (i = 0; i < unitISO.atoms().length(); ++i)
	{
		for (j = 0; j < unitISO.atoms()[i].length(); ++j)
		{
			
			// Get new position
			newPos = _atomMap[i][j]->fractional();
			newPos *= _positionConversion;
			newPos += _translations[i][j];
			ISO::moveIntoCell(newPos);
			unitISO.atoms()[i][j].fractional(newPos);
		}
	}
}



/* void ReduceISO::expandForces(OList<Vector3D >& totalForces, OList<Vector3D > primForces,
 *		const ISO& unitISO, const ISO& primISO)
 *
 * Convert forces from primitive cell
 */

void ReduceISO::expandForces(OList<Vector3D >& totalForces, OList<Vector3D > primForces, \
	const ISO& unitISO, const ISO& primISO)
{
	
	// Loop over primitive cell forces and convert them to unit cell
	int i;
	for (i = 0; i < primForces.length(); ++i)
	{
		primISO.basis().toFractional(primForces[i]);
		primForces[i] *= _positionConversion;
	}
	
	// Loop over atoms and save the forces
	int j;
	for (i = 0; i < _atomMap.length(); ++i)
	{
		for (j = 0; j < _atomMap[i].length(); ++j)
			totalForces[unitISO.atoms()[i][j].atomNumber()] += primForces[_atomMap[i][j]->atomNumber()];
	}
}
