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



#include "multi.h"
#include "num.h"
#include "ewald.h"
#include "language.h"
#include "output.h"
#include <cstdlib>
#include <cmath>

/**
 *
 * Set Ewald parameters. Expected format of line #1
 *	ewald <element #1 symbol> <element #1 charge> <element #2 symbol> <element #2 charge> <...>
 * 
 * Subsequent options:
 *	accuracy [number] - Accuracy of summation
 *  permittivity [value] - Value of permittivity constant in eV*Angstrom
 *  relative [value] - Permittivity relative to epsilon_0
 */
void Ewald::set(const Text& input)
{
	
	// Finished if empty
	if (!input.length())
		return;
	
	// Output
	Output::newline();
	Output::print("Coulomb potential using Ewald summation");
	Output::increase();
	
	// Get elements and charges from first line
	int i;
	for (i = 1; i < input[0].length(); ++i)
	{
		
		// Break if a comment
		if (Language::isComment(input[0][i]))
			break;
		
		// Found element
		if (Element::isElement(input[0][i], false))
			_elements += Element::find(input[0][i], false);
		
		// Found number
		else if (Language::isNumber(input[0][i]))
			_charges += atof(input[0][i].array());
		
		// Found something else
		else
			readError(input[0]);
	}
	
	// Charges and elements do not match
	if (_elements.length() != _charges.length())
	{
		Output::newline(ERROR);
		Output::print("Number of elements and charges must be equal in Ewald potential");
		Output::quit();
	}
	
	// Set general options
	setEwaldOptions(input, false);
	
	// Print elements and charges
	for (i = 0; i < _elements.length(); ++i)
	{
		Output::newline();
		Output::print("Charge of ");
		Output::print(_elements[i].symbol());
		Output::print(": ");
		Output::print(_charges[i]);
	}
	
	// Print other parameters
	Output::newline();
	Output::print("Relative permitivity: ");
	Output::print(_perm / Constants::eps0);
	Output::newline();
	Output::print("Permitivity: ");
	Output::printSci(_perm);
	Output::print(" eV*Ang");
	Output::newline();
	Output::print("Accuracy: ");
	Output::printSci(_accuracy);
	
	// Output
	Output::decrease();
}

/**
 * Set options for the Ewald summation, given text input.
 * 
 * See documentation of set for those options
 * @param input Options for ewald potential, as read from file
 * @param forgiving Whether to crash on unrecognized option
 */
void Ewald::setEwaldOptions(const Text& input, bool forgiving) {
	// Look for other values
	for (int i = 1; i < input.length(); ++i) {

		// Line is empty
		if (!input[i].length())
			continue;

		// Found a comment
		if (Language::isComment(input[i][0]))
			continue;

		// Line is too short
		if (input[i].length() < 2)
			readError(input[i]);

		// Found accuracy
		if (input[i][0].equal("accuracy", false, 3)) {
			if (Language::isNumber(input[i][1]))
				_accuracy = atof(input[i][1].array());
			else
				readError(input[i]);
		}
			// Found permittivity
		else if (input[i][0].equal("permittivity", false, 4)) {
			{
				if (Language::isNumber(input[i][1]))
					_perm = atof(input[i][1].array());
				else
					readError(input[i]);
			}
		}
			// Found relative permittivity
		else if (input[i][0].equal("relative", false, 3)) {
			{
				if (Language::isNumber(input[i][1]))
					_perm = Constants::eps0 * atof(input[i][1].array());
				else
					readError(input[i]);
			}
		}
			// Anything else
		else if (! forgiving)
			readError(input[i]);
	}
}



/**
 * Compute the Ewald energy
 * @param iso [in] System to evaluate
 * @param totalEnergy [out] Total energy, Ewald energy will be added to this value
 * @param totalForces [out] For on each atom, Ewald force will be added to this value
 */
void Ewald::evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces) const
{
	
	// Initialize the calculation
	initialize(iso, iso.numAtoms());
	
	// Loop over elements
	int i, j;
	int count = 0;
	double localReal = 0;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Loop over atoms of current element and add energy
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			if ((++count + Multi::rank()) % Multi::worldSize() == 0)
			{
				if (totalEnergy)
					localReal += realEnergy(iso, &iso.atoms()[i][j], true);
			}
		}
	}
	
	// Send energy between processors
	if (totalEnergy)
	{
		double temp;
		for (i = 0; i < Multi::worldSize(); ++i)
		{
			temp = localReal;
			Multi::broadcast(temp, i);
			*totalEnergy += temp;
		}
		*totalEnergy += recipEnergy(iso) - selfEnergy(iso) - chargedEnergy(iso);
	}
	
	if (totalForces) {
		computeForces(iso, totalForces);	
	}
}

/** 
 * Compute the Ewald energy using symmetry
 */
void Ewald::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, \
	OList<Vector3D >* totalForces) const
{
	
	// Do not use symmetry unless if reduces the number of atoms by at least a factor of two
	if (iso.numAtoms() / symmetry.orbits().length() < 2) {
		evaluate(iso, totalEnergy, totalForces);
		return;
	}
	
	// Initialize the calculation
	initialize(iso, symmetry.orbits().length());
	
	// Loop over unique atoms and add energy of each
	int i;
	int count = 0;
	double localReal = 0;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		if ((++count + Multi::rank()) % Multi::worldSize() == 0)
		{
			if (totalEnergy)
				localReal += symmetry.orbits()[i].atoms().length() * \
					realEnergy(iso, symmetry.orbits()[i].atoms()[0], false);
		}
	}
	
	// Send energy between processors
	if (totalEnergy)
	{
		double temp;
		for (i = 0; i < Multi::worldSize(); ++i)
		{
			temp = localReal;
			Multi::broadcast(temp, i);
			*totalEnergy += temp/2;
		}
		*totalEnergy += recipEnergy(iso) - selfEnergy(iso) - chargedEnergy(iso);
	}
	
	if (totalForces) {
		computeForces(iso, totalForces);
	}
}

/**
 * Compute the forces on each atom using 
 * @param iso [in] Structure being evaluated
 * @param totalForces [in/out] Forces on each atom, will be added to
 */
void Ewald::computeForces(const ISO& iso, OList<Vector3D>* totalForces) const {
	for (int e=0; e<iso.atoms().length(); e++) {
		for (int i=0; i<iso.atoms()[e].length(); i++) {
			int id = iso.atoms()[e][i].atomNumber();
			(*totalForces)[id] += realForce(iso, &iso.atoms()[e][i]) 
					+ recipForce(iso, &iso.atoms()[e][i]);
		}
	}
}

/**
 * Precompute parameters needed for the Ewald sum, such as:
 *  - Determining optimal mixing parameter (alpha)
 *  - Computing cutoffs in real/reciprocal space for desired accuracy
 *  - Creating an iterator to move over images of atoms in real space
 *  - Determining which reciprocal space vectors to use 
 *  - Computing prefactors used for reciprocal space summation
 */
void Ewald::initialize(const ISO& iso, int numUniqueAtoms) const
{
	
	// Set alpha value
	double w = 0.05;
	if (numUniqueAtoms > 400)
		w = 0.2;
	else if (numUniqueAtoms > 100)
		w = .2 + .15*(numUniqueAtoms - 400)/300;
	_alpha = sqrt(Constants::pi) * pow(w * iso.numAtoms(), 1.0/6.0) / pow(iso.basis().volume(), 1.0/3.0);
	
	// Set cutoffs
	double sqrtLogAcc = sqrt(-log(_accuracy));
    double realCut = sqrtLogAcc / _alpha;
    double recipCut = 2 * _alpha * sqrtLogAcc;
	
	// Set the real space image iterator
	_realIterator.setCell(iso.basis(), realCut);
	
	// Save the reciprocal space lattice vectors
	Matrix3D recipVecs = iso.basis().inverseTranspose() * (2 * Constants::pi);
	Basis recipBasis(recipVecs, false);
	
	// Set the reciprocal space lattice image iterator
	ImageIterator recipIterator;
	recipIterator.setCell(recipBasis, recipCut);
	recipIterator.reset(Vector3D(0.0), Vector3D(0.0));	
	
	// Save the reciprocal space lattice vectors with non-zero length
	_recipVectors.clear();
	while (!recipIterator.finished())
	{
		if (++recipIterator > 1e-8)
			_recipVectors.add(recipIterator.cellVector() * (2 * Constants::pi));
	}
	
	// Save the prefactors needed in reciprocal space energy evaluation
	int i = 0;
	double magSquared;
	Vector3D tempVec;
	_recipFactors.length(_recipVectors.length());
	for (Linked<Vector3D >::iterator it = _recipVectors.begin(); it != _recipVectors.end(); ++it, ++i)
	{
		tempVec = iso.basis().inverse() * *it;
		magSquared = tempVec * tempVec;
		_recipFactors[i] = exp(-magSquared / (4 * _alpha*_alpha)) / (_perm * iso.basis().volume() * magSquared);
	}
}



/**
 * Compute the real-space contribution to the Ewald energy from a single atom
 * @param iso [in] Structure being evaluated
 * @param atom [in] Pointer to atom being considered
 * @param skipLowerAtoms [in] Whether to compute contributions from interactions 
 *  with atoms with a lower index (false when symmetry is used)
 */
double Ewald::realEnergy(const ISO& iso, Atom* atom, bool skipLowerAtoms) const
{
	
	// Get the charge of the current atom
	double atomCharge = getCharge(atom->element());
	if (atomCharge == 0)
		return 0;
	
	// Evaluate the real space term
	int i, j;
	double real = 0;
	double curCharge;
	double tempEnergy;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Get the charge of the current element
		curCharge = getCharge(iso.atoms()[i][0].element());
		if (curCharge == 0)
			continue;
		
		// Loop over atoms of current element
		tempEnergy = 0;
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Skip if atom has already been counted
			if (skipLowerAtoms)
			{
				if (iso.atoms()[i][j].atomNumber() < atom->atomNumber())
					continue;
			}
			
			// Set the image iterator for current pair
			_realIterator.reset(atom->fractional(), iso.atoms()[i][j].fractional());
			
			// Getting energy of atom with itself and not using symmetry
			if ((atom->atomNumber() == iso.atoms()[i][j].atomNumber()) && (skipLowerAtoms))
			{
				while (!_realIterator.finished())
				{
					if (++_realIterator > 1e-8)
						tempEnergy += realEnergy(_realIterator.distance()) / 2;
				}
			}
			
			// Getting energy of different atoms or using symmetry
			else
			{
				while (!_realIterator.finished())
				{
					if (++_realIterator > 1e-8)
						tempEnergy += realEnergy(_realIterator.distance());
				}
			}
		}
		
		// Add energy for previous set of atoms
		real += atomCharge * curCharge * tempEnergy;
	}
	
	// Return the real energy
	return real / (4 * Constants::pi * _perm);
}

/**
 * Compute the real-space contribution to force on a single atom
 * @param iso [in] Structure being evaluated
 * @param atom [in] Atom on which force is active
 * @return Force in each directory
 */
Vector3D Ewald::realForce(const ISO& iso, Atom* atom) const {
	Vector3D force(0.0);
	// Get the charge of the current atom
	double atomCharge = getCharge(atom->element());
	if (atomCharge == 0)
		return force;
	
	// Compute useful prefactors
	double twoAoverRootPi = 2 * _alpha / sqrt(Constants::pi);
	double alphaSquared = _alpha * _alpha;
	
	// Loop over every element type
	for (int e=0; e < iso.atoms().length(); e++) { 
		// Get the charge of this atom
		double curCharge = getCharge(iso.atoms()[e][0].element());
		if (curCharge == 0)
			continue;
		
		// Loop over all atoms of that element
		for (int j=0; j < iso.atoms()[e].length(); j++) {
			
			_realIterator.reset(atom->fractional(), iso.atoms()[e][j].fractional());
			// Loop over all images of that atom
			while (!_realIterator.finished()) {
				if (++_realIterator < 1e-8) continue; 
				
				double distance = _realIterator.distance();
				double mag = erfc(_alpha * distance) / distance
					+ twoAoverRootPi * exp(-1 * alphaSquared * distance * distance);
				force += _realIterator.cartVector() * mag * curCharge / distance / distance;
			}
		}
	}
	
	// Return the real energy
	return force * (atomCharge / 4 / Constants::pi / _perm);
}

/**
 * Compute the reciprocal space term of the Ewald summation
 * 
 */
double Ewald::recipEnergy(const ISO& iso) const {
	
	// Loop over reciprocal space lattice vectors
	int i, j, k;
	int count = 0;
	double dot;
	double charge;
	double cosTerm;
	double sinTerm;
	double localRecip = 0;
	Linked<Vector3D >::iterator itVector = _recipVectors.begin();
	for (i = 0; itVector != _recipVectors.end(); ++itVector, ++i)
	{
		
		// Check if running on current processor
		if ((++count + Multi::rank()) % Multi::worldSize() != 0)
			continue;
		
		// Loop over elements
		cosTerm = 0;
		sinTerm = 0;
		for (j = 0; j < iso.atoms().length(); ++j)
		{
			
			// Get charge of current atom
			charge = getCharge(iso.atoms()[j][0].element());
			if (charge == 0)
				continue;
			
			// Loop over atoms of current element
			for (k = 0; k < iso.atoms()[j].length(); ++k)
			{
				
				// Evaluate current atom contribution
				dot = *itVector * iso.atoms()[j][k].fractional();
				cosTerm += charge * cos(dot);
				sinTerm += charge * sin(dot);
			}
		}
		
		// Add energy
		localRecip += _recipFactors[i] * (cosTerm*cosTerm + sinTerm*sinTerm);
	}
	
	// Send energy between processors
	double temp;
	double recip = 0;
	for (i = 0; i < Multi::worldSize(); ++i)
	{
		temp = localRecip;
		Multi::broadcast(temp, i);
		recip += temp;
	}

	// Return energy
	return recip / 2;
}

/**
 * Compute the reciprocal space contribution to Ewald force on an atom
 * @param iso [in] Structure being evaluated
 * @param atom [in] Atom on which force is acting
 * @return Force acting on that atom
 */
Vector3D Ewald::recipForce(const ISO& iso, Atom* atom) const {
	Vector3D force(0.0);
	// Get the charge of the current atom
	double atomCharge = getCharge(atom->element());
	if (atomCharge == 0)
		return force;
	
	// Loop over reciprocal space vectors
	Linked<Vector3D >::iterator itVector = _recipVectors.begin();
	for (int i = 0; itVector != _recipVectors.end(); ++itVector, ++i) {
		// Get dot product between this reciprocal lattice factor and the atom
		double atomDot = *itVector * atom->fractional();
		
		// Get the current reciprocal space vector in Cartesian units
		Vector3D recipVector = iso.basis().inverse() * *itVector;
		
		double cosTerm = 0, sinTerm = 0;
		// Loop over all atom types
		for (int e=0; e < iso.atoms().length(); e++) {
			double curCharge = getCharge(iso.atoms()[e][0].element());
			
			// Loop over all atoms of that type
			for (int j=0; j < iso.atoms()[e].length(); j++) {
				double dot = *itVector * iso.atoms()[e][j].fractional();
				
				cosTerm += curCharge * cos(dot);
				sinTerm += curCharge * sin(dot);
			}
		}
		
		// Add contribution from this k-point to total forces
		
		force += recipVector * (atomCharge * _recipFactors[i] * (sin(atomDot) * cosTerm +
				cos(atomDot) * sinTerm));
	}
	
	return force;
}



/**
 * Compute the self energy term of the Ewald summation
 * @param iso [in] Structure from which to compute energies
 * @return Self energy of system
 */
double Ewald::selfEnergy(const ISO& iso) const
{
	
	// Get the total of all squared charges
	int i, j;
	double totalChargeSquared = 0;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Look for current element
		for (j = 0; j < _elements.length(); ++j)
		{
			
			// Found element
			if (iso.atoms()[i][0].element() == _elements[j])
			{
				totalChargeSquared += iso.atoms()[i].length() * _charges[j] * _charges[j];
				break;
			}
		}
	}
	
	// Return the self energy
	return totalChargeSquared * _alpha / (4 * pow(Constants::pi, 3.0/2.0) * _perm);
}



/* double Ewald::chargedEnergy(const ISO& iso) const
 *
 * Return the charged cell energy correction
 */

double Ewald::chargedEnergy(const ISO& iso) const
{
	
	// Get the total of all charges
	int i, j;
	double totalCharge = 0;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Look for current element
		for (j = 0; j < _elements.length(); ++j)
		{
			
			// Found element
			if (iso.atoms()[i][0].element() == _elements[j])
			{
				totalCharge += iso.atoms()[i].length() * _charges[j];
				break;
			}
		}
	}
	
	// Return the charged cell energy correction
	return totalCharge * totalCharge / (8 * _perm * iso.basis().volume() * _alpha*_alpha);
}



/* double Ewald::getCharge(const Element& element) const
 *
 * Get the charge for an element
 */

double Ewald::getCharge(const Element& element) const
{
	for (int i = 0; i < _elements.length(); ++i)
	{
		if (_elements[i] == element)
			return _charges[i];
	}
	return 0;
}
