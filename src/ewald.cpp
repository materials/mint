/* ewald.cpp -- Ewald summation
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "multi.h"
#include "num.h"
#include "ewald.h"
#include "language.h"
#include "output.h"
#include <cstdlib>



/* void Ewald::set(const Text& input)
 *
 * Set Ewald parameters
 */

void Ewald::set(const Text& input)
{
	
	// Finished if emtpy
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
	
	// Look for other values
	for (i = 1; i < input.length(); ++i)
	{
		
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
		if (input[i][0].equal("accuracy", false, 3))
		{
			if (Language::isNumber(input[i][1]))
				_accuracy = atof(input[i][1].array());
			else
				readError(input[i]);
		}
		
		// Found permitivity
		else if (input[i][0].equal("permitivity", false, 4))
		{
			{
				if (Language::isNumber(input[i][1]))
					_perm = atof(input[i][1].array());
				else
					readError(input[i]);
			}
		}
		
		// Found relative permitivity
		else if (input[i][0].equal("relative", false, 3))
		{
			{
				if (Language::isNumber(input[i][1]))
					_perm = Constants::eps0 * atof(input[i][1].array());
				else
					readError(input[i]);
			}
		}
		
		// Anything else
		else
			readError(input[i]);
	}
	
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



/* void Ewald::evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces) const
 *
 * Return the Ewald energy
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
}



/* void Ewald::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy,
 *		OList<Vector3D >* totalForces) const
 *
 * Return the ewald energy using symmetry
 */

void Ewald::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, \
	OList<Vector3D >* totalForces) const
{
	
	// Do not use symmetry unless if reduces the number of atoms by at least a factor of two
	if (iso.numAtoms() / symmetry.orbits().length() < 2)
	{
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
}



/* void Ewald::initialize(const ISO& iso, int numUniqueAtoms) const
 *
 * Set the cells to iterate over in real and reciprocal space for Ewald potential
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



/* double Ewald::realEnergy(const ISO& iso, Atom* atom, bool skipLowerAtoms) const
 *
 * Return the energy of an atom in Ewald summation
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



/* double Ewald::recipEnergy(const ISO& iso) const
 *
 * Return the reciprocal space energy
 */

double Ewald::recipEnergy(const ISO& iso) const
{
	
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



/* double Ewald::selfEnergy(const ISO& iso) const
 *
 * Return the self energy
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
