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



#include "multi.h"
#include "pairPotential.h"
#include "language.h"
#include "output.h"
#include <cstdlib>



/* void PairPotential::evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces) const
 *
 * Return the energy of a structure
 */

void PairPotential::evaluate(const ISO& iso, double* totalEnergy, OList<Vector3D >* totalForces) const
{
	
	// Set image iterator
	_images.setCell(iso.basis(), _cutoff);
	
	// Set which elements should be updated
	OList<Element>::D2 elements(1);
	elements[0] += _element1;
	elements[0] += _element2;
	if (_element1 != _element2)
	{
		elements.length(2);
		elements[1] += _element2;
		elements[1] += _element1;
	}
	
	// Variable to store forces
	OList<Vector3D > localForces;
	if (totalForces)
	{
		localForces.length(totalForces->length());
		localForces.fill(0.0);
	}
	
	// Loop over element pairs
	int i, j, k;
	int count = 0;
	double localTotal = 0;
	for (i = 0; i < elements.length(); ++i)
	{
		
		// Loop over elements
		for (j = 0; j < iso.atoms().length(); ++j)
		{

			// Found element
			if (iso.atoms()[j][0].element() == elements[i][0])
			{

				// Loop over atoms of current element and add energy and get force
				for (k = 0; k < iso.atoms()[j].length(); ++k)
				{
					if ((++count + Multi::rank()) % Multi::worldSize() == 0)
					{
						if (totalEnergy)
							localTotal += energy(iso, &iso.atoms()[j][k], elements[i][0], elements[i][1], true);
						if (totalForces)
							localForces[iso.atoms()[j][k].atomNumber()] = \
								force(iso, &iso.atoms()[j][k], elements[i][0], elements[i][1]);
					}
				}

				// Add tail if needed
				if ((_addTail) && (totalEnergy))
					localTotal += iso.atoms()[j].length() * tail(iso, elements[i][1]) / Multi::worldSize();

				// Shift energy if needed
				if ((_shift) && (!_addTail) && (totalEnergy))
					localTotal -= iso.atoms()[j].length() * pairEnergy(_cutoff) / (2 * Multi::worldSize());

				// Finished since element was found
				break;
			}
		}
	}
	
	// Send energy between processors
	if (totalEnergy)
	{
		double temp;
		for (i = 0; i < Multi::worldSize(); ++i)
		{
			temp = localTotal;
			Multi::broadcast(temp, i);
			*totalEnergy += temp;
		}
	}
	
	// Send forces between processors
	if (totalForces)
	{
		Vector3D temp;
		for (i = 0; i < localForces.length(); ++i)
		{
			for (j = 0; j < Multi::worldSize(); ++j)
			{
				temp = localForces[i];
				Multi::broadcast(temp, j);
				(*totalForces)[i] += temp;
			}
		}
	}
}



/* void PairPotential::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy,
 *		OList<Vector3D >* totalForces) const
 *
 * Return the energy of a structure
 */

void PairPotential::evaluate(const ISO& iso, const Symmetry& symmetry, double* totalEnergy, \
	OList<Vector3D >* totalForces) const
{
	
	// Do not use symmetry unless if reduces the number of atoms by at least a factor of two
	if (iso.numAtoms() / symmetry.orbits().length() < 2)
	{
		evaluate(iso, totalEnergy, totalForces);
		return;
	}
	
	// Set image iterator
	_images.setCell(iso.basis(), _cutoff);
	
	// Set which elements should be updated
	OList<Element>::D2 elements(1);
	elements[0] += _element1;
	elements[0] += _element2;
	if (_element1 != _element2)
	{
		elements.length(2);
		elements[1] += _element2;
		elements[1] += _element1;
	}
	
	// Variable to store forces
	OList<Vector3D > localForces;
	if (totalForces)
	{
		localForces.length(totalForces->length());
		localForces.fill(0.0);
	}
	
	// Loop over unique atoms
	int i, j;
	int count = 0;
	double localTotal = 0;
	for (i = 0; i < elements.length(); ++i)
	{
		for (j = 0; j < symmetry.orbits().length(); ++j)
		{

			// Skip if element is not corrent
			if (symmetry.orbits()[j].atoms()[0]->element() != elements[i][0])
				continue;
			
			// Check if adding on current processor
			if ((++count + Multi::rank()) % Multi::worldSize() != 0)
				continue;
			
			// Add energy
			if (totalEnergy)
				localTotal += symmetry.orbits()[j].atoms().length() * \
					energy(iso, symmetry.orbits()[j].atoms()[0], elements[i][0], elements[i][1], false) / 2;
			
			// Add tail if needed
			if ((_addTail) && (totalEnergy))
				localTotal += symmetry.orbits()[j].atoms().length() * tail(iso, elements[i][1]);

			// Shift energy if needed
			if ((_shift) && (!_addTail) && (totalEnergy))
				localTotal -= symmetry.orbits()[j].atoms().length() * pairEnergy(_cutoff) / 2;
			
			// Get force
			if (totalForces)
				localForces[symmetry.orbits()[j].atoms()[0]->atomNumber()] = \
					symmetry.orbits()[j].specialPositions()[0].rotation() * \
						force(iso, symmetry.orbits()[j].atoms()[0], elements[i][0], elements[i][1]);
		}
	}
	
	// Send energy between processors
	if (totalEnergy)
	{
		double temp;
		for (i = 0; i < Multi::worldSize(); ++i)
		{
			temp = localTotal;
			Multi::broadcast(temp, i);
			*totalEnergy += temp;
		}
	}
	
	// Send forces between processors
	if (totalForces)
	{
		
		// Send forces
		Vector3D temp;
		for (i = 0; i < localForces.length(); ++i)
		{
			for (j = 0; j < Multi::worldSize(); ++j)
			{
				temp = localForces[i];
				Multi::broadcast(temp, j);
				(*totalForces)[i] += temp;
			}
		}
		
		// Apply symmetry operations to generate forces on equivalent atoms
		for (i = 0; i < symmetry.orbits().length(); ++i)
		{
			for (j = 1; j < symmetry.orbits()[i].atoms().length(); ++j)
				(*totalForces)[symmetry.orbits()[i].atoms()[j]->atomNumber()] = \
					symmetry.orbits()[i].generators()[j].rotation() * \
					(*totalForces)[symmetry.orbits()[i].atoms()[0]->atomNumber()];
		}
	}
}



/* double PairPotential::energy(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2, 
 *		bool skipLowerAtoms) const
 *
 * Return the energy of an atom in a structure
 */

double PairPotential::energy(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2, \
	bool skipLowerAtoms) const
{
	
	// Return if atom is not of correct type
	if (atom->element() != elem1)
		return 0;
	
	// Loop over elements in the structure
	int i, j;
	double res = 0;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Found second element
		if (iso.atoms()[i][0].element() == elem2)
		{
			
			// Loop over atoms of second element
			for (j = 0; j < iso.atoms()[i].length(); ++j)
			{
				
				// Skip if atom has already been counted
				if (skipLowerAtoms)
				{
					if (iso.atoms()[i][j].atomNumber() < atom->atomNumber())
						continue;
				}
				
				// Set image iterator for current atoms
				_images.reset(atom->fractional(), iso.atoms()[i][j].fractional());
				
				// Getting energy of atom with itself and not using symmetry
				if ((atom->atomNumber() == iso.atoms()[i][j].atomNumber()) && (skipLowerAtoms))
				{
					while (!_images.finished())
					{
						if (++_images > 1e-8)
						 	res += pairEnergy(_images.distance()) / 2;
					}
				}
				
				// Getting energy of different atoms or using symmetry
				else
				{
					while (!_images.finished())
					{
						if (++_images > 1e-8)
							res += pairEnergy(_images.distance());
					}
				}
			}
			
			// Break since element was found
			break;
		}
	}
	
	// Return energy
	return res;
}



/* Vector3D PairPotential::force(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2) const
 *
 * Return the force on an atom
 */

Vector3D PairPotential::force(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2) const
{
	
	// Return if atom is not of correct type
	if (atom->element() != elem1)
		return Vector3D(0.0);
	
	// Loop over elements in the structure
	int i, j;
	Vector3D temp;
	Vector3D res(0.0);
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Found second element
		if (iso.atoms()[i][0].element() == elem2)
		{
			
			// Loop over atoms of second element
			for (j = 0; j < iso.atoms()[i].length(); ++j)
			{
				
				// Set image iterator for current atoms
				_images.reset(atom->fractional(), iso.atoms()[i][j].fractional());
				
				// Get force
				while (!_images.finished())
				{
					if (++_images > 1e-8)
					{
						temp = _images.cartVector();
						temp /= _images.distance();
						temp *= pairForce(_images.distance());
						iso.basis().toFractional(temp);
						res -= temp;
					}
				}
			}
			
			// Break since element was found
			break;
		}
	}
	
	// Return energy
	return res;
}



/* void PairPotential::set(const Text& input)
 *
 * Set pair potential values
 */

void PairPotential::set(const Text& input)
{
	
	// Check that first line contains enough data
	int i;
	if (!input.length())
		return;
	if (input[0].length() < 3)
		readError(input[0]);
	
	// Get elements
	_element1 = Element::find(input[0][1], true, true);
	_element2 = Element::find(input[0][2], true, true);
	
	// Loop over lines
	for (i = 1; i < input.length(); ++i)
	{
		
		// Skip if empty
		if (!input[i].length())
			continue;
		
		// Found a comment
		if (Language::isComment(input[i][0]))
			continue;
		
		// Line is too short
		if (input[i].length() < 2)
			readError(input[i]);
		
		// Found add tail
		if (input[i][0].equal("tail", false))
		{
			if (input[i][1].equal("true", false, 1))
				_addTail = true;
			else if (input[i][1].equal("false", false, 1))
				_addTail = false;
			else
				readError(input[i]);
		}
		
		// Found shift
		else if ((input[i][0].equal("shift", false)) || (input[i][0].equal("zero", false)))
		{
			if (input[i][1].equal("true", false, 1))
				_shift = true;
			else if (input[i][1].equal("false", false, 1))
				_shift = false;
			else
				readError(input[i]);
		}
		
		// Found cutoff
		else if (input[i][0].equal("cutoff", false, 3))
		{
			if (Language::isNumber(input[i][1]))
				_cutoff = atof(input[i][1].array());
			else
				readError(input[i]);
		}
		
		// Anything else
		else
			readError(input[i]);
	}
}



/* void LennardJones::set(const Text& input)
 *
 * Set Lennard jones function
 */

void LennardJones::set(const Text& input)
{
	
	// No data
	if (!input.length())
		return;

	// Get standard data
	PairPotential::set(input);
	
	// Output
	Output::newline();
	Output::print("Lennard-Jones potential between ");
	Output::print(_element1.symbol());
	Output::print(" and ");
	Output::print(_element2.symbol());
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 5)
		readError(input[0]);
	
	// Get epsilon
	if (!Language::isNumber(input[0][3]))
		readError(input[0]);
	_eps = atof(input[0][3].array());
	
	// Get sigma
	if (!Language::isNumber(input[0][4]))
		readError(input[1]);
	_sig = atof(input[0][4].array());
	
	// Set cutoff
	if (_cutoff == -1)
		_cutoff = 2.5 * _sig;
	
	// Output
	Output::newline();
	Output::print("Epsilon: ");
	Output::print(_eps);
	Output::print(" eV");
	Output::newline();
	Output::print("Sigma: ");
	Output::print(_sig);
	Output::print(" Ang");
	PairPotential::print();
	
	// Output
	Output::decrease();
}



/* void Buckingham::set(const Text& input)
 *
 * Set Buckingham function
 */

void Buckingham::set(const Text& input)
{
	
	// No data
	if (!input.length())
		return;

	// Get standard data
	PairPotential::set(input);
	
	// Output
	Output::newline();
	Output::print("Buckingham potential between ");
	Output::print(_element1.symbol());
	Output::print(" and ");
	Output::print(_element2.symbol());
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 6)
		readError(input[0]);
	
	// Get A
	if (!Language::isNumber(input[0][3]))
		readError(input[0]);
	_A = atof(input[0][3].array());
	
	// Get rho
	if (!Language::isNumber(input[0][4]))
		readError(input[1]);
	_rho = atof(input[0][4].array());
	
	// Get C
	if (!Language::isNumber(input[0][5]))
		readError(input[1]);
	_C = atof(input[0][5].array());
	
	// Set cutoff
	if (_cutoff == -1)
		_cutoff = 2.5 / _rho;
	
	// Output
	Output::newline();
	Output::print("A: ");
	Output::print(_A);
	Output::print(" eV");
	Output::newline();
	Output::print("Rho: ");
	Output::print(_rho);
	Output::print(" Ang");
	Output::newline();
	Output::print("C: ");
	Output::print(_C);
	Output::print(" eV*Ang^6");
	PairPotential::print();
	
	// Output
	Output::decrease();
}



/* void Power::set(const Text& input)
 *
 * Set Power function
 */

void Power::set(const Text& input)
{
	
	// No data
	if (!input.length())
		return;

	// Get standard data
	PairPotential::set(input);
	
	// Output
	Output::newline();
	Output::print("Power potential between ");
	Output::print(_element1.symbol());
	Output::print(" and ");
	Output::print(_element2.symbol());
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 6)
		readError(input[0]);
	
	// Get power
	if (!Language::isInteger(input[0][3]))
		readError(input[0]);
	_power = atoi(input[0][3].array());
	
	// Get epsilson
	if (!Language::isNumber(input[0][4]))
		readError(input[1]);
	_eps = atof(input[0][4].array());
	
	// Get sigma
	if (!Language::isNumber(input[0][5]))
		readError(input[1]);
	_sig = atof(input[0][5].array());
	
	// Set cutoff
	if (_cutoff == -1)
		_cutoff = 2.5 * _sig;
	
	// Output
	Output::newline();
	Output::print("Power: ");
	Output::print(_power);
	Output::newline();
	Output::print("Epsilon: ");
	Output::print(_eps);
	Output::print(" eV");
	Output::newline();
	Output::print("Sigma: ");
	Output::print(_sig);
	Output::print(" Ang");
	PairPotential::print();
	
	// Output
	Output::decrease();
}



/* void Exponential::set(const Text& input)
 *
 * Set Exponential function
 */

void Exponential::set(const Text& input)
{
	
	// No data
	if (!input.length())
		return;

	// Get standard data
	PairPotential::set(input);
	
	// Output
	Output::newline();
	Output::print("Exponential potential between ");
	Output::print(_element1.symbol());
	Output::print(" and ");
	Output::print(_element2.symbol());
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 5)
		readError(input[0]);
	
	// Get epsilon
	if (!Language::isNumber(input[0][3]))
		readError(input[0]);
	_eps = atof(input[0][3].array());
	
	// Get rho
	if (!Language::isNumber(input[0][4]))
		readError(input[1]);
	_rho = atof(input[0][4].array());
	
	// Set cutoff
	if (_cutoff == -1)
		_cutoff = 2.5 / _rho;
	
	// Output
	Output::newline();
	Output::print("Epsilon: ");
	Output::print(_eps);
	Output::print(" eV");
	Output::newline();
	Output::print("Rho: ");
	Output::print(_rho);
	Output::print(" Ang");
	PairPotential::print();
	
	// Output
	Output::decrease();
}



/* void Covalent::set(const Text& input)
 *
 * Set Convalent function
 */

void Covalent::set(const Text& input)
{
	
	// No data
	if (!input.length())
		return;

	// Get standard data
	PairPotential::set(input);
	
	// Output
	Output::newline();
	Output::print("Covalent potential between ");
	Output::print(_element1.symbol());
	Output::print(" and ");
	Output::print(_element2.symbol());
	Output::increase();
	
	// Line is too short
	if (input[0].length() < 5)
		readError(input[0]);
	
	// Get epsilon
	if (!Language::isNumber(input[0][3]))
		readError(input[0]);
	_eps = atof(input[0][3].array());
	
	// Get sigma
	if (!Language::isNumber(input[0][4]))
		readError(input[1]);
	_sig = atof(input[0][4].array());
	
	// Set cutoff
	if (_cutoff == -1)
		_cutoff = 2.5 / _sig;
	
	// Set ideal distance
	_idealDistance = _element1.radius() + _element2.radius();
	
	// Output
	Output::newline();
	Output::print("Epsilon: ");
	Output::print(_eps);
	Output::print(" eV");
	Output::newline();
	Output::print("Sigma: ");
	Output::print(_sig);
	Output::print(" Ang");
	PairPotential::print();
	
	// Output
	Output::decrease();
}



/* void PairPotential::print()
 *
 * Print pair potential data
 */

void PairPotential::print()
{
	Output::newline();
	Output::print("Cutoff: ");
	Output::print(_cutoff);
	Output::print(" Ang");
	Output::newline();
	Output::print("Add tail: ");
	if (_addTail)
		Output::print("True");
	else
		Output::print("False");
	Output::newline();
	Output::print("Shift to zero: ");
	if (_shift)
		Output::print("True");
	else
		Output::print("False");
}
