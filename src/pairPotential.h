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



#ifndef PAIRPOTENTIAL_H
#define PAIRPOTENTIAL_H



#include "locPotential.h"
#include "iso.h"
#include "elements.h"
#include "symmetry.h"
#include "text.h"
#include "constants.h"
#include "num.h"
#include <cmath>



// Pair potential object
class PairPotential : public SingleLocalPotential
{
	
protected:
	
	// Variables
	bool _addTail;
	bool _shift;
    double _cutoff;
	Element _element1;
	Element _element2;
	mutable ImageIterator _images;
	
	// Functions
	double energy(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2, bool skipLowerAtoms) const;
	Vector3D force(const ISO& iso, Atom* atom, const Element& elem1, const Element& elem2) const;
	double density(const ISO& iso, const Element& elem2) const;
	void print();
	
	// Virtual functions
	/** Compute energy as a function of distance */
	virtual double pairEnergy(double distance) const = 0;
	/** Compute force as a function of distance */
	virtual double pairForce(double distance) const = 0;
	/** Compute "tail" energy based on density of element within a solid */
	virtual double tail(const ISO& iso, const Element& elem2) const = 0;
	
public:
	
	// Constructor
	PairPotential()	{ _addTail = false; _shift = true; _cutoff = -1; }
	
	// Setup by file input
	virtual void set(const Text& input);
	
	virtual void setElementOne(const Element& input) { _element1 = input; }
	virtual void setElementTwo(const Element& input) { _element2 = input; }
	
	// Evaluation functions
	void evaluate(const ISO& iso, double* totalEnergy = 0, OList<Vector3D >* totalForces = 0) const;
	void evaluate(const ISO& iso, const Symmetry& symmetry, double* energy = 0, OList<Vector3D>* forces = 0) const;
};



/**
 * Lennard jones potential
 * 
 * E(r) = 4 * eps * ((sig / R)^12 - (sig / R)^6)
 * 
 * Input format: eps sigma
 */
class LennardJones : public PairPotential
{
	
	// Variables
	double _eps;
	double _sig;
	
public:
	
	// Constructor
	LennardJones() : PairPotential() {}
	
	// Setup from file
	void set(const Text& input);
	
	// Evaluation functions
	double pairEnergy(double distance) const
		{ return 4 * _eps * (pow(_sig / distance, 12) - pow(_sig / distance, 6)); }
	double pairForce(double distance) const
		{ return 24 * _eps * (2 * pow(_sig / distance, 12) - pow(_sig / distance, 6)) / distance; }
	double tail(const ISO& iso, const Element& elem2) const
		{
			return 8 * Constants::pi * _eps * density(iso, elem2) * \
				(pow(_sig, 12) - 3 * pow(_sig * _cutoff, 6)) / (9 * pow(_cutoff, 9));
		}
};



/**
 * Buckingham potential
 * 
 * E(r) = A * exp(-R/rho) - C * R ^-6
 * 
 * Input format: A rho C
 */
class Buckingham : public PairPotential
{
	
	// Variables
	double _A;
	double _rho;
	double _C;

public:
	
	// Constructor
	Buckingham() : PairPotential() {}

	// Setup from file
	void set(const Text& input);

	// Evaluation functions
	double pairEnergy(double distance) const
		{ return _A * exp(-distance / _rho) - _C / pow(distance, 6); }
	double pairForce(double distance) const
		{ return _A * exp(-distance / _rho) / _rho - 6 * _C / pow(distance, 7); }
	double tail(const ISO& iso, const Element& elem2) const
		{
			return -2 * _C * Constants::pi * density(iso, elem2) / (3 * pow(_cutoff, 3)) + 2 * _A * \
				exp(-_cutoff / _rho) * Constants::pi * density(iso, elem2) * (2 + _cutoff * pow(_rho, 3) * \
				(2 + _cutoff / _rho) / _rho);
		}		
};



/**
 * Power potential. E(r) = eps * (sigma / distance) ^ power
 * 
 * Format in input file: 
 *	&lt;element #1&gt; &lt;element #2&gt; &lt;power&gt; &lt;epsilon&gt; &lt;sigma&gt;
 *  {options consistant with pair potential}
 * 
 */
class Power : public PairPotential
{
	
	// Variables
	double _power;
	double _eps;
	double _sig;

public:
	
	// Constructor
	Power() : PairPotential() {}

	// Setup from file
	void set(const Text& input);

	// Evaluation functions
	double pairEnergy(double distance) const
		{ return _eps * pow(_sig / distance, _power); }
	double pairForce(double distance) const
		{ return _eps * pow(_sig / distance, _power) / distance; }
	double tail(const ISO& iso, const Element& elem2) const
		{
			return 2 * Constants::pi * pow(_cutoff, 3 - _power) * _eps * density(iso, elem2) * \
				pow(_sig, _power) / (_power - 3);
		}
};



/**
 *  Exponential potential
 * 
 * E(R) = A * exp(-R/rho)
 * 
 * Input format: A rho
 */
class Exponential : public PairPotential
{
	
	// Variables
	double _eps;
	double _rho;

public:
	
	// Constructor
	Exponential() : PairPotential() {}

	// Setup from file
	void set(const Text& input);

	// Evaluation functions
	double pairEnergy(double distance) const
		{ return _eps * exp(-distance / _rho); }
	double pairForce(double distance) const
		{ return _eps * exp(-distance / _rho) / _rho; }
	double tail(const ISO& iso, const Element& elem2) const
		{
			return 2 * Constants::pi * _eps * density(iso, elem2) * _rho * exp(-_cutoff / _rho) * \
				(_cutoff * _cutoff + 2 * _cutoff * _rho + 2 * _rho * _rho);
		}
};



/**
 * Covalent potential
 * 
 * E(r) = -A * exp(-B * (r - R*)^2/2R)
 * 
 * Input Format: A B
 */
class Covalent : public PairPotential
{
	
	// Variables
	double _eps;
	double _sig;
	double _idealDistance;

public:
	
	// Constructor
	Covalent() : PairPotential() {}
	
	// Setup from file	
	void set(const Text& input);
	
	// Evaluation functions
	double pairEnergy(double distance) const
		{ return -_eps * exp(-_sig * (distance - _idealDistance) * (distance - _idealDistance) / (2 * distance)); }
	double pairForce(double distance) const
		{ 
			return -_eps * _sig * (distance + _idealDistance) * (distance - _idealDistance) * \
				exp(-_sig * (distance - _idealDistance) * (distance - _idealDistance) / (2 * distance)) /\
				(2 * distance * distance);
		}
	double tail(const ISO& iso, const Element& elem2) const
		{ return 0; }
};


/**
 * Potential used to describe hard-sphere interaction
 * 
 * E(R) = R &lt; R_HS ? a * (exp((R - R_HS) ^ 2) - 1) : 0
 * 
 * Input format: a R_HS
 * 
 */
class HardSphere : public PairPotential {
	
	double _force;
	double _radius;
	
public:
	
	virtual void set(const Text& input);
	
	void setForce(double input) { _force = input; }
	
	void setRadius(double input) { _radius = input; _cutoff = input; }
	
	double pairEnergy(double distance) const {
		return distance < _radius ? _force * (exp((distance - _radius) 
				* (distance - _radius)) - 1) : 0.0;
	}
	
	double pairForce(double distance) const {
		return distance < _radius ? 2 * _force * (distance - _radius) * 
				exp((distance - _radius) * (distance - _radius)) : 0.0;
	}

	double tail(const ISO& iso, const Element& elem2) const {
		return 0.0;
	}

};



// =====================================================================================================================
// Pair potential
// =====================================================================================================================

/**
 * inline double PairPotential::density(const ISO& iso, const Element& elem2) const
 *
 * Return the density of atoms of second element
 */
inline double PairPotential::density(const ISO& iso, const Element& elem2) const
{
	for (int i = 0; i < iso.atoms().length(); ++i)
	{
		if (iso.atoms()[i][0].element() == elem2)
			return iso.atoms()[i].length() / iso.basis().volume();
	}
	return 0;
}



#endif
