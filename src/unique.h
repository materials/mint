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



#ifndef UNIQUE_H
#define UNIQUE_H



#include "iso.h"
#include "symmetry.h"
#include "text.h"
#include "list.h"
#include "num.h"



// Class to find unique groups of atoms
class Unique
{
	
	// Class to store a single unique group
	class Group;

	// Variables
	OList<Group>::D2 _groups;
	
	// Settings variables
	bool _getEquivalent;
	bool _areTransitions;
	bool _useInterstitials;
	double _maxDistance;
	double _transPathRadius;
	double _voronoiTolerance;
	
	// Functions
	void initialize();
	void generateGroup(Group curGroup, const ISO& iso, const Symmetry& symmetry, const OList<Element>& elements, \
		int startGroup, int startIndex, double tol);
	void saveGroup(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol);
	bool isTransKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol);
	bool isTransGood (const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol);
	bool isGroupKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol);
	void filterTransitions(const ISO& iso, const Symmetry& symmetry, double tol);
	void generateEquivalent(const ISO& iso, const Symmetry& symmetry, double tol);
	void saveGroup(const Group& curGroup, int set);
	bool isEquivTransKnown(const Group& curGroup, int set);
	bool isEquivGroupKnown(const Group& curGroup, int set);
	void sortGroups();
	void sortGroup(int num, int left, int right);
	
public:
	
	// Constructor
	Unique();
	
	// Setup functions
	void search(const ISO& iso, const Symmetry& symmetry, OList<Element> elements, double tol);
	
	// Clear data
	void clear();
	
	// Settings
	void getEquivalent   	(bool input)	{ _getEquivalent    = input; }
	void areTransitions 	(bool input)	{ _areTransitions   = input; }
	void useInterstitials	(bool input)	{ _useInterstitials = input; }
	void maxDistance     	(double input)	{ _maxDistance      = input; }
	
	// Other functions
	void print() const;
	Word makeISO(ISO& iso, int group) const;
	
	// Access functions
	const OList<Group>::D2& groups() const	{ return _groups; }
};



// Class to store a single unique group
class Unique::Group
{
	
	// Variables
	List<Atom*> _atoms;
	List<double> _distances;
	OList<Vector3D> _vectors;
	
	// Functions
	bool isLess(const Group& rhs, int curNum) const;
	
public:
	
	// Constructors
	Group()					{}
	Group(const Group& rhs)	{ *this = rhs; }
	
	// Functions
	Group& operator= (const Group& rhs);
	bool operator< (const Group& rhs) const;
	
	// Access functions
	List<Atom*>& atoms()					{ return _atoms; }
	const List<Atom*>& atoms() const		{ return _atoms; }
	List<double>& distances()				{ return _distances; }
	const List<double>& distances() const	{ return _distances; }
	OList<Vector3D>& vectors()				{ return _vectors; }
	const OList<Vector3D>& vectors() const	{ return _vectors; }
};



/* inline Unique::Unique()
 *
 * Constructor
 */

inline Unique::Unique()
{
	initialize();
}



/* inline void Unique::initialize()
 *
 * Set defaults
 */

inline void Unique::initialize()
{
	_getEquivalent = false;
	_areTransitions = false;
	_useInterstitials = false;
	_maxDistance = -1;
	_transPathRadius = 0.25;
	_voronoiTolerance = 0.05;
}



/* inline void Unique::clear()
 *
 * Clear data
 */

inline void Unique::clear()
{
	_groups.clear();
}



/* inline void Unique::sortGroups()
 *
 * Sort groups
 */

inline void Unique::sortGroups()
{
	for (int i = 0; i < _groups.length(); ++i)
		sortGroup(i, 0, _groups[i].length() - 1);
}



/* inline void Unique::sortGroup(int num, int left, int right)
 * 
 * Quick sort of current group
 */

inline void Unique::sortGroup(int num, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	Group pivot = _groups[num][pivotIndex];
	
	// Move pivot to end
	_groups[num].swap(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (_groups[num][i] < pivot)
		{
			_groups[num].swap(i, newPivotIndex);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	_groups[num].swap(newPivotIndex, right);
	
	// Recursive calls to next sorts
	sortGroup(num, left, newPivotIndex - 1);
	sortGroup(num, newPivotIndex + 1, right);
}



/* inline Unique::Group& Unique::Group::operator= (const Unique::Group& rhs)
 *
 * Assignment operator for unique group
 */

inline Unique::Group& Unique::Group::operator= (const Unique::Group& rhs)
{
	if (this != &rhs)
	{
		_atoms = rhs._atoms;
		_vectors = rhs._vectors;
		_distances = rhs._distances;
	}
	return *this;
}



/* inline bool Unique::Group::operator< (const Unique::Group& rhs) const
 *
 * Test whether current group is less than another
 */

inline bool Unique::Group::operator< (const Unique::Group& rhs) const
{
	return isLess(rhs, 0);
}



/* inline bool Unique::Group::isLess(const Unique::Group& rhs, int curNum) const
 *
 * Return whether current value in group is less than another
 */

inline bool Unique::Group::isLess(const Unique::Group& rhs, int curNum) const
{
	if ((curNum >= _atoms.length()) || (curNum >= rhs._atoms.length()))
		return false;
	if (_atoms[curNum]->atomNumber() < rhs._atoms[curNum]->atomNumber())
		return true;
	if (_atoms[curNum]->atomNumber() > rhs._atoms[curNum]->atomNumber())
		return false;
	return isLess(rhs, curNum+1);
}



#endif
