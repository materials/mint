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



#include "unique.h"
#include "language.h"
#include "output.h"



/* void Unique::search(const ISO& iso, const Symmetry& symmetry, OList<Element> elements, double tol)
 *
 * Search for unique groups
 */

void Unique::search(const ISO& iso, const Symmetry& symmetry, OList<Element> elements, double tol)
{
	
	// Clear data
	clear();
	
	// Return if there are no elements
	if (!elements.length())
		return;
	
	// Make sure that structure contains all elements
	int i, j;
	bool found;
	for (i = 0; i < elements.length(); ++i)
	{
		
		// Loop over elements in the structure
		found = false;
		for (j = 0; j < iso.atoms().length(); ++j)
		{
			if (iso.atoms()[j][0].element() == elements[i])
			{
				found = true;
				break;
			}
		}
		
		// Did not find current element
		if (!found)
		{
			Output::newline(WARNING);
			Output::print(elements[i].symbol());
			Output::print(" is not an element in the system so no groups exist");
			return;
		}
	}
	
	// Set the max distance
	if (_maxDistance < 0)
		_maxDistance = (_areTransitions == true) ? 7.0 : 2.0;
	
	// If these are transitions, make sure that there are at least two atoms
	if ((_areTransitions) && (elements.length() == 1))
		elements += elements[0];
	
	// Output
	Output::newline();
	Output::print("Searching for unique ");
	if (_areTransitions)
		Output::print("transitions between ");
	else if (elements.length() > 1)
		Output::print("groups of ");
	if ((_areTransitions) && (elements.length() == 2))
	{
		Output::print(elements[0].symbol());
		Output::print(" ");
	}
	else
	{
		for (i = 0; i < elements.length(); ++i)
		{
			Output::print(elements[i].symbol());
			if ((elements.length() > 2) && (i != elements.length() - 1))
				Output::print(",");
			if (i == elements.length() - 2)
				Output::print(" and");
			Output::print(" ");
		}
	}
	if (_areTransitions)
		Output::print("sites");
	else
		Output::print("atoms");
	if ((elements.length() > 1) && (!_areTransitions))
	{
		Output::print(" within ");
		Output::print(_maxDistance);
		Output::print(" ang");
	}
	Output::increase();
	
	// Loop over unique atoms of first element
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		
		// Found current element
		if (symmetry.orbits()[i].atoms()[0]->element() == elements[0])
		{
			
			// Loop over atoms in group
			for (j = 0; j < symmetry.orbits()[i].atoms().length(); ++j)
			{
				
				// Skip if not interstitial match
				if (symmetry.orbits()[i].atoms()[j]->isInterstitial() != _useInterstitials)
					continue;
				
				// Create new group
				Group curGroup;
				curGroup.atoms() += symmetry.orbits()[i].atoms()[0];
				curGroup.distances() += 0;
				curGroup.vectors() += Vector3D(0.0);
			
				// Generate full group
				generateGroup(curGroup, iso, symmetry, elements, i, j, tol);
				break;
			}
		}
	}
	
	// Filter transitions if needed
	if (_areTransitions)
		filterTransitions(iso, symmetry, tol);
	
	// Generate equivalent groups if needed
	if (_getEquivalent)
	{
		generateEquivalent(iso, symmetry, tol);
		sortGroups();
	}
	
	// Print results
	Output::newline();
	Output::print("Found ");
	Output::print(_groups.length());
	Output::print(" unique group");
	if (_groups.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _groups.length(); ++i)
	{
		Output::newline();
		Output::print("Unique group ");
		Output::print(i+1);
		Output::print(": ");
		for (j = 0; j < _groups[i][0].atoms().length(); ++j)
		{
			Output::print(_groups[i][0].atoms()[j]->atomNumber() + 1);
			Output::print(" ");
		}
	}
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* void Unique::generateGroup(Group curGroup, const ISO& iso, const Symmetry& symmetry,
 *		const OList<Element>& elements, int startGroup, int startIndex, double tol)
 *
 * Recursively generate remaining atoms in group
 */

void Unique::generateGroup(Group curGroup, const ISO& iso, const Symmetry& symmetry, \
	const OList<Element>& elements, int startGroup, int startIndex, double tol)
{
	
	// Finished if the set number of atoms has been reached
	if (curGroup.atoms().length() == elements.length())
	{
		saveGroup(curGroup, iso, symmetry, tol);
		return;
	}
	
	// If finding transitions then reset start indices
	if ((_areTransitions) || (elements[curGroup.atoms().length()] != elements[curGroup.atoms().length()-1]))
		startGroup = startIndex = 0;
	
	// Allocate space in group
	curGroup.atoms().add();
	curGroup.vectors().add();
	curGroup.distances().add();
	
	// Loop over atoms in the structure of the current element
	int i, j, k;
	double curDis;
	Vector3D curVector;
	ImageIterator images;
	images.setCell(iso.basis(), _maxDistance);
	for (i = startGroup; i < symmetry.orbits().length(); ++i)
	{
		if (elements[curGroup.atoms().length()-1] == symmetry.orbits()[i].atoms()[0]->element())
		{
			
			// Loop over atoms in current orbit
			for (j = (i == startGroup) ? startIndex : 0; j < symmetry.orbits()[i].atoms().length(); ++j)
			{
				
				// Skip if not interstitial match
				if (symmetry.orbits()[i].atoms()[j]->isInterstitial() != _useInterstitials)
					continue;
				
				// Iterate over images
				images.reset(curGroup.atoms()[0]->fractional(), symmetry.orbits()[i].atoms()[j]->fractional());
				while (!images.finished())
				{
					
					// Skip if the same position
					curDis = ++images;
					if (curDis < 1e-8)
						continue;
					
					// Make sure that atom has not already been used
					curVector = images.fracVector();
					for (k = _areTransitions ? curGroup.atoms().length() - 2 : 0; k < curGroup.atoms().length()-1; ++k)
					{
						if ((symmetry.orbits()[i].atoms()[j] == curGroup.atoms()[k]) && \
							(curVector == curGroup.vectors()[k]))
							break;
					}
					if (k != curGroup.atoms().length() - 1)
						continue;
					
					// Save current atom image in group and call next addition
					curGroup.atoms().last() = symmetry.orbits()[i].atoms()[j];
					curGroup.vectors().last() = curVector;
					curGroup.distances().last() = curDis;
					generateGroup(curGroup, iso, symmetry, elements, i, j, tol);
				}
			}
		}
	}
}



/* void Unique::saveGroup(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Save current group of atoms
 */

void Unique::saveGroup(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// If transitions, then check if known
	if (_areTransitions)
	{
		if (isTransKnown(curGroup, iso, symmetry, tol))
			return;
		if (!isTransGood(curGroup, iso, symmetry, tol))
			return;
	}
	
	// If not transitions, then check if known
	else if (isGroupKnown(curGroup, iso, symmetry, tol))
		return;
	
	// Save new group
	_groups.add();
	_groups.last() += curGroup;
}



/* bool Unique::isTransKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Return whether a symmetrically equivalent transition is already known
 */

bool Unique::isTransKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// Loop over groups that are saved
	int i, j, k;
	Vector3D rot;
	Vector3D transPos;
	for (i = 0; i < _groups.length(); ++i)
	{
		
		// Atoms must be in equivalent orbits and have similar distances from first atom
		for (j = 0; j < curGroup.atoms().length(); ++j)
		{
			if (symmetry.orbitNumbers()[curGroup.atoms()[j]->atomNumber()] != \
				symmetry.orbitNumbers()[_groups[i][0].atoms()[j]->atomNumber()])
				break;
			if (Num<double>::neq(curGroup.distances()[j], _groups[i][0].distances()[j], tol))
				break;
		}
		if (j < curGroup.atoms().length())
			continue;
		
		// Loop over symmetry operations
		for (j = 0; j < symmetry.operations().length(); ++j)
		{
			
			// Operation must map first atoms
			rot = symmetry.operations()[j].rotation() * curGroup.atoms()[0]->fractional();
			for (k = 0; k < symmetry.operations()[j].translations().length(); ++k)
			{
				transPos = rot + symmetry.operations()[j].translations()[k];
				ISO::moveIntoCell(transPos);
				if (iso.basis().distance(transPos, FRACTIONAL, _groups[i][0].atoms()[0]->fractional(), \
					FRACTIONAL) < tol)
					break;
			}
			
			// Current operation does not map first atoms
			if (k >= symmetry.operations()[j].translations().length())
				continue;
			
			// Make sure that operation maps all vectors
			for (k = 1; k < _groups[i][0].atoms().length(); ++k)
			{
				if (iso.basis().absoluteDistance(_groups[i][0].vectors()[k], FRACTIONAL, \
					symmetry.operations()[j].rotation() * curGroup.vectors()[k], FRACTIONAL) > tol)
					break;
			}
			
			// Found a match
			if (k >= _groups[i][0].atoms().length())
				return true;
		}
	}
	
	// Return that transition is not known if at this point
	return false;
}



/* bool Unique::isTransGood(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Return whether a transition passes too close to atoms in the structure
 */

bool Unique::isTransGood(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// Get farthest point an atom can be from either end point
	double range = sqrt(curGroup.distances()[1]*curGroup.distances()[1] + \
		(_transPathRadius + tol)*(_transPathRadius + tol));
	
	// Save cartesian postions
	Vector3D cartPos[2];
	cartPos[0] = curGroup.atoms()[0]->cartesian();
	cartPos[1] = curGroup.atoms()[0]->fractional();
	cartPos[1] += curGroup.vectors()[1];
	iso.basis().toCartesian(cartPos[1]);
	
	// Loop over atoms in the structure
	int i, j, k;
	Vector3D cartPoint;
	ImageIterator images(iso.basis(), range);
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Loop over start and end points
			for (k = 0; k <= 1; ++k)
			{
				
				// Loop over images of atom
				images.reset(curGroup.atoms()[k]->fractional(), iso.atoms()[i][j].fractional());
				while (!images.finished())
				{
					if (++images > 1e-4)
					{
						
						// Get cartesian position of point
						cartPoint  = curGroup.atoms()[0]->fractional();
						cartPoint += curGroup.vectors()[k];
						cartPoint += images.fracVector();
						iso.basis().toCartesian(cartPoint);
						
						// Skip if points are the same
						if (cartPoint.distanceToPoint(cartPos[0]) < tol)
							continue;
						if (cartPoint.distanceToPoint(cartPos[1]) < tol)
							continue;
						
						// Check if within allowed distance
						if (cartPoint.distanceToLineSegment(cartPos[0], cartPos[1]) < _transPathRadius - tol)
							return false;
					}
				}
			}
		}
	}
	
	// Return that transition is allowed
	return true;
}



/* bool Unique::isGroupKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Return whether a symmetrically equivalent group is already known
 */

bool Unique::isGroupKnown(const Group& curGroup, const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// Loop over groups that are already known
	int i, j, k, m;
	int ref;
	Vector3D rot;
	Vector3D transPos;
	List<bool> matched(curGroup.atoms().length());
	List<double> distances(curGroup.atoms().length());
	OList<Vector3D> vectors(curGroup.atoms().length());
	for (i = 0; i < _groups.length(); ++i)
	{
		
		// Loop over symmetry operations
		for (j = 0; j < symmetry.operations().length(); ++j)
		{
			
			// Loop over atoms in the second group
			matched.fill(false);
			rot = symmetry.operations()[j].rotation() * curGroup.atoms()[0]->fractional();
			for (k = 0; k < _groups[i][0].atoms().length(); ++k)
			{
				
				// Skip if atoms are not in the same orbit
				if (symmetry.orbitNumbers()[curGroup.atoms()[0]->atomNumber()] != \
					symmetry.orbitNumbers()[_groups[i][0].atoms()[k]->atomNumber()])
					continue;
				
				// Loop over translations and check if atoms map
				for (m = 0; m < symmetry.operations()[j].translations().length(); ++m)
				{
					transPos = rot + symmetry.operations()[j].translations()[m];
					ISO::moveIntoCell(transPos);
					if (iso.basis().distance(transPos, FRACTIONAL, _groups[i][0].atoms()[k]->fractional(), \
						FRACTIONAL) < tol)
					{
						ref = k;
						matched[k] = true;
						break;
					}
				}
				
				// Matched an atom
				if (m < symmetry.operations()[j].translations().length())
					break;
			}
			
			// Not a possible match
			if (k >= _groups[i][0].atoms().length())
				continue;
			
			// Save vectors for new reference atom
			for (k = 0; k < _groups[i][0].atoms().length(); ++k)
			{
				vectors[k] = _groups[i][0].vectors()[k] - _groups[i][0].vectors()[ref];
				distances[k] = iso.basis().getCartesian(vectors[k]).magnitude();
			}
			
			// Loop over pairs of atoms
			for (k = 1; k < curGroup.atoms().length(); ++k)
			{
				rot = symmetry.operations()[j].rotation() * curGroup.vectors()[k];
				for (m = 0; m < _groups[i][0].atoms().length(); ++m)
				{
					
					// Skip if already matched
					if (matched[m])
						continue;
					
					// Skip if not in the same orbit
					if (symmetry.orbitNumbers()[curGroup.atoms()[k]->atomNumber()] != \
						symmetry.orbitNumbers()[_groups[i][0].atoms()[m]->atomNumber()])
						continue;
					
					// Skip if distances do not match
					if (Num<double>::neq(curGroup.distances()[k], distances[m], tol))
						continue;
					
					// Found a match
					if (iso.basis().absoluteDistance(rot, FRACTIONAL, vectors[m], FRACTIONAL) < tol)
					{
						matched[m] = true;
						break;
					}
				}
				
				// Did not find a match
				if (m >= _groups[i][0].atoms().length())
					break;
			}
			
			// Check if all atoms were matched
			for (k = 0; k < matched.length(); ++k)
			{
				if (!matched[k])
					break;
			}
			if (k >= matched.length())
				return true;
		}
	}
	
	// Return that group is new if at this point
	return false;
}



/* void Unique::filterTransitions(const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Save only the transitions that define Voronoi volume about each atom
 */

void Unique::filterTransitions(const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// Loop over transitions
	int i, j, k;
	int start;
	int orbitNum;
	double curDis;
	double minDis;
	Atom* nearAtom;
	Vector3D newPos;
	Vector3D tempPos;
	Vector3D curJump;
	Vector3D curCell;
	Vector3D nearCell;
	List<int> groupsToRemove;
	OList<Vector3D> norms;
	OList<Vector3D> points;
	OList<Vector3D> centers;
	for (i = 0; i < _groups.length();)
	{
		
		// Loop over groups that have the same starting atom
		start = i;
		points.length(0);
		for (; i < _groups.length(); ++i)
		{
			
			// Break if no longer on the same starting atom
			if (_groups[i][0].atoms()[0]->atomNumber() != _groups[start][0].atoms()[0]->atomNumber())
				break;
			
			// Loop over point group operations for current atom to generate equivalent jumps
			curJump = _groups[i][0].vectors()[1];
			orbitNum = symmetry.orbitNumbers()[_groups[i][0].atoms()[0]->atomNumber()];
			const OList<Matrix3D>& rotations = symmetry.orbits()[orbitNum].specialPositions()[0].rotations();
			orbitNum = symmetry.orbitNumbers()[_groups[i][0].atoms()[1]->atomNumber()];
			const Orbit& orbit = symmetry.orbits()[orbitNum];
			for (j = 0; j < rotations.length(); ++j)
			{
				
				// Save ending position for rotated jump
				newPos  = curJump;
				newPos *= rotations[j];
				newPos += _groups[i][0].atoms()[0]->fractional();
				tempPos = newPos;
				ISO::moveIntoCell(tempPos);
				
				// Find which atom is nearest to rotated position
				for (k = 0; k < orbit.atoms().length(); ++k)
				{
					curDis = iso.basis().distance(tempPos, FRACTIONAL, orbit.atoms()[k]->fractional(), \
						FRACTIONAL, &curCell);
					if ((!k) || (curDis < minDis))
					{
						minDis = curDis;
						nearCell = curCell;
						nearAtom = orbit.atoms()[k];
					}
				}
				
				// Save position of atom
				newPos -= tempPos;	
				newPos += nearAtom->fractional();
				newPos += nearCell;
				
				// Save if point is new
				for (k = 0; k < points.length(); ++k)
				{
					if (iso.basis().absoluteDistance(newPos, FRACTIONAL, points[k], FRACTIONAL) < 1e-4)
						break;
				}
				if (k >= points.length())
					points += newPos;
			}
		}
		
		// Convert all points to cartesian coordinates
		for (j = 0; j < points.length(); ++j)
			iso.basis().toCartesian(points[j]);
		
		// Get Voronoi volume for current set of atoms
		_groups[start][0].atoms()[0]->cartesian().voronoi(points, tol, 0, &centers, &norms, 0);
		
		// Loop over jumps and skip those that do not contribute to Voronoi volume
		for (i = start; i < _groups.length(); ++i)
		{
			
			// Break if no longer on the same starting atom
			if (_groups[i][0].atoms()[0]->atomNumber() != _groups[start][0].atoms()[0]->atomNumber())
				break;
			
			// Loop over planes and check if current intersection occurs within half of length
			curJump = iso.basis().getCartesian(_groups[i][0].vectors()[1]);
			minDis = curJump.magnitude() / 2 - _voronoiTolerance - tol;
			for (j = 0; j < centers.length(); ++j)
			{
				
				// Get intersection
				curDis = _groups[start][0].atoms()[0]->cartesian().nearestPointOnPlane(curJump, centers[j], norms[j]);
				if ((curDis > 0) && (Num<double>::abs(curDis) < minDis))
				{
					groupsToRemove += i;
					break;
				}
			}
		}
	}
	
	// Remove groups
	for (i = groupsToRemove.length() - 1; i >= 0; --i)
		_groups.remove(groupsToRemove[i]);
}



/* void Unique::generateEquivalent(const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Generate groups that are equivalent to unique ones
 */

void Unique::generateEquivalent(const ISO& iso, const Symmetry& symmetry, double tol)
{
	
	// Return if there are no groups
	if (_groups.length() == 0)
		return;
	
	// Set break tolerance
	tol /= 10;
	
	// Loop over unique groups
	int i, j, k, m, n;
	int nearIndex;
	double curDis;
	double nearDis;
	const Orbit* orbit;
	Group curGroup;
	Vector3D pos;
	Vector3D rotPos;
	Vector3D transPos;
	for (i = 0; i < _groups.length(); ++i)
	{
		
		// Save copy of current group
		curGroup = _groups[i][0];
		
		// Loop over symmetry operations
		for (j = 0; j < symmetry.operations().length(); ++j)
		{
			
			// Get rotated vectors
			for (k = 1; k < curGroup.vectors().length(); ++k)
				curGroup.vectors()[k] = symmetry.operations()[j].rotation() * _groups[i][0].vectors()[k];
			
			// Loop over translations
			rotPos = symmetry.operations()[j].rotation() * _groups[i][0].atoms()[0]->fractional();
			for (k = 0; k < symmetry.operations()[j].translations().length(); ++k)
			{
				
				// Save information for current first atom
				transPos = rotPos + symmetry.operations()[j].translations()[k];
				ISO::moveIntoCell(transPos);
				
				// Loop over atoms in the group
				for (m = 0; m < curGroup.atoms().length(); ++m)
				{
					
					// Generate position of atom
					pos = transPos + curGroup.vectors()[m];
					ISO::moveIntoCell(pos);
					
					// Save first atom in orbit as nearest
					nearIndex = 0;
					orbit = &symmetry.orbits()[symmetry.orbitNumbers()[_groups[i][0].atoms()[m]->atomNumber()]];
					nearDis = iso.basis().distance(orbit->atoms()[0]->fractional(), FRACTIONAL, pos, FRACTIONAL);

					// Loop over atoms in orbit and check for match
					for (n = 1; n < orbit->atoms().length(); ++n)
					{

						// Current distance is within tolerance
						if (nearDis < tol)
							break;

						// Check if a new best distance
						curDis = iso.basis().distance(orbit->atoms()[n]->fractional(), FRACTIONAL, pos, FRACTIONAL);
						if (curDis < nearDis)
						{
							nearDis = curDis;
							nearIndex = n;
						}
					}
					
					// Save nearest atom
					curGroup.atoms()[m] = orbit->atoms()[nearIndex];
				}
				
				// Save group
				saveGroup(curGroup, i);
			}
		}
	}
}



/* void Unique::saveGroup(const Group& curGroup, int set)
 *
 * Save equivalent group
 */

void Unique::saveGroup(const Group& curGroup, int set)
{
	
	// Check if already known
	if (_areTransitions)
	{
		if (isEquivTransKnown(curGroup, set))
			return;
	}
	else if (isEquivGroupKnown(curGroup, set))
		return;
	
	// Save new group
	_groups[set] += curGroup;
}



/* bool Unique::isEquivTransKnown(const Group& curGroup, int set)
 *
 * Check if equivalent translation is already known in group
 */

bool Unique::isEquivTransKnown(const Group& curGroup, int set)
{
	
	// Loop over groups already in set
	int i, j;
	for (i = 0; i < _groups[set].length(); ++i)
	{
		
		// Check if there is a 1:1 map
		for (j = 0; j < curGroup.atoms().length(); ++j)
		{
			if (_groups[set][i].atoms()[j] != curGroup.atoms()[j])
				break;
		}
		
		// Found a match
		if (j >= curGroup.atoms().length())
			return true;
	}
	
	// Return that new if at this point
	return false;
}



/* bool Unique::isEquivGroupKnown(const Group& curGroup, int set)
 * 
 * Check if equivalent group is already known in group
 */

bool Unique::isEquivGroupKnown(const Group& curGroup, int set)
{
	
	// Loop over groups already in set
	int i, j, k;
	List<bool> matched(curGroup.atoms().length());
	for (i = 0; i < _groups[set].length(); ++i)
	{
		
		// Loop over pairs of atoms
		matched.fill(false);
		for (j = 0; j < curGroup.atoms().length(); ++j)
		{
			for (k = 0; k < _groups[set][i].atoms().length(); ++k)
			{
				
				// Check if atoms are the same
				if (matched[k])
					continue;
				if (curGroup.atoms()[j] == _groups[set][i].atoms()[k])
				{
					matched[k] = true;
					break;
				}
			}
			
			// Atom was not matched
			if (k >= _groups[set][i].atoms().length())
				break;
		}
		
		// Found a match
		if (j >= curGroup.atoms().length())
			return true;
	}
	
	// Return that new if at this point
	return false;
}



/* void Unique::print() const
 *
 * Print unique groups
 */

void Unique::print() const
{
	
	// Set output method
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print number of unique groups
	int i;
	Output::newline();
	Output::print("Found ");
	Output::print(_groups.length());
	Output::print(" unique ");
	if (_areTransitions)
	{
		Output::print("transition");
		if (_groups.length() != 1)
			Output::print("s");
	}
	else
	{
		if (_groups.length() > 0)
		{
			if (_groups[0][0].atoms().length() == 1)
			{
				Output::print(_groups[0][0].atoms()[0]->element().symbol());
				Output::print(" atom");
			}
			else
				Output::print("group");
			if (_groups.length() != 1)
				Output::print("s");
			if (_groups[0][0].atoms().length() != 1)
			{
				Output::print(" of ");
				for (i = 0; i < _groups[0][0].atoms().length(); ++i)
				{
					Output::print(_groups[0][0].atoms()[i]->element().symbol());
					if ((_groups[0][0].atoms().length() > 2) && (i != _groups[0][0].atoms().length() - 1))
						Output::print(",");
					if (i == _groups[0][0].atoms().length() - 2)
						Output::print(" and");
					Output::print(" ");
				}
				Output::print("atoms");
			}
		}
		else
			Output::print("groups");
	}
	
	// Return if no groups
	if (!_groups.length())
	{
		Output::method(origMethod);
		return;
	}
	
	// Set title to use
	Word title;
	int numAtoms = _groups[0][0].atoms().length();
	if (_areTransitions)
		title = "Transition ";
	else if (numAtoms == 1)
		title = "Atom ";
	else
		title = "Group ";
	
	// Loop over groups
	int j, k;
	int cutLen = 75;
	Word curLine;
	char connector = _areTransitions ? '>' : ',';
	for (i = 0; i < _groups.length(); ++i)
	{
		
		// Print group number
		Output::newline();
		if (_getEquivalent)
			Output::newline();
		Output::tab();
		Output::print(title);
		Output::print(i + 1);
		Output::print(": ");
		
		// Printing only a single instance
		if (!_getEquivalent)
		{
			Output::print('[');
			for (j = 0; j < numAtoms; ++j)
			{
				Output::print(_groups[i][0].atoms()[j]->atomNumber() + 1);
				if (j != numAtoms - 1)
					Output::print(connector);
			}
			Output::print(']');
		}
		
		// Printing all equivalent groups
		else
		{
			
			// Print number in group
			Output::print(_groups[i].length());
			Output::print(" equivalent instance");
			if (_groups[i].length() != 1)
				Output::print("s");

			// Print equivalent groups
			curLine = "        ";
			for (j = 0; j < _groups[i].length(); ++j)
			{

				// Add current group
				curLine += '[';
				for (k = 0; k < numAtoms; ++k)
				{
					curLine += Language::numberToWord(_groups[i][j].atoms()[k]->atomNumber() + 1);
					if (k != numAtoms - 1)
						curLine += connector;
				}
				curLine += ']';
				if (j != _groups[i].length() - 1)
					curLine += ' ';

				// Print if over length
				if (curLine.length() > cutLen)
				{
					Output::newline();
					Output::print(curLine);
					curLine = "        ";
				}
			}
			if (curLine.length() > 8)
			{
				Output::newline();
				Output::print(curLine);
			}
		}
	}
	
	// Reset output
	Output::method(origMethod);
}



/* Word Unique::makeISO(ISO& iso, int group) const
 *
 * Turn unique group into a structure
 */

Word Unique::makeISO(ISO& iso, int group) const
{
	
	// Return if no atoms in group
	if (group >= _groups.length())
		return Word("");
	if (_groups[group][0].atoms().length() == 0)
		return Word("");
	
	// Output
	Output::newline();
	Output::print("Creating structure from unique group ");
	Output::print(group + 1);
	Output::increase();
	
	// Make list of atoms
	int i;
	List<int> atoms;
	atoms.length(_groups[group][0].atoms().length());
	for (i = 0; i < atoms.length(); ++i)
		atoms[i] = _groups[group][0].atoms()[i]->atomNumber();
	atoms.sort();
	
	// Loop over atoms in reverse and remove
	for (i = atoms.length() - 1; i >= 0; --i)
		iso.removeAtom(atoms[i]);
	
	// Make description to return
	Word des = "Removed atom";
	if (atoms.length() != 1)
		des += 's';
	for (i = 0; i < atoms.length(); ++i)
	{
		des += ' ';
		des += Language::numberToWord(atoms[i] + 1);
		if ((i != atoms.length() - 1) && (atoms.length() > 2))
			des += ',';
		if (i == atoms.length() - 2)
			des += " and";
	}
	
	// Output
	Output::decrease();
	
	// Return description
	return des;
}
