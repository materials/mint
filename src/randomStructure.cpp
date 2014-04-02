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



#include "randomStructure.h"
#include "output.h"
#include "constants.h"
#include <cmath>



// Maximum number of loops to make when trying to place atoms at random positions
int RandomStructure::_maxTrialLoops = 100;
double RandomStructure::_minBondFraction = 0.5;

// Wyckoff site biasing weights
double RandomStructure::_weights[143] = {0.0, 0.853446103005, 0.000751486636607, 0.00161032850702, 0.00139095042055,\
	0.00229646847957, 0.00432688268407, 0.0125605623547, 0.000266054275072, 0.0013022656622, 1.86704754437e-05,\
	0.0330747472484, 0.000135360946967, 0.000504102836979, 0.00414017792963, 0.000172701897854, 0.000872844726991,\
	0.000574117119892, 0.000140028565827, 4.20085697482e-05, 4.66761886091e-06, 0.071479915236, 0.0,\
	3.26733320264e-05, 1.86704754437e-05, 2.80057131655e-05, 0.000191372373297, 0.000242716180768, 0.0,\
	3.26733320264e-05, 0.0, 0.00101287329282, 0.0, 2.80057131655e-05, 0.000242716180768, 0.0, 7.93495206355e-05,\
	1.40028565827e-05, 0.0, 0.0, 0.0, 0.00549845501816, 0.0, 0.0, 0.0, 4.66761886091e-06, 0.0, 0.000172701897854,\
	0.0, 3.26733320264e-05, 0.0, 0.000336068557986, 0.0, 0.0, 9.33523772183e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
	0.00125092185473, 0.0, 2.33380943046e-05, 0.0, 0.0, 0.0, 0.000130693328106, 0.0, 9.33523772183e-05, 0.0,\
	7.93495206355e-05, 5.13438074701e-05, 0.0, 0.000126025709245, 0.0, 3.73409508873e-05, 0.0, 5.6011426331e-05,\
	0.0, 0.0, 0.000214710467602, 0.0, 0.0, 9.33523772183e-06, 0.0, 2.80057131655e-05, 9.33523772183e-06, 0.0,\
	1.86704754437e-05, 0.0, 0.0, 9.33523772183e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000770157112051, 0.0,\
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.66761886091e-06, 0.0};



/* void RandomStructure::perturbBasis(ISO& iso, double min, double max, Random& random)
 *
 * Randomly adjust basis vectors
 */

void RandomStructure::perturbBasis(ISO& iso, double min, double max, Random& random)
{
	
	// Check if any elements of basis are fixed
	for (int i = 0; i < 3; ++i)
	{
		if ((iso.basis().lengthFixed()[i]) || (iso.basis().angleFixed()[i]))
			return;
	}
	
	// Swap min and max if needed
	if (max < min)
		Num<double>::swap(min, max);
	
	// Output
	Output::newline();
	Output::print("Perturbing basis vectors by ");
	if (min == max)
		Output::print(min);
	else
	{
		Output::print(" an amount between ");
		Output::print(min);
		Output::print(" and ");
		Output::print(max);
	}
	Output::print(" Ang");
	Output::increase();
	
	// Save new basis
	iso.basis(generateBasis(iso.basis().vectors(), min, max, random), true);
	
	// Output
	Output::decrease();
}



/* void RandomStructure::perturbAtoms(ISO& iso, double min, double max, Random& random)
 *
 * Randomly perturb atoms
 */

void RandomStructure::perturbAtoms(ISO& iso, double min, double max, Random& random)
{
	
	// Swap min and max if needed
	if (max < min)
		Num<double>::swap(min, max);
	
	// Output
	Output::newline();
	Output::print("Perturbing positions of atoms by ");
	if (min == max)
		Output::print(min);
	else
	{
		Output::print(" an amount between ");
		Output::print(min);
		Output::print(" and ");
		Output::print(max);
	}
	Output::print(" Ang");
	Output::increase();
	
	// Loop over atoms
	int i, j, k;
	double mag;
	double norm;
	Vector3D vector;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Skip atom if fixed
			if (iso.atoms()[i][j].anyFixed())
				continue;
			
			// Generator magnitude
			mag = (min == max) ? min : random.decimal(min, max);

			// Generate random direction
			for (k = 0; k < 3; ++k)
				vector[k] = random.decimal(-1, 1);

			// Normalize vector to magnitude
			norm = vector.magnitude();
			vector *= mag / norm;
			
			// Set new position
			iso.atoms()[i][j].cartesian(iso.atoms()[i][j].cartesian() + vector);
			
			// Output
			Output::newline();
			Output::print("Position of ");
			Output::print(iso.atoms()[i][j].element().symbol());
			Output::print(" ");
			Output::print(j + 1);
			Output::print(": ");
			Output::print(iso.atoms()[i][j].fractional(), 8, false);
		}
	}
	
	// Output
    Output::decrease();
}



/* void RandomStructure::perturbBasis(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random)
 *
 * Perturb basis while preserving symmetry
 */

void RandomStructure::perturbBasis(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random)
{
	
	// Check if any elements of basis are fixed
	int i;
	for (i = 0; i < 3; ++i)
	{
		if ((iso.basis().lengthFixed()[i]) || (iso.basis().angleFixed()[i]))
			return;
	}
	
	// Swap min and max if needed
	if (max < min)
		Num<double>::swap(min, max);
	
	// Output
	Output::newline();
	Output::print("Perturbing basis vectors by ");
	if (min == max)
		Output::print(min);
	else
	{
		Output::print("an amount between ");
		Output::print(min);
		Output::print(" and ");
		Output::print(max);
	}
	Output::print(" Ang while preserving symmetry");
	Output::increase();
	
	// Generate new basis
	Matrix3D newVectors = generateBasis(iso.basis().vectors(), min, max, random, &symmetry);
	
	// Confine basis
	Vector3D newLengths = Basis::lengths(newVectors);
	Vector3D newAngles = Basis::angles(newVectors);
	Symmetry::confineBasis(newLengths, newAngles, iso.basis().latticeSystem());
	
	// Save new basis
	iso.basis(newLengths, newAngles, true);
	
	// Output
	Output::decrease();
}



/* void RandomStructure::perturbAtoms(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random)
 *
 * Perturb atoms while preserving symmetry
 */

void RandomStructure::perturbAtoms(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random)
{
	
	// Swap min and max if needed
	if (max < min)
		Num<double>::swap(min, max);
	
	// Output
	Output::newline();
	Output::print("Perturbing positions of atoms by ");
	if (min == max)
		Output::print(min);
	else
	{
		Output::print(" an amount between ");
		Output::print(min);
		Output::print(" and ");
		Output::print(max);
	}
	Output::print(" Ang while preserving symmetry");
	Output::increase();
	
	// Loop over unique atoms
	int i, j;
	double mag;
	double norm;
	Vector3D vector;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		
		// Check if any atom in orbit is fixed
		if (symmetry.orbits()[i].anyAtomsFixed())
			break;
		
		// Generator magnitude
		mag = (min == max) ? min : random.decimal(min, max);

		// Generate random direction
		for (j = 0; j < 3; ++j)
			vector[j] = random.decimal(-1, 1);
		
		// Constrain vector with site symmetry and turn into cartesian coordinate
		vector *= symmetry.orbits()[i].specialPositions()[0].rotation();
		iso.basis().toCartesian(vector);
		
		// Get norm of vector and skip if zero
		norm = vector.magnitude();
		if (!norm)
			continue;
		
		// Normalize vector and set new position
		vector *= max / norm;
		symmetry.orbits()[i].atoms()[0]->cartesian(symmetry.orbits()[i].atoms()[0]->cartesian() + vector);
		
		// Output
		Output::newline();
		Output::print("Position of atom ");
		Output::print(symmetry.orbits()[i].atoms()[0]->atomNumber() + 1);
		Output::print(": ");
		Output::print(symmetry.orbits()[i].atoms()[0]->fractional(), 8, false);
		
		// Loop over equivalent atoms and save positions
		for (j = 1; j < symmetry.orbits()[i].atoms().length(); ++j)
		{
			
			// Set position
			symmetry.orbits()[i].atoms()[j]->fractional(\
				symmetry.orbits()[i].generators()[j].rotation() * symmetry.orbits()[i].atoms()[0]->fractional() + \
				symmetry.orbits()[i].generators()[j].translations()[0]);
			
			// Output
			Output::newline();
			Output::print("Position of atom ");
			Output::print(symmetry.orbits()[i].atoms()[j]->atomNumber() + 1);
			Output::print(": ");
			Output::print(symmetry.orbits()[i].atoms()[j]->fractional(), 8, false);
		}
	}
	
	// Output
    Output::decrease();
}



/* void RandomStructure::generate(ISO& iso, Random& random, double bias, Symmetry* symmetry)
 *
 * Generate random structure
 */

void RandomStructure::generate(ISO& iso, Random& random, double bias, Symmetry* symmetry)
{
	
	// Calculate the total target volume
	double targetVolume = 0;
	for (int i = 0; i < iso.atoms().length(); ++i)
		targetVolume += iso.atoms()[i].length() * pow(2*iso.atoms()[i][0].element().radius(), 3.0);
	
	// Space group has been set
	if (iso.spaceGroup().length() > 0)
		generateWithSymmetry(iso, targetVolume, random, bias, symmetry);
	
	// Space group has not been set
	else
		generateWithoutSymmetry(iso, targetVolume, random, symmetry);
}



/* void RandomStructure::generateWithoutSymmetry(ISO& iso, double targetVolume, Random& random, Symmetry* symmetry)
 *
 * Generate a random structure
 */

void RandomStructure::generateWithoutSymmetry(ISO& iso, double targetVolume, Random& random, Symmetry* symmetry)
{
	
	// Output
	Output::newline();
	Output::print("Generating random structure");
	Output::increase();
	
	// Generate the basis
	if (iso.basis().volume() <= 0)
		iso.basis(generateBasis(targetVolume, random, iso.basis().latticeSystem()), true);
	
	// Generate the positions
	if (iso.anyUnset())
		generatePositions(iso, random);
	
	// Set symmetry to P1 if needed
	if (symmetry)
	{
		symmetry->clear();
		symmetry->setToP1(iso);
	}
	
	// Output
	Output::decrease();
}



/* void RandomStructure::generateWithSymmetry(ISO& iso, double targetVolume, Random& random, double bias,
 *		Symmetry* symmetry)
 *
 * Generate a random structure within a set space group
 */

void RandomStructure::generateWithSymmetry(ISO& iso, double targetVolume, Random& random, double bias, \
	Symmetry* symmetry)
{
	
	// Output
	Output::newline();
	Output::print("Generating random structure in space group ");
	Output::print(iso.spaceGroup());
	Output::increase();
	
	// Get the space group
	SpaceGroup spaceGroup(iso.spaceGroup(), true, true);
	iso.basis().latticeSystem(spaceGroup.system());
	
	// Generate the basis
	if (iso.basis().volume() <= 0)
		iso.basis(generateBasis(targetVolume, random, spaceGroup.system()), true);
	
	// Set general symmetry information
	if (symmetry)
	{
		symmetry->clear();
		symmetry->operations(spaceGroup.symmetry());
	}
	setSiteSymmetry(iso, spaceGroup, symmetry);
	
	// Generate the positions
	if (iso.anyUnset())
		generatePositions(iso, random, spaceGroup, bias, symmetry);
	
	// Output
	Output::decrease();
}



/* Matrix3D RandomStructure::generateBasis(const Matrix3D origMat, double min, double max,
 *		Random& random, const Symmetry* symmetry)
 * 
 * Generate matrix to change basis
 */

Matrix3D RandomStructure::generateBasis(const Matrix3D origMat, double min, double max, \
	Random& random, const Symmetry* symmetry)
{
	
	// Loop until a good basis is found
	int i, j;
	double mag;
	Vector3D vector;
	Matrix3D newBasis;
	do
	{
		
		// Loop over directions
		for (i = 0; i < 3; ++i)
		{

			// Generator magnitude
			mag = (min == max) ? min : random.decimal(min, max);

			// Generate random direction
			for (j = 0; j < 3; ++j)
				vector[j] = random.decimal(-1, 1);

			// Normalize vector to magnitude
			vector *= mag / vector.magnitude();

			// Set new basis row
			for (j = 0; j < 3; ++j)
				newBasis(i, j) = origMat(i, j) + vector[j];
		}
		
		// Confine basis
		if (symmetry)
			symmetry->refineBasis(newBasis);
		
	} while (newBasis.determinant() <= 0);
	
	// Return new matrix
	return newBasis;
}



/* Matrix3D RandomStructure::generateBasis(double targetVolume, Random& random, LatticeSystem latticeSystem)
 *
 * Generate random basis within set lattice system
 */

Matrix3D RandomStructure::generateBasis(double targetVolume, Random& random, LatticeSystem latticeSystem)
{
	
	// Generate volume
	double volume = random.decimal(0.75*targetVolume, 1.25*targetVolume);
	
	// Generate in cubic system
	double avgLength;
	Vector3D lengths;
	Vector3D angles;
	if (latticeSystem == LS_CUBIC)
	{
		lengths[0] = lengths[1] = lengths[2] = pow(volume, 1.0/3.0);
		angles[0] = angles[1] = angles[2] = Constants::pi / 2.0;
	}
	
	// Generate in hexagonal system
	else if (latticeSystem == LS_HEXAGONAL)
	{
		avgLength = pow(2 * volume / sqrt(3), 1.0/3.0);
		angles[0] = angles[1] = Constants::pi / 2.0;
		angles[2] = 2 * Constants::pi / 3.0;
		lengths[0] = lengths[1] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		lengths[2] = 2 * volume / (sqrt(3) * lengths[0] * lengths[0]);
	}
	
	// Generate in rhombohedral system
	else if (latticeSystem == LS_RHOMBOHEDRAL)
	{
		angles[0] = angles[1] = angles[2] = random.decimal(Constants::pi/3, 2*Constants::pi/3);
		lengths[0] = lengths[1] = lengths[2] = \
			pow(volume / sqrt(1 - 3*pow(cos(angles[0]), 2) + 2*pow(cos(angles[0]), 3)), 1.0/3.0);
	}
	
	// Generate in tetragonal system
	else if (latticeSystem == LS_TETRAGONAL)
	{
		avgLength = pow(volume, 1.0/3.0);
		angles[0] = angles[1] = angles[2] = Constants::pi / 2.0;
		lengths[0] = lengths[1] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		lengths[2] = volume / (lengths[0] * lengths[0]);
	}
	
	// Generate in orthorhombic system
	else if (latticeSystem == LS_ORTHORHOMBIC)
	{
		avgLength = pow(volume, 1.0/3.0);
		angles[0] = angles[1] = angles[2] = Constants::pi / 2.0;
		lengths[0] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		avgLength = sqrt(volume / lengths[0]);
		lengths[1] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		lengths[2] = volume / (lengths[0]*lengths[1]);
	}
	
	// Generate in monoclinic system
	else if (latticeSystem == LS_MONOCLINIC)
	{
		angles[0] = angles[2] = Constants::pi / 2.0;
		angles[1] = random.decimal(Constants::pi/3, 2*Constants::pi/3);
		avgLength = pow(volume / sin(angles[1]), 1.0/3.0);
		lengths[0] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		avgLength = sqrt(volume / (lengths[0] * sin(angles[1])));
		lengths[1] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		lengths[2] = volume / (lengths[0]*lengths[1]*sin(angles[1]));
	}
	
	// Generate in triclinic system
	else
	{
		angles[0] = random.decimal(Constants::pi/3, 2*Constants::pi/3);
		angles[1] = random.decimal(Constants::pi/3, 2*Constants::pi/3);
		angles[2] = random.decimal(Constants::pi/3, 2*Constants::pi/3);
		avgLength = pow(volume / sqrt(1 - pow(cos(angles[0]), 2) - pow(cos(angles[1]), 2) - pow(cos(angles[2]), 2) + \
			2 * cos(angles[0]) * cos(angles[1]) * cos(angles[2])), 1.0/3.0);
		lengths[0] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		avgLength = sqrt(volume / sqrt(1 - pow(cos(angles[0]), 2) - pow(cos(angles[1]), 2) - pow(cos(angles[2]), 2) + \
			2 * cos(angles[0]) * cos(angles[1]) * cos(angles[2])) / lengths[0]);
		lengths[1] = avgLength * pow(10, random.decimalOnNormal(0.0, 0.15, 0.45));
		lengths[2] = volume / (sqrt(1 - pow(cos(angles[0]), 2) - pow(cos(angles[1]), 2) - pow(cos(angles[2]), 2) + \
			2 * cos(angles[0]) * cos(angles[1]) * cos(angles[2]))) / (lengths[0] * lengths[1]);
	}
	
	// Return new basis
	return Basis::vectors(lengths, angles);
}



/* void RandomStructure::generatePositions(ISO& iso, Random& random)
 *
 * Generate random positions within the structure
 */

void RandomStructure::generatePositions(ISO& iso, Random& random)
{
	
	// Get the number of atoms of each element that need to be assigned
	List<Atom*>::D2 atomsToAssign = getAtomsToAssign(iso);
	
	// Output
	Output::increase();
	
	// Loop over atoms in the structure
	int i, j, k, m, n;
	double distance;
	double curWorstDist;
	double minBondLength;
	double bestWorstDistance;
	Vector3D position;
	Vector3D bestPosition;
	for (i = 0; i < atomsToAssign.length(); ++i)
	{
		for (j = 0; j < atomsToAssign[i].length(); ++j)
		{
			
			// Loop until position is found outside atom radius or _maxTrialLoops loops is reached
			for (k = 0; k < _maxTrialLoops; ++k)
			{
				
				// Generate random position
				for (m = 0; m < 3; ++m)
					position[m] = random.decimal(0, 1);
				
				// Loop over atoms and check distance to those that are defined
				curWorstDist = 0;
				for (m = 0; m < iso.atoms().length(); ++m)
				{
					minBondLength = (atomsToAssign[i][j]->element().radius() + iso.atoms()[m][0].element().radius())/2;
					for (n = 0; n < iso.atoms()[m].length(); ++n)
					{
						
						// Skip if atom is not assigned
						if (!iso.atoms()[m][n].assigned())
							continue;
						
						// Get difference between bond length and target
						distance = iso.basis().distance(position, FRACTIONAL, iso.atoms()[m][n].fractional(), \
							FRACTIONAL) - minBondLength;
						
						// Get difference between bond length and target
						if (distance < curWorstDist)
							curWorstDist = distance;
					}
				}
				
				// Found a position outside atomic distances
				if (curWorstDist >= 0)
				{
					bestPosition = position;
					break;
				}
				
				// Save distance if new best or first checked
				if ((!k) || (curWorstDist > bestWorstDistance))
				{
					bestWorstDistance = curWorstDist;
					bestPosition = position;
				}
			}
			
			// Assign atom to best position
			atomsToAssign[i][j]->fractional(bestPosition);
			
			// Output
			Output::newline();
			Output::print("Placing atom ");
			Output::print(atomsToAssign[i][j]->atomNumber() + 1);
			Output::print(" at ");
			Output::print(atomsToAssign[i][j]->fractional(), 8, false);
		}
	}
	
	// Ouput
	Output::decrease();
}



/* void RandomStructure::generatePositions(ISO& iso, Random& random, const SpaceGroup& spaceGroup, double bias,
 *		Symmetry* symmetry)
 *
 * Generate random positions within space group
 */

void RandomStructure::generatePositions(ISO& iso, Random& random, const SpaceGroup& spaceGroup, double bias, \
	Symmetry* symmetry)
{
	
	// Get the number of atoms of each element that need to be assigned
	List<Atom*>::D2 atomsToAssign = getAtomsToAssign(iso);
	
	// Output
	Output::increase();
	
	// Get the minimum radius of all atoms
	int i;
	double minRadius = 0;
	if (iso.atoms().length() > 0)
		minRadius = iso.atoms()[0][0].element().radius();
	for (i = 1; i < iso.atoms().length(); ++i)
	{
		if (iso.atoms()[i][0].element().radius() < minRadius)
			minRadius = iso.atoms()[i][0].element().radius();
	}
	
	// Get list of Wyckoff sites grouped by multiplicity
	List<Wyckoff*>::D2 wyckoffGroups = RandomStructure::getWyckoffGroups(iso, spaceGroup, minRadius);
	
	// Get list of possible combinations of Wyckoff sites
	Linked<List<int>::D2> allowedGroups; 
	getAllowedWyckoff(allowedGroups, atomsToAssign, wyckoffGroups);
	
	// Variable to store list of positions generated at each trial
	List<Vector3D >::D2 curPositions(atomsToAssign.length());
	for (i = 0; i < atomsToAssign.length(); ++i)
		curPositions[i].length(atomsToAssign[i].length());
	
	// Count the number of sites occupied in each set
	int j, k;
	double minSites = 1;
	double maxSites = 1;
	Linked<List<int>::D2>::iterator itSet;
	List<double> ratios(allowedGroups.length());
	for (itSet = allowedGroups.begin(), i = 0; itSet != allowedGroups.end(); ++itSet, ++i)
	{
		ratios[i] = 0;
		for (j = 0; j < (*itSet).length(); ++j)
		{
			for (k = 0; k < (*itSet)[j].length(); ++k)
				ratios[i] += (*itSet)[j][k];
		}
		if ((ratios[i] < minSites) || (i == 0))
			minSites = ratios[i];
		if ((ratios[i] > maxSites) || (i == 0))
			maxSites = ratios[i];
	}
	
	// Print results so far
	Output::newline();
	Output::print("Number of combinations of Wyckoff site multiplicities: ");
	Output::print(allowedGroups.length());
	Output::newline();
	Output::print("Minimum number of Wyckoff sites that could be occupied: ");
	Output::print(minSites);
	Output::newline();
	Output::print("Maximum number of Wyckoff sites that could be occupied: ");
	Output::print(maxSites);
	
	// Turn counts into ratios
	for (i = 0; i < ratios.length(); ++i)
		ratios[i] /= minSites;
	
	// Organize sites by bin
	int curBin;
	Linked<int> binNumbers;
	Linked<Linked<List<int>::D2*> > bins;
	Linked<int>::iterator itBN;
	Linked<Linked<List<int>::D2*> >::iterator itB;
	for (i = 0, itSet = allowedGroups.begin(); i < ratios.length(); ++i, ++itSet)
	{
		curBin = weightBin(ratios[i]);
		for (itBN = binNumbers.begin(), itB = bins.begin(); itB != bins.end(); ++itBN, ++itB)
		{
			if (curBin == *itBN)
			{
				*itB += &(*itSet);
				break;
			}
		}
		if (itB == bins.end())
		{
			binNumbers += curBin;
			bins.add();
			*(bins.last()) += &(*itSet);
		}
	}
	
	// Set weights for each bin
	double avgWeight = 0;
	List<double> weights(bins.length());
	for (i = 0, itBN = binNumbers.begin(); i < bins.length(); ++i, ++itBN)
	{
		weights[i] = weight(*itBN);
		avgWeight += weights[i];
	}
	
	// Set biasing
	avgWeight /= weights.length();
	for (i = 0; i < weights.length(); ++i)
	{
		weights[i] = (1 - bias)*avgWeight + bias*weights[i];
		if (i > 0)
			weights[i] += weights[i-1];
	}
		
	// Loop until a set of sites is generated in which all atoms are outside cutoff or until max loops is reached
	int m, n, p;
	int index;
	int site;
	double rand;
	double distance;
	double curWorstDist;
	double minBondLength;
	double bestWorstDistance;
	Vector3D position;
	List<Wyckoff*>::D2 curWycks(atomsToAssign.length());
	List<Wyckoff*>::D2 bestWycks;
	OList<Vector3D >::D2 curUniquePos(atomsToAssign.length());
	OList<Vector3D >::D2 bestUniquePos;
	List<Wyckoff*>::D2 remainingGroups;
	Linked<Linked<List<int>::D2*> >::iterator itBin;
	for (i = 0; i < _maxTrialLoops; ++i)
	{
		
		// Pick random bin
		rand = random.decimal(0, weights.last());
		for (site = 0; site < weights.length(); ++site)
		{
			if (rand < weights[site])
				break;
		}
		if (site >= weights.length())
			continue;
	
		// Pick random set within bin
		itBin = bins.begin() + site;
		List<int>::D2& set = **((*itBin).begin() + random.integer(0, (*itBin).length() - 1));
		
		// Loop 10 times over current trial
		for (j = 0; (j < 10) && (i < _maxTrialLoops); ++i, ++j)
		{
			
			// Reset list of groups that can be selected
			remainingGroups = wyckoffGroups;

			// Clear space
			for (k = 0; k < curUniquePos.length(); ++k)
			{
				curWycks[k].length(0);
				curUniquePos[k].length(0);
			}
			
			// Loop over unique elements of atoms to be added
			for (k = 0; k < set.length(); ++k)
			{

				// Loop over wyckoff sites
				index = -1;
				for (m = 0; m < set[k].length(); ++m)
				{

					// Loop over number of occurances to add of current wyckoff site
					for (n = 0; n < set[k][m]; ++n)
					{

						// Pick site at random from current group of wyckoff sites
						site = random.integer(0, remainingGroups[m].length() - 1);

						// Generate random position in cell
						for (p = 0; p < 3; ++p)
							position[p] = random.decimal(0, 1);

						// Save position and wyckoff site that was chosen
						curWycks[k] += remainingGroups[m][site];
						curUniquePos[k] += position;

						// Loop over orbit of wyckoff site
						for (p = 0; p < remainingGroups[m][site]->multiplicity(); ++p)
						{

							// Generate position and move into cell
							curPositions[k][++index] = remainingGroups[m][site]->rotations()[p] * position + \
								remainingGroups[m][site]->translations()[p];
							ISO::moveIntoCell(curPositions[k][index]);
						}

						// If site is limited then remove
						if (remainingGroups[m][site]->rank() == 0)
							remainingGroups[m].remove(site);
					}
				}
			}
			
			// Loop over atoms that are going to be added
			curWorstDist = 0;
			for (k = 0; k < curPositions.length(); ++k)
			{
				for (m = 0; m < curPositions[k].length(); ++m)
				{

					// Loop over assigned atoms in structure
					for (n = 0; n < iso.atoms().length(); ++n)
					{
						minBondLength = _minBondFraction * (atomsToAssign[k][m]->element().radius() + \
							iso.atoms()[n][0].element().radius());
						for (p = 0; p < iso.atoms()[n].length(); ++p)
						{
							if (!iso.atoms()[n][p].assigned())
								continue;
							distance = iso.basis().distance(curPositions[k][m], FRACTIONAL, \
								iso.atoms()[n][p].fractional(), FRACTIONAL) - minBondLength;
							if (distance < curWorstDist)
								curWorstDist = distance;
						}
					}

					// Loop over new atoms
					for (n = 0; n < curPositions.length(); ++n)
					{
						minBondLength = _minBondFraction * (atomsToAssign[k][m]->element().radius() + \
							atomsToAssign[n][0]->element().radius());
						for (p = 0; p < curPositions[n].length(); ++p)
						{
							if ((k == n) && (m == p))
								continue;
							distance = iso.basis().distance(curPositions[k][m], FRACTIONAL, curPositions[n][p], \
								FRACTIONAL) - minBondLength;
							if (distance < curWorstDist)
								curWorstDist = distance;
						}
					}
				}
			}

			// Found a position that preserves atomic distances
			if (curWorstDist >= 0)
			{
				bestWycks = curWycks;
				bestUniquePos = curUniquePos;
				break;
			}

			// Save distance if new best or first checked
			if (((!i) && (!j)) || (curWorstDist > bestWorstDistance))
			{
				bestWycks = curWycks;
				bestWorstDistance = curWorstDist;
				bestUniquePos = curUniquePos;
			}
		}
		
		// Break if finished
		if (curWorstDist >= 0)
			break;
	}
	
	// Assign positions
	bool found;
	Vector3D tempTranslation;
	for (i = 0; i < bestUniquePos.length(); ++i)
	{
		
		// Loop over unique atoms of current element
		index = -1;
		for (j = 0; j < bestUniquePos[i].length(); ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Placing atoms on Wyckoff site ");
			Output::print(bestWycks[i][j]->name());
			Output::increase();
			
			// Loop over orbit of position
			Orbit curOrbit;
			for (k = 0; k < bestWycks[i][j]->multiplicity(); ++k)
			{
				
				// Save position
				atomsToAssign[i][++index]->fractional(bestWycks[i][j]->rotations()[k] * bestUniquePos[i][j] + \
					bestWycks[i][j]->translations()[k]);
				
				// Output
				Output::newline();
				Output::print("Atom ");
				Output::print(atomsToAssign[i][index]->atomNumber() + 1);
				Output::print(" at ");
				Output::print(atomsToAssign[i][index]->fractional(), 8, false);
				
				// Set first atom in orbit
				if (!k)
					curOrbit.set(atomsToAssign[i][index], bestWycks[i][j]->rotations()[k], \
						bestWycks[i][j]->translations()[k]);
				
				// Add atom to orbit
				else
				{
					
					// Loop over symmetry operations to find which maps special positions
					found = false;
					for (m = 0; m < spaceGroup.symmetry().length(); ++m)
					{
						
						// Skip if rotations do not match
						if (!(spaceGroup.symmetry()[m].rotation()*bestWycks[i][j]->rotations()[0] == \
							bestWycks[i][j]->rotations()[k]))
							continue;
						
						// Loop over translations
						for (n = 0; n < spaceGroup.symmetry()[m].translations().length(); ++n)
						{
							
							// Generate current translation
							tempTranslation = spaceGroup.symmetry()[m].rotation()*bestWycks[i][j]->translations()[0] +\
								spaceGroup.symmetry()[m].translations()[n];
							ISO::moveIntoCell(tempTranslation);
							
							// Check if translation is a match
							if (!(tempTranslation == bestWycks[i][j]->translations()[k]))
								continue;
							
							// Found match
							curOrbit.add(atomsToAssign[i][index], spaceGroup.symmetry()[m].rotation(), \
								spaceGroup.symmetry()[m].translations()[n]);
							found = true;
							break;
						}
						if (found)
							break;
					}
					
					// Error if not matched
					if (!found)
					{
						Output::newline(ERROR);
						Output::print("Internal error: could not match special position");
						Output::quit();
					}
				}
			}
			
			// Add orbit to symmetry
			if (symmetry)
				symmetry->addOrbit(curOrbit);
			
			// Output
			Output::decrease();
		}
	}
	
	// Output
	Output::decrease();
}



/* List<Atom*>::D2 RandomStructure::getAtomsToAssign(ISO& iso)
 *
 * Get list of atoms to assign to random positions
 */

List<Atom*>::D2 RandomStructure::getAtomsToAssign(ISO& iso)
{
	
	// Loop over elements
	int i, j;
	bool found;
	List<Atom*>::D2 res;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Loop over atoms of current element
		found = false;
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Current atom has not been assigned
			if (!iso.atoms()[i][j].assigned())
			{
				
				// First time this element has been found
				if (!found)
				{
					res.add();
					found = true;
				}
				
				// Add atom to list
				res.last() += &iso.atoms()[i][j];
			}
		}
	}
	
	// Output
	Output::newline();
	Output::print("Assigning random position");
	if ((res.length() > 1) || (res[0].length() > 1))
		Output::print("s");
	Output::print(" to ");
	for (i = 0; i < res.length(); ++i)
	{
		Output::print(res[i].length());
		Output::print(" ");
		Output::print(res[i][0]->element().symbol());
		Output::print(" atom");
		if (res[i].length() != 1)
			Output::print("s");
		if ((res.length() > 2) && (i != res.length() - 1))
			Output::print(",");
		if (i == res.length() - 2)
			Output::print(" and");
		if (i != res.length() - 1)
			Output::print(" ");
	}
	
	// Return list of atoms
	return res;
}



/* List<Wyckoff*>::D2 RandomStructure::getWyckoffGroups(const ISO& iso, const SpaceGroup& spaceGroup, double minRadius)
 *
 * Return list of Wyckoff sites in space group that is grouped by multiplicity
 */

List<Wyckoff*>::D2 RandomStructure::getWyckoffGroups(const ISO& iso, const SpaceGroup& spaceGroup, double minRadius)
{
	
	// Loop over wyckoff sites in space group
	int i, j;
	bool found;
	List<Wyckoff*>::D2 res;
	for (i = 0; i < spaceGroup.wyckoff().length(); ++i)
	{
		
		// Check if multiplicity is saved
		found = false;
		for (j = 0; j < res.length(); ++j)
		{
			
			// Found group with same multiplicity
			if (spaceGroup.wyckoff()[i].multiplicity() == res[j][0]->multiplicity())
			{
				found = true;
				res[j] += &spaceGroup.wyckoff()[i];
				break;
			}
		}
		
		// Found a new multiplicity
		if (!found)
		{
			res.add();
			res.last() += &spaceGroup.wyckoff()[i];
		}
	}
	
	// Loop over groups and remove any that are limited and already filled by an atom in the structure
	int k, m, n;
	for (i = res.length() - 1; i >= 0; --i)
	{
		
		// Loop over sites in current group
		for (j = res[i].length() - 1; j >= 0; --j)
		{
			
			// Skip if not limited
			if (res[i][j]->rank() != 0)
				continue;
				
			// Loop over positions in orbit
			found = false;
			for (k = 0; k < res[i][j]->multiplicity(); ++k)
			{
				
				// Loop over atoms in the structure
				for (m = 0; m < iso.atoms().length(); ++m)
				{
					for (n = 0; n < iso.atoms()[m].length(); ++n)
					{

						// Skip if atom is not assigned
						if (!iso.atoms()[m][n].assigned())
							continue;

						// Atom is within range of site
						if (iso.basis().distance(iso.atoms()[m][n].fractional(), FRACTIONAL, \
							res[i][j]->translations()[k], FRACTIONAL) < minRadius)
						{
							res[i].remove(j);
							found = true;
							break;
						}
					}
					if (found)
						break;
				}
				if (found)
					break;
			}
		}
		
		// Remove group if empty
		if (!res[i].length())
			res.remove(i);
	}
	
	// Sort groups
	sortWyckoffGroups(res, 0, res.length() - 1);
	
	// Return sorted list
	return res;
}



/* void RandomStructure::sortWyckoffGroups(List<Wyckoff*>::D2& groups, int left, int right)
 *
 * Sort Wyckoff groups by multiplicity
 */

void RandomStructure::sortWyckoffGroups(List<Wyckoff*>::D2& groups, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	int pivot = groups[pivotIndex][0]->multiplicity();
	
	// Move pivot to end
	groups.swap(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (groups[i][0]->multiplicity() < pivot)
		{
			groups.swap(i, newPivotIndex);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	groups.swap(newPivotIndex, right);
	
	// Recursive calls to next sorts
	sortWyckoffGroups(groups, left, newPivotIndex - 1);
	sortWyckoffGroups(groups, newPivotIndex + 1, right);
}



/* void RandomStructure::getAllowedWyckoff(Linked<List<int>::D2>& allowedGroups, const List<Atom*>::D2& atomsToAssign,
 *		const List<Wyckoff*>::D2 wyckoffGroups)
 *
 * Return groups of wyckoff sites that can be used to populate atoms
 */

void RandomStructure::getAllowedWyckoff(Linked<List<int>::D2>& allowedGroups, const List<Atom*>::D2& atomsToAssign, \
	const List<Wyckoff*>::D2 wyckoffGroups)
{
	
	// Check which groups of wyckoff sites are completely limited
	int i, j;
	List<bool> limited(wyckoffGroups.length());
	limited.fill(true);
	for (i = 0; i < wyckoffGroups.length(); ++i)
	{
		for (j = 0; j < wyckoffGroups[i].length(); ++j)
		{
			if (wyckoffGroups[i][j]->rank() > 0)
			{
				limited[i] = false;
				break;
			}
		}
	}
	
	// Recursively generate combinations
	List<int> numUsed(limited.length());
	numUsed.fill(0);
	List<int> tempList(wyckoffGroups.length());
	tempList.fill(0);
	List<int>::D2 startList(atomsToAssign.length());
	startList.fill(tempList);
	recurseAllowedWyckoff(allowedGroups, atomsToAssign, wyckoffGroups, limited, numUsed, startList, 0, 0);
	
	// Did not find any allowed combinations
	if (!allowedGroups.length())
	{
		Output::newline(ERROR);
		Output::print("The number of atoms to be assigned is not commensurate with space group");
		Output::quit();
	}
}



/* void RandomStructure::recurseAllowedWyckoff(Linked<List<int>::D2>& res, const List<Atom*>::D2& atomsToAssign,
 *		const List<Wyckoff*>::D2 wyckoffGroups, const List<bool>& limited, const List<int>& numUsed, 
 *		const List<int>& curList, int curAtom, int curGroup)
 *
 * Recursively generate list of allowed combinations of wyckoff sites (by multiplicity)
 */

void RandomStructure::recurseAllowedWyckoff(Linked<List<int>::D2>& res, const List<Atom*>::D2& atomsToAssign, \
	const List<Wyckoff*>::D2 wyckoffGroups, const List<bool>& limited, const List<int>& numUsed, \
	const List<int>::D2& curList, int curAtom, int curGroup)
{
	
	// Count total number of atoms used for current element
	int i;
	int total = 0;
	for (i = 0; i < wyckoffGroups.length(); ++i)
		total += curList[curAtom][i] * wyckoffGroups[i][0]->multiplicity();
	
	// Return if out of range (as a safe guard)
	if (total > atomsToAssign[curAtom].length())
		return;
	
	// Have the correct number of atoms
	int j;
	int newAtom = curAtom;
	if (total == atomsToAssign[curAtom].length())
	{
		
		// On last atom
		if (curAtom == atomsToAssign.length() - 1)
		{
			res.add(curList);
			return;
		}
		
		// Go to next atom
		++newAtom;
		curGroup = 0;
		total = 0;
	}
	
	// Loop over remaining wyckoff groups
	int max;
	List<int>::D2 newList;
	List<int> newNumUsed;
	for (i = curGroup; i < wyckoffGroups.length(); ++i)
	{
		
		// Figure out the maximum number of current type that can be added
		max = (atomsToAssign[newAtom].length() - total) / wyckoffGroups[i][0]->multiplicity();
		if ((limited[i]) && (wyckoffGroups[i].length() - numUsed[i] < max))
			max = wyckoffGroups[i].length() - numUsed[i];
		
		// Loop over possible values
		for (j = i == curGroup ? 1 : 0; j <= max; ++j)
		{
			
			// Set new list
			newList = curList;
			newList[newAtom][i] = j;
			newNumUsed = numUsed;
			newNumUsed[i] += j;
			
			// Call next iteration
			recurseAllowedWyckoff(res, atomsToAssign, wyckoffGroups, limited, newNumUsed, newList, newAtom, i + 1);
		}
	}
}



/* void RandomStructure::setSiteSymmetry(const ISO& iso, const SpaceGroup& spaceGroup, Symmetry* symmetry)
 *
 * Set symmetry of sites that are already in structure
 */

void RandomStructure::setSiteSymmetry(const ISO& iso, const SpaceGroup& spaceGroup, Symmetry* symmetry)
{
	
	// Loop over elements
	int i, j, k;
	bool found;
	int count;
	double tol = 1e-4;
	double curDistance;
	Atom curAtom;
	Vector3D origin(0.0);
	Vector3D curPos;
	Vector3D curSpecialTranslation;
	Matrix3D curSpecialRotation;
	Linked<bool> used;
	Linked<bool>::iterator itUsed;
	Linked<int> curGens;
	Linked<int>::iterator itGens;
	Linked<int> curTrans;
	Linked<int>::iterator itTrans;
	Linked<Atom*> atoms;
	Linked<Atom*>::iterator itAtom;
	Linked<Atom*> curOrbit;
	Linked<Atom*>::iterator itOrbit;
	Linked<double> distances;
	Linked<double>::iterator itDist;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Create list of atoms of current element
		used.clear();
		atoms.clear();
		distances.clear();
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			if (!iso.atoms()[i][j].assigned())
				continue;
			atoms.add(&(iso.atoms()[i][j]));
			distances.add(iso.basis().distance(iso.atoms()[i][j].fractional(), FRACTIONAL, origin, FRACTIONAL));
			used.add(false);
		}
		
		// Loop until no atoms are left
		while (atoms.length())
		{
					
			// Set current atom
			curAtom  = **(atoms.begin());
			(*used.begin()) = true;
			
			// Loop over symmetry operations
			count = 0;
			curSpecialRotation = 0.0;
			curSpecialTranslation = 0.0;
			curGens.clear();
			curTrans.clear();
			curOrbit.clear();
			for (j = 0; j < spaceGroup.symmetry().length(); ++j)
			{
				for (k = 0; k < spaceGroup.symmetry()[j].translations().length(); ++k)
				{
					
					// Generate new position
					curPos = spaceGroup.symmetry()[j].rotation() * (*atoms.begin())->fractional();
					curAtom.fractional(curPos + spaceGroup.symmetry()[j].translations()[k]);
					curDistance = iso.basis().distance(curAtom.fractional(), FRACTIONAL, origin, FRACTIONAL);
					
					// Loop over atoms of current element
					found = false;
					itUsed = used.begin();
					itAtom = atoms.begin();
					itDist = distances.begin();
					for (; itAtom != atoms.end(); ++itUsed, ++itAtom, ++itDist)
					{
						
						// Check if atoms are the same
						if (Num<double>::abs(curDistance - *itDist) > tol)
							continue;
						if (!curAtom.equal(**itAtom, tol))
							continue;
						
						// Returned same atom
						if (itAtom == atoms.begin())
						{
							++count;
							curSpecialRotation += spaceGroup.symmetry()[j].rotation();
							curSpecialTranslation += (*atoms.begin())->fractional() - curPos;
						}
						
						// Found a new atom
						else if (!(*itUsed))
						{
							curGens.add(j);
							curTrans.add(k);
							curOrbit.add(*itAtom);
							*itUsed = true;
						}
						
						// Break since atom was found
						found = true;
						break;
					}
					
					// No match
					if (!found)
					{
						Output::newline(ERROR);
						Output::print("Atoms in structure do not fit space group symmetry");
						Output::quit();
					}
				}
			}
			
			// Save orbit if symmetry was passed
			if (symmetry)
			{
				
				// Add new orbit
				Orbit orbit;
				curSpecialRotation /= count;
				curSpecialTranslation /= count;
				orbit.set(*(atoms.begin()), curSpecialRotation, curSpecialTranslation);

				// Save atoms in orbit
				itGens = curGens.begin();
				itTrans = curTrans.begin();
				for (itOrbit = curOrbit.begin(); itOrbit != curOrbit.end(); ++itOrbit, ++itGens, ++itTrans)
					orbit.add(*itOrbit, spaceGroup.symmetry()[*itGens].rotation(), \
						spaceGroup.symmetry()[*itGens].translations()[*itTrans]);

				// Add orbit
				symmetry->addOrbit(orbit);
			}
			
			// Loop until no atoms are removed
			do
			{
				found = false;
				itUsed = used.begin();
				itAtom = atoms.begin();
				itDist = distances.begin();
				for (; itAtom != atoms.end(); ++itUsed, ++itAtom, ++itDist)
				{
					if (*itUsed)
					{
						used.remove(itUsed);
						atoms.remove(itAtom);
						distances.remove(itDist);
						found = true;
						break;
					}
				}
			} while (found);
		}
	}
}
