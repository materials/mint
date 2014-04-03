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



#include "kmc.h"
#include "symmetry.h"
#include "unique.h"
#include "language.h"
#include <cmath>
#include <cstdlib>



/* void KMC::generateJumps(const ISO& iso, const Element& element, bool useInterstitials, double minImageDis,
 *		double tol, StructureFormat format, double maxJumpDistance)
 * 
 * Generate jumps file for kmc simulation
 */

void KMC::generateJumps(const ISO& iso, const Element& element, bool useInterstitials, double minImageDis, \
	double tol, StructureFormat format, double maxJumpDistance)
{
	
	// Output
	Output::newline();
	Output::print("Generating nodes and jumps for kmc simulation of diffusion");
	Output::increase();
	
	// Output
	Output::newline();
	Output::print("Searching for unique jumps between sites in the primitive cell");
	Output::increase();
	
	// Get the transformation to reduced primitive cell
	Matrix3D unitToPrim = iso.primitiveTransformation(tol, true);
	Matrix3D unitToRedPrim = Basis::reducedTransformation(unitToPrim * iso.basis().vectors()) * unitToPrim;
	
	// Get the primitive cell and symmetry information
	ISO prim = iso;
	Symmetry primSymmetry;
	prim.transform(unitToRedPrim, tol, true);
	primSymmetry.set(prim, tol, false);
	
	// Get the unique jumps in the primitive cell
	Unique jumps;
	jumps.getEquivalent(false);
	jumps.areTransitions(true);
	jumps.useInterstitials(useInterstitials);
	jumps.maxDistance(maxJumpDistance);
	jumps.search(prim, primSymmetry, OList<Element>(element), tol);
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Generating the ideal cell and searching for all jumps within it");
	Output::increase();
	
	// Get the max jump length that was generated
	int i;
	double maxLength = -1;
	for (i = 0; i < jumps.groups().length(); ++i)
	{
		if (jumps.groups()[i][0].distances()[1] > maxLength)
			maxLength = jumps.groups()[i][0].distances()[1];
	}
	
	// The minimimum image distance should be the greater of minImageDis and maxLength + 0.5
	maxLength += 0.5;
	if (minImageDis > maxLength)
		maxLength = minImageDis;
	
	// Generate ideal cell and get its symmetry
	ISO ideal = iso;
	Matrix3D unitToIdeal = Symmetry::idealTransformation(ideal, maxLength, tol);
	ideal.transform(unitToIdeal, tol);
	ideal.orderAtomNumbers();
	Symmetry idealSymmetry (ideal, tol);
	
	// Get the unique sites in the ideal cell
	Unique sites;
	sites.getEquivalent(true);
	sites.areTransitions(false);
	sites.useInterstitials(useInterstitials);
	sites.search(ideal, idealSymmetry, OList<Element>(element), tol);
	
	// Get the unique jumps in the ideal cell
	jumps.clear();
	jumps.getEquivalent(true);
	jumps.search(ideal, idealSymmetry, OList<Element>(element), tol);
	
	// No jumps were found
	if (jumps.groups().length() == 0)
	{
		Output::newline(WARNING);
		Output::print("No jumps were found for KMC simulation");
	}
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Saving the conversion from ideal to conventional cell");
	Output::increase();
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Creating KMC jumps file and generating structures");
	Output::increase();
	
	// Loop over jumps and figure out which ones to use to minimize the number of calculations
	int j, k;
	bool found;
	List<int> save;
	List<int>::D2 jumpTypes;
	jumpTypes.length(jumps.groups().length());
	for (i = 0; i < jumps.groups().length(); ++i)
	{
		jumpTypes[i] += i;
		jumpTypes[i] += save.length();
		jumpTypes[i] += 1;
		for (j = 0, found = false; j < save.length(); ++j)
		{
			for (k = 0; k < jumps.groups()[save[j]].length(); ++k)
			{
				
				// Found reverse jump already saved
				if ((jumps.groups()[i][0].atoms()[0] == jumps.groups()[save[j]][k].atoms()[1]) && \
					(jumps.groups()[i][0].atoms()[1] == jumps.groups()[save[j]][k].atoms()[0]))
				{
					jumpTypes[i][0] = save[j];
					jumpTypes[i][1] = jumpTypes[save[j]][1];
					jumpTypes[i][2] = -1;
					found = true;
					break;
				}
			}
			if (found)
				break;
		}
		if (!found)
			save += i;
	}
	
	// Print unique sites and create structures
	Word curDir;
	Output message;
	message.addLines(sites.groups().length());
	List<Atom*> atomsToRemove;
	for (i = 0; i < sites.groups().length(); ++i)
	{
		
		// Save site to output message
		message.addLine();
		message.add("SITE");
		message.add(i+1);
		message.add(" ");
		message.add("[[initial structure file]]");
		message.add("[[relaxed structure file]]");
		
		// Create the structure
		ISO tempISO = ideal;
		if (!useInterstitials)
			tempISO.removeAtom(sites.groups()[i][0].atoms()[0]->atomNumber(), false);
		else
		{
			
			// Get list of atoms to remove
			atomsToRemove.length(0);
			for (j = 0; j < sites.groups().length(); ++j)
			{
				for (k = 0; k < sites.groups()[j].length(); ++k)
				{
					if ((j != i) || (k != 0))
						atomsToRemove += tempISO.atom(sites.groups()[j][k].atoms()[0]->atomNumber());
				}
			}
			
			// Remove atoms
			for (j = 0; j < atomsToRemove.length(); ++j)
				tempISO.removeAtom(atomsToRemove[j]->atomNumber(), false);
		}
		
		// Create directory
		curDir = "site_";
		curDir += Language::numberToWord(i+1);
		Directory::create(curDir, true);
		
		// Write structure to file
		StructureIO::write(Directory::makePath(curDir, "structure"), tempISO, format, FRACTIONAL);
	}
	
	// Print the pure structure
	ISO idealCopy = ideal;
	if (useInterstitials)
	{
		for (i = idealCopy.atoms().length() - 1; i >= 0; --i)
		{
			for (j = idealCopy.atoms()[i].length() - 1; j >= 0; --j)
			{
				if ((idealCopy.atoms()[i][j].element() == element) && (idealCopy.atoms()[i][j].isInterstitial()))
					idealCopy.removeAtom(idealCopy.atoms()[i][j].atomNumber(), false);
			}
		}
	}
	Directory::create("supercell", true);
	StructureIO::write(Directory::makePath("supercell", "structure"), idealCopy, format, FRACTIONAL);
	
	// Open file to print
	int origStream = Output::streamID();
	int kmcStream = Output::addStream("kmc.setup");
	Output::setStream(kmcStream);
	Output::newline();
	Output::print(message, RIGHT);
	Output::setStream(origStream);
	
	// Create unique jumps lines
	int m, n;
	int orbit;
	Matrix3D rot;
	Vector3D cent;
	Vector3D trans;
	message.clear();
	message.addLines(save.length());
	for (i = 0; i < save.length(); ++i)
	{
		
		// Initialize line
		message.addLine();
		message.addWords(30);
		message.add("JUMP");
		message.add(i+1);
		
		// Save the vector needed to center starting point
		cent = 0.5;
		cent -= jumps.groups()[save[i]][0].atoms()[0]->fractional() + (jumps.groups()[save[i]][0].vectors()[1] / 2);
		
		// Save information about the jump
		for (j = 0; j < 2; ++j)
		{
			
			// Find orbit of corresponding site
			orbit = idealSymmetry.orbitNumbers()[jumps.groups()[save[i]][0].atoms()[j]->atomNumber()];
			for (k = 0; k < sites.groups().length(); ++k)
			{
				if (orbit == idealSymmetry.orbitNumbers()[sites.groups()[k][0].atoms()[0]->atomNumber()])
				{
					
					// Save the site index
					message.add(k+1);
					
					// Loop over atoms in orbit to find site
					for (m = 0; m < idealSymmetry.orbits()[orbit].atoms().length(); ++m)
					{
						if (idealSymmetry.orbits()[orbit].atoms()[m] == sites.groups()[k][0].atoms()[0])
						{
							rot = idealSymmetry.orbits()[orbit].generators()[m].rotation().inverse();
							trans = idealSymmetry.orbits()[orbit].generators()[m].translations()[0] * -1;
							break;
						}
					}
					if (m >= idealSymmetry.orbits()[orbit].atoms().length())
					{
						Output::newline(ERROR);
						Output::print("Internal error - failed to match atom in orbit");
						Output::quit();
					}
					
					// Loop over atoms in orbit to find jump atom
					for (m = 0; m < idealSymmetry.orbits()[orbit].atoms().length(); ++m)
					{
						if (idealSymmetry.orbits()[orbit].atoms()[m] == jumps.groups()[save[i]][0].atoms()[j])
						{
							trans += rot * idealSymmetry.orbits()[orbit].generators()[m].translations()[0];
							rot = rot * idealSymmetry.orbits()[orbit].generators()[m].rotation();
							break;
						}
					}
					if (m >= idealSymmetry.orbits()[orbit].atoms().length())
					{
						Output::newline(ERROR);
						Output::print("Internal error - failed to match jump atom in orbit");
						Output::quit();
					}
					
					// Save translation
					trans += cent;
					ISO::moveIntoCell(trans);
					
					// Print operation
					for (m = 0; m < 3; ++m)
					{
						for (n = 0; n < 3; ++n)
							message.add(Language::numberToFraction(rot(m, n)));
					}
					for (m = 0; m < 3; ++m)
						message.add(trans[m], 10);
					
					// Break since match was found
					break;
				}
			}
		}
	}
	
	// Print jumps
	Output::setStream(kmcStream);
	Output::newline();
	Output::newline();
	Output::print(message, RIGHT);
	Output::setStream(origStream);
	
	// Print type of simulation
	Output::setStream(kmcStream);
	Output::newline();
	Output::newline();
	if (!useInterstitials)
		Output::print("complement");
	Output::setStream(origStream);
	
	// Print basis
	Output::setStream(kmcStream);
	Output::newline();
	for (i = 0; i < 3; ++i)
	{
		Output::print(ideal.basis().lengths()[i]);
		Output::print(" ");
	}
	for (i = 0; i < 3; ++i)
	{
		Output::print(Num<double>::toDegrees(ideal.basis().angles()[i]));
		if (i != 2)
			Output::print(" ");
	}
	Output::setStream(origStream);
	
	// Print conversion from ideal to unit cell
	Matrix3D idealToUnit = unitToIdeal.inverse();
	Output::setStream(kmcStream);
	Output::newline();
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			Output::print(Language::numberToFraction(idealToUnit(i, j)));
			if ((i != 2) || (j != 2))
				Output::print(" ");
		}
	}
	Output::setStream(origStream);
	
	// Print conversion from unit cell to reduced primitive cell
	Output::setStream(kmcStream);
	Output::newline();
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			Output::print(Language::numberToFraction(unitToRedPrim(i, j)));
			if ((i != 2) || (j != 2))
				Output::print(" ");
		}
	}
	Output::setStream(origStream);
	
	// Print all of the jumps that are allowed in the ideal cell
	Output::setStream(kmcStream);
	for (i = 0; i < sites.groups().length(); ++i)
	{
		for (j = 0; j < sites.groups()[i].length(); ++j)
		{
			Output::newline();
			Output::print(i+1);
			Output::print(" ");
			Output::print(sites.groups()[i][j].atoms()[0]->fractional(), 14, false);
			message.clear();
			for (k = 0; k < jumps.groups().length(); ++k)
			{
				for (m = 0; m < jumps.groups()[k].length(); ++m)
				{
					if (sites.groups()[i][j].atoms()[0] == jumps.groups()[k][m].atoms()[0])
					{
						message.addLine();
						message.addWords(5);
						message.add(jumpTypes[k][1]+1);
						message.add(jumpTypes[k][2]);
						for (n = 0; n < 3; ++n)
							message.add(jumps.groups()[k][m].vectors()[1][n], 14);
					}
				}
			}
			Output::newline();
			Output::print(message, RIGHT);
		}
	}
	Output::setStream(origStream);
	
	// Close kmc file
	Output::removeStream(kmcStream);
	
	// Output
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* bool KMC::set(const Text& input)
 *
 * Read data from KMC file and return whether simulation is ready to run
 */

bool KMC::set(const Text& input)
{
	
	// Output
	Output::newline();
	Output::print("Setting properties of KMC simulation");
	Output::increase();
	
	// Look for a jump line in the file
	for (int i = 0; i < input.length(); ++i)
	{
		if (input[i].length() > 0)
		{
			if (input[i][0].equal("jump", false))
			{
				
				// This is a setup file so generate jump structures
				if (input[i].length() > 20)
				{
					setJumps(input);
					Output::decrease();
					return false;
				}
				
				// This is a full simulation file so save properties
				setSimulation(input);
				Output::decrease();
				return true;
			}
		}
	}
	
	// Did not find a jump line
	Output::newline(ERROR);
	Output::print("Did not find a jumps line in file");
	Output::print("    Cannot continue with KMC setup");
	Output::quit();
	Output::decrease();
	return false;
}



/* void KMC::setJumps(const Text& input)
 *
 * Intermediate step of KMC setup
 * Generate structures to get unique barriers
 */

void KMC::setJumps(const Text& input)
{
	
	// Output
	Output::newline();
	Output::print("Generating structures to describe unique jumps");
	Output::increase();
	
	// Get the transformations between cells
	getTransformations(input);
	
	// Output
	Output::newline();
	Output::print("Saving files containing structures");
	Output::increase();
	
	// Loop over lines in file to get site structure files
	int i;
	int index;
	Output message;
	OList<Word>::D2 siteFiles;
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Check if a site line
		if (input[i][0].equal("site", false))
		{
			
			// Make sure that followed by three words
			if (input[i].length() < 4)
			{
				Output::newline(ERROR);
				Output::print("Line ");
				Output::print(i+1);
				Output::print(" of KMC file does not contain enough information");
				Output::quit();
			}
			
			// Check that second character is an integer
			if (!Language::isInteger(input[i][1]))
			{
				Output::newline(ERROR);
				Output::print("Expecting an integer at second place of line ");
				Output::print(i+1);
				Output::print(" of KMC file (found ");
				Output::print(input[i][1]);
				Output::print(")");
				Output::quit();
			}
			
			// Save index
			index = atoi(input[i][1].array()) - 1;
			
			// Save files
			if (index >= siteFiles.length())
				siteFiles.length(index+1);
			siteFiles[index].length(2);
			siteFiles[index][0] = input[i][2];
			siteFiles[index][1] = input[i][3];
			
			// Save line for file
			message.addLine();
			message.add("SITE");
			message.add(index + 1);
			message.add(" ");
			message.add("[[Formation energy]]");
			message.add("[[Total energy]]");
			message.add("[[Name of file with phonon frequencies - or remove]]");
			
			// Output
			Output::newline();
			Output::print("Site ");
			Output::print(index + 1);
			Output::print(": ");
			Output::print(siteFiles[index][0]);
			Output::print(" and ");
			Output::print(siteFiles[index][1]);
		}
	}
	
	// Open file to print
	int origStream = Output::streamID();
	int kmcStream = Output::addStream("kmc.in");
	Output::setStream(kmcStream);
	Output::newline();
	Output::print(message, RIGHT);
	Output::setStream(origStream);
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Reading initial and final structures from file");
	Output::increase();
	
	// Save structures for each unique site
	Text content;
	StructureFormat format = SF_UNKNOWN;
	OList<ISO>::D2 siteISO(siteFiles.length());
	for (i = 0; i < siteFiles.length(); ++i)
	{
		if (siteFiles[i].length() < 2)
			continue;
		siteISO[i].length(2);
		Output::newline();
		Output::print("Reading initial structure for site ");
		Output::print(i+1);
		Output::increase();
		siteISO[i][0] = StructureIO::read(siteFiles[i][0]);
		Output::decrease();
		Output::newline();
		Output::print("Reading relaxed structure for site ");
		Output::print(i+1);
		Output::increase();
		content = Read::text(siteFiles[i][1]);
		if (format == SF_UNKNOWN)
			format = StructureIO::getFormat(content);
		siteISO[i][1] = StructureIO::read(content);
		Output::decrease();
	}
	
	// Output
	Output::decrease();
	
	// Output
	Output::newline();
	Output::print("Creating structures to describe jumps and writing to input file");
	Output::increase();
	
	// Loop over lines in file to get jumps
	int j, k, m;
	int startIndex;
	int finalIndex;
	Matrix3D rotStart;
	Matrix3D rotFinal;
	Vector3D transStart;
	Vector3D transFinal;
	ISO start;
	ISO final;
	Word curDir;
	message.clear();
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Check if a jump line
		if (input[i][0].equal("jump", false))
		{
			
			// Line is too short
			if (input[i].length() < 28)
			{
				Output::newline(ERROR);
				Output::print("Line ");
				Output::print(i+1);
				Output::print(" of KMC file does not contain enough information");
				Output::quit();
			}
			
			// Save jump index
			index = atoi(input[i][1].array()) - 1;
			
			// Save site indexes
			startIndex = atoi(input[i][2].array()) - 1;
			finalIndex = atoi(input[i][15].array()) - 1;
			
			// Skip if site was not saved
			if ((startIndex >= siteISO.length()) || (finalIndex >= siteISO.length()))
				continue;
			if ((siteISO[startIndex].length() == 0) || (siteISO[finalIndex].length() == 0))
				continue;
			
			// Output
			Output::newline();
			Output::print("Setting information for jump ");
			Output::print(index + 1);
			Output::increase();
			
			// Build conversion operations
			for (j = 0, m = 3; j < 3; ++j)
			{
				for (k = 0; k < 3; ++k, ++m)
					rotStart(j, k) = Language::fractionToNumber(input[i][m]);
			}
			for (j = 0, k = 12; j < 3; ++j, ++k)
				transStart[j] = Language::fractionToNumber(input[i][k]);
			for (j = 0, m = 16; j < 3; ++j)
			{
				for (k = 0; k < 3; ++k, ++m)
					rotFinal(j, k) = Language::fractionToNumber(input[i][m]);
			}
			for (j = 0, k = 25; j < 3; ++j, ++k)
				transFinal[j] = Language::fractionToNumber(input[i][k]);
			
			// Create structures
			makeStructures(start, final, siteISO[startIndex][0], siteISO[startIndex][1], siteISO[finalIndex][0], \
				siteISO[finalIndex][1], rotStart, transStart, rotFinal, transFinal);
			
			// Create directory to store results
			curDir = "jump_";
			curDir += Language::numberToWord(index + 1);
			Directory::create(curDir, true);
			
			// Write structures
			StructureIO::write(Directory::makePath(curDir, "initial"), start, format, FRACTIONAL);
			StructureIO::write(Directory::makePath(curDir, "final"), final, format, FRACTIONAL);
			
			// Save line
			message.addLine();
			message.add("JUMP");
			message.add(index + 1);
			message.add(" ");
			message.add(startIndex + 1);
			message.add(finalIndex + 1);
			message.add("[[Total energy at transition state]]");
			message.add("[[Name of file with transition state phonon frequencies - or remove]]");
			
			// Output
			Output::decrease();
		}
	}
	
	// Print to file
	Output::setStream(kmcStream);
	Output::newline();
	Output::newline();
	Output::print(message, RIGHT);
	Output::setStream(origStream);
	
	// Print all of the network lines
	Output::setStream(kmcStream);
	Output::newline();
	message.clear();
	for (i = 0; i < input.length(); ++i)
	{
		if (input[i].length() == 0)
			continue;
		if (((input[i][0][0] >= '0') && (input[i][0][0] <= '9')) || (input[i][0][0] == '-'))
		{
			if (input[i].length() == 5)
			{
				message.addLine();
				for (j = 0; j < 5; ++j)
					message.add(input[i][j]);
			}
			else
			{
				if (message.numLines() != 0)
				{
					Output::newline();
					Output::print(message, RIGHT);
					message.clear();
				}
				Output::newline();
				for (j = 0; j < input[i].length(); ++j)
				{
					Output::print(input[i][j]);
					if (j != input[i].length() - 1)
						Output::print(" ");
				}
			}
		}
		else if (input[i][0].equal("com", false, 3))
		{
			Output::newline();
			Output::print(input[i][0]);
		}
	}
	if (message.numLines() != 0)
	{
		Output::newline();
		Output::print(message, RIGHT);
		message.clear();
	}
	Output::setStream(origStream);
	
	// Close kmc file
	Output::removeStream(kmcStream);
	
	// Output
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* void KMC::getTransformations(const Text& input)
 *
 * Get information about the conversion from ideal to unit cell
 */

void KMC::getTransformations(const Text& input)
{
	int i, j, k;
	bool unitSet = false;
	for (i = 0; i < input.length(); ++i)
	{
		if (input[i].length() == 9)
		{
			if (((input[i][0][0] >= '0') && (input[i][0][0] <= '9')) || (input[i][0][0] == '-'))
			{
				unitSet = true;
				for (j = 0; j < 3; ++j)
				{
					for (k = 0; k < 3; ++k)
						_idealToUnit(j, k) = Language::fractionToNumber(input[i][3*j+k]);
				}
				break;
			}
		}
	}
	if (!unitSet)
	{
		Output::newline(ERROR);
		Output::print("Did not find conversion matrix to unit cell in KMC input file");
		Output::quit();
	}
}



/* void KMC::makeStructures(ISO& start, ISO& final, ISO startInitial, const ISO& startRelaxed, ISO finalInitial, 
 *		const ISO& finalRelaxed, const Matrix3D& rotStart, const Vector3D& transStart, const Matrix3D& rotFinal,
 *		const Vector3D& transFinal)
 *
 * Create structures for starting and ending points of jumps
 */

void KMC::makeStructures(ISO& start, ISO& final, ISO startInitial, const ISO& startRelaxed, ISO finalInitial, \
	const ISO& finalRelaxed, const Matrix3D& rotStart, const Vector3D& transStart, const Matrix3D& rotFinal, \
	const Vector3D& transFinal)
{
	
	// Make copies of transformed relaxed structures
	start = startRelaxed;
	final = finalRelaxed;
	
	// Apply symmetry operations to structures
	transformStructure(startInitial, rotStart, transStart);
	transformStructure(start, rotStart, transStart);
	transformStructure(finalInitial, rotFinal, transFinal);
	transformStructure(final, rotFinal, transFinal);
	
	// Loop over atoms in the unrelaxed structures to line up atom orders
	int i, j, k;
	int numFound = 0;
	double tol = 1e-2;
	for (i = 0; i < startInitial.atoms().length(); ++i)
	{
		List<bool> matched(startInitial.atoms()[i].length(), false);
		List<int> positions(finalInitial.atoms()[i].length(), -1);
		for (j = 0; j < startInitial.atoms()[i].length(); ++j)
		{
			
			// Loop over atoms in the final structure
			for (k = 0; k < finalInitial.atoms()[i].length(); ++k)
			{
				
				// Atoms are at the same positions
				if (start.basis().distance(startInitial.atoms()[i][j].fractional(), FRACTIONAL, \
					finalInitial.atoms()[i][k].fractional(), FRACTIONAL) < tol)
				{
					matched[j] = true;
					positions[k] = j;
					break;
				}
			}
			
			// Did not find overlapping atoms
			if (k >= finalInitial.atoms()[i].length())
				++numFound;
		}
		
		// Figure out if any atoms were not matched
		for (j = 0, k = 0; j < positions.length(); ++j)
		{
			if (positions[j] == -1)
			{
				for (; k < matched.length(); ++k)
				{
					if (!matched[k])
					{
						positions[j] = k;
						break;
					}
				}
			}
		}
		
		// Save current atom numbers
		List<int> atomNums(finalInitial.atoms()[i].length(), 0);
		for (j = 0; j < atomNums.length(); ++j)
			atomNums[j] = finalInitial.atoms()[i][j].atomNumber();
		
		// Order atom numbers
		for (j = 0; j < atomNums.length(); ++j)
		{
			for (k = 0; k < atomNums.length(); ++k)
			{
				if (positions[k] == j)
				{
					final.setAtomInElement(atomNums[k], j);
					break;
				}
			}
		}
	}
	final.orderAtomNumbers();
	
	// Print a warning if all atoms overlapped
	if (numFound == 0)
	{
		Output::newline(WARNING);
		Output::print("It appears that no atoms are moving between the initial and final structures");
	}
	
	// Print a warning if multiple atoms were found to not overlap
	if (numFound > 1)
	{
		Output::newline(WARNING);
		Output::print("Found ");
		Output::print(numFound);
		Output::print(" atoms that appear to be moving when expecting exactly one");
		Output::newline(WARNING);
		Output::print("    The results of the simulation cannot be trusted");
	}
}



/* void KMC::transformStructure(ISO& iso, const Matrix3D& rotation, const Vector3D& translation, const Matrix3D& P,
 *		const Matrix3D& Q)
 *
 * Apply symmetry operation to structure
 */

void KMC::transformStructure(ISO& iso, const Matrix3D& rotation, const Vector3D& translation)
{
	
	// Loop over atoms in the structure
	int i, j;
	Vector3D newPos;
	Vector3D shiftVec;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			newPos = rotation * iso.atoms()[i][j].fractional();
			newPos += translation;
			ISO::moveIntoCell(newPos);
			iso.atoms()[i][j].fractional(newPos);
		}
	}
}



/* void KMC::setSimulation(const Text& input)
 *
 * Save properties for KMC simulation
 */

void KMC::setSimulation(const Text& input)
{
	
	// Clear information
	clear();
	
	// Output
	Output::newline();
	Output::print("Saving information for KMC simulation");
	Output::increase();
	
	// Loop over lines in file to get number of unique sites and jumps
	int i;
	int maxSite = 0;
	int maxJump = 0;
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found complement line
		if (input[i][0].equal("com", false, 3))
			_runComplement = true;
		
		// Found a site line
		else if (input[i][0].equal("site", false))
		{
			if (atoi(input[i][1].array()) > maxSite)
				maxSite = atoi(input[i][1].array());
		}
		
		// Found a jump line
		else if (input[i][0].equal("jump", false))
		{
			if (atoi(input[i][1].array()) > maxJump)
				maxJump = atoi(input[i][1].array());
			if (atoi(input[i][2].array()) > maxSite)
				maxSite = atoi(input[i][2].array());
			if (atoi(input[i][3].array()) > maxSite)
				maxSite = atoi(input[i][3].array());
		}
		
		// Found a network line
		else if ((input[i][0][0] >= '0') && (input[i][0][0] <= '9'))
		{
			if (input[i].length() == 4)
			{
				if (atoi(input[i][0].array()) > maxSite)
					maxSite = atoi(input[i][0].array());
			}
			else if (input[i].length() == 5)
			{
				if (atoi(input[i][0].array()) > maxJump)
					maxJump = atoi(input[i][0].array());
			}
		}
	}
	
	// Get information about the unique sites
	int curSite;
	List<bool> siteSet (maxSite, false);
	List<double> siteFormationEnergy (maxSite);
	List<double> siteTotalEnergy (maxSite);
	OList<Word> sitePhononFile (maxSite);
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found a site line
		if (input[i][0].equal("site", false))
		{
			
			// Not enough information
			if (input[i].length() < 4)
			{
				Output::newline(ERROR);
				Output::print("Not enough information on line ");
				Output::print(i+1);
				Output::print(" of kmc input file");
				Output::quit();
			}
			
			// Save energies
			curSite = atoi(input[i][1].array()) - 1;
			siteSet[curSite] = true;
			siteFormationEnergy[curSite] = atof(input[i][2].array());
			siteTotalEnergy[curSite] = atof(input[i][3].array());
			
			// Save phonon file if available
			if (input[i].length() > 4)
			{
				if (!Language::isComment(input[i][4]))
					sitePhononFile[curSite] = input[i][4];
			}
		}
	}
	
	// Get information about unique jumps
	int curJump;
	List<bool> jumpSet (maxJump, false);
	List<int> jumpStartSite (maxJump);
	List<int> jumpEndSite (maxJump);
	List<double> jumpTotalEnergy (maxJump);
	OList<Word> jumpPhononFile (maxJump);
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found a jump line
		if (input[i][0].equal("jump", false))
		{
			
			// Not enough information
			if (input[i].length() < 5)
			{
				Output::newline(ERROR);
				Output::print("Not enough information on line ");
				Output::print(i+1);
				Output::print(" of kmc input file");
				Output::quit();
			}
			
			// Save information about jump
			curJump = atoi(input[i][1].array()) - 1;
			jumpSet[curJump] = true;
			jumpStartSite[curJump] = atoi(input[i][2].array()) - 1;
			jumpEndSite[curJump] = atoi(input[i][3].array()) - 1;
			jumpTotalEnergy[curJump] = atof(input[i][4].array());
			
			// Save phonon file if available
			if (input[i].length() > 5)
			{
				if (!Language::isComment(input[i][5]))
					jumpPhononFile[curJump] = input[i][5];
			}
		}
	}
	
	// Save attempt values
	double curBarrier;
	_attempts.length(jumpSet.length());
	for (i = 0; i < jumpSet.length(); ++i)
	{
		
		// Skip if jump was not set
		if (!jumpSet[i])
			continue;
		
		// Skip if sites were not set
		if ((!siteSet[jumpStartSite[i]]) || (!siteSet[jumpEndSite[i]]))
		{
			jumpSet[i] = false;
			continue;
		}
		
		// Save forward jump information
		_attempts[i].length(2);
		curBarrier = jumpTotalEnergy[i] - siteTotalEnergy[jumpStartSite[i]];
		if ((sitePhononFile[jumpStartSite[i]].length() > 0) && (jumpPhononFile[i].length() > 0))
			_attempts[i][0].set(curBarrier, sitePhononFile[jumpStartSite[i]], jumpPhononFile[i]);
		else
			_attempts[i][0].set(curBarrier, 1e13);
		
		// Save reverse jump information
		curBarrier = jumpTotalEnergy[i] - siteTotalEnergy[jumpEndSite[i]];
		if ((sitePhononFile[jumpEndSite[i]].length() > 0) && (jumpPhononFile[i].length() > 0))
			_attempts[i][1].set(curBarrier, sitePhononFile[jumpEndSite[i]], jumpPhononFile[i]);
		else
			_attempts[i][1].set(curBarrier, 1e13);
	}
	
	// Loop over file and save information about network nodes
	int j, k;
	StorageNode* curNode = 0;
	Vector3D curVector;
	for (i = 0; i < input.length(); ++i)
	{
		
		// Skip if line is empty or commented
		if (input[i].length() == 0)
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Found a possible network line
		else if ((input[i][0][0] >= '0') && (input[i][0][0] <= '9'))
		{
			
			// Found a starting node line
			if (input[i].length() == 4)
			{
				curSite = atoi(input[i][0].array()) - 1;
				if (siteSet[curSite])
				{
					for (j = 0, k = 1; j < 3; ++j, ++k)
						curVector[j] = atof(input[i][k].array());
					_nodes.add();
					_nodes.last().position(curVector);
					_nodes.last().formationEnergy(siteFormationEnergy[curSite]);
					curNode = &_nodes.last();
				}
				else
					curNode = 0;
			}
			
			// Found a jump line
			else if ((input[i].length() == 5) && (curNode != 0))
			{
				curJump = atoi(input[i][0].array()) - 1;
				if (jumpSet[curJump])
				{
					for (j = 0, k = 2; j < 3; ++j, ++k)
						curVector[j] = atof(input[i][k].array());
					if (atoi(input[i][1].array()) > 0)
						curNode->add(curVector, &_attempts[curJump][0]);
					else
						curNode->add(curVector, &_attempts[curJump][1]);
				}
			}
		}
	}
	
	// Get information about the basis
	bool basisSet = false;
	for (i = 0; i < input.length(); ++i)
	{
		if (input[i].length() == 6)
		{
			if ((input[i][0][0] >= '0') && (input[i][0][0] <= '9'))
			{
				basisSet = true;
				Vector3D lengths;
				Vector3D angles;
				for (j = 0, k = 3; j < 3; ++j, ++k)
				{
					lengths[j] = atof(input[i][j].array());
					angles[j] = Num<double>::toRadians(atof(input[i][k].array()));
				}
				_basis.set(lengths, angles, false);
			}
		}
	}
	if (!basisSet)
	{
		Output::newline(ERROR);
		Output::print("Did not find a basis in the KMC input file");
		Output::quit();
	}
	
	// Get information about transformations between cells
	getTransformations(input);
	
	// Get the number of sites that were set
	int numSites = 0;
	for (i = 0; i < siteSet.length(); ++i)
		numSites += siteSet[i] == true ? 1 : 0;
	if (numSites == 0)
	{
		Output::newline(ERROR);
		Output::print("No site information was set in kmc input file");
		Output::quit();
	}
	Output::newline();
	Output::print("Information for ");
	Output::print(numSites);
	Output::print(" of ");
	Output::print(maxSite);
	Output::print(" site");
	if (maxSite != 1)
		Output::print("s");
	Output::print(" has been set");
	
	// Get the number of jumps that were set
	int numJumps = 0;
	for (i = 0; i < jumpSet.length(); ++i)
		numJumps += jumpSet[i] == true ? 1 : 0;
	if (numJumps == 0)
	{
		Output::newline(ERROR);
		Output::print("No jump information was set in kmc input file");
		Output::quit();
	}
	Output::newline();
	Output::print("Information for ");
	Output::print(numJumps);
	Output::print(" of ");
	Output::print(maxJump);
	Output::print(" jump");
	if (numJumps != 1)
		Output::print("s");
	Output::print(" has been set");
	
	// Output
	Output::decrease();
}



/* void KMC::run(Random& random)
 *
 * Run KMC simulation
 */

void KMC::run(Random& random)
{
	if (_runComplement)
		runComplement(random);
	else
		runPrimary(random);
}



/* void KMC::runPrimary(Random& random)
 *
 * Run KMC simulation tracking single defect
 */

void KMC::runPrimary(Random& random)
{
	
	// Make the list of temperatures to run at
	List<double> sizes;
	List<double> invTemps;
	const int numSteps = 7;
	for (double cur = 1.0/1000.0; cur <= 1.0/300.0 + 1e-5; cur += (1.0/300.0 - 1.0/1000.0)/numSteps)
	{
		sizes += 1;
		invTemps += cur;
	}
	
	// Run the simulation
	runSimulations(random, sizes, invTemps);
}



/* void KMC::runComplement(Random& random)
 *
 * Run KMC simulation tracking atoms other than the defect
 */

void KMC::runComplement(Random& random)
{
	
	// Based on how many sites are in the cell, figure out the minimum number of cells to use
	double startSize = Num<double>::ceil (pow(  500.0 / _nodes.length(), 1.0/3.0));
	double endSize   = Num<double>::floor(pow(25000.0 / _nodes.length(), 1.0/3.0));
	
	// Set the temperatures at which simulations will be performed
	List<double> sizes;
	List<double> invTemps;
	Functor<KMC> concFun(this, &KMC::concentrationPerCell);
	for (double curSize = startSize; curSize < endSize + 0.1; ++curSize)
	{
		sizes += curSize;
		invTemps.add();
		Solve<KMC>::findRoot(concFun, 1.0/_nodes.length()/pow(curSize, 3), 1e-10, 1e-10, 1e-4, 1.0, invTemps.last());
	}
	
	// Run the simulations
	runSimulations(random, sizes, invTemps);
}



/* void KMC::runSimulations(Random& random, const List<double>& sizes, const List<double>& invTemps)
 *
 * Run KMC simulations at set temperatures and cell sizes
 */

void KMC::runSimulations(Random& random, const List<double>& sizes, const List<double>& invTemps)
{
	
	// Output
	Output::newline();
	Output::print("Running ");
	Output::print(invTemps.length());
	Output::print(" KMC simulation");
	if (invTemps.length() == 1)
	{
		Output::print(" at ");
		Output::print(1.0/invTemps[0]);
		Output::print(" K (");
		Output::print(invTemps[0]);
		Output::print(" 1/K)");
	}
	else
	{
		Output::print("s");
		Output::print(" from ");
		Output::print(1.0/invTemps[0], 4);
		Output::print(" K (");
		Output::print(invTemps[0], 8);
		Output::print(" 1/K) to ");
		Output::print(1.0/invTemps.last(), 4);
		Output::print(" K (");
		Output::print(invTemps.last(), 8);
		Output::print(" 1/K)");
	}
	Output::increase();
	
	// Save conversion of vector from ideal cell to reference cell
	Matrix3D idealVecToUnit = _idealToUnit.inverse().transpose();
	
	// Save the reference basis
	Basis unitBasis (_idealToUnit * _basis.vectors(), false);
	
	// Loop over simulations to run
	int i, j;
	int start;
	Basis simCellBasis;
	OList<Vector3D> D (invTemps.length());
	for (i = 0; i < invTemps.length(); ++i)
	{
		
		// Output
		Output::newline();
		Output::print("Performing simulation at ");
		Output::print(1.0/invTemps[i], 4);
		Output::print(" K (");
		Output::print(invTemps[i], 8);
		Output::print(" 1/K)");
		Output::increase();
		
		// Create the simulation cell
		OList<Node> simNodes;
		createSimulationNodes(simNodes, sizes[i], simCellBasis);
		
		// Pick the start node
		start = pickStartNode(simNodes, invTemps[i], random);
		
		// Initialize tracers in the cell
		OList<Tracer> tracers;
		initializeTracers(simNodes, tracers, start);
		
		// Set temperature for simulation
		for (j = 0; j < _attempts.length(); ++j)
		{
			_attempts[j][0].setRate(1.0/invTemps[i]);
			_attempts[j][1].setRate(1.0/invTemps[i]);
		}
		for (j = 0; j < simNodes.length(); ++j)
			simNodes[j].setRates();
		
		// Run simulation at current temperature
		D[i] = runSingleSimulation(random, simNodes.length(), &simNodes[start], tracers);
		
		// Convert diffusivity into reference basis
		simCellBasis.toFractional(D[i]);
		D[i] *= idealVecToUnit;
		D[i] *= sizes[i];
		
		// Output
		Output::decrease();
	}
	
	// Output
	Output::decrease();
	
	// Analyze the results
	analyzeResults(unitBasis, invTemps, D);
}



/* Vector3D KMC::runSingleSimulation(Random& random, int numNodes, Node* curNode, const OList<Tracer>& tracers)
 *
 * Run simulation with current cell and settings
 */

Vector3D KMC::runSingleSimulation(Random& random, int numNodes, Node* curNode, const OList<Tracer>& tracers)
{
	
	// Output
	Output::newline();
	Output::print("Running simulation");
	Output::increase();
	
	// Variables to store information during the run
	double avgD;
	double avgLocalD;
	double avgDsquared;
	double totalTime;
	Vector3D D;
	Vector3D localD;
	Vector3D totalVecSquared;
	
	// Total variables
	double netAvgD = 0;
	double netAvgDsquared = 0;
	Vector3D netD (0.0);
	
	// Loop until converged
	int i;
	double ran;
	int numJumps;
	int numSimulations = 0;
	Tracer* curTracer;
	Vector3D prevVec;
	const int jumpsBeforeCheck = _jumpsPerAtom * numNodes;
	while (1)
	{
		
		// Reset variables for new simulation
		avgD = 0;
		avgLocalD = 0;
		avgDsquared = 0;
		totalTime = 0;
		D = 0.0;
		localD = 0.0;
		totalVecSquared = 0.0;
		numJumps = 0;
		prevVec = 0.0;
		
		// Loop until current simulation is converged
		while (1)
		{
			
			// Make a jump
			curTracer = 0;
			ran = random.decimal(0, curNode->rates().last());
			for (i = 0; i < curNode->rates().length(); ++i)
			{
				if (ran < curNode->rates()[i])
				{
				
					// Save the tracer motion
					if (curNode->tracer())
					{
						curTracer = curNode->tracer();
						prevVec = curTracer->vector();
						curTracer->add(curNode->vectors()[i]);
						curNode->endNodes()[i]->tracer() = curTracer;
						curNode->tracer() = 0;
					}
					else
					{
						curTracer = curNode->endNodes()[i]->tracer();
						prevVec = curTracer->vector();
						curTracer->subtract(curNode->vectors()[i]);
						curNode->tracer() = curTracer;
						curNode->endNodes()[i]->tracer() = 0;
					}
				
					// Update time
					while ((ran = random.decimal(0, 1)) == 0.0) {}
					totalTime += -log(ran) / curNode->rates().last();
				
					// Update node
					curNode = curNode->endNodes()[i];
					break;
				}
			}
		
			// Skip if something went wrong
			if (curTracer == 0)
				continue;
		
			// Update the cumulative total squared displacement
			for (i = 0; i < 3; ++i)
			{
				totalVecSquared[i] -= prevVec[i]*prevVec[i];
				totalVecSquared[i] += curTracer->vector()[i]*curTracer->vector()[i];
			}
		
			// Get the diffusivity tensor for current time step
			localD  = totalVecSquared;
			localD /= 2.0;
			localD /= totalTime;
			localD /= tracers.length();
			for (i = 0; i < 3; ++i)
				localD[i] = Num<double>::abs(localD[i]);
		
			// Update the total diffusivity
			for (i = 0; i < 3; ++i)
				D[i] = (numJumps * D[i] + localD[i]) / (numJumps + 1);
			avgD = (D[0] + D[1] + D[2]) / 3.0;
		
			// Update the squared diffusivity
			avgLocalD = (localD[0] + localD[1] + localD[2]) / 3.0;
			avgDsquared = (numJumps * avgDsquared + avgLocalD * avgLocalD) / (numJumps + 1);
		
			// Check for convergence if needed
			if (++numJumps >= jumpsBeforeCheck)
			{
				if (sqrt((avgDsquared - avgD*avgD) / numJumps) / avgD < _convergence)
					break;
			}
		}
		
		// Update variables
		for (i = 0; i < 3; ++i)
			netD[i] = (numSimulations * netD[i] + D[i]) / (numSimulations + 1);
		netAvgD = (netD[0] + netD[1] + netD[2]) / 3.0;
		avgLocalD = (D[0] + D[1] + D[2]) / 3.0;
		netAvgDsquared = (numSimulations * netAvgDsquared + avgLocalD * avgLocalD) / (numSimulations + 1);
		
		// Check for convergence if needed
		if (++numSimulations >= 1000)
		{
			if (sqrt((netAvgDsquared - netAvgD*netAvgD) / numSimulations) / netAvgD < _convergence)
				break;
		}
	}
	
	// Output
	Output::decrease();
	
	// Convert from A^2 to cm^2
	for (i = 0; i < 3; ++i)
		netD[i] /= 1e16;
	netAvgD = (netD[0] + netD[1] + netD[2]) / 3.0;
	
	// Print the result
	Output::newline();
	Output::print("Average diffusivity: ");
	Output::printSci(netAvgD, 4);
	Output::print(" cm^2/s");
	
	// Return the diffusivity
	return netD;
}



/* void KMC::createSimulationNodes(OList<Node>& simNodes, double cellSize, Basis& simCellBasis)
 *
 * Create list of nodes in the KMC simulation
 */

void KMC::createSimulationNodes(OList<Node>& simNodes, double cellSize, Basis& simCellBasis)
{
	
	// Output
	Output::newline();
	Output::print("Creating simulation cell");
	Output::increase();
	
	// Save the new basis
	Matrix3D conv(cellSize, 0, 0, 0, cellSize, 0, 0, 0, cellSize);
	simCellBasis.set(conv * _basis.vectors(), false);
	
	// Set the number of nodes that will be in the simulation
	int numCells = (int)Num<double>::round(cellSize*cellSize*cellSize, 1.0);
	simNodes.length(numCells * _nodes.length());
	
	// Variable to store the nodes grouped by type
	int i;
	List<Node*>::D2 nodesByType (_nodes.length());
	for (i = 0; i < _nodes.length(); ++i)
		nodesByType[i].length(numCells);
	
	// Variable to store the end node types
	List<int>::D3 endTypes;
	endTypes.length(_nodes.length());
	for (i = 0; i < _nodes.length(); ++i)
		endTypes[i].length(numCells);
	
	// Variable to store the distance of a site from the origin
	Vector3D origin(0.0);
	List<double>::D2 distances;
	distances.length(_nodes.length());
	for (i = 0; i < _nodes.length(); ++i)
		distances[i].length(numCells);
	
	// Loop over cells to add
	int j, k;
	bool found;
	int curNode = 0;
	int curImage;
	Vector3D curVec;
	Vector3D curPos;
	Vector3D endPos;
	List<int> curEndTypes;
	const double tol = 1e-2;
	for (i = 0; i < _nodes.length(); ++i)
	{
		
		// Get the end node types
		curEndTypes.length(_nodes[i].vectors().length());
		for (j = 0; j < _nodes[i].vectors().length(); ++j)
		{
			found = false;
			endPos = _nodes[i].position();
			endPos += _nodes[i].vectors()[j];
			for (k = 0; k < _nodes.length(); ++k)
			{
				if (simCellBasis.distance(endPos, FRACTIONAL, _nodes[k].position(), FRACTIONAL) < tol)
				{
					found = true;
					curEndTypes[j] = k;
					break;
				}
			}
			if (!found)
			{
				Output::newline(ERROR);
				Output::print("Internal error: Failed to match kmc jump to node during setup");
				Output::quit();
			}
		}
		
		// Loop over images
		curImage = 0;
		for (curVec[0] = 1; curVec[0] < cellSize+0.1; ++curVec[0])
		{
			for (curVec[1] = 1; curVec[1] < cellSize+0.1; ++curVec[1])
			{
				for (curVec[2] = 1; curVec[2] < cellSize+0.1; ++curVec[2])
				{
					
					// Save current position
					for (j = 0; j < 3; ++j)
						curPos[j] = (_nodes[i].position()[j] + curVec[j] - 1) / cellSize;
					distances[i][curImage] = simCellBasis.distance(origin, FRACTIONAL, curPos, FRACTIONAL);
					
					// Save information about the node
					simNodes[curNode].position(curPos);
					simNodes[curNode].formationEnergy(_nodes[i].formationEnergy());
					for (j = 0; j < _nodes[i].vectors().length(); ++j)
					{
						simNodes[curNode].addVector(_nodes[i].vectors()[j] / cellSize);
						simNodes[curNode].addAttempt(_nodes[i].attempts()[j]);
					}
					nodesByType[i][curImage] = &simNodes[curNode];
					endTypes[i][curImage] = curEndTypes;
					
					// Increment counters
					++curNode;
					++curImage;
				}
			}
		}
	}
	
	// Set the end nodes for all jumps
	int m;
	double curDis;
	for (i = 0; i < nodesByType.length(); ++i)
	{
		for (j = 0; j < nodesByType[i].length(); ++j)
		{
			for (k = 0; k < nodesByType[i][j]->vectors().length(); ++k)
			{
				
				// Set the current end position
				endPos = nodesByType[i][j]->position();
				endPos += nodesByType[i][j]->vectors()[k];
				ISO::moveIntoCell(endPos);
				curDis = simCellBasis.distance(origin, FRACTIONAL, endPos, FRACTIONAL);
				
				// Find the end node
				found = false;
				for (m = 0; m < nodesByType[endTypes[i][j][k]].length(); ++m)
				{
					if (Num<double>::eq(curDis, distances[endTypes[i][j][k]][m], tol))
					{
						if (simCellBasis.distance(endPos, FRACTIONAL, nodesByType[endTypes[i][j][k]][m]->position(), \
							FRACTIONAL) < tol)
						{
							found = true;
							nodesByType[i][j]->addEndNode(nodesByType[endTypes[i][j][k]][m]);
							break;
						}
					}
				}
				if (!found)
				{
					Output::newline(ERROR);
					Output::print("Internal error: Failed to match kmc jump to node");
					Output::quit();
				}
			}
		}
	}
	
	// Convert everything to cartesian coordinates
	for (i = 0; i < simNodes.length(); ++i)
	{
		simNodes[i].position(simCellBasis.getCartesian(simNodes[i].position()));
		for (j = 0; j < simNodes[i].vectors().length(); ++j)
			simNodes[i].vectors()[j] = simCellBasis.getCartesian(simNodes[i].vectors()[j]);
	}
	
	// Output
	Output::decrease();
}



/* void KMC::initializeTracers(OList<Node>& simNodes, OList<Tracer>& tracers, int start)
 *
 * Set tracers at start of the simulation
 */

void KMC::initializeTracers(OList<Node>& simNodes, OList<Tracer>& tracers, int start)
{
	
	// Running in complement mode
	if (_runComplement)
	{
		tracers.length(simNodes.length() - 1);
		for (int i = 0, j = 0; i < simNodes.length(); ++i)
		{
			if (i == start)
				simNodes[i].tracer() = 0;
			else
				simNodes[i].tracer() = &tracers[j++];
		}
	}
	
	// Running in primary mode
	else
	{
		tracers.length(1);
		for (int i = 0; i < simNodes.length(); ++i)
			simNodes[i].tracer() = 0;
		simNodes[start].tracer() = &tracers[0];
	}
}



/* void KMC::analyzeResults(Basis& unitBasis, const List<double>& invTemps, const OList<Vector3D>& Dfrac)
 *
 * Analayze results from kmc simulation
 */

void KMC::analyzeResults(Basis& unitBasis, const List<double>& invTemps, const OList<Vector3D>& Dfrac)
{
	
	// Output
	Output::newline();
	Output::print("Analyzing results");
	Output::increase();
	
	// Convert diffusivities to cartesian coordinates
	int i;
	OList<Vector3D> Dcart (Dfrac.length());
	for (i = 0; i < Dfrac.length(); ++i)
		Dcart[i] = unitBasis.getCartesian(Dfrac[i]);
	
	// Initialize fitting list
	List<double>::D2 fitList;
	fitList.length(invTemps.length());
	for (i = 0; i < invTemps.length(); ++i)
	{
		fitList[i].length(2);
		fitList[i][0] = invTemps[i];
	}
	
	// Save results of diffusivities along each direction
	int j;
	Vector polyFit;
	for (i = 0; i < 3; ++i)
	{
		
		// Fit linear function to data
		for (j = 0; j < fitList.length(); ++j)
			fitList[j][1] = log(Num<double>::abs(Dfrac[j][i] * unitBasis.lengths()[i]));
		polyFit = Fit::polynomial(fitList, 0, 2);
		
		// Save results
		_activationEnergyVector[i] = -Constants::kb * polyFit[1];
		_D0vector[i] = exp(polyFit[0]);
	}
	
	// Save results of average diffusivity
	for (i = 0; i < fitList.length(); ++i)
		fitList[i][1] = log((Dcart[i][0] + Dcart[i][1] + Dcart[i][2]) / 3.0);
	polyFit = Fit::polynomial(fitList, 0, 2);
	_activationEnergy = -Constants::kb * polyFit[1];
	_D0 = exp(polyFit[0]);
	
	// Print results to output
	Output::newline();
	Output::print("Activation energy: ");
	Output::print(_activationEnergy, 4);
	Output::print(" eV (");
	Output::print(_activationEnergyVector, 4, true);
	Output::print(")");
	Output::newline();
	Output::print("D0: ");
	Output::printSci(_D0, 4);
	Output::print(" cm^2/s (");
	for (i = 0; i < 3; ++i)
	{
		Output::printSci(_D0vector[i], 4);
		if (i != 2)
			Output::print(", ");
	}
	Output::print(")");
	
	// Print results to file
	print("kmc.out");
	
	// Print raw results
	int newStreamID;
	int origStreamID = Output::streamID();
	Output message;
	for (i = 0; i < 4; ++i)
	{
		
		// Add data
		message.clear();
		message.addLines(invTemps.length());
		for (j = 0; j < invTemps.length(); ++j)
		{
			message.addLine();
			message.addSci(invTemps[j], 6);
			if (i < 3)
				message.addSci(Num<double>::abs(Dfrac[j][i] * unitBasis.lengths()[i]), 6);
			else
				message.addSci((Dcart[j][0] + Dcart[j][1] + Dcart[j][2]) / 3.0, 6);
		}
		
		// Set the stream to file
		     if (i == 0) newStreamID = Output::addStream("Da.out");
		else if (i == 1) newStreamID = Output::addStream("Db.out");
		else if (i == 2) newStreamID = Output::addStream("Dc.out");
		else             newStreamID = Output::addStream("Davg.out");
		Output::setStream(newStreamID);
		
		// Print data
		Output::newline();
		Output::print(message, RIGHT);
		
		// Reset stream
		Output::setStream(origStreamID);
		Output::removeStream(newStreamID);
	}
	
	// Output
	Output::decrease();
}



/* void KMC::print(const Word& file) const
 *
 * Print results of KMC simulation
 */

void KMC::print(const Word& file) const
{
	
	// Setup output if file was set
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file.length() > 0)
	{
		
		// Open file for writing if needed
		if (file != "stdout")
			Output::setStream(Output::addStream(file));

		// Set output method
		Output::method(STANDARD);
	}
	
	// Print results
	Output::newline();
	Output::print("Activation energy");
	Output::newline();
	Output::print("    Avg: ");
	Output::print(_activationEnergy, 4);
	Output::print(" eV");
	Output::newline();
	Output::print("      a: ");
	Output::print(_activationEnergyVector[0], 4);
	Output::print(" eV");
	Output::newline();
	Output::print("      b: ");
	Output::print(_activationEnergyVector[1], 4);
	Output::print(" eV");
	Output::newline();
	Output::print("      c: ");
	Output::print(_activationEnergyVector[2], 4);
	Output::print(" eV");
	Output::newline();
	Output::newline();
	Output::print("Diffusivity prefactor (D0)");
	Output::newline();
	Output::print("    Avg: ");
	Output::printSci(_D0, 4);
	Output::print(" cm^2/s");
	Output::newline();
	Output::print("      a: ");
	Output::printSci(_D0vector[0], 4);
	Output::print(" cm^2/s");
	Output::newline();
	Output::print("      b: ");
	Output::printSci(_D0vector[1], 4);
	Output::print(" cm^2/s");
	Output::newline();
	Output::print("      c: ");
	Output::printSci(_D0vector[2], 4);
	Output::print(" cm^2/s");
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}



/* void KMC::Attempt::set(double barrier, const Word& gsModesFile, const Word& tsModesFile)
 *
 * Set barrier and phonon frequencies along kmc transition
 */

void KMC::Attempt::set(double barrier, const Word& gsModesFile, const Word& tsModesFile)
{
	
	// Save barrier
	set(barrier);
	
	// Read frequencies from phonon files
	Text gsContent = Read::text(gsModesFile);
	Text tsContent = Read::text(tsModesFile);
	
	// Save ground state modes
	int i, j;
	Linked<double> gsModes;
	for (i = 0; i < gsContent.length(); ++i)
	{
		for (j = 0; j < gsContent[i].length(); ++j)
		{
			if (!Language::isNumber(gsContent[i][j]))
			{
				Output::newline(ERROR);
				Output::print("Found a non-numerical value (");
				Output::print(gsContent[i][j]);
				Output::print(") in phonon frequency file ");
				Output::print(gsModesFile);
				Output::quit();
			}
			gsModes += atof(gsContent[i][j].array());
		}
	}
	
	// Save transition state modes
	Linked<double> tsModes;
	for (i = 0; i < tsContent.length(); ++i)
	{
		for (j = 0; j < tsContent[i].length(); ++j)
		{
			if (!Language::isNumber(tsContent[i][j]))
			{
				Output::newline(ERROR);
				Output::print("Found a non-numerical value (");
				Output::print(tsContent[i][j]);
				Output::print(") in phonon frequency file ");
				Output::print(tsModesFile);
				Output::quit();
			}
			tsModes += atof(tsContent[i][j].array());
		}
	}
	
	// Modes are not the same length
	if (gsModes.length() != tsModes.length())
	{
		Output::newline(ERROR);
		Output::print("Found ");
		Output::print(gsModes.length());
		Output::print(" mode");
		if (gsModes.length() != 1)
			Output::print("s");
		Output::print(" in ");
		Output::print(gsModesFile);
		Output::print(" but ");
		Output::print(tsModes.length());
		Output::print(" mode");
		if (tsModes.length() != 1)
			Output::print("s");
		Output::print(" in ");
		Output::print(tsModesFile);
		Output::quit();
	}
	
	// Save modes in units of 1/s
	_gsModes.length(gsModes.length());
	_tsModes.length(tsModes.length());
	Linked<double>::iterator itGS = gsModes.begin();
	Linked<double>::iterator itTS = tsModes.begin();
	for (i = 0; i < gsModes.length(); ++i, ++itGS, ++itTS)
	{
		_gsModes[i] = Num<double>::frequencyToHz(*itGS);
		_tsModes[i] = Num<double>::frequencyToHz(*itTS);
	}
	
	// Make sure that lists are sorted
	_gsModes.sort();
	_tsModes.sort();
}



/* double KMC::Attempt::setRate(double temperature)
 *
 * Set the transition rate at a temperature
 */

double KMC::Attempt::setRate(double temperature)
{
	
	// Set the frequency if needed
	if (_gsModes.length() > 0)
	{
		
		// Constant scaling factor
		double scale = 2 * Constants::kb * temperature / Constants::h;
	
		// Loop over modes and calculate frequency
		_frequency = scale;
		for (int i = 4; i < _gsModes.length(); ++i)
			_frequency *= sinh(_gsModes[i] / scale) / sinh(_tsModes[i] / scale);
		_frequency *= sinh(_gsModes[3] / scale);
	}
	
	// Set rate
	_rate = _frequency * exp(-_barrier / (Constants::kb * temperature));
	return _rate;
}
