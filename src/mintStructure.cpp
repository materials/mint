/* mintStructure.cpp -- Structure for mint files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "num.h"
#include "mintStructure.h"
#include "elements.h"
#include "symmetry.h"
#include "spaceGroup.h"
#include "language.h"
#include "output.h"
#include <cstdlib>



/* ISO MintStructure::read(const Text& content, double tol, double clusterTol)
 * 
 * Read structure from file
 */

ISO MintStructure::read(const Text& content, double tol, double clusterTol)
{
	
	// Output
	Output::newline();
	Output::print("Converting mint structure file to an internal structure object");
	Output::increase();
	
	// Variable to store result
	ISO iso;
	
	// Look for basis
	int i, j, k;
	bool sysFound = false;
	LatticeSystem system;
	List<bool> basisFixed;
	List<double> numbers;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if empty
		if (!content[i].length())
			continue;
		
		// Found basis
		if (getLabel(content[i][0]) == BASIS)
		{
			
			// Loop over lines
			for (j = i; j < content.length(); ++j)
			{
				
				// Skip if line is empty
				if (!content[j].length())
					continue;
				
				// Break if a label is found
				if (j != i)
				{
					if (getLabel(content[j][0]) != UNKNOWN)
						break;
				}
				
				// Loop over words on line
				for (k = (j == i) ? 1 : 0; k < content[j].length(); ++k)
				{

					// Found a comment
					if (Language::isComment(content[j][k]))
						break;

					// Found a system
					system = ISO::system(content[j][k]);
					if (system != LS_UNKNOWN)
					{
						sysFound = true;
						iso.basis().latticeSystem(system);
						break;
					}

					// Found a number
					else if (Language::isNumber(content[j][k]))
						numbers += atof(content[j][k].array());

					// Found t/f
					else if (content[j][k].equal("t", false, 1))
						basisFixed += true;
					else if (content[j][k].equal("f", false, 1))
						basisFixed += false;

					// Anything else
					else
					{
						Output::newline(ERROR);
						Output::print("Did not recognize basis setting (");
						Output::print(content[j][k]);
						Output::print(")");
						Output::quit();
					}
				}
				
				// Break if system was set
				if (sysFound)
					break;
			}
		}
	}
	
	// Found at least 1 number
	bool basisSet = false;
	if (numbers.length() > 0)
	{
		
		// System was also found
		if (sysFound)
		{
			Output::newline(ERROR);
			Output::print("Found numbers and lattice system describing basis in mint structure file");
			Output::quit();
		}
		
		// Found 6 numbers
		if (numbers.length() == 6)
		{
			basisSet = true;
			iso.basis(Vector3D(numbers[0], numbers[1], numbers[2]), \
				Vector3D(Num<double>::toRadians(numbers[3]), Num<double>::toRadians(numbers[4]), \
				Num<double>::toRadians(numbers[5])), true);
		}
		
		// Found 9 numbers
		else if (numbers.length() == 9)
		{
			basisSet = true;
			iso.basis(Matrix3D(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], \
				numbers[6], numbers[7], numbers[8]), true);
		}
		
		// Found anything else
		else
		{
			Output::newline(ERROR);
			Output::print("Found ");
			Output::print(numbers.length());
			Output::print(" number");
			if (numbers.length() != 1)
				Output::print("s");
			Output::print(" when setting basis but expecting 0, 6, or 9");
			Output::quit();
		}
	}
	
	// Values of fixed basis have been set
	if (basisFixed.length() > 0)
	{
		
		// Basis was not set
		if (!basisSet)
		{
			Output::newline(ERROR);
			Output::print("Cannot fixed basis in mint structure file when basis is not provided");
			Output::quit();
		}
		
		// Set which values are fixed
		int max = (basisFixed.length() < 3) ? basisFixed.length() : 3;
		for (i = 0; i < max; ++i)
			iso.basis().lengthFixed(i, basisFixed[i]);
		max = (basisFixed.length() < 6) ? basisFixed.length() : 6;
		for (i = 3; i < max; ++i)
			iso.basis().angleFixed(i-3, basisFixed[i]);
	}
	
	// Check if anything is fixed
	bool anythingFixed = false;
	for (i = 0; i < basisFixed.length(); ++i)
	{
		if (basisFixed[i])
		{
			anythingFixed = true;
			break;
		}
	}
	
	// Output
	Output::newline();
	Output::print("Fixed basis metrics: ");
	if (anythingFixed)
	{
		if (iso.basis().lengthFixed()[0])
			Output::print("a ");
		if (iso.basis().lengthFixed()[1])
			Output::print("b ");
		if (iso.basis().lengthFixed()[2])
			Output::print("c ");
		if (iso.basis().angleFixed()[0])
			Output::print("alpha ");
		if (iso.basis().angleFixed()[1])
			Output::print("beta ");
		if (iso.basis().angleFixed()[2])
			Output::print("gamma ");
	}
	else
		Output::print("none");
	
	// Look for space group line
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if empty
		if (!content[i].length())
			continue;
		
		// Found space group
		if (getLabel(content[i][0]) == SPACEGROUP)
		{
			
			// More data on line
			Word spaceGroup;
			if (content[i].length() > 1)
			{
				
				// Skip next word if group
				j = 1;
				if (content[i][1].equal("group", false, 1))
					++j;
				
				// If next word is not a comment then save space group
				if (!Language::isComment(content[i][j]))
				{
					for (; j < content[i].length(); ++j)
					{
						if (Language::isComment(content[i][j]))
							break;
						spaceGroup += content[i][j];
					}
					iso.spaceGroup(spaceGroup);
					break;
				}
			}
			
			// Go to next line
			if (++i < content.length())
			{
				
				// Save space group
				for (j = 0; j < content[i].length(); ++j)
				{
					if (Language::isComment(content[i][j]))
						break;
					spaceGroup += content[i][j];
				}
				iso.spaceGroup(spaceGroup);
				break;
			}
		}
	}
	
	// System and space group were set
	if ((iso.spaceGroup().length()) && (sysFound))
	{
		Output::newline(ERROR);
		Output::print("Space group and lattice system were both set in mint structure file");
		Output::quit();
	}
	
	// Print space group if set
	if (iso.spaceGroup().length())
	{
		Output::newline();
		Output::print("Setting structure in space group ");
		Output::print(iso.spaceGroup());
	}
	
	// Print lattice system if set
	else if (sysFound)
	{
		Output::newline();
		Output::print("Setting basis in ");
		Output::print(ISO::system(iso.basis().latticeSystem()).tolower());
		Output::print(" lattice system");
	}
	
	// Output
	Output::newline();
	Output::print("Adding atoms");
	Output::increase();
	
	// Temporary variable to store atom data
	bool setElement;
	bool setFractional;
	bool setCartesian;
	bool setFixed;
	bool setMult;
	bool interstitial;
	bool expand;
	bool fixed[3];
	int numToAdd;
	int numValues;
	List<int> steps;
	Element element;
	Vector3D position;
	
	// Variable to store atoms that will be added based on symmetry of structure
	OList<Element> expandElement;
	List<CoordinateType> expandPositionType;
	OList<Vector3D > expandPosition;
	List<bool> expandFixedSet;
	OList<List<bool> > expandFixed;
	List<bool> expandInterstitial;
	
	// Look for atom definition line
	int m;
	Label curLabel;
	List<Label> labels;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is empty
		if (!content[i].length())
			continue;
		
		// Current line does not define atoms
		curLabel = getLabel(content[i][0]);
		if ((curLabel == UNKNOWN) || (curLabel == BASIS) || (curLabel == SPACEGROUP))
			continue;
		
		// Get labels
		labels.clear();
		for (j = 0; j < content[i].length(); ++j)
		{
			
			// Break if a comment
			if (Language::isComment(content[i][j]))
				break;
			
			// Label not recognized
			curLabel = getLabel(content[i][j]);
			if ((curLabel == UNKNOWN) || (curLabel == BASIS) || (curLabel == SPACEGROUP))
			{
				Output::newline(ERROR);
				Output::print("Unknown label to describe atom in mint structure file (");
				Output::print(content[i][j]);
				Output::print(")");
				Output::quit();
			}
			
			// Save label
			labels += curLabel;
		}
		
		// Check which labels are being set and count number of values that should be on each line
		setElement = false;
		setFractional = false;
		setCartesian = false;
		setFixed = false;
		setMult = false;
		numValues = 0;
		steps.clear();
		for (j = 0; j < labels.length(); ++j)
		{
			if (labels[j] == ELEMENT)
			{
				setElement = true;
				steps += 1;
				numValues += 1;
			}
			else if (labels[j] == FRACPOS)
			{
				setFractional = true;
				steps += 3;
				numValues += 3;
			}
			else if (labels[j] == CARTPOS)
			{
				setCartesian = true;
				steps += 3;
				numValues += 3;
			}
			else if (labels[j] == FIXED)
			{
				setFixed = true;
				steps += 3;
				numValues += 3;
			}
			else if (labels[j] == MULTIPLICITY)
			{
				setMult = true;
				steps += 1;
				numValues += 1;
			}
			else if (labels[j] == EXPAND)
			{
				steps += 1;
				numValues += 1;
			}
			else if (labels[j] == INTERSTITIAL)
			{
				steps += 1;
				numValues += 1;
			}
			else
			{
				Output::newline(ERROR);
				Output::print("Internal error: unknown label id in mint structure file");
				Output::quit();
			}
		}
		
		// Make sure that element is a label
		if (!setElement)
		{
			Output::newline(ERROR);
			Output::print("Element must be defined in mint structure file");
			Output::quit();
		}
		
		// Both fractional and cartesian are set
		if ((setFractional) && (setCartesian))
		{
			Output::newline(ERROR);
			Output::print("Setting both fractional and cartesian coordinates in mint structure file");
			Output::quit();
		}
		
		// Multiplicity is set with anything other than element
		if ((setMult) && ((setFractional) || (setCartesian) || (setFixed)))
		{
			Output::newline(ERROR);
			Output::print("Multiplicity of atom can only be paired with element");
			Output::quit();
		}
		
		// Did not set fractional or cartesian coordinates
		if ((!setMult) && (!setFractional) && (!setCartesian))
		{
			Output::newline(ERROR);
			Output::print("Must set coordinates of atom in mint structure file");
			Output::quit();
		}
		
		// Setting coordinates when basis has not been set
		if ((!basisSet) && ((setFractional) || (setCartesian)))
		{
			Output::newline(ERROR);
			Output::print("Cannot set atom positions without a basis");
			Output::quit();
		}
		
		// Loop over rest of lines or until a label is found
		for (++i; i < content.length(); ++i)
		{
			
			// Skip if a blank line
			if (!content[i].length())
				continue;
			
			// Skip if a comment
			if (Language::isComment(content[i][0]))
				continue;
			
			// Found a label
			if (getLabel(content[i][0]) != UNKNOWN)
			{
				--i;
				break;
			}
			
			// Line is too short
			if (content[i].length() < numValues)
			{
				Output::newline(ERROR);
				Output::print("Not enough data for atom definition in mint structure file");
				Output::quit();
			}
			
			// Get element
			for (j = 0, k = 0; j < labels.length(); ++j)
			{
				
				// At element label
				if (labels[j] == ELEMENT)
				{
					element = Element::find(content[i][k], true, true);
					break;
				}
				
				// Move word index
				k += steps[j];
			}
			
			// Using multiplicity
			if (setMult)
			{
				
				// Loop over labels
				for (j = 0, k = 0; j < labels.length(); ++j)
				{
					
					// At multiplicity label
					if (labels[j] == MULTIPLICITY)
					{
						
						// Not an integer
						if (!Language::isInteger(content[i][k]))
						{
							Output::newline(ERROR);
							Output::print("Did not recognize multiplicity (");
							Output::print(content[i][k]);
							Output::print(")");
							Output::quit();
						}
						
						// Save multiplicity
						numToAdd = atoi(content[i][k].array());
						
						// Output
						Output::newline();
						Output::print("Adding ");
						Output::print(numToAdd);
						Output::print(" ");
						Output::print(element.symbol());
						Output::print(" atom");
						if (numToAdd != 1)
							Output::print("s");
						Output::print(" that will be assigned ");
						if (numToAdd == 1)
							Output::print("a ");
						Output::print("random position");
						if (numToAdd != 1)
							Output::print("s");
						
						// Add atoms
						for (m = 0; m < numToAdd; ++m)
							iso.addAtom(element);
						
						// Break since label was found
						break;
					}

					// Move word index
					k += steps[j];
				}
			}
			
			// Not using multiplicity
			else
			{
				
				// Loop over labels to get data
				expand = false;
				interstitial = false;
				fixed[0] = fixed[1] = fixed[2] = false;
				for (j = 0, k = 0; j < labels.length(); ++j)
				{
					
					// Get coordinates
					if ((labels[j] == FRACPOS) || (labels[j] == CARTPOS))
					{
						for (m = 0; m < 3; ++m)
						{
							if (Language::isNumber(content[i][k+m]))
								position[m] = atof(content[i][k+m].array());
							else
							{
								Output::newline(ERROR);
								Output::print("Did not recognize coordinate as a number (");
								Output::print(content[i][k+m]);
								Output::print(")");
								Output::quit();
							}
						}
					}
					
					// Get whether atom is fixed
					else if (labels[j] == FIXED)
					{
						for (m = 0; m < 3; ++m)
						{
							if ((content[i][k+m][0] == 't') || (content[i][k+m][0] == 'T'))
								fixed[m] = true;
							else if ((content[i][k+m][0] == 'f') || (content[i][k+m][0] == 'F'))
								fixed[m] = false;
							else
							{
								Output::newline(ERROR);
								Output::print("Did not recognize fixed value as true or false (");
								Output::print(content[i][k+m]);
								Output::print(")");
								Output::quit();
							}
						}
					}
					
					// Get whether an interstitial
					else if (labels[j] == INTERSTITIAL)
					{
						if ((content[i][k][0] == 't') || (content[i][k][0] == 'T'))
							interstitial = true;
						else if ((content[i][k][0] == 'f') || (content[i][k][0] == 'F'))
							interstitial = false;
						else
						{
							Output::newline(ERROR);
							Output::print("Did not recognize interstitial value as true or false (");
							Output::print(content[i][k]);
							Output::print(")");
							Output::quit();
						}
					}
					
					// Check if expanding
					else if (labels[j] == EXPAND)
					{
						if ((content[i][k][0] == 't') || (content[i][k][0] == 'T'))
							expand = true;
						else if ((content[i][k][0] == 'f') || (content[i][k][0] == 'F'))
							expand = false;
						else
						{
							Output::newline(ERROR);
							Output::print("Did not recognize expand value as true or false (");
							Output::print(content[i][k]);
							Output::print(")");
							Output::quit();
						}
					}
					
					// Move word index
					k += steps[j];
				}
				
				// Expanding current atom
				if (expand)
				{
					expandElement += element;
					if (setFractional)
						expandPositionType += FRACTIONAL;
					else
						expandPositionType += CARTESIAN;
					expandPosition += position;
					expandFixedSet += setFixed;
					expandFixed += List<bool>(3, fixed);
					expandInterstitial += interstitial;
				}
				
				// Not expanding current atom so add
				else
					addAtom(iso, element, position, setFractional, setFixed, fixed, interstitial);
			}
		}
	}
	
	// Output
	Output::decrease();
	
	// Run this section only if there are atoms that need to be expanded
	if (expandElement.length() > 0)
	{
		
		// Output
		Output::newline();
		Output::print("Generating atomic positions from symmetry");
		Output::increase();
		
		// Space group is set
		Symmetry symmetry;
		SpaceGroup spaceGroup;
		if (iso.spaceGroup().length() > 0)
			spaceGroup.set(iso.spaceGroup());
		
		// Get symmetry of structure
		else
			symmetry.set(iso, tol);
		
		// Save operations list
		const OList<SymmetryOperation>& operations = (iso.spaceGroup().length() > 0) ? spaceGroup.symmetry() : \
			symmetry.operations();
		
		// Loop over atoms to add
		Atom tempAtom;
		for (i = 0; i < expandElement.length(); ++i)
		{
			
			// Save atom properties
			tempAtom.basis(&iso.basis());
			tempAtom.element(expandElement[i]);
			if (expandPositionType[i] == FRACTIONAL)
				tempAtom.fractional(expandPosition[i]);
			else
				tempAtom.cartesian(expandPosition[i]);
			tempAtom.fixed(expandFixed[i].array());
			tempAtom.isInterstitial(expandInterstitial[i]);
			
			// Output
			Output::newline();
			Output::print("Expanding ");
			Output::print(tempAtom.element().symbol());
			Output::print(" at ");
			Output::print(tempAtom.fractional(), 8);
			if (expandFixedSet[i])
			{
				Output::print(" that is ");
				for (j = 0; j < 3; ++j)
				{
					if (tempAtom.fixed()[j])
						Output::print("fixed");
					else
						Output::print("free");
					if (j != 2)
						Output::print(", ");
					if (j == 1)
						Output::print("and ");
				}
			}
			Output::increase();
			
			// Add atom
			Symmetry::addAtom(iso, tempAtom, operations, clusterTol);
			
			// Output
			Output::decrease();
		}
		
		// Output
		Output::decrease();
	}
	
	// Output
	Output::decrease();
	
	// Return structure
	return iso;
}



/* void MintStructure::addAtom(ISO& iso, const Element& element, const Vector3D& position, bool useFractional,
 *		bool setFixed, const bool* fixed, bool interstitial)
 *
 * Add atom to structure
 */

void MintStructure::addAtom(ISO& iso, const Element& element, const Vector3D& position, bool useFractional, \
	bool setFixed, const bool* fixed, bool interstitial)
{
	
	// Output
	Output::newline();
	Output::print("Adding atom ");
	Output::print(iso.numAtoms() + 1);
	Output::print(" (");
	Output::print(element.symbol());
	Output::print(")");
	Output::increase();

	// Add atom
	Atom* atom = iso.addAtom(element);
	if (useFractional)
		atom->fractional(position);
	else
		atom->cartesian(position);
	if (setFixed)
		atom->fixed(fixed);
	atom->isInterstitial(interstitial);
	
	// Print coordinates
	Output::newline();
	Output::print("Fractional: ");
	Output::print(atom->fractional(), 8, false);
	Output::newline();
	Output::print("Cartesian: ");
	Output::print(atom->cartesian(), 8, false);

	// Print whether fixed if being set
	if (setFixed)
	{
		Output::newline();
		Output::print("Fixed: ");
		for (int i = 0; i < 3; ++i)
		{
			if (atom->fixed()[i])
				Output::print("true");
			else
				Output::print("false");
			if (i != 2)
				Output::print(", ");
		}
	}

	// Output
	Output::decrease();
}



/* void MintStructure::write(const Word& file, const ISO& iso, CoordinateType coordinates)
 *
 * Write structure to mint file
 */

void MintStructure::write(const Word& file, const ISO& iso, CoordinateType coordinates)
{
	
	// Precision for printing numbers
	int prec = 14;
	
	// Error if any sites are partially occupied
	if (iso.anyPartiallyOccupied())
	{
		Output::newline(ERROR);
		Output::print("Cannot print to Mint structure format with partially occupancies");
		Output::quit();
	}
	
	// Setup output if file was set
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::setStream(Output::addStream(file));
	}
	
	// Check if any basis metrics are fixed
	int i;
	bool printFixed = false;
	for (i = 0; i < 3; ++i)
	{
		if (iso.basis().lengthFixed()[i])
			printFixed = true;
	}
	for (i = 0; i < 3; ++i)
	{
		if (iso.basis().angleFixed()[i])
			printFixed = true;
	}
	
	// Set output method
	if (file.length() > 0)
		Output::method(STANDARD);
	
	// Do not print basis if it has not been set
	int j;
	Output message;
	if (iso.basis().volume() > 1e-8)
	{
		
		// Print the basis when using fractional coordinates
		if (coordinates == FRACTIONAL)
		{
			message.addLine();
			for (i = 0; i < 3; ++i)
				message.add(iso.basis().lengths()[i], prec);
			if (printFixed)
			{
				for (i = 0; i < 3; ++i)
				{
					if (iso.basis().lengthFixed()[i])
						message.add("t");
					else
						message.add("f");
				}
			}
			message.addLine();
			for (i = 0; i < 3; ++i)
				message.add(Num<double>::toDegrees(iso.basis().angles()[i]), prec);
			if (printFixed)
			{
				for (i = 0; i < 3; ++i)
				{
					if (iso.basis().angleFixed()[i])
						message.add("t");
					else
						message.add("f");
				}
			}
			Output::newline();
			Output::print("Basis");
			Output::newline();
			Output::print(message, RIGHT);
		}
		
		// Print the basis when using cartesian coordinates
		else
		{
			for (i = 0; i < 3; ++i)
			{
				message.addLine();
				for (j = 0; j < 3; ++j)
					message.add(iso.basis().vectors()(i, j), prec);
			}
			Output::newline();
			Output::print("Basis");
			Output::newline();
			Output::print(message, RIGHT);
			if (printFixed)
			{
				Output::newline();
				for (i = 0; i < 3; ++i)
				{
					if (iso.basis().lengthFixed()[i])
						Output::print("t ");
					else
						Output::print("f ");
				}
				for (i = 0; i < 3; ++i)
				{
					if (iso.basis().angleFixed()[i])
						Output::print("t ");
					else
						Output::print("f ");
				}
			}
		}
	}
	
	// Make list of atoms that have been assigned and those that have not
	int numAssigned = 0;
	int numUnassigned = 0;
	bool printInt = false;
	OList<Linked<Atom*> > assigned(iso.atoms().length());
	OList<Linked<Atom*> > unassigned(iso.atoms().length());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			if (iso.atoms()[i][j].assigned())
			{
				assigned[i] += &iso.atoms()[i][j];
				++numAssigned;
				if (iso.atoms()[i][j].isInterstitial())
					printInt = true;
			}
			else
			{
				unassigned[i] += &iso.atoms()[i][j];
				++numUnassigned;
			}
		}
	}
	
	// Print fixed atoms if present
	Linked<Atom*>::iterator it;
	if (numAssigned > 0)
	{
		
		// Figure out which things to print
		bool printFractional = (coordinates == FRACTIONAL) ? true : false;
		printFixed = iso.anyFixed();

		// Print line with information about atoms
		Output::newline();
		Output::print("Element");
		if (printFractional)
			Output::print(" Fractional");
		else
			Output::print(" Cartesian");
		if (printFixed)
			Output::print(" Fixed");
		if (printInt)
			Output::print(" Interstitial");

		// Create alignment list
		List<PrintAlign> align;
		align += LEFT;
		for (i = 0; i < 3; ++i)
			align += RIGHT;
		if (printFixed)
		{
			for (i = 0; i < 3; ++i)
				align += LEFT;
		}
		if (printInt)
			align += LEFT;

		// Loop over atoms
		message.clear();
		message.addLines(numAssigned);
		for (i = 0; i < assigned.length(); ++i)
		{
			for (it = assigned[i].begin(); it != assigned[i].end(); ++it)
			{

				// Add element
				message.addLine();
				message.add((**it).element().symbol());

				// Add position
				if (printFractional)
				{
					for (j = 0; j < 3; ++j)
						message.add((**it).fractional()[j], prec);
				}
				else
				{
					for (j = 0; j < 3; ++j)
						message.add((**it).cartesian()[j], prec);
				}

				// Add fixed
				if (printFixed)
				{
					for (j = 0; j < 3; ++j)
					{
						if ((**it).fixed()[j])
							message.add("t");
						else
							message.add("f");
					}
				}
				
				// Add interstitial
				if (printInt)
				{
					if ((**it).isInterstitial())
						message.add("t");
					else
						message.add("f");
				}
			}
		}

		// Print atoms
		Output::newline();
		Output::print(message, align);
	}
	
	// Print atoms that are not fixed if present
	if (numUnassigned > 0)
	{
		
		// Print line with information about atoms
		Output::newline();
		Output::print("Element Multiplicity");

		// Loop over elements
		message.clear();
		for (i = 0; i < unassigned.length(); ++i)
		{
			
			// Skip if no atoms of current element
			if (!unassigned[i].length())
				continue;
			
			// Add element
			message.addLine();
			message.add((**unassigned[i].begin()).element().symbol());
			message.add(unassigned[i].length());
		}
		
		// Set alignment
		List<PrintAlign> align (2);
		align[0] = LEFT;
		align[1] = RIGHT;
		
		// Print atoms
		Output::newline();
		Output::print(message, align);
	}
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}



/* bool MintStructure::isFormat(const Text& content)
 *
 * Return whether file contains mint structure
 */

bool MintStructure::isFormat(const Text& content)
{
	
	// Loop over lines in the file
	Label curLabel;
	for (int i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is too short
		if (content[i].length() < 2)
			continue;
		
		// First word is not a label or is basis/spacegroup
		curLabel = getLabel(content[i][0]);
		if ((curLabel == UNKNOWN) || (curLabel == BASIS) || (curLabel == SPACEGROUP))
			continue;
		
		// Second word is not a label
		if (getLabel(content[i][1]) == UNKNOWN)
			return false;
		
		// First two words are labels
		return true;
	}
	
	// Return that file is not mint structure format if at this point
	return false;
}



/* MintStructure::Label MintStructure::getLabel(const Word& word)
 *
 * Return the label cooresponding to a word
 */

MintStructure::Label MintStructure::getLabel(const Word& word)
{
	if (word.equal("space", false, 5))
		return SPACEGROUP;
	if (word.equal("basis", false, 3))
		return BASIS;
	if (word.equal("element", false, 4))
		return ELEMENT;
	if (word.equal("fractional", false, 4))
		return FRACPOS;
	if (word.equal("cartesian", false, 4))
		return CARTPOS;
	if (word.equal("fixed", false, 3))
		return FIXED;
	if (word.equal("multiplicity", false, 4))
		return MULTIPLICITY;
	if (word.equal("expand", false, 3))
		return EXPAND;
	if (word.equal("interstitial", false, 3))
		return INTERSTITIAL;
	return UNKNOWN;
}
