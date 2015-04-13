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



#include "cif.h"
#include "elements.h"
#include "symmetry.h"
#include "spaceGroup.h"
#include "language.h"
#include "output.h"
#include "num.h"
#include <cstdlib>



/* ISO CIF::read(const Text& content, double clusterTol)
 *
 * Convert CIF file to structure
 */

ISO CIF::read(const Text& content, double clusterTol)
{
	
	// Output
	Output::newline();
	Output::print("Reading CIF file");
	Output::increase();
	
	// Variables to store cell properties
	int Z = 1;
	Word spaceGroupName;
	Vector3D lengths(0.0);
	Vector3D angles(0.0);
	Vector3D position;
	Linked<Atom> atoms;
	OList<SymmetryOperation> symmetry;
	OList<Element> elements;
	List<double> composition;
	
	// Loop over contents of file
	int i, j, k;
	bool saveAny;
	bool addAtom;
	bool unknown;
	int index;
	int elemEnd;
	Tag tag;
	List<Tag> tags;
	Words curLine;
	Words curData;
	Atom tempAtom;
	Matrix3D rotation;
	Vector3D translation;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is empty
		if (!content[i].length())
			continue;
		
		// Found cell parameter
		tag = getTag(content[i][0]);
		checkTag(tag, content[i], lengths,angles, spaceGroupName, elements, composition, Z);
		
		// Found start of loop
		if (content[i][0] == "loop_")
		{
			
			// Loop over lines to get tags
			tags.length(0);
			for (++i; i < content.length(); ++i)
			{
				if (!content[i].length())
					continue;
				if (content[i][0][0] != '_')
					break;
				tags += getTag(content[i][0]);
				if (checkTag(tags.last(), content[i], lengths,angles, spaceGroupName, elements, composition, Z))
					tags.remove(tags.length() - 1);
			}
			
			// Check if any tags are needed
			saveAny = false;
			for (j = 0; j < tags.length(); ++j)
			{
				if (tags[j] != tag_unknown)
				{
					saveAny = true;
					break;
				}
			}
			
			// Loop until end of data
			index = 0;
			addAtom = false;
			tempAtom.clear();
			position = 0.0;
			curData.length(tags.length());
			for (; i < content.length(); ++i)
			{
				
				// Check if line is empty or loop is over
				if (!content[i].length())
					break;
				if (Language::isComment(content[i][0]))
					break;
				if ((content[i][0][0] == '_') || (content[i][0] == "loop_"))
				{
					--i;
					break;
				}
				
				// Skip line if not saving any data
				if (!saveAny)
					continue;
				
				// Group data on line
				curLine = getLine(content[i]);
				
				// Add contents
				for (j = 0; j < curLine.length(); ++j)
				{
					
					// Add word
					curData[index++] = curLine[j];
					
					// Finished with current set of data
					if (index == curData.length())
					{
						
						// Analyze data
						unknown = false;
						for (k = 0; k < tags.length(); ++k)
						{
							
							// Check if not known
							if (curData[k][0] == '?')
								unknown = true;

							// Save symmetry operation
							if ((tags[k] == space_group_symop_operation_xyz) || (tags[k] == symmetry_equiv_pos_as_xyz))
							{
								JonesFaithful::fromString(rotation, translation, curData[k]);
								symmetry.add();
								symmetry.last().setRotation(rotation);
								symmetry.last().addTranslation(translation);
							}

							// Setting atom property
							if ((tags[k] == atom_site_label) || (tags[k] == atom_site_fract_x) || \
								(tags[k] == atom_site_fract_y) || (tags[k] == atom_site_fract_z))
								addAtom = true;

							// Get element
							if ((!tempAtom.element().number()) && \
								((tags[k] == atom_site_label) || (tags[k] == atom_site_type_symbol)))
							{
								elemEnd = 0;
								if (curData[k].length() >= 3)
								{
									if (Element::isElement(Word(3, curData[j]), false))
										elemEnd = 3;
								}
								if ((!elemEnd) && (curData[k].length() >= 2))
								{
									if (Element::isElement(Word(2, curData[k]), false))
										elemEnd = 2;
								}
								if (!elemEnd)
								{
									if (Element::isElement(Word(1, curData[k]), false))
										elemEnd = 1;
								}
								if (elemEnd > 0)
									tempAtom.element(Element::find(Word(elemEnd, curData[k])));
							}

							// Save position
							else if (tags[k] == atom_site_fract_x)
								position[0] = atof(curData[k].array());
							else if (tags[k] == atom_site_fract_y)
								position[1] = atof(curData[k].array());
							else if (tags[k] == atom_site_fract_z)
								position[2] = atof(curData[k].array());

							// Save occupancy
							else if (tags[k] == atom_site_occupancy)
								tempAtom.occupancy(atof(curData[k].array()));
						}
						
						// Add atom if needed
						if ((addAtom) && (!unknown))
						{

							// Element was not set
							if (!tempAtom.element().number())
							{
								Output::newline(ERROR);
								Output::print("Saving atom when the element has not been set");
								Output::quit();
							}

							// Add atom
							tempAtom.fractional(position);
							atoms += tempAtom;
						}
						
						// Reset variables
						index = 0;
						addAtom = false;
						tempAtom.clear();
						position = 0.0;
					}
				}
			}
			
			// Number of items found does not match number of tags
			if (index != 0)
			{
				Output::newline(ERROR);
				Output::print("Did not find the expected amount of data (");
				Output::print(tags.length());
				Output::print(")");
				Output::quit();
			}
			
			// Go to next line
			continue;
		}
	}
	
	// Lattice parameters were not set
	for (i = 0; i < 3; ++i)
	{
		if ((lengths[i] == 0) || (angles[i] == 0))
		{
			Output::newline(ERROR);
			Output::print("Basis not properly set");
			Output::quit();
		}
	}
	
	// Figure out symmetry to use
	bool setBySpaceGroup = false;
	SpaceGroup spaceGroup;
	if ((!symmetry.length()) && (!spaceGroupName.length()))
	{
		Output::newline();
		Output::print("No symmetry information was found so assuming P1 space group");
		symmetry.add();
		symmetry.last().setRotation(Matrix3D::identity());
		symmetry.last().addTranslation(Vector3D(0.0));
	}
	else if (symmetry.length())
	{
		Output::newline();
		Output::print("Found ");
		Output::print(symmetry.length());
		Output::print(" symmetry operation");
		if (symmetry.length() != 1)
			Output::print("s");
		Output::increase();
		for (i = 0; i < symmetry.length(); ++i)
		{
			Output::newline();
			Output::print(JonesFaithful::toString(symmetry[i].rotation(), &symmetry[i].translations()[0]), true, false);
		}
		Output::decrease();
	}
	else
	{
		setBySpaceGroup = true;
		spaceGroup.set(spaceGroupName);
	}
	
	// Save symmetry to use
	const OList<SymmetryOperation>& operations = (setBySpaceGroup) ? spaceGroup.symmetry() : symmetry;
	
	// Save space group
	ISO iso;
	if (spaceGroupName.length() > 0)
		iso.spaceGroup(spaceGroupName);
	
	// Save basis
	iso.basis(lengths, angles);
	
	// Output
	Output::newline();
	Output::print("Adding atoms");
	Output::increase();
	
	// Loop over atoms to add
	for (Linked<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
	{
		
		// Output
		Output::newline();
		Output::print("Adding unique ");
		Output::print((*it).element().symbol());
		Output::print(" atom at ");
		Output::print((*it).fractional(), 8);
		if ((*it).occupancy() != 1)
		{
			Output::print(" with occupancy ");
			Output::print((*it).occupancy());
		}
		Output::increase();
		
		// Add atom
		Symmetry::addAtom(iso, *it, operations, clusterTol);
		
		// Output
		Output::decrease();
	}
	
	// Only run this section if composition was set
	if (elements.length() > 0)
	{
		
		// Loop over elements that were defined in the composition
		double numToAdd;
		for (i = 0; i < elements.length(); ++i)
		{
			
			// Loop over elements in the structure
			numToAdd = Z * composition[i];
			for (j = 0; j < iso.atoms().length(); ++j)
			{
				
				// Elements are the same
				if (elements[i] == iso.atoms()[j][0].element())
				{
					
					// Set the number to add
					numToAdd -= iso.atoms()[j].length();
					break;
				}
			}
			
			// Not an integer number of atoms
			if (Num<double>::abs(numToAdd - Num<double>::round(numToAdd, 1)) > 1e-6)
			{
				Output::newline(ERROR);
				Output::print("Composition sets ");
				Output::print(numToAdd);
				Output::print(" ");
				Output::print(elements[i].symbol());
				Output::print(" atom");
				if (numToAdd != 1)
					Output::print("s");
				Output::print(" (expecting an integer)");
				Output::quit();
			}
			
			// Add elements
			numToAdd = Num<double>::round(numToAdd, 1);
			if (numToAdd > 0)
			{
				
				// Output
				Output::newline();
				Output::print("Adding ");
				Output::print(numToAdd);
				Output::print(" ");
				Output::print(elements[i].symbol());
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
				for (j = 0; j < numToAdd; ++j)
					iso.addAtom(elements[i]);
			}
			
			// Composition was wrong
			else if (numToAdd < 0)
			{
				Output::newline(ERROR);
				Output::print("Generated ");
				Output::print(iso.atoms()[j].length());
				Output::print(" ");
				Output::print(iso.atoms()[j][0].element().symbol());
				Output::print(" atom");
				if (iso.atoms()[j].length() != 1)
					Output::print("s");
				Output::print(" when expecting ");
				Output::print((int)(Z*composition[i]));
				Output::quit();
			}
		}
	}
	
	// Output
	Output::decrease();
	
	// Output
	Output::decrease();
	
	// Return structure
	return iso;
}



/* bool CIF::checkTag(Tag tag, OList<Word>& curLine, Vector3D& lengths, Vector3D& angles, \
 *		Word& spaceGroupName, OList<Element>& elements, List<double>& composition, int& Z)
 *
 * Check whether one in a list of tags is found on line
 */

bool CIF::checkTag(Tag tag, OList<Word>& curLine, Vector3D& lengths, Vector3D& angles, \
	Word& spaceGroupName, OList<Element>& elements, List<double>& composition, int& Z)
{
	
	// Check for cell parameter
	if ((tag == cell_length_a) && (curLine.length() > 1))
	{
		lengths[0] = atof(curLine[1].array());
		return true;
	}
	if ((tag == cell_length_b) && (curLine.length() > 1))
	{
		lengths[1] = atof(curLine[1].array());
		return true;
	}
	if ((tag == cell_length_c) && (curLine.length() > 1))
	{
		lengths[2] = atof(curLine[1].array());
		return true;
	}
	if ((tag == cell_angle_alpha) && (curLine.length() > 1))
	{
		angles[0] = Num<double>::toRadians(atof(curLine[1].array()));
		return true;
	}
	if ((tag == cell_angle_beta) && (curLine.length() > 1))
	{
		angles[1] = Num<double>::toRadians(atof(curLine[1].array()));
		return true;
	}
	if ((tag == cell_angle_gamma) && (curLine.length() > 1))
	{
		angles[2] = Num<double>::toRadians(atof(curLine[1].array()));
		return true;
	}
	
	// Found space group
	if ((tag == symmetry_Int_Tables_number) || (tag == space_group_IT_number) 
			|| (tag == space_group_name_H_M) || (tag == space_group_name_Hall) )
	{
		Words line = getLine(curLine);
		if (line.length() > 1)
		{
			spaceGroupName = getLine(curLine)[1];
			return true;
		}
	}
	
	// Found composition
	if (tag == chemical_formula_sum)
	{
		
		// Loop until end of line
		Word element;
		Word number;
		for (int i = 1; i < curLine.length(); ++i)
		{
			
			// Remove ' symbols
			curLine[i].strip('\'');
			
			// Get element and multiplicity
			getElementAndNumber(element, number, curLine[i]);
			elements += Element::find(element, false, true);
			composition += atof(number.array());
		}
		return true;
	}
	
	// Found number of formula units in cell
	if ((tag == cell_formula_units_Z) && (curLine.length() >= 1))
	{
		curLine[1].strip('\'');
		Z = atoi(curLine[1].array());
		return true;
	}
	
	// Return that tag was not found
	return false;
}



/**
 * Write CIF file
 * @param file [in] Filename
 * @param iso [in] Structure to be written
 * @param tol [in] Tolerance used when computing symmetry of structure
 */
void CIF::write(const Word& file, const ISO& iso, double tol)
{
	
	// Precision for printing numbers
	int prec = 14;
	
	// Get symmetry of the structure
	Symmetry symmetry(iso, tol);

	// Get the space group of the structure
	SpaceGroup spg(iso, tol);
	
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
	
	// Get the composition
	int i;
	int composition[iso.atoms().length()];
	for (i = 0; i < iso.atoms().length(); ++i)
		composition[i] = iso.atoms()[i].length();
	
	// Set the formula unit
	int Z = Num<int>::gcf(iso.atoms().length(), composition);
	for (i = 0; i < iso.atoms().length(); ++i)
		composition[i] /= Z;
	
	// Print top line
	Output::newline();
	Output::print("data_");
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		Output::print(iso.atoms()[i][0].element().symbol());
		if (iso.atoms()[i].length() != 1)
			Output::print(iso.atoms()[i].length());
	}
	Output::print("_mint");
	
	// Print the chemical formula
	Output::newline();
	Output::print("_chemical_formula_sum \'");
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		Output::print(iso.atoms()[i][0].element().symbol());
		if (composition[i] != 1)
			Output::print(composition[i]);
		if (i != iso.atoms().length() - 1)
			Output::print(" ");
	}
	Output::print("\'");
	
	// Print the number of formula units
	Output::newline();
	Output::print("_cell_formula_units_Z ");
	Output::print(Z);
	
	// Print the basis
	Output::newline();
	Output::print("_cell_length_a ");
	Output::print(iso.basis().lengths()[0], prec);
	Output::newline();
	Output::print("_cell_length_b ");
	Output::print(iso.basis().lengths()[1], prec);
	Output::newline();
	Output::print("_cell_length_c ");
	Output::print(iso.basis().lengths()[2], prec);
	Output::newline();
	Output::print("_cell_angle_alpha ");
	Output::print(Num<double>::toDegrees(iso.basis().angles()[0]), prec);
	Output::newline();
	Output::print("_cell_angle_beta ");
	Output::print(Num<double>::toDegrees(iso.basis().angles()[1]), prec);
	Output::newline();
	Output::print("_cell_angle_gamma ");
	Output::print(Num<double>::toDegrees(iso.basis().angles()[2]), prec);

	// Print space group
	Output::newline();
	Output::print("_space_group_name_Hall '");
	Output::print(spg.hall());
	Output::print("'");
	
	// Print symmetry operations
	int j;
	Output::newline();
	Output::print("loop_");
	Output::newline();
	Output::print("_space_group_symop_operation_xyz");
	for (i = 0; i < symmetry.operations().length(); ++i)
	{
		for (j = 0; j < symmetry.operations()[i].translations().length(); ++j)
		{
			Output::newline();
			Output::print("'");
			Output::print(symmetry.operations()[i].getString(j), true, false);
			Output::print("'");
		}
	}
	
	// Allocate space for atoms to print
	Output message;
	message.addLines(symmetry.orbits().length());
	
	// Decide whether to print occupancies
	bool anyPartiallyOccupied = iso.anyPartiallyOccupied();
	
	// Save unique atoms
	int elemCount = 0;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		if (i)
		{
			if (symmetry.orbits()[i].atoms()[0]->element() != symmetry.orbits()[i-1].atoms()[0]->element())
				elemCount = 0;
		}
		++elemCount;
		message.addLine();
		message.addWords(5);
		message.add(symmetry.orbits()[i].atoms()[0]->element().symbol() + Language::numberToWord(elemCount));
		message.add(symmetry.orbits()[i].atoms()[0]->element().symbol());
		for (j = 0; j < 3; ++j)
			message.add(symmetry.orbits()[i].atoms()[0]->fractional()[j], prec);
		if (anyPartiallyOccupied)
			message.add(symmetry.orbits()[i].atoms()[0]->occupancy());
	}
	
	// Set alignment
	List<PrintAlign> align(5);
	align[0] = align[1] = LEFT;
	align[2] = align[3] = align[4] = RIGHT;
	if (anyPartiallyOccupied)
	{
		align.length(6);
		align[5] = LEFT;
	}
	
	// Print atoms
	Output::newline();
	Output::print("loop_");
	Output::newline();
	Output::print("_atom_site_label");
	Output::newline();
	Output::print("_atom_site_type_symbol");
	Output::newline();
    Output::print("_atom_site_fract_x");
	Output::newline();
	Output::print("_atom_site_fract_y");
	Output::newline();
	Output::print("_atom_site_fract_z");
	if (anyPartiallyOccupied)
	{
		Output::newline();
		Output::print("_atom_site_occupancy");
	}
	Output::newline();
	Output::print(message, align);
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}



/* bool CIF::isFormat(const Text& content)
 *
 * Return whether file is in CIF format
 */

bool CIF::isFormat(const Text& content)
{
	
	// Loop over lines in the file and look for _atom_site... tag
	for (int i = 0; i < content.length(); ++i)
	{
		
		// Line is empty
		if (!content[i].length())
			continue;
		
		// Found tag
		if (content[i][0].equal("_atom_site", false, 10))
			return true;
		if (content[i][0].equal("_chemical", false, 9))
			return true;
	}
	
	// Return that file is not in correct format if at this point
	return false;
}



/* Words CIF::getLine(const OList<Word>& line)
 *
 * Get words from line
 */

Words CIF::getLine(const OList<Word>& line)
{
	
	// Loop over words on line
	Words res;
	for (int i = 0; i < line.length(); ++i)
	{
		
		// Add current word
		res += line[i];
		
		// If current word starts with ' then add words until closing ' is found
		if ((line[i][0] == '\'') && (line[i].last() != '\''))
		{
			for (++i; i < line.length(); ++i)
			{
				res.last() += ' ';
				res.last() += line[i];
				if (line[i].last() == '\'')
					break;
			}
		}
		
		// Strip any ' that appear in word
		res.last().strip('\'');
	}
	
	// Return result
	return res;
}



/* CIF::Tag CIF::getTag(const Word& word)
 *
 * Get tag from word
 */

CIF::Tag CIF::getTag(const Word& word)
{
	if (word == "_atom_site_label")					return atom_site_label;
	if (word == "_atom_site_type_symbol")			return atom_site_type_symbol;
	if (word == "_atom_site_fract_x")				return atom_site_fract_x;
	if (word == "_atom_site_fract_y")				return atom_site_fract_y;
	if (word == "_atom_site_fract_z")				return atom_site_fract_z;
	if (word == "_atom_site_occupancy")				return atom_site_occupancy;
	if (word == "_cell_length_a")					return cell_length_a;
	if (word == "_cell_length_b")					return cell_length_b;
	if (word == "_cell_length_c")					return cell_length_c;
	if (word == "_cell_angle_alpha")				return cell_angle_alpha;
	if (word == "_cell_angle_beta")					return cell_angle_beta;
	if (word == "_cell_angle_gamma")				return cell_angle_gamma;
	if (word == "_space_group_IT_number")			return space_group_IT_number;
	if (word == "_space_group_symop_operation_xyz")	return space_group_symop_operation_xyz;
	if (word == "_symmetry_space_group_name_H-M")	return space_group_name_H_M;
	if (word == "_space_group_name_Hall")			return space_group_name_Hall;
	if (word == "_symmetry_Int_Tables_number")		return symmetry_Int_Tables_number;
	if (word == "_symmetry_equiv_pos_as_xyz")		return symmetry_equiv_pos_as_xyz;
	if (word == "_chemical_formula_sum")			return chemical_formula_sum;
	if (word == "_cell_formula_units_Z")			return cell_formula_units_Z;
	return tag_unknown;
}



/* void CIF::getElementAndNumber(Word& element, Word& number, const Word& word)
 *
 * Split a word into element and number
 */

void CIF::getElementAndNumber(Word& element, Word& number, const Word& word)
{
	
	// Clear space
	element.clear();
	number.clear();
	
	// Loop over characters in word
	for (int i = 0; i < word.length(); ++i)
	{
		
		// Found a number
		if ((word[i] == '0') || (word[i] == '1') || (word[i] == '2') || (word[i] == '3') || (word[i] == '4') || \
			(word[i] == '5') || (word[i] == '6') || (word[i] == '7') || (word[i] == '8') || (word[i] == '9') || \
			(word[i] == '.'))
			number += word[i];
		
		// Found anything else
		else
			element += word[i];
	}
	
	// If no number was found then set to 1
	if (!number.length())
		number += '1';
}
