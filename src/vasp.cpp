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
#include "vasp.h"
#include "kpoints.h"
#include "structureIO.h"
#include "language.h"
#include "output.h"
#include <fstream>
#include <cstdlib>
#include <unistd.h>
using namespace std;



// Name of directory to run in
Word Vasp::Job::_dirName = "VaspRun";

// Name of previous run directory
bool Vasp::Job::_completed = false;
Word Vasp::Job::_prevDir;


/* Vasp::Accuracy Vasp::accuracy(const Word& word)
 *
 * Return accuracy from word
 */

Vasp::Accuracy Vasp::accuracy(const Word& word)
{
	if (word.equal("high", false, 4))
		return VA_HIGH;
	if ((word.equal("medium", false, 3)) || (word.equal("normal", false, 4)))
		return VA_NORMAL;
	if (word.equal("low", false, 3))
		return VA_LOW;
	return VA_UNKNOWN;
}



/* Word Vasp::accuracy(Accuracy accuracy)
 *
 * Return word from accuracy
 */

Word Vasp::accuracy(Accuracy accuracy)
{
	switch (accuracy)
	{
		case VA_HIGH:
			return Word("High");
		case VA_NORMAL:
			return Word("Normal");
		case VA_LOW:
			return Word("Low");
		default:
			return Word("Unknown");
	}
}



/* bool Vasp::Job::submit(const ISO& iso, bool restart, bool quitIfError)
 *
 * Submit vasp job
 */

bool Vasp::Job::submit(const ISO& iso, bool restart, bool quitIfError)
{
	
	// Clear storage
	_energy.clear();
	_convergence.clear();
	_forces.clear();
	_structure.clear();
	
	// Save current directory
	Word origDir;
	Directory::current(origDir);
	
	// Decide whether to restart
	bool allowRestart = false;
	if (_prevDir.length())
	{
		if (Directory::exists(_prevDir))
			allowRestart = ((restart) && (_completed));
	}
	
	// Create directory if it does not already exist
	Word runDir;
	if (allowRestart)
		runDir = _prevDir;
	else
	{
		if (_saveFiles)
		{
			for (int i = 1;; ++i)
			{
				runDir = Directory::makePath(origDir, _dirName + Word("_") + Language::numberToWord(i));
				if (!Directory::exists(runDir))
					break;
			}
		}
		else
			runDir = Directory::makePath(origDir, _dirName);
		Directory::create(runDir, true);
	}
	
	// Switch to run directory
	Directory::change(runDir);
	
	// Write files
	writeFiles(iso, allowRestart);
	
	// Run vasp
	Multi::external(_executable, 0, "stdout");
	
	// Get job results
	_completed = false;
	gatherResults();

	// Error if job did not complete
	if (!_completed)
	{
		if (quitIfError)
			Output::newline(ERROR);
		else
			Output::newline(WARNING);
		Output::print("VASP did not exit properly");
		if (quitIfError)
			Output::quit();
	}
	
	// Update previous directory
	_prevDir = runDir;

	// Switch back to run directory
	Directory::change(origDir);
	
	// Return whether run was successful
	return _completed;
}



/* bool Vasp::Job::submitNEB(const OList<ISO>& iso, bool restart, bool quitIfError)
 *
 * Submit Vasp NEB job
 */

bool Vasp::Job::submitNEB(const OList<ISO>& iso, bool restart, bool quitIfError)
{
	
	// Clear storage
	_energy.clear();
	_convergence.clear();
	_forces.clear();
	_structure.clear();
	
	// Save current directory
	Word origDir;
	Directory::current(origDir);
	
	// Decide whether to restart
	bool allowRestart = false;
	if (_prevDir.length())
	{
		if (Directory::exists(_prevDir))
			allowRestart = ((restart) && (_completed));
	}
	
	// Create directory if it does not already exist
	Word runDir;
	if (allowRestart)
		runDir = _prevDir;
	else
	{
		if (_saveFiles)
		{
			for (int i = 1;; ++i)
			{
				runDir = Directory::makePath(origDir, _dirName + Word("_") + Language::numberToWord(i));
				if (!Directory::exists(runDir))
					break;
			}
		}
		else
			runDir = Directory::makePath(origDir, _dirName);
		Directory::create(runDir, true);
	}
	
	// Switch to run directory
	Directory::change(runDir);
	
	// Write files
	writeNEBFiles(iso, allowRestart);
	
	// Run vasp
	Multi::external(_executable, 0, "stdout");
	
	// Get job results
	_completed = false;
	gatherNEBResults(iso.length());

	// Error if job did not complete
	if (!_completed)
	{
		if (quitIfError)
			Output::newline(ERROR);
		else
			Output::newline(WARNING);
		Output::print("VASP did not exit properly");
		if (quitIfError)
			Output::quit();
	}
	
	// Update previous directory
	_prevDir = runDir;

	// Switch back to run directory
	Directory::change(origDir);
	
	// Return whether run was successful
	return _completed;
}



/* void Vasp::Job::writeFiles(const ISO& iso, bool restart)
 *
 * Write files for job
 */

void Vasp::Job::writeFiles(const ISO& iso, bool restart)
{
	
	// Write potential
	List<double> encut;
	Potential::write("POTCAR", _elements, _potentials, iso, &encut);
	
	// Write settings
	Settings settings;
	settings.defaults(_calculation, _relaxation, _material, _accuracy, &encut);
	if (restart)
	{
		settings.add("ISTART", "1", true, false);
		settings.add("ICHARG", "0", true, false);
	}
	else
	{
		settings.add("ISTART", "0", true, false);
		settings.add("ICHARG", "1", true, false);
	}
	for (int i = 0; i < _tags.length(); ++i)
		settings.add(_tags[i][0], _tags[i][1], true, true);
	settings.write("INCAR");
	
	// Write structure
	StructureIO::write("POSCAR", iso, SF_VASP4);
	
	// Write kpoints file
	KPoints kpoints;
	if (_material == VM_METAL)
	{
		if (_accuracy == VA_HIGH)
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraMetalHigh());
		else if (_accuracy == VA_NORMAL)
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraMetalNormal());
		else
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraMetalLow());
	}
	else
	{
		if (_accuracy == VA_HIGH)
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraInsulatorHigh());
		else if (_accuracy == VA_NORMAL)
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraInsulatorNormal());
		else
			kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), KPoints::kppraInsulatorLow());
	}
	kpoints.write("KPOINTS", KPoints::KP_VASP);
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void Vasp::Job::writeNEBFiles(const OList<ISO>& iso, bool restart)
 *
 * Write files to be used during NEB run
 */

void Vasp::Job::writeNEBFiles(const OList<ISO>& iso, bool restart)
{
	
	// Write potential
	List<double> encut;
	Potential::write("POTCAR", _elements, _potentials, iso[0], &encut);
	
	// Write settings
	int i;
	Settings settings;
	settings.defaults(_calculation, _relaxation, _material, _accuracy, &encut);
	if (restart)
	{
		settings.add("ISTART", "1", true, false);
		settings.add("ICHARG", "0", true, false);
	}
	else
	{
		settings.add("ISTART", "0", true, false);
		settings.add("ICHARG", "1", true, false);
	}
	for (i = 0; i < _tags.length(); ++i)
		settings.add(_tags[i][0], _tags[i][1], true, true);
	settings.write("INCAR");
	
	// Write kpoints file
	KPoints kpoints;
	if (_material == VM_METAL)
	{
		if (_accuracy == VA_HIGH)
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraMetalHigh());
		else if (_accuracy == VA_NORMAL)
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraMetalNormal());
		else
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraMetalLow());
	}
	else
	{
		if (_accuracy == VA_HIGH)
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraInsulatorHigh());
		else if (_accuracy == VA_NORMAL)
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraInsulatorNormal());
		else
			kpoints.mesh(iso[0].numAtoms(), iso[0].basis().lengths(), KPoints::kppraInsulatorLow());
	}
	kpoints.write("KPOINTS", KPoints::KP_VASP);
	
	// Save current directory
	Word origDir;
	Directory::current(origDir);
	
	// Create run directories and copy files
	Word dir;
	for (i = 0; i < iso.length(); ++i)
	{
		
		// Create directory
		dir  = (i < 10) ? "0" : "";
		dir += Language::numberToWord(i);
		dir  = Directory::makePath(origDir, dir);
		Directory::create(dir, true);
		
		// Move to directory and write file
		Directory::change(dir);
		StructureIO::write("POSCAR", iso[i], SF_VASP4);
		Directory::change(origDir);
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void Vasp::Job::gatherResults()
 *
 * Save the results from Vasp job
 */

void Vasp::Job::gatherResults()
{
	
	// Open vasprun.xml file to check for completion
	ifstream input("vasprun.xml");
	if (!input.is_open())
		return;
	
	// Process vasprun.xml file
	bool save = false;
	int count;
	char temp[1000];
	while (input >> temp)
	{
		
		// Check for end of run
		if (Text::equal(temp, "</modeling>"))
		{
			_completed = true;
			continue;
		}
		
		// Check for start of forces
		if (Text::contains(temp, "name=\"forces\""))
		{
			save = true;
			_forces.add();
			continue;
		}
		
		// Check for end of forces
		if (Text::equal(temp, "</varray>"))
		{
			save = false;
			continue;
		}
		
		// Save force
		if (save)
		{
			
			// Add new force if on <v>
			if (Text::equal(temp, "<v>"))
			{
				count = 0;
				_forces.last().add();
				_forces.last().last() = 0.0;
				continue;
			}
			
			// Skip if junk
			if ((Text::contains(temp, "<")) || (Text::contains(temp, ">")))
				continue;
			
			// Save force
			_forces.last().last()[count++] = atof(temp);
		}
	}
	
	// Close file
	input.close();
	
	// Return if run did not exit properly
	if (!_completed)
		return;
	
	// Save final structure
	Output::quietOn();
	if (!File::exists("CONTCAR"))
	{
		_completed = false;
		return;
	}
	Text contcar = Read::text("CONTCAR");
	if (Vasp::Structure::isVersion4(contcar))
		_structure = Vasp::Structure::read(contcar, false);
	else if (Vasp::Structure::isVersion5(contcar))
		_structure = Vasp::Structure::read(contcar, true);
	else
	{
		_completed = false;
		Output::quietOff();
		return;
	}
	Output::quietOff();
	
	// Save data in OSZICAR file
	int i;
	Text oszicar = Read::text("OSZICAR");
	for (i = 0; i < oszicar.length(); ++i)
	{
		
		// Skip if line is too short
		if (oszicar[i].length() < 3)
			continue;
		
		// Found energy
		if (oszicar[i][1] == "F=")
		{
			_energy += atof(oszicar[i][2].array());
			continue;
		}
		
		// Found start of new electronic convergence loop
		if (oszicar[i][1] == "1")
			_convergence.add();
		
		// Found electronic convergence step
		if (Language::isNumber(oszicar[i][1]))
			_convergence.last() += atof(oszicar[i][2].array());
	}
	
	// Check if the last electronic convergence loop converged
	if (_convergence.last().length() == 60)
		_completed = false;
}



/* void Vasp::Job::gatherNEBResults(int numDirs)
 *
 * Save the results from Vasp NEB job
 */

void Vasp::Job::gatherNEBResults(int numDirs)
{
	
	// Make sure there is enough room in structures
	_structures.length(numDirs);
	
	// Save current directory
	Word origDir;
	Directory::current(origDir);
	
	// Loop over directories and make sure that each completed and get final energy for each image
	int i;
	bool save;
	double energy;
	double maxEnergy;
	char temp[1000];
	Word dir;
	Word tsDir;
	Text data;
	ifstream input;
	for (i = 1; i < numDirs - 1; ++i)
	{
		
		// Move into current directory
		dir  = (i < 10) ? "0" : "";
		dir += Language::numberToWord(i);
		dir  = Directory::makePath(origDir, dir);
		Directory::change(dir);
		
		// Open vasprun.xml file to check for completion
		input.open("vasprun.xml");
		if (!input.is_open())
			return;
		
		// Process file
		while (input >> temp)
		{
			if (Text::equal(temp, "</modeling>"))
			{
				_completed = true;
				break;
			}
		}
		
		// Close vasprun.xml
		input.close();
		
		// Return if run did not exit properly
		if (!_completed)
			return;
		
		// Save final structure
		Output::quietOn();
		if (!File::exists("CONTCAR"))
		{
			_completed = false;
			return;
		}
		data = Read::text("CONTCAR");
		if (Vasp::Structure::isVersion4(data))
			_structures[i] = Vasp::Structure::read(data, false);
		else if (Vasp::Structure::isVersion5(data))
			_structures[i] = Vasp::Structure::read(data, true);
		else
		{
			_completed = false;
			return;
		}
		Output::quietOff();
		
		// Get last energy for current structure
		input.open("OSZICAR");
		if (!input.is_open())
			return;
		
		// Read file and look for energies
		save = false;
		while (input >> temp)
		{
			if (save)
			{
				save = false;
				energy = atof(temp);
			}
			else if (Text::equal(temp, "F="))
				save = true;
		}
		
		// Close OSZICAR
		input.close();
		
		// Save max energy
		if ((i == 1) || (energy > maxEnergy))
		{
			tsDir = dir;
			maxEnergy = energy;
			_structure = _structures[i];
		}
		
		// Return to original directory
		Directory::change(origDir);
	}
	
	// Save transition state energy
	Directory::change(tsDir);
	data = Read::text("OSZICAR");
	for (i = 0; i < data.length(); ++i)
	{
		
		// Skip if line is too short
		if (data.length() < 3)
			continue;
		
		// Found energy line
		if (data[i][1] == "F=")
		{
			_energy += atof(data[i][2].array());
			continue;
		}
		
		// Found start of new electronic convergence loop
		if (data[i][1] == "1")
			_convergence.add();
		
		// Found electronic convergence step
		if (Language::isNumber(data[i][1]))
			_convergence.last() += atof(data[i][2].array());
	}
	
	// Return to original directory
	Directory::change(origDir);
	
	// Check if the last electronic convergence loop converged
	if (_convergence.last().length() == 60)
		_completed = false;
}



/* void Vasp::Job::updateStructure(ISO& iso) const
 *
 * Get structure from vasp job
 */

void Vasp::Job::updateStructure(ISO& iso) const
{
	
	// Not set
	if (!_completed)
		return;
	
	// Assign new basis
	iso.basis(_structure.basis().vectors(), false);
	
	// Loop over atoms and save new positions
	int i, j;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
			iso.atoms()[i][j].fractional(_structure.atoms()[i][j].fractional());
	}
}



/* void Vasp::Job::updateStructures(OList<ISO>& isos) const
 * 
 * Update structures from NEB run
 */

void Vasp::Job::updateStructures(OList<ISO>& isos) const
{
	
	// Not set
	if (!_completed)
		return;
	
	// Loop over structures
	int i, j, k;
	for (i = 0; i < _structures.length(); ++i)
	{
		
		// Assign new basis
		isos[i].basis(_structures[i].basis().vectors(), false);

		// Loop over atoms and save new positions
		for (j = 0; j < isos[i].atoms().length(); ++j)
		{
			for (k = 0; k < isos[i].atoms()[j].length(); ++k)
				isos[i].atoms()[j][k].fractional(_structures[i].atoms()[j][k].fractional());
		}
	}
}



/* ISO Vasp::Structure::read(const Text& content, bool useVasp5)
 *
 * Read structure from file
 */

ISO Vasp::Structure::read(const Text& content, bool useVasp5)
{
	
	// Variable to store result
	ISO iso;
    
    // Output
    Output::newline();
    Output::print("Converting structure in VASP ");
    if (useVasp5)
        Output::print("5");
    else
        Output::print("4");
    Output::print(" file to an internal structure object");
    Output::increase();

	// File is empty
	if (content.length() < 1)
	{
		Output::newline(ERROR);
		Output::print("Trying to read VASP structure from an emtpy file");
		Output::quit();
	}
    
    // Get the elements in the system and comment for vasp 4
	int i, j;
	bool found;
	Words elementNames;
	Words comment;
    if (!useVasp5)
    {
		
		// Loop over words on first line
		found = false;
		for (i = 0; i < content[0].length(); ++i)
		{
			
			// Look for : in word
			for (j = 0; j < content[0][i].length(); ++j)
			{
				if (content[0][i][j] == ':')
				{
					
					// Add current word
					if (j)
					{
						elementNames.add();
						elementNames.last().set(j, content[0][i]);
					}
					
					// Add next word if possible
					if (j < content[0][i].length() - 1)
					{
						comment.add();
						comment.last().set(j+1, content[0][i].length() - j - 1, content[0][i]);
					}
					
					// Break since found
					found = true;
					break;
				}
			}
			
			// Found :
			if (found)
				break;
			
			// Add element
			elementNames += content[0][i];
		}
	
		// Add rest of words as comment
		for (++i; i < content[0].length(); ++i)
			comment += content[0][i];
    }
    
    // Get the elements in the system and comment for vasp 5
    else
    {
	
		// Not enough lines in file
		if (content.length() < 6)
		{
			Output::newline(ERROR);
			Output::print("Not enough lines in VASP structure file");
			Output::quit();
		}
	
		// Get the comment
		for (i = 0; i < content[0].length(); ++i)
			comment += content[0][i];
	
		// Get the elements
		for (i = 0; i < content[5].length(); ++i)
			elementNames += content[5][i];
    }
    
    // Set the comment
	iso.comment(comment);
	
	// No elements
	if (!elementNames.length())
	{
		Output::newline(ERROR);
		Output::print("No elements were found in VASP structure file");
		Output::quit();
	}
	
	// Get the basis multiplier
	if (!content[1].length())
	{
		Output::newline(ERROR);
		Output::print("Could not find the basis multiplier in VASP structure file");
		Output::quit();
	}
	double multiplier = atof(content[1][0].array());
	
	// Make sure basis is present
	for (i = 2; i < 5; ++i)
	{
		if (content[i].length() < 3)
		{
			Output::newline(ERROR);
			Output::print("Could not find basis vectors in VASP structure file");
			Output::quit();
		}
	}
	
	// Get the basis vectors
	Matrix3D basisMat;
	for (i = 2; i < 5; ++i)
	{
		for (j = 0; j < 3; ++j)
			basisMat(i-2, j) = multiplier * atof(content[i][j].array());
	}
    
    // Set the basis vectors
    iso.basis(basisMat, true);
    
	// Make sure there are enough lines
	int curLine = 5;
	if (useVasp5)
		curLine = 6;
	if (content.length() < curLine + 1)
	{
		Output::newline(ERROR);
		Output::print("Not enough lines in VASP structure file");
		Output::quit();
	}
	
	// Number of atom types do not match elements
	if (elementNames.length() != content[curLine].length())
	{
		Output::newline(ERROR);
        if (!useVasp5)
            Output::print("The number of elements defined on the first line does not match sixth line");
        else
            Output::print("The number of elements defined on the sixth line does not match seventh line");
        Output::quit();
	}
	
	// Remove potential information from element name
	for (i = 0; i < elementNames.length(); ++i)
	{
		if (elementNames[i][1] == '_')
			elementNames[i].cutAt(0);
		else if (elementNames[i][2] == '_')
			elementNames[i].cutAt(1);
	}
	
	// Get the elements
	OList<Element> elements (elementNames.length());
	for (i = 0; i < elementNames.length(); ++i)
		elements[i] = Element::find(elementNames[i], true, true);
	
	// Save the number of atoms of each type of element
	List<int> numPerElement (elementNames.length());
	for (i = 0; i < content[curLine].length(); ++i)
		numPerElement[i] = atoi(content[curLine][i].array());
	
	// Print the elements
	Output::newline();
	Output::print("Atom");
	if ((elements.length() != 1) || (numPerElement[0] != 1))
		Output::print("s");
	Output::print(" in the system: ");
	for (i = 0; i < elements.length(); ++i)
	{
		if (!numPerElement[i])
			continue;
		Output::print(numPerElement[i]);
		Output::print(" ");
		Output::print(elements[i].symbol());
		if ((elements.length() > 2) && (i != elements.length() - 1))
			Output::print(",");
		if (i == elements.length() - 2)
			Output::print(" and");
		if (i != elements.length() - 1)
			Output::print(" ");
	}
	
	// Check if using selective dynamics
	++curLine;
	bool useSelective = false;
	if (content[curLine].length())
	{
		if ((content[curLine][0][0] == 's') || (content[curLine][0][0] == 'S'))
		{
			useSelective = true;
			++curLine;
		}
	}
	
	// Check if using fractional or cartesian coordinates
	CoordinateType coordinates = FRACTIONAL;
	if ((content[curLine][0][0] == 'c') || (content[curLine][0][0] == 'C') || \
		(content[curLine][0][0] == 'k') || (content[curLine][0][0] == 'K'))
		coordinates = CARTESIAN;
	
	// Print the coordinate type
	Output::newline();
	Output::print("Using ");
	if (coordinates == FRACTIONAL)
		Output::print(" direct (fractional) ");
	else
		Output::print(" cartesian ");
	Output::print("coordinates");
	
	// Count the total number of atoms
	int totalAtoms = 0;
	for (i = 0; i < numPerElement.length(); ++i)
		totalAtoms += numPerElement[i];
	
	// Print the total number of atoms
	Output::newline();
	Output::print("Adding ");
	Output::print(totalAtoms);
	Output::print(" atom");
	if (totalAtoms != 1)
		Output::print("s");
	Output::increase();
	
	// Get the atomic coordinates
	int curElement = 0;
	int curInElement = 0;
	Atom* curAtom;
	Vector3D position;
	for (i = 0; i < totalAtoms; ++i)
	{
		
		// Go to next line in file
		curLine++;
		
		// Go to next type if needed
        while (curInElement >= numPerElement[curElement])
        {
            curInElement = 0;
            curElement++;
        }
		curInElement++;
		
		// Make sure that file length has not been exceded
		if (curLine >= content.length())
		{
			Output::newline(ERROR);
			Output::print("Expecting ");
			Output::print(totalAtoms);
			Output::print(" atom");
			if (totalAtoms != 1)
				Output::print("s");
			Output::print(" but found ");
			Output::print(i);
			Output::quit();
		}
		
		// Make sure line has content
		if (content[curLine].length() < 3)
		{
			Output::newline(ERROR);
			Output::print("Not enough data on position line for atom ");
			Output::print(i + 1);
			Output::quit();
		}
		
		// Get position of atom
		for (j = 0; j < 3; ++j)
			position[j] = atof(content[curLine][j].array());
		
		// Output
		Output::newline();
		Output::print("Setting atom ");
		Output::print(i + 1);
		Output::print(" (");
		Output::print(elements[curElement].symbol());
		Output::print(" ");
		Output::print(curInElement);
		Output::print(")");
		Output::increase();
		
		// Add atom
		curAtom = iso.addAtom(elements[curElement]);
		if (coordinates == FRACTIONAL)
			curAtom->fractional(position);
		else
			curAtom->cartesian(position * multiplier);
		
		// Get fixed or free
		if (useSelective)
		{
			for (j = 3; j < content[curLine].length(); ++j)
			{
				if ((content[curLine][j][0] == 'f') || (content[curLine][j][0] == 'F'))
					curAtom->fixed(j-3, true);
				if (j == 5)
					break;
			}
		}
		
		// Print properties
		Output::newline();
		Output::print("Fractional: ");
		Output::print(curAtom->fractional(), 8, false);
		Output::newline();
		Output::print("Cartesian: ");
		Output::print(curAtom->cartesian(), 8, false);
		if (useSelective)
		{
			Output::newline();
			Output::print("Fixed: ");
			for (j = 0; j < 3; ++j)
			{
				if (curAtom->fixed()[j])
					Output::print("true");
				else
					Output::print("false");
				if (j != 2)
					Output::print(", ");
			}
		}
		
		// Output
		Output::decrease();
	}
	
	// Output
	Output::decrease();
	Output::decrease();
	
	// Return structure
	return iso;
}



/* void Vasp::Structure::write(const Word& file, const ISO& iso, CoordinateType coordinates, bool useVasp5)
 *
 * Write structure to file
 */

void Vasp::Structure::write(const Word& file, const ISO& iso, CoordinateType coordinates, bool useVasp5)
{
	
	// Precision for printing positions
	int prec = 14;
	
	// Error if any sites are partially occupied
	if (iso.anyPartiallyOccupied())
	{
		Output::newline(ERROR);
		Output::print("Cannot print to VASP format with partially occupancies");
		Output::quit();
	}
	
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
	
	// Print elements if in vasp 4
	int i;
	Output::newline();
    if (!useVasp5)
    {
        for (i = 0; i < iso.atoms().length(); ++i)
        {
			Output::print(iso.atoms()[i][0].element().symbol());
			Output::print(" ");
        }
        if (iso.comment().length())
			Output::print(": ");
    }
    
    // Print comment
	for (i = 0; i < iso.comment().length(); ++i)
	{
		Output::print(iso.comment()[i]);
		Output::print(" ");
	}
	if ((!iso.comment().length()) && (useVasp5))
		Output::print("comment");

	// Print lattice multiplier
	Output::newline();
    Output::print("1.0");

	// Form basis vectors
	int j;
	Output message;
	message.addLines(3);
    for (i = 0; i < 3; ++i)
    {
		message.addLine();
		message.addWords(4);
		message.add(" ");
		for (j = 0; j < 3; ++j)
			message.add(iso.basis().vectors()(i, j), prec);
    }

	// Print basis vectors
	Output::newline();
	Output::print(message, RIGHT);

	// Set the elements if in vasp 5
	message.clear();
    if (useVasp5)
    {
		message.addLine();
		message.addWords(iso.atoms().length());
        for (i = 0; i < iso.atoms().length(); ++i)
			message.add(iso.atoms()[i][0].element().symbol());
    }

	// Set the number of atoms
	message.addLine();
	message.addWords(iso.atoms().length());
	for (i = 0; i < iso.atoms().length(); ++i)
		message.add(iso.atoms()[i].length());
	
	// Print number of atoms
	Output::newline();
	Output::print(message, RIGHT);
	
	// Check if selective dynamics should be used
	bool useSelectiveDynamics = iso.anyFixed();
	if (useSelectiveDynamics)
	{
		Output::newline();
		Output::print("Selective Dynamics");
    }

	// Write cartesian or direct
	Output::newline();
    if (coordinates == CARTESIAN)
    	Output::print("Cartesian");
    else
		Output::print("Direct");
	
	// Write coordinates
	int k;
	message.clear();
	if (iso.atoms().length())
		message.addLines(iso.atoms().last().last().atomNumber() + 1);
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			message.addLine();
			if (useSelectiveDynamics)
				message.addWords(7);
			else
				message.addWords(4);
			message.add(" ");
			for (k = 0; k < 3; ++k)
			{
				if (coordinates == FRACTIONAL)
					message.add(iso.atoms()[i][j].fractional()[k], prec);
				else
					message.add(iso.atoms()[i][j].cartesian()[k], prec);
			}
			if (useSelectiveDynamics)
			{
				for (k = 0; k < 3; ++k)
				{
					if (iso.atoms()[i][j].fixed()[k])
						message.add("F");
					else
						message.add("T");
				}
			}
		}
	}
	Output::newline();
	Output::print(message, RIGHT);
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}



/* bool Vasp::Structure::isVersion4(const Text& content)
 *
 * Test whether file is a VASP version 4 structure
 */

bool Vasp::Structure::isVersion4(const Text& content)
{
	
	// File is too short
    if (content.length() < 7)
        return false;
    
	// Check that multiplier is present
	if (!content[1].length())
		return false;
	if (!Language::isNumber(content[1][0]))
		return false;
	
	// Check that basis is present
	int i, j;
	for (i = 2; i <= 4; ++i)
	{
		if (content[i].length() < 3)
			return false;
		for (j = 0; j < 3; ++j)
		{
			if (!Language::isNumber(content[i][j]))
				return false;
		}
	}

	// Check if selective dynamics is turned on
	int curLine = 6;
	if (content[curLine].length())
	{
		if ((content[curLine][0][0] == 's') || (content[curLine][0][0] == 'S'))
			curLine++;
	}
    
    // File is too short
    if (content.length() <= curLine)
        return false;
    
	// Check for cartesian or direct
	if (!content[curLine].length())
		return false;
	if ((content[curLine][0][0] != 'd') && (content[curLine][0][0] != 'D') && (content[curLine][0][0] != 'c') && \
		(content[curLine][0][0] != 'C') && (content[curLine][0][0] != 'k') && (content[curLine][0][0] != 'K'))
		return false;
	
	// Loop over rest of lines and make sure that coordinates are present
	for (++curLine; curLine < content.length(); ++curLine)
	{
		
		// Break if line is empty
		if (!content[curLine].length())
			break;
		
		// Line is too short
		if (content[curLine].length() < 3)
			return false;
		
		// Does not contain numbers
		if ((!Language::isNumber(content[curLine][0])) || (!Language::isNumber(content[curLine][1])) || \
			(!Language::isNumber(content[curLine][2])))
			return false;
	}
    
    // Return that file type was found
    return true;
}



/* bool Vasp::Structure::isVersion5(const Text& content)
 *
 * Test whether file is a VASP version 5 structure
 */

bool Vasp::Structure::isVersion5(const Text& content)
{
	
	// File is too short
    if (content.length() < 8)
        return false;
    
	// Check that multiplier is present
	if (!content[1].length())
		return false;
	if (!Language::isNumber(content[1][0]))
		return false;
	
	// Check that basis is present
	int i, j;
	for (i = 2; i <= 4; ++i)
	{
		if (content[i].length() < 3)
			return false;
		for (j = 0; j < 3; ++j)
		{
			if (!Language::isNumber(content[i][j]))
				return false;
		}
	}

	// Check if selective dynamics is turned on
	int curLine = 7;
	if (content[curLine].length())
	{
		if ((content[curLine][0][0] == 's') || (content[curLine][0][0] == 'S'))
			curLine++;
	}
    
    // File is too short
    if (content.length() <= curLine)
        return false;
    
	// Check for cartesian or direct
	if (!content[curLine].length())
		return false;
	if ((content[curLine][0][0] != 'd') && (content[curLine][0][0] != 'D') && (content[curLine][0][0] != 'c') && \
		(content[curLine][0][0] != 'C') && (content[curLine][0][0] != 'k') && (content[curLine][0][0] != 'K'))
		return false;
	
	// Loop over rest of lines and make sure that coordinates are present
	for (++curLine; curLine < content.length(); ++curLine)
	{
		
		// Break if line is empty
		if (!content[curLine].length())
			break;
		
		// Line is too short
		if (content[curLine].length() < 3)
			return false;
		
		// Does not contain numbers
		if ((!Language::isNumber(content[curLine][0])) || (!Language::isNumber(content[curLine][1])) || \
			(!Language::isNumber(content[curLine][2])))
			return false;
	}
    
    // Return that file type was found
    return true;
}



/* void Vasp::Settings::defaults(Calculation calculation, Relaxation relaxation, Material material, Accuracy accuracy,
 *		const List<double>* encut)
 *
 * Set defaults
 */

void Vasp::Settings::defaults(Calculation calculation, Relaxation relaxation, Material material, Accuracy accuracy, \
	const List<double>* encut)
{
	
	// Get the max encut value
	double maxEncut = 700;
	if (encut)
	{
		if (encut->length() > 0)
			maxEncut = (*encut)[0];
		for (int i = 1; i < encut->length(); ++i)
		{
			if ((*encut)[i] > maxEncut)
				maxEncut = (*encut)[i];
		}
	}
	
	// Set general defaults
	add("NELM",   "60",      false, false);
	add("NELMIN", "3",       false, false);
	add("ALGO",   "Fast",    false, false);
	add("LCHARG", ".TRUE.",  true,  false);
	add("LWAVE",  ".TRUE.",  true,  false);
	add("EDIFF",  "1E-04",   true,  false);
	add("EDIFFG", "-2E-02",  true,  false);
	add("LREAL",  "Auto",    true,  false);
	add("ISMEAR", "0",       true,  false);
	add("SIGMA",  "0.2",     true,  false);
	
	// Set accuracy
	if (accuracy == VA_HIGH)
	{
		add("PREC",   "Accurate", true, false);
		add("EDIFF",  "1e-05",    true, false);
		add("EDIFFG", "-1E-02",   true, false);
		add("LREAL",  ".FALSE.",  true, false);
	}
	else if (accuracy == VA_LOW)
		add("PREC", "Low", true, false);
	else
		add("PREC", "Normal", true, false);
	
	// Set ENCUT
	if (calculation == VC_RELAX_ALL)
		add("ENCUT", Language::numberToWord((int)(1.3*maxEncut)), true, false);
	else
		remove("ENCUT");
	
	// Set defaults based on calculation, relaxation, material, and accuracy types
	if (calculation == VC_SINGLE)
	{
		add("NSW", "0", true, false);
		add("IBRION", "-1", true, false);
	}
	else if ((calculation == VC_RELAX_ALL) || (calculation == VC_RELAX_POSITIONS))
	{
		add("NSW", "100", true, false);
		if (relaxation == VR_QUASI_NEWTON)
		{
			add("IBRION", "1", true, false);
			add("NELMIN", "4", true, false);
			add("POTIM", "0.5", true, false);
		}
		else if (relaxation == VR_DAMPED_MD)
		{
			add("IBRION", "3", true, false);
			add("SMASS", "0.5", true, false);
			add("POTIM", "0.5", true, false);
		}
		else
		{
			add("IBRION", "2", true, false);
			add("POTIM", "0.5", true, false);
		}
		if (calculation == VC_RELAX_ALL)
			add("ISIF", "3", true, false);
		else
			add("ISIF", "2", true, false);
	}
	else if (calculation == VC_MD)
	{
		add("NELMIN", "4",        true, false);
		add("IBRION", "0",        true, false);
		add("ISYM",   "0",        true, false);
		add("TEBEG",  "300",      true, false);
		add("TEEND",  "300",      true, false);
		add("NSW",    "10000",    true, false);
		add("POTIM",  "2.0",      true, false);
		add("ALGO",   "VeryFast", true, false);
	}
	else if ((calculation == VC_NEB) || (calculation == VC_CINEB))
	{
		add("NSW",    "100",     true, false);
		add("ICHAIN", "0",       true, false);
		add("SPRING", "-5",      true, false);
		add("IMAGES", "5",       true, false);
		add("IOPT",   "2",       true, false);
		add("IBRION", "3",       true, false);
		add("POTIM",  "0",       true, false);
		if (calculation == VC_NEB)
			add("LCLIMB", ".FALSE.", true, false);
		else
			add("LCLIMB", ".TRUE.", true, false);
	}
}



/* bool Vasp::Settings::add(const Word& name, const Word& value, bool replaceExisting, bool forced)
 *
 * Add setting
 */

bool Vasp::Settings::add(const Word& name, const Word& value, bool replaceExisting, bool forced)
{
	
	// Check if a setting
	IncarGroup tempGroup;
	IncarTag tempTag;
	if (!isSetting(name, &tempGroup, &tempTag))
		return false;
	
	// Loop over settings to see if it already exists
	for (int i = 0; i < _settings.length(); ++i)
	{
		if (tempTag == _settings[i].tag)
		{
			if (forced)
			{
				_settings[i].value = value;
				_settings[i].forced = true;
			}
			else if ((!_settings[i].forced) && (replaceExisting))
				_settings[i].value = value;
			return true;
		}
	}
	
	// Save new tag
	_settings.add();
	_settings.last().forced = forced;
	_settings.last().group = tempGroup;
	_settings.last().tag = tempTag;
	_settings.last().name = name.toupper();
	_settings.last().value = value;
	
	// Return that a setting
	return true;
}



/* void Vasp::Settings::remove(const Word& name)
 *
 * Remove setting
 */

void Vasp::Settings::remove(const Word& name)
{
	
	// Check if a setting
	IncarGroup tempGroup;
	IncarTag tempTag;
	if (!isSetting(name, &tempGroup, &tempTag))
		return;
	
	// Loop over settings to find
	for (int i = 0; i < _settings.length(); ++i)
	{
		if (tempTag == _settings[i].tag)
		{
			if (!_settings[i].forced)
				_settings.remove(i);
			return;
		}
	}
}



/* bool Vasp::Settings::isSetting(const Word& name, IncarGroup* group, IncarTag* tag)
 *
 * Return whether a word is a setting
 */

bool Vasp::Settings::isSetting(const Word& name, IncarGroup* group, IncarTag* tag)
{
	
	// Check if a tag
	if (match(name, group, tag, "SYSTEM", IG_GENERAL, IT_SYSTEM))
		return true;
	if (match(name, group, tag, "LCOMPAT", IG_GENERAL, IT_LCOMPAT))
		return true;
	if (match(name, group, tag, "PREC", IG_ELECTRONIC, IT_PREC))
	    return true;
	if (match(name, group, tag, "ENMAX", IG_ELECTRONIC, IT_ENMAX))
	    return true;
	if (match(name, group, tag, "ENAUG", IG_ELECTRONIC, IT_ENAUG))
	    return true;
	if (match(name, group, tag, "EDIFF", IG_ELECTRONIC, IT_EDIFF))
	    return true;
	if (match(name, group, tag, "ALGO", IG_ELECTRONIC, IT_ALGO))
	    return true;
	if (match(name, group, tag, "IALGO", IG_ELECTRONIC, IT_IALGO))
	    return true;
	if (match(name, group, tag, "IWAVPR", IG_ELECTRONIC, IT_IWAVPR))
	    return true;
	if (match(name, group, tag, "NBANDS", IG_ELECTRONIC, IT_NBANDS))
	    return true;
	if (match(name, group, tag, "NELECT", IG_ELECTRONIC, IT_NELECT))
	    return true;
	if (match(name, group, tag, "ISMEAR", IG_ELECTRONIC_SMEARING, IT_ISMEAR))
	    return true;
	if (match(name, group, tag, "SIGMA", IG_ELECTRONIC_SMEARING, IT_SIGMA))
	    return true;
	if (match(name, group, tag, "LREA", IG_ELECTRONIC_PROJECTORS, IT_LREA))
	    return true;
	if (match(name, group, tag, "ROPT", IG_ELECTRONIC_PROJECTORS, IT_ROPT))
	    return true;
	if (match(name, group, tag, "LMAXPAW", IG_ELECTRONIC_PROJECTORS, IT_LMAXPAW))
	    return true;
	if (match(name, group, tag, "LMAXMIX", IG_ELECTRONIC_PROJECTORS, IT_LMAXMIX))
	    return true;
	if (match(name, group, tag, "ISTART", IG_ELECTRONIC_STARTUP, IT_ISTART))
	    return true;
	if (match(name, group, tag, "ICHARG", IG_ELECTRONIC_STARTUP, IT_ICHARG))
	    return true;
	if (match(name, group, tag, "INIWAV", IG_ELECTRONIC_STARTUP, IT_INIWAV))
	    return true;
	if (match(name, group, tag, "ISPIN", IG_ELECTRONIC_SPIN, IT_ISPIN))
	    return true;
	if (match(name, group, tag, "LNONCOLLINEAR", IG_ELECTRONIC_SPIN, IT_LNONCOLLINEAR))
	    return true;
	if (match(name, group, tag, "MAGMOM", IG_ELECTRONIC_SPIN, IT_MAGMOM))
	    return true;
	if (match(name, group, tag, "NUPDOWN", IG_ELECTRONIC_SPIN, IT_NUPDOWN))
	    return true;
	if (match(name, group, tag, "LSORBIT", IG_ELECTRONIC_SPIN, IT_LSORBIT))
	    return true;
	if (match(name, group, tag, "SAXIS", IG_ELECTRONIC_SPIN, IT_SAXIS))
	    return true;
	if (match(name, group, tag, "LSPIRAL", IG_ELECTRONIC_SPIN, IT_LSPIRAL))
	    return true;
	if (match(name, group, tag, "QSPIRAL", IG_ELECTRONIC_SPIN, IT_QSPIRAL))
	    return true;
	if (match(name, group, tag, "LZEROZ", IG_ELECTRONIC_SPIN, IT_LZEROZ))
	    return true;
	if (match(name, group, tag, "GGA", IG_ELECTRONIC_EXCHANGE_CORRELATION, IT_GGA))
	    return true;
	if (match(name, group, tag, "VOSKOWN", IG_ELECTRONIC_EXCHANGE_CORRELATION, IT_VOSKOWN))
	    return true;
	if (match(name, group, tag, "LASPH", IG_ELECTRONIC_EXCHANGE_CORRELATION, IT_LASPH))
	    return true;
	if (match(name, group, tag, "LMETAGGA", IG_ELECTRONIC_EXCHANGE_CORRELATION, IT_LMETAGGA))
	    return true;
	if (match(name, group, tag, "GGA2", IG_ELECTRONIC_EXCHANGE_CORRELATION, IT_GGA2))
	    return true;
	if (match(name, group, tag, "NELM", IG_ELECTRONIC_CONVERGENCE, IT_NELM))
	    return true;
	if (match(name, group, tag, "NELMDL", IG_ELECTRONIC_CONVERGENCE, IT_NELMDL))
	    return true;
	if (match(name, group, tag, "NELMIN", IG_ELECTRONIC_CONVERGENCE, IT_NELMIN))
	    return true;
	if (match(name, group, tag, "ENINI", IG_ELECTRONIC_CONVERGENCE, IT_ENINI))
	    return true;
	if (match(name, group, tag, "LDIAG", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_LDIAG))
	    return true;
	if (match(name, group, tag, "WEIMIN", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_WEIMIN))
	    return true;
	if (match(name, group, tag, "EBREAK", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_EBREAK))
	    return true;
	if (match(name, group, tag, "DEPER", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_DEPER))
	    return true;
	if (match(name, group, tag, "NRMM", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_NRMM))
	    return true;
	if (match(name, group, tag, "TIME", IG_ELECTRONIC_CONVERGENCE_DETAIL, IT_TIME))
	    return true;
	if (match(name, group, tag, "AMIX", IG_ELECTRONIC_MIXER, IT_AMIX))
	    return true;
	if (match(name, group, tag, "BMIX", IG_ELECTRONIC_MIXER, IT_BMIX))
	    return true;
	if (match(name, group, tag, "AMIN", IG_ELECTRONIC_MIXER, IT_AMIN))
	    return true;
	if (match(name, group, tag, "AMIX_MAG", IG_ELECTRONIC_MIXER, IT_AMIX_MAG))
	    return true;
	if (match(name, group, tag, "BMIX_MAG", IG_ELECTRONIC_MIXER, IT_BMIX_MAG))
	    return true;
	if (match(name, group, tag, "IMIX", IG_ELECTRONIC_MIXER_DETAILS, IT_IMIX))
	    return true;
	if (match(name, group, tag, "MAXMIX", IG_ELECTRONIC_MIXER_DETAILS, IT_MAXMIX))
	    return true;
	if (match(name, group, tag, "WC", IG_ELECTRONIC_MIXER_DETAILS, IT_WC))
	    return true;
	if (match(name, group, tag, "INIMIX", IG_ELECTRONIC_MIXER_DETAILS, IT_INIMIX))
	    return true;
	if (match(name, group, tag, "MIXPRE", IG_ELECTRONIC_MIXER_DETAILS, IT_MIXPRE))
	    return true;
	if (match(name, group, tag, "MREMOVE", IG_ELECTRONIC_MIXER_DETAILS, IT_MREMOVE))
	    return true;
	if (match(name, group, tag, "IDIPOL", IG_ELECTRONIC_DIPOLE_CORRECTION, IT_IDIPOL))
	    return true;
	if (match(name, group, tag, "LDIPOL", IG_ELECTRONIC_DIPOLE_CORRECTION, IT_LDIPOL))
	    return true;
	if (match(name, group, tag, "DIPOL", IG_ELECTRONIC_DIPOLE_CORRECTION, IT_DIPOL))
	    return true;
	if (match(name, group, tag, "EFIELD", IG_ELECTRONIC_DIPOLE_CORRECTION, IT_EFIELD))
	    return true;
	if (match(name, group, tag, "NGX", IG_GRIDS, IT_NGX))
	    return true;
	if (match(name, group, tag, "NGY", IG_GRIDS, IT_NGY))
	    return true;
	if (match(name, group, tag, "NGZ", IG_GRIDS, IT_NGZ))
	    return true;
	if (match(name, group, tag, "NGXF", IG_GRIDS, IT_NGXF))
	    return true;
	if (match(name, group, tag, "NGYF", IG_GRIDS, IT_NGYF))
	    return true;
	if (match(name, group, tag, "NGZF", IG_GRIDS, IT_NGZF))
	    return true;
	if (match(name, group, tag, "ADDGRID", IG_GRIDS, IT_ADDGRID))
	    return true;
	if (match(name, group, tag, "NSW", IG_IONIC, IT_NSW))
	    return true;
	if (match(name, group, tag, "IBRION", IG_IONIC, IT_IBRION))
	    return true;
	if (match(name, group, tag, "ISIF", IG_IONIC, IT_ISIF))
	    return true;
	if (match(name, group, tag, "PSTRESS", IG_IONIC, IT_PSTRESS))
	    return true;
	if (match(name, group, tag, "EDIFFG", IG_IONIC, IT_EDIFFG))
	    return true;
	if (match(name, group, tag, "NFREE", IG_IONIC, IT_NFREE))
	    return true;
	if (match(name, group, tag, "POTIM", IG_IONIC, IT_POTIM))
	    return true;
	if (match(name, group, tag, "SMASS", IG_IONIC, IT_SMASS))
	    return true;
	if (match(name, group, tag, "TEBEG", IG_IONIC_MD, IT_TEBEG))
	    return true;
	if (match(name, group, tag, "TEEND", IG_IONIC_MD, IT_TEEND))
	    return true;
	if (match(name, group, tag, "NBLOCK", IG_IONIC_MD, IT_NBLOCK))
	    return true;
	if (match(name, group, tag, "KBLOCK", IG_IONIC_MD, IT_KBLOCK))
	    return true;
	if (match(name, group, tag, "NPACO", IG_IONIC_MD, IT_NPACO))
	    return true;
	if (match(name, group, tag, "APACO", IG_IONIC_MD, IT_APACO))
	    return true;
	if (match(name, group, tag, "ISYM", IG_SYMMETRY, IT_ISYM))
	    return true;
	if (match(name, group, tag, "SYMPREC", IG_SYMMETRY, IT_SYMPREC))
	    return true;
	if (match(name, group, tag, "LORBIT", IG_DOS, IT_LORBIT))
	    return true;
	if (match(name, group, tag, "RWIGS", IG_DOS, IT_RWIGS))
	    return true;
	if (match(name, group, tag, "NEDOS", IG_DOS, IT_NEDOS))
	    return true;
	if (match(name, group, tag, "EMIN", IG_DOS, IT_EMIN))
	    return true;
	if (match(name, group, tag, "EMAX", IG_DOS, IT_EMAX))
	    return true;
	if (match(name, group, tag, "NWRITE", IG_WRITING, IT_NWRITE))
	    return true;
	if (match(name, group, tag, "LWAVE", IG_WRITING, IT_LWAVE))
	    return true;
	if (match(name, group, tag, "LCHARG", IG_WRITING, IT_LCHARG))
	    return true;
	if (match(name, group, tag, "LPARD", IG_WRITING, IT_LPARD))
	    return true;
	if (match(name, group, tag, "LVTOT", IG_WRITING, IT_LVTOT))
	    return true;
	if (match(name, group, tag, "LELF", IG_WRITING, IT_LELF))
	    return true;
	if (match(name, group, tag, "LOPTICS", IG_WRITING, IT_LOPTICS))
	    return true;
	if (match(name, group, tag, "STM", IG_WRITING, IT_STM))
	    return true;
	if (match(name, group, tag, "NPAR", IG_PERFORMANCE, IT_NPAR))
	    return true;
	if (match(name, group, tag, "NSIM", IG_PERFORMANCE, IT_NSIM))
	    return true;
	if (match(name, group, tag, "NBLK", IG_PERFORMANCE, IT_NBLK))
	    return true;
	if (match(name, group, tag, "LPLANE", IG_PERFORMANCE, IT_LPLANE))
	    return true;
	if (match(name, group, tag, "LCRITICAL_MEM", IG_PERFORMANCE, IT_LCRITICAL_MEM))
	    return true;
	if (match(name, group, tag, "LSCALAPACK", IG_PERFORMANCE, IT_LSCALAPACK))
	    return true;
	if (match(name, group, tag, "LSCALU", IG_PERFORMANCE, IT_LSCALU))
	    return true;
	if (match(name, group, tag, "LASYNC", IG_PERFORMANCE, IT_LASYNC))
	    return true;
	if (match(name, group, tag, "IDIOT", IG_MISCELLANEOUS, IT_IDIOT))
	    return true;
	if (match(name, group, tag, "LMUSIC", IG_MISCELLANEOUS, IT_LMUSIC))
	    return true;
	if (match(name, group, tag, "POMASS", IG_MISCELLANEOUS, IT_POMASS))
	    return true;
	if (match(name, group, tag, "LCORR", IG_MISCELLANEOUS, IT_LCORR))
	    return true;
	if (match(name, group, tag, "LREAL", IG_OTHER, IT_LREAL))
	    return true;
	if (match(name, group, tag, "LREAL_COMPAT", IG_OTHER, IT_LREAL_COMPAT))
	    return true;
	if (match(name, group, tag, "GGA_COMPAT", IG_OTHER, IT_GGA_COMPAT))
	    return true;
	if (match(name, group, tag, "LBERRY", IG_OTHER, IT_LBERRY))
	    return true;
	if (match(name, group, tag, "ICORELEVEL", IG_OTHER, IT_ICORELEVEL))
	    return true;
	if (match(name, group, tag, "LDAU", IG_OTHER, IT_LDAU))
	    return true;
	if (match(name, group, tag, "I_CONSTRAINED_M", IG_OTHER, IT_I_CONSTRAINED_M))
	    return true;
	if (match(name, group, tag, "ICHAIN", IG_VARIATIONAL_TRANSITION_STATE, IT_ICHAIN))
	    return true;
	if (match(name, group, tag, "IMAGES", IG_VARIATIONAL_TRANSITION_STATE, IT_IMAGES))
	    return true;
	if (match(name, group, tag, "SPRING", IG_VARIATIONAL_TRANSITION_STATE, IT_SPRING))
	    return true;
	if (match(name, group, tag, "LCLIMB", IG_VARIATIONAL_TRANSITION_STATE, IT_LCLIMB))
	    return true;
	if (match(name, group, tag, "LTANGENTOLD", IG_VARIATIONAL_TRANSITION_STATE, IT_LTANGENTOLD))
	    return true;
	if (match(name, group, tag, "LDNEB", IG_VARIATIONAL_TRANSITION_STATE, IT_LDNEB))
	    return true;
	if (match(name, group, tag, "NEBCELL", IG_VARIATIONAL_TRANSITION_STATE, IT_NEBCELL))
	    return true;
	if (match(name, group, tag, "IOPT", IG_VARIATIONAL_TRANSITION_STATE, IT_IOPT))
		return true;
    
    // Return that the value was not recognized as a tag
    return false;
}



/* bool Vasp::Settings::match(const Word& name, IncarGroup* group, IncarTag* tag, const char* comp,
 *		IncarGroup compGroup, IncarTag compTag)
 *
 * Compare name and setting
 */

bool Vasp::Settings::match(const Word& name, IncarGroup* group, IncarTag* tag, const char* comp, \
	IncarGroup compGroup, IncarTag compTag)
{
	if (name.equal(comp, false))
	{
		if (group)
			*group = compGroup;
		if (tag)
			*tag = compTag;
		return true;
	}
	return false;
}



/* void Vasp::Settings::read(const Text& content)
 *
 * Read vasp settings file
 */

void Vasp::Settings::read(const Text& content)
{
	
	// Clear space
	clear();
	
	// Copy text and remove = signs
	Text copy(content);
	copy.split('=');
	
	// Loop over lines and add settings
	int i, j;
	Word curValue;
	for (i = 0; i < copy.length(); ++i)
	{
		
		// Skip if line is too short
		if (copy[i].length() < 2)
			continue;
		
		// Skip if a comment
		if (Language::isComment(copy[i][0]))
			continue;
		
		// Save value
		curValue.clear();
		for (j = 1; j < copy[i].length(); ++j)
		{
			if (Language::isComment(copy[i][j]))
				break;
			if (j > 1)
				curValue += ' ';
			curValue += copy[i][j];
		}
		
		// Save setting if a value was found
		if (curValue.length() > 0)
			add(copy[i][0], curValue, true, false);
	}
}



/* bool Vasp::Settings::isFormat(const Text& content)
 *
 * Return whether file is in vasp settings format
 */

bool Vasp::Settings::isFormat(const Text& content)
{
	
	// Copy text and remove = signs
	Text copy(content);
	copy.split('=');
	
	// Loop over lines and count the number that contain a tag
	int tags = 0;
	int total = 0;
	for (int i = 0; i < copy.length(); ++i)
	{
		
		// Skip if blank or a comment
		if (!copy[i].length())
			continue;
		if (Language::isComment(copy[i][0]))
			continue;
		
		// Check if a tag
		++total;
		if (isSetting(copy[i][0]))
			++tags;
	}
	
	// Not a settings file if empty
	if (!total)
		return false;
	
	// Check whether a settings file
	if ((double)tags / total >= 0.5)
		return true;
	return false; 
}



/* void Vasp::Settings::write(const Word& file) const
 *
 * Write vasp settings file
 */

void Vasp::Settings::write(const Word& file) const
{
	
	// Open file for writing if needed
	int origStream = Output::streamID();
	if (file != "stdout")
		Output::setStream(Output::addStream(file));
	
	// Print settings
	write(IG_GENERAL, "# General");
	write(IG_ELECTRONIC, "# Electronic");
	write(IG_ELECTRONIC_SMEARING, "# Electronic smearing");
	write(IG_ELECTRONIC_PROJECTORS, "# Electronic projectors");
	write(IG_ELECTRONIC_STARTUP, "# Electronic startup");
	write(IG_ELECTRONIC_SPIN, "# Electronic spin");
	write(IG_ELECTRONIC_EXCHANGE_CORRELATION, "# Electronic exchange-correlation");
	write(IG_ELECTRONIC_CONVERGENCE, "# Electronic convergence");
	write(IG_ELECTRONIC_CONVERGENCE_DETAIL, "# Electronic convergence details");
	write(IG_ELECTRONIC_MIXER, "# Electronic mixer");
	write(IG_ELECTRONIC_MIXER_DETAILS, "# Electronic mixer details");
	write(IG_ELECTRONIC_DIPOLE_CORRECTION, "# Electronic dipole correction");
	write(IG_GRIDS, "# Grids");
	write(IG_IONIC, "# Ionic");
	write(IG_IONIC_MD, "# Ionic MD");
	write(IG_SYMMETRY, "# Symmetry");
	write(IG_DOS, "# DOS");
	write(IG_WRITING, "# Writing");
	write(IG_PERFORMANCE, "# Performance");
	write(IG_MISCELLANEOUS, "# Miscellaneous");
	write(IG_OTHER, "# Other");
	write(IG_VARIATIONAL_TRANSITION_STATE, "# Variational transition state");
	
	// Reset output
	if (file != "stdout")
		Output::removeStream(Output::streamID());
	Output::setStream(origStream);
}



/* void Vasp::Settings::write(IncarGroup group, const char* header) const
 *
 * Write vasp settings
 */

void Vasp::Settings::write(IncarGroup group, const char* header) const
{
	
	// Make list of settings to print
	int i;
	List<Setting*> toPrint;
	for (i = 0; i < _settings.length(); ++i)
	{
		if (group == _settings[i].group)
			toPrint += &_settings[i];
	}
	
	// Return if no settings are in current group
	if (!toPrint.length())
		return;
	
	// Set the method
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print the header
	Output::newline();
	Output::print(header);
	
	// Loop over settings and write
	for (i = 0; i < toPrint.length(); ++i)
	{
		Output::newline();
		Output::print("    ");
		Output::print(toPrint[i]->name);
		Output::print(" = ");
		Output::print(toPrint[i]->value);
	}
	
	// Reset the method
	Output::method(origMethod);
}



/* void Vasp::Potential::write(const Word& file, const OList<Element>& elements, const OList<Word>& files,
 *		const ISO& iso, List<double>* encut)
 *
 * Write potential file
 */

void Vasp::Potential::write(const Word& file, const OList<Element>& elements, const OList<Word>& files, \
	const ISO& iso, List<double>* encut)
{	
	
	// Only run on root processor
	int i, j;
	if (Multi::rank() == 0)
	{
		
		// Open file for writing
	    ofstream outfile (file.array());
	    if (!outfile.is_open())
	    {
	        Output::newline(ERROR);
	        Output::print("Could not open ");
	        Output::print(file);
	        Output::print(" for writing");
	        Output::quit();
	    }

		// Loop over elements in the structure
		bool found;
		int lineLength = 500;
		char line[lineLength];
		ifstream infile;
		for (i = 0; i < iso.atoms().length(); ++i)
		{

			// Loop over elements that were supplied
			found = false;
			for (j = 0; j < elements.length(); ++j)
			{

				// Found element
				if (iso.atoms()[i][0].element() == elements[j])
				{

					// Open file
					infile.open(files[j].array());
					if (!infile.is_open())
					{
						Output::newline(ERROR);
						Output::print("Error opening file: ");
						Output::print(files[j]);
						Output::quit();
					}

					// Copy contents
					while (infile.getline(line, lineLength))
						outfile << line << endl;

					// Close file
					infile.close();

					// Save that found
					found = true;
					break;
				}
			}

			// Element was not found
			if (!found)
			{
				Output::newline(ERROR);
				Output::print("VASP potential was not supplied for ");
				Output::print(iso.atoms()[i][0].element().symbol());
				Output::quit();
			}
		}
		
		// Close output file
		outfile.close();
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
	
	// If encut was not passed then return
	if (!encut)
		return;
	encut->clear();
	
	// Read potential file that was written
	Text content = Read::text(file);
	content.split(';');
	content.split('=');
	
	// Find encut values
	for (i = 0; i < content.length(); ++i)
	{
		for (j = 0; j < content[i].length(); ++j)
		{
			if (content[i][j].equal("ENMAX", false))
			{
				if (content[i].length() > j + 1)
					*encut += atof(content[i][j+1].array());
			}
		}
	}
}
