/* espresso.cpp -- Quantum Espresso functionality
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "multi.h"
#include "espresso.h"
#include "structureIO.h"
#include "language.h"
#include "output.h"
#include <cstdlib>
#include <fstream>
using namespace std;



// Name of directory to run in
Word Espresso::Job::_dirName = "QERun";

// Name of files to use
Word Espresso::Job::_inputFileName = "input.qe";
Word Espresso::Job::_outputFileName = "output.qe";



/* Espresso::Accuracy Espresso::accuracy(const Word& word)
 *
 * Return the accuracy level as an enumerated value
 */

Espresso::Accuracy Espresso::accuracy(const Word& word)
{
	if (word.equal("high", false, 4))
		return EA_HIGH;
	if ((word.equal("medium", false, 3)) || (word.equal("normal", false, 4)))
		return EA_NORMAL;
	if (word.equal("low", false, 3))
		return EA_LOW;
	return EA_UNKNOWN;
}



/* Word Espresso::accuracy(Accuracy accuracy)
 *
 * Return the accuracy level as a word
 */

Word Espresso::accuracy(Accuracy accuracy)
{
	switch (accuracy)
	{
		case EA_HIGH:
			return Word("High");
		case EA_NORMAL:
			return Word("Normal");
		case EA_LOW:
			return Word("Low");
		default:
			return Word("Unknown");
	}
}



/* bool Espresso::Job::submit(const ISO& iso, bool quitIfError)
 *
 * Submit quantum espresso job
 */

bool Espresso::Job::submit(const ISO& iso, bool quitIfError)
{

	// Clear storage
	_completed = false;
	_energy.clear();
	_forces.clear();
	_structure.clear();
	
	// Switch to run directory
	Word origDir;
	Directory::current(origDir);
	Word runDir;
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
	Directory::change(runDir);
	
	// Create kpoints object
	KPoints kpoints;
	int kppra = 0;
	if (_material == EM_INSULATOR)
	{
		if (_accuracy == EA_HIGH)
			kppra = KPoints::kppraInsulatorHigh();
		else if (_accuracy == EA_LOW)
			kppra = KPoints::kppraInsulatorLow();
		else
			kppra = KPoints::kppraInsulatorNormal();
	}
	else
	{
		if (_accuracy == EA_HIGH)
			kppra = KPoints::kppraMetalHigh();
		else if (_accuracy == EA_LOW)
			kppra = KPoints::kppraMetalLow();
		else
			kppra = KPoints::kppraMetalNormal();
	}
	kpoints.mesh(iso.numAtoms(), iso.basis().lengths(), kppra);

	// Write files
	Espresso::Files::write(_inputFileName, iso, FRACTIONAL, this, &kpoints);

	// Run job
	Multi::external(_executable, (Word("-inp ")+_inputFileName).array(), _outputFileName.array());
	
	// Get job results
	gatherResults(&iso.basis());
	
	// Error if job did not complete
	if (!_completed)
	{
		if (quitIfError)
			Output::newline(ERROR);
		else
			Output::newline(WARNING);
		Output::print("Quantum Espresso did not exit properly");
		if (quitIfError)
			Output::quit();
	}
	
	// Switch back to run directory
	Directory::change(origDir);
	
	// Return runtime directory
	if (!_saveFiles)
		Directory::remove(runDir);
	
	// Return whether run was successful
	return _completed;
}



/* void Espresso::Job::gatherResults(const Basis* basis)
 *
 * Get results from Quantum Espresso job
 */

void Espresso::Job::gatherResults(const Basis* basis)
{
	
	// Open file with results
	if (!File::exists(_outputFileName))
		return;
	Text results = Read::text(_outputFileName);
	
	// Loop over lines in file
	int i, j;
	Text lastStructure;
	for (i = 0; i < results.length(); ++i)
	{
		
		// Check for an energy line
		if (results[i].length() >= 5)
		{
			if ((results[i][0] == "total") && (results[i][1] == "energy") && (results[i][2] == "="))
				_energy += atof(results[i][3].array());
			else if ((results[i][1] == "total") && (results[i][2] == "energy") && (results[i][3] == "="))
				_energy += atof(results[i][4].array());
		}
		
		// Check for forces
		if (results[i].length() > 4)
		{
			if ((results[i][0] == "Forces") && (results[i][1] == "acting") && \
				(results[i][2] == "on") && (results[i][3] == "atoms"))
			{
				
				// Add new set of forces
				_forces.add();
				
				// Loop until end of forces
				for (++i; i < results.length(); ++i)
				{
					
					// Skip if line is empty
					if (!results[i].length())
						continue;
					
					// Save force if line starts with "atom"
					if ((results[i][0].equal("atom", false)) && (results[i].length() >= 9))
					{
						_forces.last().add();
						for (j = 0; j < 3; ++j)
							_forces.last().last()[j] = Num<double>::forceFromAU(atof(results[i][j+6].array()));
					}
					
					// Break if line contains something else
					else
					{
						--i;
						break;
					}
				}
			}
		}
		
		// Check for start of structure
		if (results[i].length() >= 1)
		{
			if ((results[i][0] == "CELL_PARAMETERS") || (results[i][0] == "ATOMIC_POSITIONS"))
			{
				
				// Clear space
				lastStructure.clear();
				
				// Save basis lines through ATOMIC_POSITIONS
				if (results[i][0] == "CELL_PARAMETERS")
				{
					for (j = 0; (i < results.length()) && (j < 6); ++i, ++j)
					{
						lastStructure.addLine();
						lastStructure.addWords(results[i]);
					}
				}
				else
				{
					lastStructure.addLine();
					lastStructure.addWords(results[i]);
					++i;
				}
				
				// Save atom lines
				for (; i < results.length(); ++i)
				{
					
					// Break if line is empty
					if (!results[i].length())
						break;
					
					// Break if line does not start with an element
					if (!Element::isElement(results[i][0], false))
					{
						--i;
						break;
					}
					
					// Save line
					lastStructure.addLine();
					lastStructure.addWords(results[i]);
				}
			}
		}
		
		// Check for job done
		if (results[i].length() >= 2)
		{
			if ((results[i][0] == "JOB") && (results[i][1].contains("DONE")))
				_completed = true;
		}
	}

	// Save structure if job completed
	if (lastStructure.length())
		_structure = Files::read(lastStructure, false, basis);
}



/* void Espresso::Job::updateStructure(ISO& iso) const
 *
 * Get structure from quantum espresso job
 */

void Espresso::Job::updateStructure(ISO& iso) const
{
	
	// Not set
	if ((!_completed) || (!_structure.numAtoms()))
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



/* ISO Espresso::Files::read(Text content, bool readAsInput, const Basis* basis)
 *
 * Read quantum espresso file
 */

ISO Espresso::Files::read(Text content, bool readAsInput, const Basis* basis)
{
	
	// Output
	if (readAsInput)
	{
		Output::newline();
		Output::print("Converting structure in Quantum Espresso format to an internal structure object");
		Output::increase();
	}
	
	// Remove = and () signs from data
	content.split('=');
	
	// Variable to store result
	ISO iso;
	if (basis)
		iso.basis(*basis, false);
	
	// Variables to store data
	int ibrav = 0;
	bool fixed[3];
	double alat = 0;
	double celldm[6] = {-100, -100, -100, -100, -100, -100};
	double metrics[6] = {-100, -100, -100, -100, -100, -100};
	Vector3D position;
	Matrix3D cellParams = 0.0;
	ValueType cellParamType = VT_ALAT;
	ValueType posType = VT_ALAT;
	
	// Loop over lines in the file to get basis
	int i, j, k;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Loop over words on line
		for (j = 0; j < content[i].length(); ++j)
		{
			
			// Break if a comment
			if (isComment(content[i][j]))
				break;
			
			// Found ibrav
			if ((content[i][j].equal("ibrav", false)) && (content[i].length() > j+1))
				ibrav = atoi(content[i][++j].array());
			
			// Found a celldm value
			else if ((content[i][j].equal("celldm(1)", false)) && (content[i].length() > j+1))
				celldm[0] = Num<double>::distanceFromAU(atof(content[i][++j].array()));
			else if ((content[i][j].equal("celldm(2)", false)) && (content[i].length() > j+1))
				celldm[1] = Num<double>::distanceFromAU(atof(content[i][++j].array()));
			else if ((content[i][j].equal("celldm(3)", false)) && (content[i].length() > j+1))
				celldm[2] = Num<double>::distanceFromAU(atof(content[i][++j].array()));
			else if ((content[i][j].equal("celldm(4)", false)) && (content[i].length() > j+1))
				celldm[3] = atof(content[i][++j].array());
			else if ((content[i][j].equal("celldm(5)", false)) && (content[i].length() > j+1))
				celldm[4] = atof(content[i][++j].array());
			else if ((content[i][j].equal("celldm(6)", false)) && (content[i].length() > j+1))
				celldm[5] = atof(content[i][++j].array());
			
			// Found A, B, C, etc.
			else if ((content[i][j].equal("A", false)) && (content[i].length() > j+1))
				metrics[0] = atof(content[i][++j].array());
			else if ((content[i][j].equal("B", false)) && (content[i].length() > j+1))
				metrics[1] = atof(content[i][++j].array());
			else if ((content[i][j].equal("C", false)) && (content[i].length() > j+1))
				metrics[2] = atof(content[i][++j].array());
			else if ((content[i][j].equal("cosAB", false)) && (content[i].length() > j+1))
				metrics[3] = atof(content[i][++j].array());
			else if ((content[i][j].equal("cosAC", false)) && (content[i].length() > j+1))
				metrics[4] = atof(content[i][++j].array());
			else if ((content[i][j].equal("cosBC", false)) && (content[i].length() > j+1))
				metrics[5] = atof(content[i][++j].array());
		
			// Found lattice parameters
			else if (content[i][j].equal("CELL_PARAMETERS"))
			{
				
				// Check for value type
				for (k = j; k < content[i].length(); ++k)
				{
					if (isComment(content[i][k]))
						break;
					if (content[i][k].contains("alat", false))
						cellParamType = VT_ALAT;
					else if (content[i][k].contains("bohr", false))
						cellParamType = VT_BOHR;
					else if (content[i][k].contains("ang", false))
						cellParamType = VT_ANGSTROM;
				}
				
				// Check for alat value if needed
				if (cellParamType == VT_ALAT)
				{
					for (; j < content[i].length(); ++j)
					{
						content[i][j].strip(')');
						if (isComment(content[i][j]))
							break;
						if (Language::isNumber(content[i][j]))
							alat = Num<double>::distanceFromAU(atof(content[i][j].array()));
					}
				}
				
				// Loop over next three lines to get cell parameters
				for (++i, j = 0; ((i < content.length()) && (j < 3)); ++i, ++j)
				{
					if (content[i].length() < 3)
						continue;
					for (k = 0; k < 3; ++k)
						cellParams(j, k) = atof(content[i][k].array());
				}
				
				// Break since parameters were found
				break;
			}
		}
	}
	
	// ibrav was set
	if (ibrav != 0)
	{
		
		// Variables
		double a = 0;
		double b = 0;
		double c = 0;
		double cosAB = 0;
		double cosAC = 0;
		double cosBC = 0;
		
		// Using celldm
		if (celldm[0] != -100)
		{
			
			// Metrics were also set
			if (metrics[0] != -100)
			{
				Output::newline(ERROR);
				Output::print("Both celldm(1) and A have been set");
				Output::quit();
			}
			
			// Save values
			a = celldm[0];
			b = celldm[1] * celldm[0];
			c = celldm[2] * celldm[0];
			cosAB = celldm[3];
			cosAC = celldm[4];
			cosBC = celldm[5];
			if ((ibrav == 5) || (ibrav == -5))
				cosBC = celldm[3];
		}
		
		// Using metrics
		else if (metrics[0] != -100)
		{
			a = metrics[0];
			b = metrics[1];
			c = metrics[2];
			cosAB = metrics[3];
			cosAC = metrics[4];
			cosBC = metrics[5];
		}
		
		// Nothing was set
		else
		{
			Output::newline(ERROR);
			Output::print("ibrav = ");
			Output::print(ibrav);
			Output::print(" but neither celldm(1) or A have been set");
			Output::quit();
		}
		
		// Simple cubic
		if (ibrav == 1)
			setBasis(iso, a, a, a, 1, 1, 1, readAsInput);
		
		// FCC
		else if (ibrav == 2)
			setBasis(iso, -a/2, 0, a/2, 0, a/2, a/2, -a/2, a/2, 0, readAsInput);
				
		// BCC
		else if (ibrav == 3)
			setBasis(iso, a/2, a/2, a/2, -a/2, a/2, a/2, -a/2, -a/2, a/2, readAsInput);
		
		// Hexagonal or primitive trigonal
		else if (ibrav == 4)
			setBasis(iso, a, 0, 0, -a/2, a*sqrt(3.0)/2, 0, 0, 0, c, readAsInput);
		
		// Trigonal R-centered with 3fold around c(ibrav=5) or <111>(ibrav=-5)
		else if ((ibrav == 5) || (ibrav == -5))
		{
			double tx = a*sqrt((1-cosBC)/2);
			double ty = a*sqrt((1-cosBC)/6);
			double tz = a*sqrt((1+2*cosBC)/3);
			if (ibrav == 5)
				setBasis(iso, tx, -ty, tz, 0, 2*ty, tz, -tx, -ty, tz, readAsInput);
			else
			{
				double u = (tz - 2*sqrt(2.0)*ty)/sqrt(3.0);
				double v = (tz + sqrt(2.0)*ty)/sqrt(3.0);
				setBasis(iso, u, v, v, v, u, v, v, v, u, readAsInput);
			}
		}
		
		// Primitive tetragonal
		else if (ibrav == 6)
			setBasis(iso, a, a, c, 1, 1, 1, readAsInput);
		
		// Body centered tetragonal
		else if (ibrav == 7)
			setBasis(iso, a/2, -a/2, c/2, a/2, a/2, c/2, -a/2, -a/2, c/2, readAsInput);
		
		// Primitive orthorhombic
		else if (ibrav == 8)
			setBasis(iso, a, b, c, 1, 1, 1, readAsInput);
		
		// Base centered orthorhombic
		else if (ibrav == 9)
			setBasis(iso, a/2, b/2, 0, -a/2, b/2, 0, 0, 0, c, readAsInput);
		
		// Basis centered orthorhombic
		else if (ibrav == -9)
			setBasis(iso, a/2, -b/2, 0, a/2, -b/2, 0, 0, 0, c, readAsInput);
		
		// Face centered orthorhombic
		else if (ibrav == 10)
			setBasis(iso, a/2, 0, c/2, a/2, b/2, 0, 0, b/2, c/2, readAsInput);
		
		// Body centering orthorhombic
		else if (ibrav == 11)
			setBasis(iso, a/2, b/2, c/2, -a/2, b/2, c/2, -a/2, -b/2, c/2, readAsInput);
		
		// Primitive monoclinic with unique axis c
		else if (ibrav == 12)
			setBasis(iso, a, 0, 0, b*cosAB, b*sin(acos(cosAB)), 0, 0, 0, c, readAsInput);
		
		// Primitive monoclinic with unique axis b
		else if (ibrav == -12)
			setBasis(iso, a, 0, 0, 0, b, 0, c*sin(acos(cosBC)), 0, c*cosBC, readAsInput);
		
		// Base centered monoclinic
		else if (ibrav == 13)
			setBasis(iso, a/2, 0, -c/2, b*cosAB, b*sin(acos(cosAB)), 0, a/2, 0, c/2, readAsInput);
		
		// Triclinic
		else if (ibrav == 14)
			setBasis(iso, a, b, c, cosAB, cosAC, cosBC, readAsInput);
		
		// Anything else
		else
		{
			Output::newline(ERROR);
			Output::print("Unknown ibrav value (");
			Output::print(ibrav);
			Output::print(")");
			Output::quit();
		}
		
		// Set alat
		alat = iso.basis().lengths()[0];
	}
	
	// ibrav was not set or was set to zero
	else
	{
		
		// Basis was not set
		if ((cellParams == Matrix3D(0.0)) && (basis == 0))
		{
			Output::newline(ERROR);
			Output::print("Cell parameters were expected but not set properly");
			Output::quit();
		}
		
		// Set alat
		if ((alat == 0) && (cellParamType == VT_ALAT))
		{
			
			// Alat set by cell dim
			if (celldm[0] != -100)
				alat = celldm[0];
			
			// Alat set by A
			else if (metrics[0] != -100)
				alat = metrics[0];
			
			// Alat was not set
			else if (basis == 0)
			{
				Output::newline(ERROR);
				Output::print("alat was not set properly");
				Output::quit();
			}
		}
		
		// Adjust basis vectors to alat
		if (cellParamType == VT_ALAT)
			cellParams *= alat;
		
		// Turn basis vectors from Bohr to Angstrom
		else if (cellParamType == VT_BOHR)
		{
			for (i = 0; i < 3; ++i)
			{
				for (j = 0; j < 3; ++j)
					cellParams(i, j) = Num<double>::distanceFromAU(cellParams(i, j));
			}
		}
		
		// Save basis
		if (cellParams != Matrix3D(0.0))
			iso.basis(cellParams, readAsInput);
	}
	
	// Output
	if (readAsInput)
	{
		Output::newline();
		Output::print("Adding atoms");
		Output::increase();
	}
	
	// Find start of atoms definitions
	Atom* curAtom;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is empty
		if (content[i].length() == 0)
			continue;
		
		// Found ATOMIC_POSITIONS
		if (content[i][0].equal("ATOMIC_POSITIONS", false))
		{
			
			// Check for position type
			for (j = 1; j < content[i].length(); ++j)
			{
				if (isComment(content[i][j]))
					break;
				if (content[i][j].contains("alat", false))
					posType = VT_ALAT;
				else if (content[i][j].contains("bohr", false))
					posType = VT_BOHR;
				else if (content[i][j].contains("ang", false))
					posType = VT_ANGSTROM;
				else if (content[i][j].contains("cryst", false))
					posType = VT_CRYSTAL;
			}
			
			// Loop over lines
			for (++i; i < content.length(); ++i)
			{
				
				// Break if line is too short
				if (content[i].length() < 4)
					break;
				
				// Skip if line is commented
				if (isComment(content[i][0]))
					continue;
				
				// Break if line does not start with an element
				if (!Element::isElement(content[i][0]))
					break;
				
				// Break if line does not contain positions
				if ((!Language::isNumber(content[i][1])) || (!Language::isNumber(content[i][2])) || \
					(!Language::isNumber(content[i][3])))
					break;
				
				// Save position
				for (j = 1; j < 4; ++j)
					position[j-1] = atof(content[i][j].array());
				if (posType == VT_ALAT)
					position *= alat;
				else if (posType == VT_BOHR)
				{
					for (j = 0; j < 3; ++j)
						position[j] = Num<double>::distanceFromAU(position[j]);
				}
				
				// Save whether fixed
				fixed[0] = fixed[1] = fixed[2] = false;
				for (j = 4; (j < content[i].length()) && (j < 7); ++j)
				{
					if (isComment(content[i][j]))
						break;
					if (content[i][j] == "0")
						fixed[j-4] = true;
				}
				
				// Save new atom
				curAtom = iso.addAtom(content[i][0]);
				if (posType == VT_CRYSTAL)
					curAtom->fractional(position);
				else
					curAtom->cartesian(position);
				curAtom->fixed(fixed);
				
				// Output
				if (readAsInput)
				{
					Output::newline();
					Output::print("Setting atom ");
					Output::print(curAtom->atomNumber() + 1);
					Output::increase();
					Output::newline();
					Output::print("Fractional: ");
					Output::print(curAtom->fractional(), 8, false);
					Output::newline();
					Output::print("Cartesian: ");
					Output::print(curAtom->cartesian(), 8, false);
					if ((curAtom->fixed()[0]) || (curAtom->fixed()[1]) || (curAtom->fixed()[2]))
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
					Output::decrease();
				}
			}
			
			// Break since positions were found
			break;
		}
	}
	
	// Output
	if (readAsInput)
		Output::decrease();
	
	// Output
	if (readAsInput)
		Output::decrease();
	
	// Return structure
	return iso;
}



/* void Espresso::Files::write(const Word& file, const ISO& iso, CoordinateType coordinates, Job* job, KPoints* kpoints)
 *
 * Write quantum espresso file
 */

void Espresso::Files::write(const Word& file, const ISO& iso, CoordinateType coordinates, Job* job, KPoints* kpoints)
{
	
	// Precision for printing positions
	int prec = 14;
	char indent[] = "   ";
	
	// Error if any sites are partially occupied
	if (iso.anyPartiallyOccupied())
	{
		Output::newline(ERROR);
		Output::print("Cannot print to Quantum Espresso format with partially occupancies");
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
	
	// Print control
	Output message;
	List<PrintAlign> align;
	if (job)
	{
		
		// Header
		Output::newline();
		Output::print("&CONTROL");
		
		// Calculation type
		message.addLine();
		message.add(indent);
		message.add("calculation");
		message.add("=");
		if (job->calculation() == EC_RELAX_ALL)
			message.add("\'vc-relax\'");
		else if (job->calculation() == EC_RELAX_POSITIONS)
			message.add("\'relax\'");
		else
			message.add("\'scf\'");
		
		// Forces should always be shown
		message.addLine();
		message.add("tprnfor");
		message.add("=");
		message.add(".TRUE.");
		
		// Potential directory
		message.addLine();
		message.add(indent);
		message.add("pseudo_dir");
		message.add("=");
		message.add(Word("\'")+job->potDirectory()+"\'");
		
		// Details
		message.addLine();
		message.add(indent);
		message.add("disk_io");
		message.add("=");
		message.add("\'low\'");
		
		// Print section
		align.length(4);
		align.fill(LEFT);
		align[1] = RIGHT;
		Output::newline();
		Output::print(message, align);
		
		// End of section
		Output::newline();
		Output::print("/");
		Output::newline();
	}
	
	// Print system settings
	Output::newline();
	Output::print("&SYSTEM");
	
	// Basis vectors will be given
	message.clear();
	message.addLine();
	message.add(indent);
	message.add("ibrav");
	message.add("=");
	message.add("0");
	
	// Print number of atoms
	message.addLine();
	message.add(indent);
	message.add("nat");
	message.add("=");
	message.add(iso.numAtoms());
	
	// Print number of atom types
	message.addLine();
	message.add(indent);
	message.add("ntyp");
	message.add("=");
	message.add(iso.atoms().length());
	
	// Job settings
	if (job)
	{
		
		// Cutoff
		message.addLine();
		message.add(indent);
		message.add("ecutwfc");
		message.add("=");
		message.add("40");
		
		// Occupations
		message.addLine();
		message.add(indent);
		message.add("occupations");
		message.add("=");
		message.add("\'smearing\'");
		
		// Smearing
		message.addLine();
		message.add(indent);
		message.add("smearing");
		message.add("=");
		message.add("\'gaussian\'");
		
		// Smearing width
		message.addLine();
		message.add(indent);
		message.add("degauss");
		message.add("=");
		message.add(Num<double>::energyToAU(0.2));
	}
	
	// Print section
	align.length(4);
	align.fill(LEFT);
	align[1] = RIGHT;
	Output::newline();
	Output::print(message, align);
	
	// Print end of system settings
	Output::newline();
	Output::print("/");
	
	// Electron settings
	if (job)
	{
		Output::newline();
		Output::newline();
		Output::print("&ELECTRONS");
		Output::newline();
		Output::print("/");
	}
	
	// Ions settings
	if (job)
	{
		if ((job->calculation() == EC_RELAX_ALL) || (job->calculation() == EC_RELAX_POSITIONS))
		{
			Output::newline();
			Output::newline();
			Output::print("&IONS");
			Output::newline();
			Output::print("/");
		}
	}
	
	// Cell settings
	if (job)
	{
		if (job->calculation() == EC_RELAX_ALL)
		{
			Output::newline();
			Output::newline();
			Output::print("&CELL");
			Output::newline();
			Output::print("/");
		}
	}
	
	// Form basis vectors
	int i, j;
	message.clear();
	message.addLines(3);
    for (i = 0; i < 3; ++i)
    {
		message.addLine();
		message.addWords(4);
		message.add(indent);
		for (j = 0; j < 3; ++j)
			message.add(iso.basis().vectors()(i, j), prec);
    }
	
	// Print basis vectors
	Output::newline();
	Output::newline();
	Output::print("CELL_PARAMETERS angstrom");
	Output::newline();
	Output::print(message, RIGHT);
	
	// Form atom species
	bool found;
	message.clear();
	message.addLines(iso.atoms().length());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Add atom and mass
		message.addLine();
		message.addWords(3);
		message.add(indent);
		message.add(iso.atoms()[i][0].element().symbol());
		message.add(iso.atoms()[i][0].element().mass());
		
		// If job is defined then add potential
		if (job)
		{
			found = false;
			for (j = 0; j < job->elements().length(); ++j)
			{
				if (job->elements()[j] == iso.atoms()[i][0].element())
				{
					message.add(job->potentials()[j]);
					found = true;
					break;
				}
			}
			if (!found)
			{
				Output::newline(ERROR);
				Output::print("A quantum espresso potential was not supplied for ");
				Output::print(iso.atoms()[i][0].element().symbol());
				Output::quit();
			}
		}
	}
	align.length(4);
	align.fill(LEFT);
	align[2] = RIGHT;
	
	// Print atomic species
	Output::newline();
	Output::newline();
	Output::print("ATOMIC_SPECIES");
	Output::newline();
	Output::print(message, align);
	
	// Form atoms
	int k;
	message.clear();
	message.addLines(iso.numAtoms());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			message.addLine();
			message.addWords(5);
			message.add(indent);
			message.add(iso.atoms()[i][j].element().symbol());
			for (k = 0; k < 3; ++k)
			{
				if (coordinates == FRACTIONAL)
					message.add(iso.atoms()[i][j].fractional()[k], prec);
				else
					message.add(iso.atoms()[i][j].cartesian()[k], prec);
			}
			if ((iso.atoms()[i][j].fixed()[0]) || (iso.atoms()[i][j].fixed()[1]) || (iso.atoms()[i][j].fixed()[2]))
			{
				message.addWords(3);
				for (k = 0; k < 3; ++k)
				{
					if (iso.atoms()[i][j].fixed()[k])
						message.add("0");
					else
						message.add("1");
				}
			}
		}
	}
	align.length(8);
	align.fill(RIGHT);
	align[0] = align[1] = LEFT;
	
	// Print positions
	Output::newline();
	Output::newline();
	Output::print("ATOMIC_POSITIONS ");
	if (coordinates == FRACTIONAL)
		Output::print("crystal");
	else
		Output::print("angstrom");
	Output::newline();
	Output::print(message, align);
	
	// Print kpoints if needed
	if (kpoints)
	{
		
		// Print start of kpoints
		Output::newline();
		Output::newline();
		Output::print("K_POINTS");
		
		// Gamma only
		if ((kpoints->mesh()[0] == 1) && (kpoints->mesh()[1] == 1) && (kpoints->mesh()[2] == 1))
			Output::print(" gamma");
		
		// Print points
		else
		{
			Output::print(" automatic");
			Output::newline();
			Output::print(indent);
			Output::print(" ");
			for (i = 0; i < 3; ++i)
			{
				Output::print(kpoints->mesh()[i]);
				Output::print(" ");
			}
			Output::print("0 0 0");
		}
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



/* bool Espresso::Files::isFormat(const Text& content)
 *
 * Return whether file is in quantum espresso format
 */

bool Espresso::Files::isFormat(const Text& content)
{
	
	// Variables that must be found
	bool foundSystem = false;
	bool foundSpecies = false;
	bool foundPositions = false;
	
	// Loop over lines in file
	int i, j;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Loop over words on line
		for (j = 0; j < content[i].length(); ++j)
		{
			
			// Found keyword
			if (content[i][j].equal("&system", false))
				foundSystem = true;
			else if (content[i][j].equal("atomic_species", false))
				foundSpecies = true;
			else if (content[i][j].equal("atomic_positions", false))
				foundPositions = true;
		}
		
		// Break if all keywords have been found
		if ((foundSystem) && (foundSpecies) && (foundPositions))
			return true;
	}
	
	// Return that file is not in correct format if at this point
	return false;
}



/* Espresso::Files::Keyword Espresso::Files::keyword(const Word& word)
 *
 * Return the keyword corresponding to a word
 */

Espresso::Files::Keyword Espresso::Files::keyword(const Word& word)
{
	if (word.equal("&control", false))
		return EK_CONTROL;
	if (word.equal("&system", false))
		return EK_SYSTEM;
	if (word.equal("&electrons", false))
		return EK_ELECTRONS;
	if (word.equal("&ions", false))
		return EK_IONS;
	if (word.equal("&cell", false))
		return EK_CELL;
	if (word.equal("atomic_species", false))
		return EK_SPECIES;
	if (word.equal("atomic_positions", false))
		return EK_POSITIONS;
	if (word.equal("k_points", false))
		return EK_KPOINTS;
	if (word.equal("cell_parameters", false))
		return EK_PARAMETERS;
	if (word.equal("constraints", false))
		return EK_CONSTRAINTS;
	if (word.equal("occupations", false))
		return EK_OCCUPATIONS;
	return EK_UNKNOWN;
}
