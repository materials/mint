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
#include "launcher.h"
#include "settings.h"
#include "about.h"
#include "help.h"
#include "randomStructure.h"
#include "unique.h"
#include "pointGroup.h"
#include "spaceGroup.h"
#include "interstitial.h"
#include "gaPredict.h"
#include "pdf.h"
#include "fileSystem.h"
#include "language.h"
#include "timer.h"
#include "output.h"
#include <cstdlib>



/* void Storage::addISO()
 *
 * Add structure to current list
 */

void Storage::addISO()
{
	_updateSymmetry += true;
	_id += ++_curID;
	_iso.add();
	_needsGeneration += true;
	_format.add();
	_symmetry.add();
	_baseName.add();
	_history.add();
}



/* void Storage::removeISO(int index)
 *
 * Remove structure from current list
 */

void Storage::removeISO(int index)
{
	_updateSymmetry.remove(index);
	_id.remove(index);
	_iso.remove(index);
	_needsGeneration.remove(index);
	_format.remove(index);
	_symmetry.remove(index);
	_baseName.remove(index);
	_history.remove(index);
}



/* void Storage::printLabel(int index, int length)
 *
 * Print label for structure
 */

void Storage::printLabel(int index, int length)
{
	
	// Return if there is only one structure in memory
	length = (length == -1) ? _iso.length() : length;
	if (length <= 1)
		return;
	
	// Set method to standard
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Get the order of the number
	int numLength = 0;
	double temp = (double)_id[index];
	while (temp > 0.9999)
	{
		++numLength;
		temp /= 10;
	}
	int width = 12 + numLength + _history[index].length();
	
	// Print header
	int i;
	Output::newline();
	for (i = 0; i < width; ++i)
		Output::print("=");
	
	// Print id
	Output::newline();
	Output::print("Structure ");
	Output::print(_id[index]);
	Output::print(": ");
	Output::print(_history[index]);
	
	// Print footer
	Output::newline();
	for (i = 0; i < width; ++i)
		Output::print("=");
	
	// Reset output method
	Output::method(origMethod);
}



/* void Launcher::start(int argc, char** argv)
 *
 * Launch Mint functions
 */

void Launcher::start(int argc, char** argv)
{
	
	// Set output to restricted
	Output::method(RESTRICTED);
	
	// Start the timer
	Timer time;

	// Get command line arguments
	Words arguments = getArguments(argc, argv);
	
	// Turn arguments into function calls
	Functions functions = getFunctions(arguments);
	
	// Print about message and quit if no functions were found (one function should always exist)
	if (functions.length() == 1)
	{
		About::print();
		return;
	}
	
	// Check for settings console run
	if (runSettings(functions))
		return;
	
	// Check for help
	if (runHelp(functions))
		return;
	
	// Read files from input and get settings
	Storage data;
	bool keepFree = false;
	runSetup(data, functions, keepFree);
	fixCellParams(data, functions);
	
	// Loop over functions and run
	bool print = false;
	bool forcePrint = false;
	for (int i = 0; i < functions.length(); ++i)
	{
		switch (functions[i].keyword())
		{
			
			// Set output
			case KEY_OUTPUT:
				output(functions[i]);
				break;
			
			// Set tolerance
			case KEY_TOLERANCE:
				tolerance(functions[i]);
				break;
			
			// Set print
			case KEY_PRINT:
				forcePrint = true;
				break;
			
			// Remove atoms from the cell
			case KEY_REMOVE:
				generateStructure(data);
				removeAtoms(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Convert structure to reduced cell
			case KEY_REDUCED:
				generateStructure(data);
				reduced(data);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Convert structure to primitive cell
			case KEY_PRIMITIVE:
				generateStructure(data);
				primitive(data);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Convert structure to conventional cell
			case KEY_CONVENTIONAL:
				generateStructure(data);
				conventional(data);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Convert structure to most ideal cell
			case KEY_IDEAL:
				generateStructure(data);
				ideal(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Shift the origin
			case KEY_SHIFT:
				generateStructure(data);
				shift(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Transform the cell
			case KEY_TRANSFORM:
				generateStructure(data);
				transform(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Rotate positions
			case KEY_ROTATE:
				generateStructure(data);
				rotate(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Get symmetry
			case KEY_SYMMETRY:
				generateStructure(data);
				symmetry(data, functions[i]);
				print = false;
				break;
			
			// Unique groups
			case KEY_UNIQUE:
			case KEY_EQUIVALENT:
				generateStructure(data);
				unique(data, functions[i]);
				print = false;
				break;
			
			// Print details about structure
			case KEY_ABOUT:
				generateStructure(data);
				about(data, functions[i]);
				print = false;
				break;
			
			// Get point group
			case KEY_POINTGROUP:
				generateStructure(data);
				pointGroup(data, functions[i]);
				print = false;
				break;
			
			// Get space group
			case KEY_SPACEGROUP:
				generateStructure(data);
				spaceGroup(data, functions[i]);
				print = false;
				break;
			
			// Refine structure
			case KEY_REFINE:
				generateStructure(data);
				refine(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = true;
				break;
			
			// Get energy
			case KEY_ENERGY:
				generateStructure(data);
				energy(data);
				print = false;
				break;
			
			// Get forces
			case KEY_FORCES:
				generateStructure(data);
				forces(data);
				print = false;
				break;
			
			// Get diffraction diffraction pattern
			case KEY_DIFFRACTION:
				generateStructure(data);
				diffraction(data, functions[i]);
				print = false;
				break;
			
			// Perturb structure
			case KEY_PERTURB:
				generateStructure(data);
				perturb(data, functions[i]);
				print = true;
				break;
			
			// Relax structure
			case KEY_RELAX:
				generateStructure(data);
				relax(data, functions[i]);
				print = true;
				break;
			
			// Calculate phonons
			case KEY_PHONONS:
				generateStructure(data);
				phonons(data, functions[i]);
				print = false;
				break;
			
			// KMC simulation
			case KEY_KMC:
				generateStructure(data);
				kmc(data, functions[i]);
				print = false;
				break;
			
			// Optimize structure
			case KEY_OPT:
				optimize(data, functions[i]);
				forcePrint = true;
				print = true;
				break;
			
			// Get neighbors
			case KEY_NEIGHBORS:
				generateStructure(data);
				neighbors(data, functions[i]);
				print = false;
				break;
			
			// Get neighbors shells
			case KEY_SHELLS:
				generateStructure(data);
				shells(data, functions[i]);
				print = false;
				break;
			
			// Get coordination
			case KEY_COORDINATION:
				generateStructure(data);
				coordination(data, functions[i]);
				print = false;
				break;
			
			// Get interstitial sites
			case KEY_INTERSTITIAL:
				generateStructure(data);
				interstitial(data, functions[i]);
				data.updateSymmetry().fill(true);
				print = false;
				break;
			
			// Compare structures
			case KEY_COMPARE:
				generateStructure(data);
				compare(data, functions[i]);
				print = false;
				break;
			
			// Anything else
			default:
				break;
		}
	}
	
	// Print the structure
	if ((print) || (forcePrint))
	{
		generateStructure(data, keepFree);
		printStructure(data, functions);
	}
	
	// Print time if needed
	if (Settings::value<bool>(TIME_SHOW))
	{
		PrintMethod origMethod = Output::method();
		Output::method(STANDARD);
		Output::newline();
		Output::print("Elapsed time: ");
		Output::print(time.current(Settings::value<bool>(TIME_FORMAT), Settings::value<int>(TIME_PRECISION)));
		Output::method(origMethod);
	}
}



/* Words Launcher::getArguments(int argc, char** argv)
 *
 * Get arguments from command line
 */

Words Launcher::getArguments(int argc, char** argv)
{
	
	// Variable to store result
	Words arguments (argc - 1);
	
	// Save arguments
	for (int i = 1; i < argc; ++i)
		arguments[i-1] = argv[i];
	
	// Return result
	return arguments;
}



/* Launcher::Functions Launcher::getFunctions(const Words& arguments)
 *
 * Get functions from command line arguments
 */

Launcher::Functions Launcher::getFunctions(const Words& arguments)
{
	
	// Variable to store results
	Functions functions;
	
	// Create first function to store initial information
	functions.add();
	functions.last().keyword(KEY_NONE);
	
	// Loop over arguments
	Keyword curKeyword;
	for (int i = 0; i < arguments.length(); ++i)
	{
		
		// Get current keyword id
		curKeyword = getKeyword(arguments[i]);
		
		// Found a new function
		if (curKeyword != KEY_NONE)
		{
			functions.add();
			functions.last().keyword(curKeyword);
		}
		
		// Save argument to current function
		else
			functions.last().addArgument(arguments[i]);
	}
	
	// Return result
	return functions;
}



/* Keyword Launcher::getKeyword(const Word& argument)
 *
 * Return keyword id for a word
 */

Keyword Launcher::getKeyword(const Word& argument)
{
	
	// Return if word does not start with -
	if (argument[0] != '-')
		return KEY_NONE;
	
	// Return if a number
	if (argument.length() > 1)
	{
		if (((argument[1] >= '0') && (argument[1] <= '9')) || (argument[1] == '.'))
			return KEY_NONE;
	}
	
	// Settings
	if (argument.equal("-settings", false, 4))
		return KEY_SETTINGS;
	
	// Help
	if ((argument.equal("-h", false)) || (argument.equal("-help", 5)))
		return KEY_HELP;
	
	// Number of processors
	if ((argument.equal("-n", false)) || (argument.equal("-np", false)))
		return KEY_NUMPROCS;
	
	// Output
	if ((argument.equal("-output", false, 4)) || (argument.equal("-display", false, 5)))
		return KEY_OUTPUT;
	
	// Time
	if (argument.equal("-time", false, 4))
		return KEY_TIME;
	
	// Tolerance
	if (argument.equal("-tolerance", false, 4))
		return KEY_TOLERANCE;
	
	// Print
	if ((argument.equal("-print", false, 5)) || (argument.equal("-write", false, 3)))
		return KEY_PRINT;
	
	// File name
	if ((argument.equal("-name", false, 5)) || (argument.equal("-file", false, 5)))
		return KEY_NAME;
	
	// Remove atoms
	if (argument.equal("-remove", false, 4))
		return KEY_REMOVE;
	
	// Fix cell
	if (argument.equal("-fix", false, 4))
		return KEY_FIX;
	
	// Find neighbors
	if (argument.equal("-neighbors", false, 6))
		return KEY_NEIGHBORS;
	
	// Find neighbors shells
	if (argument.equal("-shells", false, 5))
		return KEY_SHELLS;
	
	// Find coordination
	if (argument.equal("-coordination", false, 6))
		return KEY_COORDINATION;
	
	// Transform cell
	if (argument.equal("-transform", false, 5))
		return KEY_TRANSFORM;
	
	// Rotate positions
	if (argument.equal("-rotate", false, 4))
		return KEY_ROTATE;
	
	// Make reduced cell
	if (argument.equal("-reduced", false, 4))
		return KEY_REDUCED;
	
	// Make primitive cell
	if (argument.equal("-primitive", false, 5))
		return KEY_PRIMITIVE;
	
	// Shift the origin
	if (argument.equal("-shift", false, 5))
		return KEY_SHIFT;
	
	// Make conventional cell
	if (argument.equal("-conventional", false, 5))
		return KEY_CONVENTIONAL;
	
	// Make ideal cell
	if (argument.equal("-ideal", false, 4))
		return KEY_IDEAL;
	
	// Symmetry
	if (argument.equal("-symmetry", false, 4))
		return KEY_SYMMETRY;
	
	// Unique
	if (argument.equal("-unique", false, 4))
		return KEY_UNIQUE;
	
	// Equivalent
	if (argument.equal("-equivalent", false, 3))
		return KEY_EQUIVALENT;
	
	// About
	if (argument.equal("-about", false, 6))
		return KEY_ABOUT;
	
	// Point group
	if (argument.equal("-point", false, 5))
		return KEY_POINTGROUP;
	
	// Space group
	if (argument.equal("-space", false, 5))
		return KEY_SPACEGROUP;
	
	// Refine
	if (argument.equal("-refine", false, 4))
		return KEY_REFINE;
	
	// Energy
	if (argument.equal("-energy", false, 5))
		return KEY_ENERGY;
	
	// Forces
	if (argument.equal("-force", false, 5))
		return KEY_FORCES;
	
	// Phonons
	if (argument.equal("-phonons", false, 5))
		return KEY_PHONONS;
	
	// KMC
	if (argument.equal("-kmc", false, 4))
		return KEY_KMC;
	
	// Diffraction
	if (argument.equal("-diffraction", false, 5))
		return KEY_DIFFRACTION;
	
	// Perturb
	if (argument.equal("-perturb", false, 5))
		return KEY_PERTURB;
	
	// Relax
	if (argument.equal("-relax", false, 4))
		return KEY_RELAX;
	
	// Optimize
	if (argument.equal("-optimize", false, 4))
		return KEY_OPT;
	
	// Interstitials
	if (argument.equal("-interstitial", false, 4))
		return KEY_INTERSTITIAL;
	
	// Compare
	if (argument.equal("-compare", false, 5))
		return KEY_COMPARE;
	
	// Did not recognize word as a keyword
	Output::newline(ERROR);
	Output::print("Did not recognize function: ");
	Output::print(argument);
	Output::quit();
	return KEY_NONE;
}



/* void Launcher::runSetup(Storage& data, Functions& functions, bool& keepFree)
 *
 * Read files from command line input
 */

void Launcher::runSetup(Storage& data, Functions& functions, bool& keepFree)
{
	
	// Loop over arguments before first function
	int i;
	OList<Word> files;
	for (i = functions[0].arguments().length() - 1; i >= 0; --i)
	{
		
		// Current word is a file
		if (File::exists(functions[0].arguments()[i]))
		{
			files += functions[0].arguments()[i];
			functions[0].removeArgument(i);
		}
		
		// File does not exist
		else
		{
			Output::newline(ERROR);
			Output::print("File ");
			Output::print(functions[0].arguments()[i]);
			Output::print(" does not exist");
			Output::quit();
		}
	}
	
	// Check for global settings file
	Text settings;
	Words settingsFiles;
	if (File::exists(Settings::globalFile()))
	{
		settings = Read::text(Settings::globalFile());
		settings.split('=');
		Settings::set(settings);
		settingsFiles += Settings::globalFile();
	}
	
	// Read all files
	OList<Text> content (files.length());
	for (i = 0; i < files.length(); ++i)
		content[i] = Read::text(files[i]);

	// Loop over files to check for settings files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Save settings file and remove
		if (Settings::isFormat(content[i]))
		{
			Settings::set(content[i]);
			settingsFiles += files[i];
			content.remove(i);
			files.remove(i);
		}
	}
	
	// Get settings from command line
	getCommandLineSettings(functions, keepFree);
	
	// Output
	Output::newline();
	Output::print("Reading data from files");
	Output::increase();
	
	// Print settings
	if (settingsFiles.length())
	{
		Output::newline();
		Output::print("Reading settings from ");
		Output::print(settingsFiles);
	}
	
	// Look for structure files
	StructureFormat structureFormat;
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Check if a structure file
		structureFormat = StructureIO::getFormat(content[i]);
		if (structureFormat != SF_UNKNOWN)
		{
			
			// Output
			Output::newline();
			Output::print("Reading structure from ");
			Output::print(files[i]);
			Output::increase();
			
			// Read file
			data.addISO();
			data.iso().last() = StructureIO::read(content[i], structureFormat, Settings::value<double>(TOLERANCE), \
				Settings::value<double>(CLUSTERTOL));
			data.baseName().last() = files[i];
			data.history().last() += "Read from ";
			data.history().last() += files[i];
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Save the format 
			data.format().last() = structureFormat;
			
			// Output
			Output::decrease();
		} 
	}
	
	// Look for pdf files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Found a PDF file
		if (PDF::isFormat(content[i]))
		{
			
			// Diffraction has already been set
			if (data.diffraction().isSet())
			{
				Output::newline(ERROR);
				Output::print("Only one PDF file may be passed as input");
				Output::quit();
			}
			
			// Output
			Output::newline();
			Output::print("Reading PDF file from ");
			Output::print(files[i]);
			Output::increase();
			
			// Read data
			if (data.iso().length() == 0)
			{
				data.addISO();
				PDF::read(data.diffraction(), content[i], &data.iso().last());
			}
			else
				PDF::read(data.diffraction(), content[i]);
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Output
			Output::decrease();
		}
	}
	
	// Look for potential files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Found a potential file
		if (Potential::isFormat(content[i]))
		{
			
			// Output
			Output::newline();
			Output::print("Reading potential from ");
			Output::print(files[i]);
			Output::increase();
			
			// Read data
			data.potential().set(content[i]);
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Output
			Output::decrease();
		}
	}
	
	// Look for KMC files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Found a kmc file
		if (KMC::isKMCFile(content[i]))
		{
			
			// KMC has already been set
			if (data.kmc().length() > 0)
			{
				Output::newline(ERROR);
				Output::print("Only one KMC file may be passed as input");
				Output::quit();
			}
			
			// Output
			Output::newline();
			Output::print("Saving KMC data from ");
			Output::print(files[i]);
			Output::increase();
			
			// Save data
			data.kmc() = content[i];
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Output
			Output::decrease();
		}
	}
	
	// Look for force constant files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Found a phonons file
		if (Phonons::isForceConstantFile(content[i]))
		{
			
			// Phonons has already been set
			if (data.phonons().isSet())
			{
				Output::newline(ERROR);
				Output::print("Only one force constant file may be passed as input");
				Output::quit();
			}
			
			// Output
			Output::newline();
			Output::print("Reading force constant data from ");
			Output::print(files[i]);
			Output::increase();
			
			// Read data
			data.phonons().set(content[i]);
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Output
			Output::decrease();
		}
	}
	
	// Look for diffraction files
	for (i = content.length() - 1; i >= 0; --i)
	{
		
		// Found a diffraction file
		if (Diffraction::isFormat(content[i]))
		{
			
			// Diffraction has already been set
			if (data.diffraction().isSet())
			{
				Output::newline(ERROR);
				Output::print("Only one diffraction file may be passed as input");
				Output::quit();
			}
			
			// Output
			Output::newline();
			Output::print("Reading diffraction data from ");
			Output::print(files[i]);
			Output::increase();
			
			// Read data
			data.diffraction().set(content[i]);
			
			// Remove file
			content.remove(i);
			files.remove(i);
			
			// Output
			Output::decrease();
		}
	}
	
	// Did not recognize file
	if (files.length())
	{
		Output::newline(ERROR);
		Output::print("Could not determine format of file ");
		Output::print(files[0]);
		Output::quit();
	}
	
	// Set clusterTol if smaller than tolerance
	if (Settings::value<double>(TOLERANCE) > Settings::value<double>(CLUSTERTOL))
		Settings::value<double>(CLUSTERTOL, Settings::value<double>(TOLERANCE));
	
	// Set random structure generation max number of trial loops
	RandomStructure::maxTrialLoops(Settings::value<int>(RANDSTR_MAXLOOPS));
	RandomStructure::minBondFraction(Settings::value<double>(RANDSTR_MINBOND));
	
	// Output
	Output::decrease();
}



/* void Launcher::getCommandLineSettings(Functions& functions, bool& keepFree)
 *
 * Get settings from command line
 */

void Launcher::getCommandLineSettings(Functions& functions, bool& keepFree)
{
	
	// Set number of processors to use for external mpi program calls
	setupJobSize(functions);
	
	// Setup output after checking for level being set on command line
	setupOutput(functions);

	// Set tolerance
	setupTolerance(functions);

	// Set the time and precision to use
	setupTime(functions);
	
	// Set the print settings
	setupPrint(functions, keepFree);
}



/* void Launcher::setupJobSize(Functions& functions)
 *
 * Setup external job size
 */

void Launcher::setupJobSize(Functions& functions)
{
	
	// Set number of processors based on current settings
	Multi::jobSize(Settings::value<int>(NUMPROCS));
	
	// Loop over functions and check for numprocs call
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_NUMPROCS)
		{
			for (int j = 0; j < functions[i].arguments().length(); ++j)
			{
				if (Language::isInteger(functions[i].arguments()[j]))
				{
					Multi::jobSize(atoi(functions[i].arguments()[j].array()));
					functions.remove(i);
					return;
				}
			}
		}
	}
}



/* void Launcher::setupOutput(Functions& functions)
 *
 * Setup output levels
 */

void Launcher::setupOutput(Functions& functions)
{
	
	// Loop over functions and check for output call
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_OUTPUT)
		{
			output(functions[i]);
			functions.remove(i);
			break;
		}
	}
	
	// Setup output
	Output::maxLevel(Settings::value<int>(OUTPUT_LEVEL));
	Output::spacesPerLevel(Settings::value<int>(OUTPUT_TAB));
}



/* void Launcher::setupTime(Functions& functions)
 *
 * Setup time output and precision
 */

void Launcher::setupTime(Functions& functions)
{
	
	// Loop over functions and check for time call
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_TIME)
		{
			
			// Check for precision being set and format
			for (int j = 0; j < functions[i].arguments().length(); ++j)
			{
				if (Language::isInteger(functions[i].arguments()[j]))
					Settings::value<int>(TIME_PRECISION, atoi(functions[i].arguments()[j].array()));
				if (functions[i].arguments()[j].equal("seconds", false, 3))
					Settings::value<bool>(TIME_FORMAT, false);
			}
			
			// Save that time should be shown and finish
			Settings::value<bool>(TIME_SHOW, true);
			functions.remove(i);
			break;
		}
	}
}



/* void Launcher::setupTolerance(Functions& functions)
 *
 * Setup tolerance
 */

void Launcher::setupTolerance(Functions& functions)
{
	
	// Loop over functions and check for tolerance call
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_TOLERANCE)
		{
			
			// Save tolerance and finish
			tolerance(functions[i]);
			functions.remove(i);
			break;
		}
	}
}



/* void Launcher::setupPrint(Functions& functions)
 *
 * Setup print functions
 */

void Launcher::setupPrint(Functions& functions, bool& keepFree)
{
	
	// Loop over functions and check for print call
	int i, j;
	for (i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_PRINT)
		{
			
			// Check for format being set
			StructureFormat tempFormat = StructureIO::structureFormat(functions[i].arguments());
			if (tempFormat != SF_UNKNOWN)
				Settings::value<StructureFormat>(STRUCTURE_FORMAT, tempFormat);
			
			// Look for keywords
			for (j = 0; j < functions[i].arguments().length(); ++j)
			{
				
				// Found print to stdout
				if ((functions[i].arguments()[j].equal("stdout", false)) || \
					(functions[i].arguments()[j].equal("screen", false)))
					Settings::value<bool>(USE_STDOUT, true);
				
				// Fractional coordinates
				else if ((functions[i].arguments()[j].equal("fractional", false, 4)) || \
					(functions[i].arguments()[j].equal("direct", false, 3)))
					Settings::value<CoordinateType>(COORDINATES, FRACTIONAL);
				
				// Cartesian coordinates
				else if (functions[i].arguments()[j].equal("cartesian", false, 4))
					Settings::value<CoordinateType>(COORDINATES, CARTESIAN);
				
				// Keep unassigned atoms free
				else if (functions[i].arguments()[j].equal("free", false, 4))
					keepFree = true;
			}
		}
	}
}



/* void Launcher::removeAtoms(Storage& data, const Function& function)
 *
 * Remove atoms from the structure
 */

void Launcher::removeAtoms(Storage& data, const Function& function) {
	
	// Return if there are no structures
	if (!data.iso().length())
		return;
	
	// Output
	Output::newline();
	Output::print("Removing atoms from the structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Loop over structures to remove atoms
	int i, j;
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Removing atoms in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get atoms to remove - they should already be sorted from getAtoms function
		List<Atom*> atoms = getAtoms(data.iso()[i], function, true, false);
		
		// Remove atoms
		for (j = atoms.length() - 1; j >= 0; --j)
			data.iso()[i].removeAtom(atoms[j]->atomNumber());
		
		// Add change to comment
		data.history()[i] += " > removed atoms";
		
		// Output
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::fixCellParams(Storage& data, Functions& functions)
 *
 * Fix any cell parameters that are set explicitly
 */

void Launcher::fixCellParams(Storage& data, Functions& functions)
{
	
	// Return if there are no structures
	if (!data.iso().length())
		return;
	
	// Look for fix command in functions
	int i, j, k, m;
	for (i = 0; i < functions.length(); ++i)
	{
		
		// Found fix
		if (functions[i].keyword() == KEY_FIX)
		{
			
			// Check for basis or atoms in arguments
			bool fixBasis = false;
			bool fixAtoms = false;
			for (j = 0; j < functions[i].arguments().length(); ++j)
			{

				// Found basis or lattice
				if ((functions[i].arguments()[j].equal("basis", false, 3)) || \
					(functions[i].arguments()[j].equal("lattice", false, 3)))
					fixBasis = true;

				// Found atoms or positions
				else if ((functions[i].arguments()[j].equal("atoms", false, 4)) || \
					(functions[i].arguments()[j].equal("positions", false, 3)))
					fixAtoms = true;
			}
			if ((!fixBasis) && (!fixAtoms))
				fixBasis = fixAtoms = true;

			// Output
			Output::newline();
			Output::print("Fixing ");
			if (fixBasis)
			{
				Output::print("basis");
				if (fixAtoms)
					Output::print(" and ");
			}
			if (fixAtoms)
				Output::print("atomic positions");
			Output::print(" where already set");
			Output::increase();
			
			// Loop over structures
			for (j = 0; j < data.iso().length(); ++j)
			{
				
				// Fix basis
				if (fixBasis)
				{
					
					// Basis has been set
					if (data.iso()[j].basis().volume() > 0)
					{
						for (k = 0; k < 3; ++k)
						{
							data.iso()[j].basis().lengthFixed(k, true);
							data.iso()[j].basis().angleFixed(k, true);
						}
					}
				}
				
				// Fix atoms
				if (fixAtoms)
				{
					
					// Loop over atoms and fix those that have been set
					for (k = 0; k < data.iso()[j].atoms().length(); ++k)
					{
						for (m = 0; m < data.iso()[j].atoms()[k].length(); ++m)
						{

							// Current atom has been set
							if (data.iso()[j].atoms()[k][m].assigned())
								data.iso()[j].atoms()[k][m].fixed(true);
						}
					}
				}
			}

			// Output
			Output::decrease();
			
			// Remove function and finish
			functions.remove(i);
			break;
		}
	}
}



/* void Launcher::output(const Function& function)
 *
 * Run output setup
 */

void Launcher::output(const Function& function)
{
	
	// Loop over arguments
	for (int i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found number
		if (Language::isInteger(function.arguments()[i]))
		{
			Output::maxLevel(atoi(function.arguments()[i].array()));
			Settings::value<int>(OUTPUT_LEVEL, atoi(function.arguments()[i].array()));
			break;
		}
		
		// Found all or everything
		if ((function.arguments()[i].equal("all", false, 3)) || (function.arguments()[i].equal("everything", false, 4)))
		{
			Output::maxLevel(1000);
			Settings::value<int>(OUTPUT_LEVEL, 1000);
			break;
		}
		
		// Found none, off, or nothing
		if ((function.arguments()[i].equal("none", false, 3)) || (function.arguments()[i].equal("off", false, 3)) || \
			(function.arguments()[i].equal("nothing", false, 3)))
		{
			Output::maxLevel(0);
			Settings::value<int>(OUTPUT_LEVEL, 0);
			break;
		}
	}
}



/* void Launcher::tolerance(const Function& function)
 *
 * Set the tolerance
 */

void Launcher::tolerance(const Function& function)
{
	
	// Loop over arguments
	for (int i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found a number so save and finish
		if (Language::isNumber(function.arguments()[i], true))
		{
			Settings::value<double>(TOLERANCE, Language::fractionToNumber(function.arguments()[i]));
			return;
		}
	}
}



/* void Launcher::generateStructure(Storage& data, bool keepFree)
 *
 * Generate random structures if needed
 */

void Launcher::generateStructure(Storage& data, bool keepFree)
{
	
	// Check if any structures are kept free and printed to mint - in this case they will not be generated
	int i;
	StructureFormat format;
	List<bool> keepISOFree (data.iso().length());
	keepISOFree.fill(false);
	if (keepFree)
	{
		for (i = 0; i < data.iso().length(); ++i)
		{

			// Skip if structure is already set or not keeping
			if (!data.needsGeneration()[i])
				continue;

			// Get format for printing structure
			format = Settings::value<StructureFormat>(STRUCTURE_FORMAT);
			if ((format == SF_AUTO) || (format == SF_UNKNOWN))
				format = data.format()[i];

			// Keep free if printing to mint format only
			if (format == SF_MINT)
				keepISOFree[i] = true;
		}
	}
	
	// Check if anything needs to be updated in structures
	bool update = false;
	for (i = 0; i < data.needsGeneration().length(); ++i)
	{
		if (data.needsGeneration()[i])
		{
			if ((!data.iso()[i].anyUnset()) && (data.iso()[i].basis().volume() > 1e-8))
				data.needsGeneration()[i] = false;
			else if (!keepISOFree[i])
				update = true;
		}
	}
	
	// Return if nothing needs to be updated 
	if (!update)
		return;
	
	// Output
	Output::newline();
	Output::print("Generating random data in structure");
	if (data.iso().length() > 1)
		Output::print("s");
	Output::increase();
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Current structure needs to be updated
		if ((data.needsGeneration()[i]) && (!keepISOFree[i]))
		{
			
			// Output
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Setting random data in structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Generate structure
			RandomStructure::generate(data.iso()[i], data.random(), Settings::value<double>(WYCKOFFBIAS));
			
			// Add change to comment
			data.history()[i] += " > randomly generated";
			
			// Output
			if (data.iso().length() > 1)
				Output::decrease();
			
			// Save that symmetry needs to be updated
			data.updateSymmetry()[i] = true;
			
			// Save that structure is set
			data.needsGeneration()[i] = false;
		}
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::printStructure(Storage& data, const Functions& functions)
 *
 * Print the structures
 */

void Launcher::printStructure(Storage& data, const Functions& functions)
{
	
	// Loop over functions and check if file name is being set
	int i, j, k;
	bool setName = false;
	for (i = 0; i < functions.length(); ++i)
	{
		
		// Found file name function
		if (functions[i].keyword() == KEY_NAME)
		{
			
			// Loop over arguments
			for (j = 0; j < functions[i].arguments().length(); ++j)
			{
				
				// Skip common words
				if ((functions[i].arguments()[j].equal("to", false)) || \
					(functions[i].arguments()[j].equal("as", false)))
					continue;
				
				// Save name to all files
				for (k = 0; k < data.baseName().length(); ++k)
					data.baseName()[k] = functions[i].arguments()[j];
				setName = true;
				break;
			}
			
			// Finish if name is set
			if (setName)
				break;
		}
	}
	
	// Set format of each structure to use
	List<StructureFormat> formats (data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		formats[i] = data.format()[i];
		if ((Settings::value<StructureFormat>(STRUCTURE_FORMAT) != SF_AUTO) && \
			(Settings::value<StructureFormat>(STRUCTURE_FORMAT) != SF_UNKNOWN))
			formats[i] = Settings::value<StructureFormat>(STRUCTURE_FORMAT);
	}
	
	// Figure out file names if needed
	int start;
	int length;
	OList<Word> fileNames (data.iso().length());
	if (!Settings::value<bool>(USE_STDOUT))
	{
		
		// Loop over files and strip extensions
		for (i = 0; i < data.iso().length(); ++i)
		{
			start = 0;
			length = data.baseName()[i].length();
			for (j = 0; j < data.baseName()[i].length(); ++j)
			{
				if (data.baseName()[i][j] == '/')
				{
					start = j+1;
					length = data.baseName()[i].length() - start;
				}
			}
			for (j = data.baseName()[i].length() - 1; j > start; --j)
			{
				if (data.baseName()[i][j] == '.')
				{
					length = j - start;
					break;
				}
			}
			fileNames[i].set(start, length, data.baseName()[i]);
		}
		
		// Check which file names appear more than once
		List<bool> multiple (fileNames.length());
		multiple.fill(false);
		for (i = 0; i < fileNames.length(); ++i)
		{
			if (multiple[i])
				continue;
			for (j = i + 1; j < fileNames.length(); ++j)
			{
				if (multiple[j])
					continue;
				if (fileNames[i] == fileNames[j])
					multiple[i] = multiple[j] = true;
			}
		}
		
		// Loop over files to set names
		int curNum;
		Word tempName;
		for (i = 0; i < fileNames.length(); ++i)
		{
			
			// Current name only appears once
			if (!multiple[i])
			{
				
				// Allowing structure overwrites so just add extension if needed, otherwise done
				if (Settings::value<bool>(OVERWRITE))
				{
					if (Settings::value<bool>(ADD_EXTENSION))
						StructureIO::addExtension(fileNames[i], formats[i]);
				}
				
				// Not allowing overwrites
				else
				{
					
					// File exists (if not then leave name as is)
					tempName = fileNames[i];
					if (Settings::value<bool>(ADD_EXTENSION))
						StructureIO::addExtension(tempName, formats[i]);
					if (File::exists(tempName))
					{
						
						// Loop until file name is found that does not exist
						curNum = 0;
						do
						{
							tempName = fileNames[i];
							tempName += '_';
							tempName += Language::numberToWord(++curNum);
							if (Settings::value<bool>(ADD_EXTENSION))
								StructureIO::addExtension(tempName, formats[i]);
						} while ((File::exists(tempName)) || (fileNameUsed(fileNames, tempName, i)));
					}
					
					// Save file name
					fileNames[i] = tempName;
				}
			}
			
			// Current name appears multiple times
			else
			{
				
				// Loop until file name is found that does not exist
				curNum = 0;
				do
				{
					tempName = fileNames[i];
					tempName += '_';
					tempName += Language::numberToWord(++curNum);
					if (Settings::value<bool>(ADD_EXTENSION))
						StructureIO::addExtension(tempName, formats[i]);
				} while ((fileNameUsed(fileNames, tempName, i)) || \
					((File::exists(tempName)) && (!Settings::value<bool>(OVERWRITE))));
				
				// Save file name
				fileNames[i] = tempName;
			}
		}
	}
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Printing to screen
		if (Settings::value<bool>(USE_STDOUT))
		{
			data.printLabel(i);
			StructureIO::write("stdout", data.iso()[i], formats[i], Settings::value<CoordinateType>(COORDINATES), \
				Settings::value<double>(TOLERANCE));
		}
		
		// Printing to file
		else
		{
			
			// Output
			Output::newline();
			Output::print("Writing structure ");
			Output::print(data.id()[i]);
			Output::print(" to ");
			Output::print(fileNames[i]);
			Output::increase();
					
			// Print
			StructureIO::write(fileNames[i], data.iso()[i], formats[i], Settings::value<CoordinateType>(COORDINATES), \
				Settings::value<double>(TOLERANCE));
			
			// Output
			Output::decrease();
		}
	}
}



/* bool Launcher::fileNameUsed(const OList<Word>& fileNames, const Word& curName, int numToCheck)
 *
 * Return whether a file name is used already
 */

bool Launcher::fileNameUsed(const OList<Word>& fileNames, const Word& curName, int numToCheck)
{
	for (int i = 0; i < numToCheck; ++i)
	{
		if (fileNames[i] == curName)
			return true;
	}
	return false;
}



/* void Launcher::reduced(Storage& data)
 *
 * Convert structure to reduced cell
 */

void Launcher::reduced(Storage& data)
{
	
	// Output
	Output::newline();
	Output::print("Transforming structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::print(" to reduced form");
	Output::increase();
	
	// Loop over structures and convert to reduced cell
	for (int i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting reduced cell for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		data.iso()[i].makeReduced();
		
		// Add change to comment
		data.history()[i] += " > reduced cell";
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::primitive(Storage& data)
 *
 * Convert structure to primitive cell
 */

void Launcher::primitive(Storage& data)
{
	
	// Output
	Output::newline();
	Output::print("Transforming structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::print(" to primitive form");
	Output::increase();
	
	// Loop over structures and convert to primitive cell
	for (int i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting primitive cell for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		data.iso()[i].makePrimitive(Settings::value<double>(TOLERANCE));
		
		// Add change to comment
		data.history()[i] += " > primitive cell";
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::conventional(Storage& data)
 *
 * Convert structure to conventional cell
 */

void Launcher::conventional(Storage& data)
{
	
	// Output
	Output::newline();
	Output::print("Transforming structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::print(" to conventional form");
	Output::increase();
	
	// Loop over structures and convert to conventional cell
	for (int i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting conventional cell for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		SpaceGroup::makeConventional(data.iso()[i], Settings::value<double>(TOLERANCE));
		
		// Add change to comment
		data.history()[i] += " > conventional cell";
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::ideal(Storage& data, const Function& function)
 *
 * Create ideal cell
 */

void Launcher::ideal(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Transforming structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::print(" to most ideal form");
	Output::increase();
	
	// Check if the number of atoms is being set
	int i;
	bool numAtomsAsMin = true;
	bool runForMinDis = false;
	int numAtoms = -1;
	double minDis = -1;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found an integer
		if (Language::isInteger(function.arguments()[i]))
		{
			numAtoms = atoi(function.arguments()[i].array());
			minDis = (double) numAtoms;
		}
		
		// Found a number
		else if (Language::isNumber(function.arguments()[i]))
			minDis = atof(function.arguments()[i].array());
		
		// Found max
		else if ((function.arguments()[i].equal("max", false)) || (function.arguments()[i].equal("most", false)))
			numAtomsAsMin = false;
		
		// Found distance
		else if (function.arguments()[i].equal("distance", false, 3))
			runForMinDis = true;
	}
	
	// Error if running for minimum distance and it has not been set
	if ((runForMinDis) && (minDis < 0))
	{
		Output::newline(ERROR);
		Output::print("Must set value of nearest image distance to make ideal cell by this method");
		Output::quit();
	}
	
	// Loop over structures and convert to conventional cell
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting ideal cell for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		if (runForMinDis)
			Symmetry::makeIdeal(data.iso()[i], minDis, Settings::value<double>(TOLERANCE));
		else
			Symmetry::makeIdeal(data.iso()[i], numAtoms == -1 ? data.iso()[i].numAtoms() : numAtoms, \
				numAtomsAsMin, Settings::value<double>(TOLERANCE));
		
		// Add change to comment
		data.history()[i] += " > ideal cell";
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::shift(Storage& data, const Function& function)
 * 
 * Shift the origin of the cell
 */

void Launcher::shift(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Shift positions of atoms in structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Loop over function arguments to get numbers
	int i;
	List<double> numbers;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (Language::isNumber(function.arguments()[i], true))
			numbers += Language::fractionToNumber(function.arguments()[i]);
	}
	
	// Create shift vector
	Vector3D shift;
	if (numbers.length() == 1)
		shift = numbers[0];
	else if (numbers.length() == 3)
	{
		for (i = 0; i < 3; ++i)
			shift[i] = numbers[i];
	}
	else
	{
		Output::newline(ERROR);
		Output::print("Expecting 1 or 3 numbers to shift positions but found ");
		Output::print(numbers.length());
		Output::quit();
	}
	
	// Create history tag
	Word tag = " > Shift by (";
	for (i = 0; i < 3; ++i)
	{
		tag += Language::numberToFraction(shift[i]);
		if (i != 2)
			tag += ",";
	}
	tag += ")";
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Shifting positions of atoms in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		data.iso()[i].shift(shift);
		
		// Add change to comment
		data.history()[i] += tag;
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::transform(Storage& data, const Function& function)
 *
 * Transform the cell
 */

void Launcher::transform(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Transforming cell");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Loop over function arguments to get numbers
	int i;
	List<double> numbers;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (Language::isNumber(function.arguments()[i], true))
			numbers += Language::fractionToNumber(function.arguments()[i]);
	}
	
	// Create transformation matrix
	int j;
	Matrix3D transformation = 0.0;
	if (numbers.length() == 3)
	{
		for (i = 0; i < 3; ++i)
			transformation(i, i) = numbers[i];
	}
	else if (numbers.length() == 9)
	{
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
				transformation(i, j) = numbers[3*i + j];
		}
	}
	else
	{
		Output::newline(ERROR);
		Output::print("Expecting 3 or 9 numbers to transform cell but found ");
		Output::print(numbers.length());
		Output::quit();
	}
	
	// Create history tag
	Word tag = " > [(";
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			tag += Language::numberToFraction(transformation(i, j));
			if (j != 2)
				tag += ',';
		}
		if (i != 2)
			tag += "),(";
	}
	tag += ")] cell";
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Applying cell transformation to structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Change the cell
		data.iso()[i].transform(transformation, Settings::value<double>(TOLERANCE));
		
		// Add change to comment
		data.history()[i] += tag;
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::rotate(Storage& data, const Function& function)
 *
 * Rotate position in the cell
 */

void Launcher::rotate(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Rotating positions in cell");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Loop over function arguments to get numbers
	int i;
	List<double> numbers;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (Language::isNumber(function.arguments()[i], true))
			numbers += Language::fractionToNumber(function.arguments()[i]);
	}
	
	// Create transformation matrix
	int j;
	Matrix3D rotation = 0.0;
	if (numbers.length() == 9)
	{
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
				rotation(i, j) = numbers[3*i + j];
		}
	}
	else
	{
		Output::newline(ERROR);
		Output::print("Expecting exactly 9 numbers to rotate positions but found ");
		Output::print(numbers.length());
		Output::quit();
	}
	
	// Create history tag
	Word tag = " > rotated by [(";
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			tag += Language::numberToFraction(rotation(i, j));
			if (j != 2)
				tag += ',';
		}
		if (i != 2)
			tag += "),(";
	}
	tag += ")]";
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Applying rotation to positions in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Rotate positions in the cell
		data.iso()[i].rotatePositions(rotation);
		
		// Add change to comment
		data.history()[i] += tag;
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::updateSymmetry(Storage& data)
 *
 * Get the symmetry of structures that need to be updated
 */

void Launcher::updateSymmetry(Storage& data)
{
	
	// Check if any structures need to be updated
	int i;
	for (i = 0; i < data.iso().length(); ++i)
	{
		if (data.updateSymmetry()[i])
			break;
	}
	if (i >= data.iso().length())
		return;
	
	// Output
	Output::newline();
	Output::print("Calculating symmetry of structure");
	if (data.iso().length() > 1)
		Output::print("s");
	Output::increase();
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Current symmetry needs to be updated
		if (data.updateSymmetry()[i])
		{
			
			// Output
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Setting symmetry for structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Get symmetry
			data.symmetry()[i].set(data.iso()[i], Settings::value<double>(TOLERANCE));
			
			// Output
			if (data.iso().length() > 1)
				Output::decrease();
			
			// Save that symmetry is done being updated
			data.updateSymmetry()[i] = false;
		}
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::symmetry(Storage& data, const Function& function)
 *
 * Get the symmetry of a structure
 */

void Launcher::symmetry(Storage& data, const Function& function)
{
	
	// Check if printing in matrix form
	int i;
	bool useJonesFaithful = true;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (function.arguments()[i].equal("matrix", false, 3))
			useJonesFaithful = false;
			
	}
	
	// Check if the symmetry needs to be updated
	updateSymmetry(data);
	
	// Loop over structures and print
	for (i = 0; i < data.symmetry().length(); ++i)
	{
		data.printLabel(i);
		data.symmetry()[i].print(useJonesFaithful);
	}
}



/* void Launcher::unique(Storage& data, const Function& function)
 *
 * Get unique groups
 */

void Launcher::unique(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Searching for unique groups");
	Output::increase();
	
	// Loop over arguments to get settings
	int i;
	bool getInt = false;
	bool areTransitions = false;
	double maxDistance = -1;
	Element tempElement;
	OList<Element> elements;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found element
		tempElement = Element::find(function.arguments()[i], false, false);
		if (tempElement.number() != 0)
			elements += tempElement;
		
		// Found jumps
		else if (function.arguments()[i].equal("jumps", false, 4))
			areTransitions = true;
		
		// Found a number
		else if (Language::isNumber(function.arguments()[i]))
			maxDistance = atof(function.arguments()[i].array());
		
		// Found interstitial
		else if (function.arguments()[i].equal("interstitial", false, 3))
			getInt = true;
	}
	
	// Variable to store unique groups
	OList<Unique> unique(data.iso().length());
	
	// Set settings for search
	for (i = 0; i < unique.length(); ++i)
	{
		
		// Save if transitions
		unique[i].areTransitions(areTransitions);
		
		// Save if using interstitials
		unique[i].useInterstitials(getInt);
		
		// Save the max distance
		if (maxDistance >= 0)
			unique[i].maxDistance(maxDistance);
		else if (areTransitions)
			unique[i].maxDistance(Settings::value<double>(MAXJUMPDISTANCE));
		
		// Save whether to generate equivalent jumps
		if (function.keyword() == KEY_EQUIVALENT)
			unique[i].getEquivalent(true);
		else
			unique[i].getEquivalent(false);
	}
	
	// Check if the symmetry needs to be updated
	updateSymmetry(data);
	
	// Loop over structures and get unique groups
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Searching for unique groups in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get unique groups
		unique[i].search(data.iso()[i], data.symmetry()[i], elements, Settings::value<double>(TOLERANCE));
		
		// Output
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Print data about unique groups
	for (i = 0; i < data.iso().length(); ++i)
	{
		data.printLabel(i);
		unique[i].print();
	}
	
	// Create unique groups if needed
	if (!areTransitions)
	{
		
		// Output
		Output::newline();
		Output::print("Creating structures for each unique group");
		Output::increase();
		
		// Loop over structures
		int j;
		int origLen = data.iso().length();
		List<int> indices;
		for (i = 0; i < origLen; ++i)
		{
			
			// Output
			if (origLen > 1)
			{
				Output::newline();
				Output::print("Creating unique groups in structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Skip if there are no groups for current structure
			if (unique[i].groups().length() == 0)
				continue;
			
			// Save that symmetry needs to be updated
			data.updateSymmetry()[i] = true;
			
			// Create new structures
			indices.length(unique[i].groups().length() - 1);
			for (j = 1; j < unique[i].groups().length(); ++j)
			{
				data.addISO();
				data.iso().last() = data.iso()[i];
				data.baseName().last() = data.baseName()[i];
				data.history().last() = data.history()[i];
				data.format().last() = data.format()[i];
				indices[j-1] = data.iso().length() - 1;
			}
			
			// Make structure for additional groups
			for (j = 0; j < indices.length(); ++j)
			{
				data.history()[indices[j]] += " > ";
				data.history()[indices[j]] += unique[i].makeISO(data.iso()[indices[j]], j+1);
			}
			
			// Create structure for first unique group
			data.history()[i] += " > ";
			data.history()[i] += unique[i].makeISO(data.iso()[i], 0);
			
			// Output
			if (origLen > 1)
				Output::decrease();
		}
		
		// Output
		Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::about(Storage& data, const Function& function)
 *
 * Print details about the structure or atoms in it
 */

void Launcher::about(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Generating information about structure");
	if (data.iso().length() > 1)
		Output::print("s");
	Output::increase();
	
	// Update symmetry information
	updateSymmetry(data);
	
	// Check if any atoms are being passed
	int i = 0;
	for (; i < function.arguments().length(); ++i)
	{
		if (function.arguments()[i].equal("atom", false, 4))
		{
			aboutAtoms(data, function);
			break;
		}
	}
	
	// Run as about structure
	if (i >= function.arguments().length())
		aboutStructure(data);
	
	// Outut
	Output::decrease();
}



/* void Launcher::aboutStructure(Storage& data)
 *
 * Print details about the structure
 */

void Launcher::aboutStructure(Storage& data)
{
	
	// Loop over structures to get space groups
	int i;
	OList<SpaceGroup> spaceGroups(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting space group for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get symmetry
		spaceGroups[i].set(data.iso()[i], Settings::value<double>(TOLERANCE));
		
		// Output
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Metric alignment
	List<PrintAlign> align(3);
	align.fill(RIGHT);
	align[2] = LEFT;
	
	// Loop over structures to print data
	int j, k;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print header
		data.printLabel(i);
		
		// Print space group
		spaceGroups[i].print();
		spaceGroups[i].pointGroup().print();
		
		// Print cell volume
		Output::method(STANDARD);
		Output::newline();
		Output::print("Unit cell volume: ");
		Output::print(data.iso()[i].basis().volume(), 8);
		Output::print(" Ang^3");
		Output::method(origMethod);
		
		// Print basis vectors
		Output::method(STANDARD);
		Output vectors;
		vectors.addLines(3);
		for (j = 0; j < 3; ++j)
		{
			vectors.addLine();
			vectors.add("   ");
			for (k = 0; k < 3; ++k)
				vectors.add(data.iso()[i].basis().vectors()(j, k), 8);
		}
		Output::newline();
		Output::print("Unit cell basis vectors");
		Output::newline();
		Output::print(vectors, RIGHT);
		Output::method(origMethod);
		
		// Print basis metrics
		Output::method(STANDARD);
		Output metrics;
		metrics.addLines(6);
		metrics.addLine();
		metrics.add("        A:");
		metrics.add(data.iso()[i].basis().lengths()[0]);
		metrics.add("Ang");
		metrics.addLine();
		metrics.add("        B:");
		metrics.add(data.iso()[i].basis().lengths()[1]);
		metrics.add("Ang");
		metrics.addLine();
		metrics.add("        C:");
		metrics.add(data.iso()[i].basis().lengths()[2]);
		metrics.add("Ang");
		metrics.addLine();
		metrics.add("    Alpha:");
		metrics.add(Num<double>::toDegrees(data.iso()[i].basis().angles()[0]));
		metrics.add("deg");
		metrics.addLine();
		metrics.add("     Beta:");
		metrics.add(Num<double>::toDegrees(data.iso()[i].basis().angles()[1]));
		metrics.add("deg");
		metrics.addLine();
		metrics.add("    Gamma:");
		metrics.add(Num<double>::toDegrees(data.iso()[i].basis().angles()[2]));
		metrics.add("deg");
		Output::newline();
		Output::print("Unit cell basis metrics");
		Output::newline();
		Output::print(metrics, align);
		Output::method(origMethod);
		
		// Print symmetry
		data.symmetry()[i].print();
		data.symmetry()[i].printSites();
	}
}



/* void Launcher::aboutAtoms(Storage& data, const Function& function)
 *
 * Print information about the atoms in a structure
 */

void Launcher::aboutAtoms(Storage& data, const Function& function)
{
	
	// Loop over structures and gather information about atoms
	int i;
	List<Atom*>::D2 atoms(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
		atoms[i] = getAtoms(data.iso()[i], function);
	
	// Loop over structures and print information for atoms in each
	int j, k, m;
	Word curLine;
	Orbit* orbit;
	Output message;
	List<int> atomNumbers;
	Words curOperation;
	SpecialPosition* special;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print header
		data.printLabel(i);
		
		// Loop over atoms
		atoms[i].sort();
		for (j = 0; j < atoms[i].length(); ++j)
		{
			
			// Print site number and label
			Output::method(STANDARD);
			Output::newline();
			Output::print("Atom ");
			Output::print(atoms[i][j]->atomNumber() + 1);
			Output::print(": ");
			Output::print(atoms[i][j]->element().symbol());
			
			// Save the coordinates
			message.clear();
			message.addLine();
			message.addTab();
			message.addTab();
			message.add("Fractional:");
			for (k = 0; k < 3; ++k)
				message.add(atoms[i][j]->fractional()[k], 8);
			message.addLine();
			message.addTab();
			message.addTab();
			message.add("Cartesian:");
			for (k = 0; k < 3; ++k)
				message.add(atoms[i][j]->cartesian()[k], 8);
			
			// Print the coordinates
			Output::newline();
			Output::tab();
			Output::print("Coordinates");
			Output::newline();
			Output::print(message, RIGHT);
			
			// Find the atom in the symmetry object of the current structure
			special = 0;
			for (k = 0; k < data.symmetry()[i].orbits().length(); ++k)
			{
				for (m = 0; m < data.symmetry()[i].orbits()[k].atoms().length(); ++m)
				{
					if (data.symmetry()[i].orbits()[k].atoms()[m] == atoms[i][j])
					{
						orbit = &(data.symmetry()[i].orbits()[k]);
						special = &(data.symmetry()[i].orbits()[k].specialPositions()[m]);
						break;
					}
				}
				if (m < data.symmetry()[i].orbits()[k].atoms().length())
					break;
			}
			
			// Generate the symmetry list
			message.clear();
			message.addLines(special->rotations().length());
			for (k = 0; k < special->rotations().length(); ++k)
			{
				curOperation = JonesFaithful::toString(special->rotations()[k], &(special->translations()[k]));
				message.addLine();
				message.addTab();
				message.addTab();
				for (m = 0; m < 3; ++m)
					message.add(curOperation[m]);
			}
			
			// Print the site symmetry
			Output::newline();
			Output::tab();
			Output::print("Site symmetry: ");
			Output::print("(");
			Output::print(special->getString(), true, false);
			Output::print(")");
			Output::newline();
			Output::print(message, LEFT);
			
			// Save list of atoms in the orbit
			atomNumbers.length(orbit->atoms().length());
			atomNumbers[0] = -1;
			for (k = 1; k < orbit->atoms().length(); ++k)
				atomNumbers[k] = orbit->atoms()[k]->atomNumber() + 1;
			atomNumbers.sort();
			atomNumbers[0] = orbit->atoms()[0]->atomNumber() + 1;
			
			// Print the equivalent atoms
			curLine = "    Atoms in orbit: ";
			for (k = 0; k < atomNumbers.length(); ++k)
			{
				curLine += Language::numberToWord(atomNumbers[k]);
				if (k != atomNumbers.length() - 1)
					curLine += ',';
				if (curLine.length() > 75)
				{
					Output::newline();
					Output::print(curLine);
					curLine = "        ";
				}
				else
					curLine += ' ';
			}
			if (curLine.length() > 8)
			{
				Output::newline();
				Output::print(curLine);
			}
			
			// Finished with atom
			Output::method(origMethod);
		}
	}
}



/* void Launcher::pointGroup(Storage& data, const Function& function)
 *
 * Get point group of structure or print information about point group
 */

void Launcher::pointGroup(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Searching point groups");
	Output::increase();
	
	// Loop over arguments and check for point group names
	int i;
	PointGroup tempPointGroup;
	OList<PointGroup> pointGroups;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (tempPointGroup.set(function.arguments()[i], false))
			pointGroups += tempPointGroup;
	}
	
	// Print point groups if found
	if (pointGroups.length())
	{
		for (i = 0; i < pointGroups.length(); ++i)
			pointGroups[i].print();
	}
	
	// No names were found
	else
	{
		
		// Loop over structures to find point groups
		pointGroups.length(data.iso().length());
		for (i = 0; i < data.iso().length(); ++i)
		{
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Determining the point group of structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Get point group information
			pointGroups[i].set(data.iso()[i], Settings::value<double>(TOLERANCE));
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
				Output::decrease();
		}
		
		// Print information about each point group
		for (i = 0; i < data.iso().length(); ++i)
		{
			data.printLabel(i);
			pointGroups[i].print();
		}
	}
	
	// Nothing was found so print all point group information
	if (!pointGroups.length())
		PointGroup::printAll();
	
	// Output
	Output::decrease();
}



/* void Launcher::spaceGroup(Storage& data, const Function& function)
 *
 * Get space group of structure of print information about space group
 */

void Launcher::spaceGroup(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Searching space groups");
	Output::increase();
	
	// Loop over arguments and check for space group names
	int i;
	SpaceGroup tempSpaceGroup;
	OList<SpaceGroup> spaceGroups;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (tempSpaceGroup.set(function.arguments()[i], true, false))
			spaceGroups += tempSpaceGroup;
	}
	
	// Print space groups if found
	if (spaceGroups.length())
	{
		for (i = 0; i < spaceGroups.length(); ++i)
			spaceGroups[i].printComplete();
	}
	
	// No names were found
	else
	{
		
		// Loop over structures to find space groups
		spaceGroups.length(data.iso().length());
		for (i = 0; i < data.iso().length(); ++i)
		{
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Determining the space group of structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Get space group information
			spaceGroups[i].set(data.iso()[i], Settings::value<double>(TOLERANCE));
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
				Output::decrease();
		}
		
		// Print information about each space group
		for (i = 0; i < data.iso().length(); ++i)
		{
			data.printLabel(i);
			spaceGroups[i].print();
		}
	}
	
	// Nothing was found so print all space group information
	if (!spaceGroups.length())
		SpaceGroup::printAll();
	
	// Output
	Output::decrease();
}



/* void Launcher::refine(Storage& data, const Function& function)
 *
 * Refine the structure
 */

void Launcher::refine(Storage& data, const Function& function)
{
	
	// Refine against diffraction pattern if available
	int i;
	if (data.diffraction().isSet())
	{
		
		// Output
		Output::newline();
		Output::print("Refining internal coordinates against diffraction pattern");
		Output::increase();
		
		// Check if the symmetry needs to be updated
		updateSymmetry(data);
		
		// Create history tag
		Word tag = " > Refined against diffraction pattern";
		
		// Loop over structures
		Diffraction tempPattern;
		for (i = 0; i < data.iso().length(); i++)
		{

			// Output if there is more than one structure
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Refining structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}

			// Refine the structure
			tempPattern.refine(data.iso()[i], data.symmetry()[i], data.diffraction());

			// Add change to comment
			data.history()[i] += tag;

			// Output if there is more than one structure
			if (data.iso().length() > 1)
				Output::decrease();
		}
		
		// Output
		Output::decrease();
		
		// Finished with function call
		return;
	}
	
	// Check for basis or atoms in arguments
	bool refineBasis = false;
	bool refineAtoms = false;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found basis or lattice
		if ((function.arguments()[i].equal("basis", false, 3)) || (function.arguments()[i].equal("lattice", false, 3)))
			refineBasis = true;
		
		// Found atoms or positions
		else if ((function.arguments()[i].equal("atoms", false, 4)) || \
			(function.arguments()[i].equal("positions", false, 3)))
			refineAtoms = true;
	}
	if ((!refineBasis) && (!refineAtoms))
		refineBasis = refineAtoms = true;
	
	// Output
	Output::newline();
	Output::print("Refining ");
	if (refineBasis)
	{
		Output::print("basis");
		if (refineAtoms)
			Output::print(" and ");
	}
	if (refineAtoms)
		Output::print("atomic positions");
	Output::increase();
	
	// Create history tag
	Word tag = " > Refined ";
	if (refineBasis)
	{
		tag += "basis";
		if (refineAtoms)
			tag += " and ";
	}
	if (refineAtoms)
		tag += "atomic positions";
	
	// Check if the symmetry needs to be updated
	updateSymmetry(data);
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Refining structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Refine the cell
		if ((refineBasis) && (refineAtoms))
			data.symmetry()[i].refine(data.iso()[i], Settings::value<double>(CLUSTERTOL));
		else if (refineBasis)
			data.symmetry()[i].refineBasis(data.iso()[i]);
		else if (refineAtoms)
			data.symmetry()[i].refineAtoms(data.iso()[i], Settings::value<double>(CLUSTERTOL));
		
		// Add change to comment
		data.history()[i] += tag;
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::energy(Storage& data)
 *
 * Get the energy of a structure
 */

void Launcher::energy(Storage& data)
{
	
	// Output
	Output::newline();
	Output::print("Calculating the energy of ");
	if (data.iso().length() == 1)
		Output::print("the structure");
	else
		Output::print(" each structure");
	Output::increase();
	
	// Check if symmetry should be used
	bool useSymmetry = data.potential().usesSymmetry();
	
	// Check if the symmetry needs to be updated
	if (useSymmetry)
		updateSymmetry(data);
	
	// Loop over structures and get energies
	int i;
	double energies[data.iso().length()];
	if (useSymmetry)
	{
		for (i = 0; i < data.iso().length(); ++i)
			data.potential().single(data.iso()[i], data.symmetry()[i], &energies[i], 0, false, true);
	}
	else
	{
		for (i = 0; i < data.iso().length(); ++i)
			data.potential().single(data.iso()[i], &energies[i], 0, false, true);
	}
	
	// Print each energy
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print the label
		data.printLabel(i);
		
		// Print the energy
		Output::method(STANDARD);
		Output::newline();
		if (data.potential().useReferences())
			Output::print("Formation");
		else
			Output::print("Total");
		Output::print(" energy: ");
		Output::print(energies[i], 8);
		Output::print(" eV");
		Output::newline();
		Output::print("Energy per atom: ");
		Output::print(energies[i] / data.iso()[i].numAtoms(), 8);
		Output::print(" eV");
		Output::method(origMethod);
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::forces(Storage& data)
 *
 * Get the forces in a structure
 */

void Launcher::forces(Storage& data)
{
	
	// Output
	Output::newline();
	Output::print("Calculating the forces in ");
	if (data.iso().length() == 1)
		Output::print("the structure");
	else
		Output::print(" each structure");
	Output::increase();
	
	// Check if symmetry should be used
	bool useSymmetry = data.potential().usesSymmetry();
	
	// Check if the symmetry needs to be updated
	if (useSymmetry)
		updateSymmetry(data);
	
	// Loop over structures and get forces
	int i;
	OList<Vector3D >::D2 forces(data.iso().length());
	if (useSymmetry)
	{
		for (i = 0; i < data.iso().length(); ++i)
			data.potential().single(data.iso()[i], data.symmetry()[i], 0, &forces[i], false, true);
	}
	else
	{
		for (i = 0; i < data.iso().length(); ++i)
			data.potential().single(data.iso()[i], 0, &forces[i], false, true);
	}
	
	// Alignment for printing forces
	List<PrintAlign> align(6);
	align.fill(RIGHT);
	align[1] = LEFT;
	
	// Print each set of forces
	int j, k;
	Word atomNumber;
	Output message;
	Vector3D tempForce;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print the label
		data.printLabel(i);
		
		// Make the forces output
		message.clear();
		message.addLines(data.iso()[i].numAtoms());
		for (j = 0; j < forces[i].length(); ++j)
		{
			message.addLine();
			message.addWords(6);
			message.add("   ");
			message.add("Atom");
			atomNumber = Language::numberToWord(j+1);
			atomNumber += ':';
			message.add(atomNumber);
			tempForce = data.iso()[i].basis().getCartesian(forces[i][j]);
			for (k = 0; k < 3; ++k)
				message.addSci(tempForce[k], 8);
		}
		
		// Print the forces
		Output::method(STANDARD);
		Output::newline();
		Output::print("Forces along cartesian axes (eV/ang):");
		Output::newline();
		Output::print(message, align);
		Output::method(origMethod);
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::diffraction(Storage& data, const Function& function)
 *
 * Calculate diffraction pattern
 */

void Launcher::diffraction(Storage& data, const Function& function)
{
	
	// Look over functions and check if there are any files that have been passed
	int i;
	OList<Text> files;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Current argument is a file
		if (File::exists(function.arguments()[i]))
		{
			
			// Current file is not an diffraction pattern
			files += Read::text(function.arguments()[i]);
			if (!Diffraction::isFormat(files.last()))
				files.length(files.length() - 1);
		}
	}
	
	// Loop over arguments and look for keywords
	bool fwhmOn = false;
	bool varOn = false;
	bool resOn = false;
	bool waveOn = false;
	bool minOn = false;
	bool maxOn = false;
	bool accOn = false;
	bool broaden = false;
	bool useGaussians = false;
	double fwhm = -1;
	double variance = -1;
	double resolution = -1;
	double wavelength = -1;
	double minTwoTheta = -1;
	double maxTwoTheta = -1;
	double accuracy = -1;
	double temp;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found fwhm
		if (function.arguments()[i].equal("fwhm", false))
		{
			fwhmOn = true;
			useGaussians = true;
		}
		
		// Found variance
		else if (function.arguments()[i].equal("variance", false, 3))
		{
			varOn = true;
			useGaussians = true;
		}
		
		// Found resolution
		else if (function.arguments()[i].equal("resolution", false, 3))
			resOn = true;
		
		// Found wavelength
		else if (function.arguments()[i].equal("wavelength", false, 4))
			waveOn = true;
		
		// Found minimum
		else if (function.arguments()[i].equal("minimum", false, 3))
			minOn = true;
		
		// Found maximum
		else if (function.arguments()[i].equal("maximum", false, 3))
			maxOn = true;
		
		// Found maximum
		else if (function.arguments()[i].equal("accuracy", false, 3))
			accOn = true;
		
		// Found broaden
		else if (function.arguments()[i].equal("broaden", false, 5))
			broaden = true;
		
		// Found a number
		else if (Language::isNumber(function.arguments()[i]))
		{
			temp = Language::fractionToNumber(function.arguments()[i]);
			if (fwhmOn)
				fwhm = temp;
			else if (varOn)
				variance = temp;
			else if (resOn)
				resolution = temp;
			else if (waveOn)
				wavelength = temp;
			else if (minOn)
				minTwoTheta = temp;
			else if (maxOn)
				maxTwoTheta = temp;
			else if (accOn)
				accuracy = temp;
			else
			{
				Output::newline(ERROR);
				Output::print("Found an unmatched number passed to diffraction function (");
				Output::print(function.arguments()[i]);
				Output::print(")");
				Output::quit();
			}
			fwhmOn = varOn = resOn = waveOn = minOn = maxOn = accOn = false;
		}
	}
	
	// Diffraction patterns have been set
	if ((data.diffraction().isSet()) && (files.length() > 0))
	{
		
		// Variable to store matches
		List<double> overlaps (files.length());
		
		// Output
		Output::newline();
		Output::print("Calculating overlap between patterns");
		Output::increase();
		
		// Loop over files and get overlaps
		Diffraction comp;
		for (i = 0; i < files.length(); ++i)
		{
			
			// Output
			if (files.length() > 1)
			{
				Output::newline();
				Output::print("Calculating overlap for pattern ");
				Output::print(i + 1);
				Output::increase();
			}
			
			// Get overlap
			comp.set(files[i]);
			overlaps[i] = comp.rFactor(data.diffraction());
			
			// Output
			if (files.length() > 1)
				Output::decrease();
		}
		
		// Output
		Output::decrease();
		
		// Print overlaps
		PrintMethod origMethod = Output::method();
		for (i = 0; i < overlaps.length(); ++i)
		{
			Output::method(STANDARD);
			Output::newline();
			Output::print("R factor");
			if (overlaps.length() > 1)
			{
				Output::print(" of pattern ");
				Output::print(i + 1);
			}
			Output::print(": ");
			Output::print(overlaps[i], 6);
			Output::method(origMethod);
		}
		
		// Finished
		return;
	}
	
	// Output
	Output::newline();
	Output::print("Calculating diffraction diffraction pattern of ");
	if (data.iso().length() == 1)
		Output::print("the structure");
	else
		Output::print(" each structure");
	Output::increase();
	
	// Check if the symmetry needs to be updated
	updateSymmetry(data);
	
	// Loop over structures and get diffraction pattern
	OList<Diffraction> patterns(data.iso().length());
	List<double> match(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Setting pattern for structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get match
		if (data.diffraction().isSet())
			match[i] = patterns[i].set(data.iso()[i], data.symmetry()[i], &data.diffraction(), true);
		else
			patterns[i].set(data.iso()[i], data.symmetry()[i]);
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Print each pattern
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print pattern
		data.printLabel(i);
		patterns[i].print("stdout", broaden);
		
		// If diffraction file was supplied, then print comparison
		if (data.diffraction().isSet())
		{
			Output::method(STANDARD);
			Output::newline();
			Output::print("R factor: ");
			Output::print(match[i], 6);
			Output::method(origMethod);
		}
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::perturb(Storage& data, const Function& function)
 *
 * Randomly perturb the structure
 */

void Launcher::perturb(Storage& data, const Function& function)
{
	
	// Loop over arguments
	int i;
	double min = -1;
	double max = -1;
	bool pertBasis = false;
	bool pertAtoms = false;
	bool preserveSymmetry = false;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found a number
		if (Language::isNumber(function.arguments()[i]))
		{
			if (min == -1)
				min = max = atof(function.arguments()[i].array());
			else if (min == max)
				max = atof(function.arguments()[i].array());
			else
			{
				Output::newline(ERROR);
				Output::print("Expecting 0, 1, or 2 numbers passed to perturb function");
				Output::quit();
			}
		}
		
		// Found basis
		else if ((function.arguments()[i].equal("basis", false, 3)) || \
			(function.arguments()[i].equal("lattice", false, 3)))
			pertBasis = true;
		
		// Found atoms
		else if ((function.arguments()[i].equal("atoms", false, 4)) || \
			(function.arguments()[i].equal("coordinates", false, 4)) || \
			(function.arguments()[i].equal("positions", false, 3)))
			pertAtoms = true;
		
		// Found symmetry
		else if (function.arguments()[i].equal("symmetry", false, 3))
			preserveSymmetry = true;
	}
	
	// Set what should be perturbed
	if ((!pertBasis) && (!pertAtoms))
		pertBasis = pertAtoms = true;
	
	// Magnitude was not set
	if (min == -1)
		min = max = 0.25;
	
	// Update the symmetry if needed 
	if (preserveSymmetry)
		updateSymmetry(data);
	else
		data.updateSymmetry().fill(true);
	
	// Output
	Output::newline();
	Output::print("Perturbing ");
	if ((!pertAtoms) && (pertBasis))
		Output::print("basis vectors ");
	else if ((pertAtoms) && (!pertBasis))
		Output::print("atomic coordinates ");
	else
		Output::print("basis vectors and atomic coordinates");
	Output::increase();
	
	// Create history tag
	Word tag = " > perturbed ";
	if (pertBasis)
	{
		tag += "basis";
		if (pertAtoms)
			tag += " and ";
	}
	if (pertAtoms)
		tag += "atomic positions";
	
	// Loop over structures
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Perturbing structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Perturb structure
		if (pertBasis)
		{
			if (preserveSymmetry)
				RandomStructure::perturbBasis(data.iso()[i], data.symmetry()[i], min, max, data.random());
			else
				RandomStructure::perturbBasis(data.iso()[i], min, max, data.random());
		}
		if (pertAtoms)
		{
			if (preserveSymmetry)
				RandomStructure::perturbAtoms(data.iso()[i], data.symmetry()[i], min, max, data.random());
			else
				RandomStructure::perturbAtoms(data.iso()[i], min, max, data.random());
		}
		
		// Add change to comment
		data.history()[i] += tag;
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::optimize(Storage& data, const Function& function)
 *
 * Optimize structure
 */

void Launcher::optimize(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Optimizing structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Setup ga
	GAPredict ga;
	ga.numScreens(Settings::value<int>(GAOPT_SCREENNUM));
	ga.numSimulations(Settings::value<int>(GAOPT_NUMSIM));
	ga.maxCrossoverLoops(Settings::value<int>(RANDSTR_MAXLOOPS));
	ga.minBondFraction(Settings::value<double>(RANDSTR_MINBOND));
	ga.wyckoffBias(Settings::value<double>(WYCKOFFBIAS));
	ga.cellMutationProbability(Settings::value<double>(GAOPT_CELLMUTPROB));
	ga.positionMutationProbability(Settings::value<double>(GAOPT_POSMUTPROB));
	ga.wyckoffMutationProbability(Settings::value<double>(GAOPT_WYCKMUTPROB));
	ga.metricToOptimize(Settings::value<GAPredictMetric>(GAOPT_METRICTOOPT));
	ga.metricToScreen(Settings::value<GAPredictMetric>(GAOPT_SCREENMETHOD));
	ga.populationSize(Settings::value<int>(GAOPT_POPSIZE));
	ga.convergeOver(Settings::value<int>(GAOPT_CONVERGEOVER));
	ga.maxGenerations(Settings::value<int>(GAOPT_MAXGENS));
	ga.numToKeep(Settings::value<int>(GAOPT_NUMTOKEEP));
	ga.selectionMethod(Settings::value<GASelectionMethod>(GAOPT_SELECTION));
	ga.energyTolerance(Settings::value<double>(GAOPT_ENERGYTOL));
	ga.diffractionTolerance(Settings::value<double>(GAOPT_DIFFRACTIONTOL));
	
	// Loop over function arguments to check for a metric to optimize
	int i;
	GAPredictMetric curMetric;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		curMetric = GAPredict::metric(function.arguments()[i]);
		if (curMetric != GAPM_UNKNOWN)
			ga.metricToOptimize(curMetric);
	}
	
	// Get current directory
	Word origDir;
	Directory::current(origDir);
	
	// Loop over structures
	Word directory;
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Optimizing structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Set directory for run
		directory = Directory::makePath(origDir, Word("OPT_STR_"));
		directory += Language::numberToWord(data.id()[i]);
		Directory::create(directory, true);
		Directory::change(directory);
		
		// Optimize structure
		ga.run(data.iso()[i], data.random(), &data.potential(), &data.diffraction());
		
		// Add change to comment
		data.history()[i] += " > optimized";
		
		// Return to original directory
		Directory::change(origDir);
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
		
		// Save that structure is set
		data.needsGeneration()[i] = false;
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::relax(Storage& data, const Function& function)
 *
 * Relax the structure
 */

void Launcher::relax(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Relaxing structure");
	if (data.iso().length() != 1)
		Output::print("s");
	Output::increase();
	
	// Check if symmetry should be used
	bool useSymmetry = data.potential().usesSymmetry();
	
	// Check if the symmetry needs to be updated
	if (useSymmetry)
		updateSymmetry(data);
	
	// Loop over structures
	int i;
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Relaxing structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Relax structure
		if (useSymmetry)
			data.potential().relax(data.iso()[i], data.symmetry()[i], 0, 0, true);
		else
			data.potential().relax(data.iso()[i], 0, 0, true);
		
		// Add change to comment
		data.history()[i] += " > relaxed";
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::phonons(Storage& data, const Function& function)
 *
 * Run phonon calculations
 */

void Launcher::phonons(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Calculating phonons");
	Output::increase();
	
	// Force constant file has already been set
	int i;
	OList<Phonons> phon;
	if (data.phonons().isSet())
		phon += data.phonons();
	
	// Get phonons data from structures
	else
	{
		
		// Update the symmetry
		updateSymmetry(data);
		
		// Loop over structures
		phon.length(data.iso().length());
		for (i = 0; i < data.iso().length(); ++i)
		{

			// Output if there is more than one structure
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Calculating force constants for structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Get phonons
			phon[i].generateForceConstants(data.iso()[i], data.symmetry()[i], data.potential(), \
				Language::numberToWord(data.id()[i]));
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
				Output::decrease();
		}
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::kmc(Storage& data, const Function& function)
 *
 * Run KMC simulation
 */

void Launcher::kmc(Storage& data, const Function& function)
{
	
	// KMC file has already been set
	if (data.kmc().length() > 0)
	{
		
		// Output
		Output::newline();
		Output::print("Preparing KMC simulation");
		Output::increase();
		
		// Pass kmc file to object and run if ready
		KMC kmc;
		if (kmc.set(data.kmc()))
		{
			
			// Set the run properties
			kmc.convergence(Settings::value<double>(KMC_CONVERGENCE));
			kmc.jumpsPerAtom(Settings::value<double>(KMC_JUMPSPERATOM));
			
			// Run the simulation
			Output::decrease();
			Output::newline();
			Output::print("Running KMC simulation");
			Output::increase();
			kmc.run(data.random());
			
			// Print results
			kmc.print("stdout");
		}
		
		// Output
		Output::decrease();
	}
	
	// Generate kmc input file
	else
	{
		
		// Output
		Output::newline();
		Output::print("Initializing KMC simulation");
		Output::increase();
		
		// Initialize the structure format
		StructureFormat format = StructureIO::structureFormat(function.arguments());
		if (format == SF_UNKNOWN)
			format = Settings::value<StructureFormat>(STRUCTURE_FORMAT);
		
		// Loop over function arguments to get settings
		int i;
		bool getInt = false;
		Element element;
		Element tempElement;
		for (i = 0; i < function.arguments().length(); ++i)
		{
			
			// Check if an element
			tempElement = Element::find(function.arguments()[i], false, false);
			if (tempElement.number() != 0)
				element = tempElement;
			
			// Check if interstitial
			if (function.arguments()[i].equal("interstitial", false, 3))
				getInt = true;
		}
		
		// Element was not passed
		if (element.number() == 0)
		{
			Output::newline(ERROR);
			Output::print("Must pass an element to KMC function");
			Output::quit();
		}
		
		// Save current directory
		Word origDir;
		Directory::current(origDir);
		
		// Loop over structures
		Word curDir;
		StructureFormat curFormat;
		for (i = 0; i < data.iso().length(); ++i)
		{

			// Output if there is more than one structure
			if (data.iso().length() > 1)
			{
				Output::newline();
				Output::print("Setting data for KMC simulation of structure ");
				Output::print(data.id()[i]);
				Output::increase();
			}
			
			// Make directory
			curDir = "KMC_";
			curDir += Language::numberToWord(i+1);
			Directory::create(curDir, true);
			Directory::change(curDir);
			
			// Setup KMC simulation
			curFormat = (format == SF_AUTO) ? data.format()[i] : format;
			KMC::generateJumps(data.iso()[i], element, getInt, Settings::value<double>(MINIMAGEDISTANCE), \
				Settings::value<double>(TOLERANCE), curFormat, Settings::value<double>(MAXJUMPDISTANCE));
			
			// Return to original directory
			Directory::change(origDir);
			
			// Output if there is more than one structure
			if (data.iso().length() > 1)
				Output::decrease();
		}
		
		// Output
		Output::decrease();
	}
}



/* void Launcher::neighbors(Storage& data, const Function& function)
 *
 * Calculate nearest neighbors
 */

void Launcher::neighbors(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Calculating nearest neighbor distances");
	Output::increase();
	
	// Loop over structures to calculate distances
	int i, j, k, m;
	int index;
	List<Atom*>::D2 atoms(data.iso().length());
	List<double>::D3 distances(data.iso().length());
	List<Word>::D3 elements(data.iso().length());
	OList<Vector3D >::D3 cells(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Searching neighbors in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get atoms
		atoms[i] = getAtoms(data.iso()[i], function);
		
		// Loop over atoms and calculate distances
		cells[i].length(atoms[i].length());
		distances[i].length(atoms[i].length());
		elements[i].length(atoms[i].length());
		for (j = 0; j < atoms[i].length(); ++j)
		{
			
			// Loop over atoms in the structure and calculate distances
			cells[i][j].length(data.iso()[i].numAtoms());
			distances[i][j].length(data.iso()[i].numAtoms());
			elements[i][j].length(data.iso()[i].numAtoms());
			for (k = 0; k < data.iso()[i].atoms().length(); ++k)
			{
				for (m = 0; m < data.iso()[i].atoms()[k].length(); ++m)
				{
					index = data.iso()[i].atoms()[k][m].atomNumber();
					elements[i][j][index] = data.iso()[i].atoms()[k][m].element().symbol();
					if (index == atoms[i][j]->atomNumber())
						distances[i][j][index] = data.iso()[i].basis().secondDistance(atoms[i][j]->fractional(), \
							FRACTIONAL, atoms[i][j]->fractional(), FRACTIONAL, &cells[i][j][index]);
					else
						distances[i][j][index] = data.iso()[i].basis().distance(atoms[i][j]->fractional(), \
							FRACTIONAL, data.iso()[i].atoms()[k][m].fractional(), FRACTIONAL, &cells[i][j][index]);
				}
			}
		}
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Variables for printing
	Output message;
	message.addSpaces(false);
	List<PrintAlign> align(12, RIGHT);
	
	// Print neighbors lists
	int minIndex;
	double minDistance;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print label
		data.printLabel(i);
		Output::method(STANDARD);
		
		// Loop over atoms
		for (j = 0; j < atoms[i].length(); ++j)
		{
			
			// Find the nearest atom
			minIndex = 0;
			minDistance = distances[i][j][0];
			for (k = 1; k < distances[i][j].length(); ++k)
			{
				if (distances[i][j][k] < minDistance)
				{
					minIndex = k;
					minDistance = distances[i][j][k];
				}
			}
			
			// Save distances as a message
			message.clear();
			message.addLines(distances[i][j].length() + 1);
			for (k = 0; k < distances[i][j].length(); ++k)
			{
				message.addLine();
				message.addWords(12);
				message.add("    ");
				message.add(Word("Atom ") + Language::numberToWord(k+1));
				message.add(Word(" (") + elements[i][j][k] + "): ");
				message.add(distances[i][j][k], 4);
				message.add(" (");
				for (m = 0; m < 3; ++m)
				{
					message.add((int)cells[i][j][k][m]);
					if (m != 2)
						message.add(", ");
				}
				message.add(")");
				if (k == minIndex)
					message.add(" < min");
			}
			
			// Print distances
			Output::newline();
			if (j > 0)
				Output::newline();
			Output::print("Distances from atom ");
			Output::print(atoms[i][j]->atomNumber() + 1);
			Output::newline();
			Output::print(message, align);
		}
		
		// Reset method
		Output::method(origMethod);
	}
	
	// Output
	Output::decrease(); 
}



/* void Launcher::shells(Storage& data, const Function& function)
 *
 * Calculate nearest neighbor shells
 */

void Launcher::shells(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Calculating nearest neighbor shells");
	Output::increase();
	
	// Get the distance out to which to calculate shells
	int i;
	bool details = false;
	double maxDistance = 3;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (Language::isDecimal(function.arguments()[i]))
			maxDistance = atof(function.arguments()[i].array());
		if (function.arguments()[i].equal("det", false, 4))
			details = true;
	}
	
	// Loop over structures to calculate distances
	int j;
	List<Atom*>::D2 atoms(data.iso().length());
	OList<Atom>::D4 shells(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Searching for shells in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get atoms
		atoms[i] = getAtoms(data.iso()[i], function, false);
		
		// Calculate shells
		shells[i].length(atoms[i].length());
		for (j = 0; j < atoms[i].length(); ++j)
			shells[i][j] = data.iso()[i].shells(atoms[i][j], maxDistance, Settings::value<double>(TOLERANCE));
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Print results
	int k, m, n;
	double avgDis;
	Output message;
	List<int> counts;
	List<double> distances;
	OList<Element> elements;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print label
		data.printLabel(i);
		Output::method(STANDARD);
		
		// Loop over atoms
		for (j = 0; j < atoms[i].length(); ++j)
		{
			
			// Output
			Output::newline();
			if (j > 0)
				Output::newline();
			Output::print("Shells for atom ");
			Output::print(atoms[i][j]->atomNumber() + 1);
			Output::print(" (");
			Output::print(atoms[i][j]->element().symbol());
			Output::print("):");
			
			// Loop over shells
			for (k = 0; k < shells[i][j].length(); ++k)
			{
				
				// Get distances to each atom in current shell
				distances.length(shells[i][j][k].length());
				for (m = 0; m < shells[i][j][k].length(); ++m)
					distances[m] = data.iso()[i].basis().absoluteDistance(atoms[i][j]->fractional(), FRACTIONAL, \
						shells[i][j][k][m].fractional(), FRACTIONAL);
				
				// Get average distance for atoms in shell
				avgDis = 0;
				for (m = 0; m < distances.length(); ++m)
					avgDis += distances[m];
				if (distances.length() > 0)
					avgDis /= distances.length();
				
				// Print standard information
				Output::newline();
				Output::tab();
				Output::print("Shell ");
				Output::print(k + 1);
				Output::print(" (");
				Output::print(avgDis, 3);
				Output::print(" Ang): ");
				
				// Get number of atoms per element
				counts.clear();
				elements.clear();
				for (m = 0; m < shells[i][j][k].length(); ++m)
				{
					for (n = 0; n < elements.length(); ++n)
					{
						if (elements[n] == shells[i][j][k][m].element())
						{
							++counts[n];
							break;
						}
					}
					if (n >= elements.length())
					{
						counts += 1;
						elements += shells[i][j][k][m].element();
					}
				}
				
				// Print elements
				for (m = 0; m < elements.length(); ++m)
				{
					if (m != 0)
						Output::print(", ");
					Output::print(counts[m]);
					Output::print(" ");
					Output::print(elements[m].symbol());
				}
				
				// Print details if needed
				if (details)
				{
					message.clear();
					message.addSpaces(false);
					message.addLines(shells[i][j][k].length());
					for (m = 0; m < shells[i][j][k].length(); ++m)
					{
						message.addLine();
						message.addTab();
						message.addTab();
						message.add("Atom ");
						message.add(shells[i][j][k][m].atomNumber() + 1);
						message.add(Word(" (") + shells[i][j][k][m].element().symbol() + ")");
						message.add(" ");
						message.add(distances[m], 3);
						message.add(":");
						for (n = 0; n < 3; ++n)
						{
							message.add(" ");
							message.add(shells[i][j][k][m].fractional()[n]);
						}
					}
					Output::newline();
					Output::print(message, RIGHT);
				}
			}
		}
		
		// Reset method
		Output::method(origMethod);
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::coordination(Storage& data, const Function& function)
 * 
 * Get the coordination of atoms in the structure
 */

void Launcher::coordination(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Calculating coordination numbers");
	Output::increase();
	
	// Loop over structures
	int i, j;
	List<Atom*>::D2 atoms(data.iso().length());
	List<Atom*>::D3 coordinations(data.iso().length());
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Calculating coordinations in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get atoms
		atoms[i] = getAtoms(data.iso()[i], function);
		
		// Loop over atoms to calculate
		coordinations[i].length(atoms[i].length());
		for (j = 0; j < atoms[i].length(); ++j)
			coordinations[i][j] = data.iso()[i].coordination(atoms[i][j]);
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Variables for printing
	Output message;
	message.addSpaces(false);
	List<PrintAlign> align(5, RIGHT);
	align.last() = LEFT;
	
	// Print coordinations
	int k;
	Word curString;
	PrintMethod origMethod = Output::method();
	for (i = 0; i < data.iso().length(); ++i)
	{
		
		// Print label
		data.printLabel(i);
		Output::method(STANDARD);
		
		// Loop over atoms
		message.clear();
		message.addLines(atoms[i].length());
		for (j = 0; j < atoms[i].length(); ++j)
		{
			
			// Save coordination
			message.addLine();
			message.add("Coordination for atom ");
			message.add(atoms[i][j]->atomNumber() + 1);
			message.add(": ");
			message.add(coordinations[i][j].length());
			if (coordinations[i][j].length())
			{
				curString = " (atom";
				if (coordinations[i][j].length() != 1)
					curString += 's';
				for (k = 0; k < coordinations[i][j].length(); ++k)
				{
					curString += " ";
					curString += Language::numberToWord(coordinations[i][j][k]->atomNumber() + 1);
					if ((coordinations[i][j].length() > 2) && (k != coordinations[i][j].length() - 1))
						curString += ",";
					if (k == coordinations[i][j].length() - 2)
						curString += " and";
				}
				curString += ")";
				message.add(curString);
			}
		}
		
		// Print coordinations
		Output::newline();
		Output::print(message, align);
		
		// Reset method
		Output::method(origMethod);
	}
	
	
	// Output
	Output::decrease();
}



/* void Launcher::interstitial(Storage& data, const Function& function)
 *
 * Search for interstitial sites
 */

void Launcher::interstitial(Storage& data, const Function& function)
{
	
	// Output
	Output::newline();
	Output::print("Searching for interstitial sites in the structure");
	if (data.iso().length() > 1)
		Output::print("s");
	Output::increase();
	
	// Get arguments passed to function
	int i;
	bool expand = false;
	int density = 25;
	double scale = 0.25;
	Element element;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found an element
		if (Element::isElement(function.arguments()[i], false))
			element = Element::find(function.arguments()[i], false);
		
		// Found an integer
		else if (Language::isInteger(function.arguments()[i]))
			density = atoi(function.arguments()[i].array());
		
		// Found a float
		else if (Language::isDecimal(function.arguments()[i]))
			scale = atof(function.arguments()[i].array());
		
		// Found expand
		else if (function.arguments()[i].equal("expand", false, 3))
			expand = true;
	}
	
	// Set element
	if (element.number() == 0)
	{
		element = Element::find(Language::numberToWord(1), true);
		Output::newline(WARNING);
		Output::print("No element was passed to interstitial function so assuming ");
		Output::print(element.symbol());
	}
	
	// Loop over structures
	OList<Interstitial> interstitial(data.iso().length());
	for (i = 0; i < data.iso().length(); i++)
	{
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
		{
			Output::newline();
			Output::print("Searching for sites in structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Get sites
		interstitial[i].evaluate(data.iso()[i], density, Settings::value<double>(CLUSTERTOL), scale);
	//	interstitial[i].voronoi(data.iso()[i], Settings::value<double>(CLUSTERTOL));
		
		// Output if there is more than one structure
		if (data.iso().length() > 1)
			Output::decrease();
	}
	
	// Output
	Output::newline();
	Output::print("Creating structures for each interstitial site");
	Output::increase();
	
	// Create each structure
	int j;
	int origLen = data.iso().length();
	Word history;
	Atom* atom;
	List<int> indices;
	for (i = 0; i < origLen; ++i)
	{
		
		// Output
		if (origLen > 0)
		{
			Output::newline();
			Output::print("Adding interstitials to structure ");
			Output::print(data.id()[i]);
			Output::increase();
		}
		
		// Skip if empty
		if (interstitial[i].sites().length() == 0)
			continue;
		
		// Create new structures
		indices.length(interstitial[i].sites().length() - 1);
		for (j = 1; j < interstitial[i].sites().length(); ++j)
		{
			data.addISO();
			data.iso().last() = data.iso()[i];
			data.baseName().last() = data.baseName()[i];
			data.history().last() = data.history()[i];
			data.format().last() = data.format()[i];
			indices[j-1] = data.iso().length() - 1;
		}
		
		// Output
		Output::newline();
		Output::print("Adding ");
		Output::print(element.symbol());
		Output::print(" at ");
		Output::print(interstitial[i].sites()[0], 8, true);
		
		// Add first position
		atom = data.iso()[i].addAtom(element);
		atom->fractional(interstitial[i].sites()[0]);
		atom->isInterstitial(true);
		if (expand)
			Symmetry::addAtom(data.iso()[i], *atom, data.symmetry()[i].operations(), \
				Settings::value<double>(CLUSTERTOL));
		data.history()[i] += Word(" > Interstitial ") + Language::numberToWord(1);
		
		// Add additional positions
		for (j = 0; j < indices.length(); ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Adding ");
			Output::print(element.symbol());
			Output::print(" at ");
			Output::print(interstitial[i].sites()[j+1], 8, true);
			
			// Add atom
			atom = data.iso()[indices[j]].addAtom(element);
			atom->fractional(interstitial[i].sites()[j+1]);
			atom->isInterstitial(true);
			if (expand)
				Symmetry::addAtom(data.iso()[indices[j]], *atom, data.symmetry()[i].operations(), \
					Settings::value<double>(CLUSTERTOL));
			data.history()[indices[j]] += Word(" > Interstitial ") + Language::numberToWord(j+2);
		}
		
		// Output
		if (origLen > 0)
			Output::decrease();
	}
	
	// Output
	Output::decrease();
	
	// Print each set of sites
	PrintMethod origMethod = Output::method();
	for (i = 0; i < origLen; ++i)
	{
		
		// Print label
		data.printLabel(i, origLen);
		Output::method(STANDARD);
		
		// Print sites
		Output::newline();
		Output::print("Found ");
		Output::print(interstitial[i].sites().length());
		Output::print(" unique site");
		if (interstitial[i].sites().length() != 1)
			Output::print("s");
		for (j = 0; j < interstitial[i].sites().length(); ++j)
		{
			Output::newline();
			Output::tab();
			Output::print(interstitial[i].sites()[j], 8, false);
		}
		
		// Reset method
		Output::method(origMethod);
	}
	
	// Output
	Output::decrease();
}



/* void Launcher::compare(Storage& data, const Function& function)
 *
 * Compare structures that have been passed
 */

void Launcher::compare(Storage& data, const Function& function)
{
	
	// Return if there is less than two structures
	if (data.iso().length() < 2)
	{
		Output::newline(WARNING);
		Output::print("Must supply at least two structures to use comparison function");
		return;
	}
	
	// Output
	Output::newline();
	Output::print("Comparing structures to determine which are equivalent");
	Output::increase();
	
	// Loop over arguments to get settings
	int i;
	bool compareVolume = false;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		if (function.arguments()[i].equal("volume", false, 3))
			compareVolume = true;
	}
	
	// Loop over pairs of structures and compare them
	int j;
	Linked<int> origID;
	Linked<int> compID;
	Linked<bool> areEqual;
	for (i = 0; i < data.iso().length(); ++i)
	{
		for (j = i + 1; j < data.iso().length(); ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Comparing structures ");
			Output::print(data.id()[i]);
			Output::print(" and ");
			Output::print(data.id()[j]);
			Output::increase();
			
			// Compare structures
			origID += data.id()[i];
			compID += data.id()[j];
			areEqual += data.iso()[i].equivalent(data.iso()[j], Settings::value<double>(TOLERANCE), compareVolume);
			
			// Output
			Output::decrease();
		}
	}
	
	// Set output method
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	// Print results
	Output message;
	Linked<int>::iterator itOrigID = origID.begin();
	Linked<int>::iterator itCompID = compID.begin();
	Linked<bool>::iterator itEqual = areEqual.begin();
	for (; itOrigID != origID.end(); ++itOrigID, ++itCompID, ++itEqual)
	{
		
		// Print whether structures are the same
		Output::newline();
		Output::print("Structures ");
		Output::print(*itOrigID);
		Output::print(" and ");
		Output::print(*itCompID);
		Output::print(" are");
		if (*itEqual == false)
			Output::print(" not");
		Output::print(" the same");
	}
	
	// Reset method
	Output::method(origMethod);
	
	// Output
	Output::decrease();
}



/* bool Launcher::runSettings(const Functions& functions)
 *
 * Check whether to run the settings console
 */

bool Launcher::runSettings(const Functions& functions)
{
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_SETTINGS)
		{
//			Settings::Console::launch();
			return true;
		}
	}
	return false;
}



/* bool Launcher::runHelp(const Functions& functions)
 *
 * Check whether to show help
 */

bool Launcher::runHelp(const Functions& functions)
{
	for (int i = 0; i < functions.length(); ++i)
	{
		if (functions[i].keyword() == KEY_HELP)
		{
			Help::run(functions[i].arguments());
			return true;
		}
	}
	return false;
}



/* List<Atom*> Launcher::getAtoms(const ISO& iso, const Function& function, bool allowPositions, bool autoPopAll)
 *
 * Return a list of atoms that are referenced in a function call
 */

List<Atom*> Launcher::getAtoms(const ISO& iso, const Function& function, bool allowPositions, bool autoPopAll)
{
	
	// Loop over arguments to get data out
	int i;
	int index;
	int posIndex = 0;
	Atom* atom;
	Element curElement;
	Element tempElement;
	Vector3D pos;
	List<Atom*> res;
	for (i = 0; i < function.arguments().length(); ++i)
	{
		
		// Found an element
		tempElement = Element::find(function.arguments()[i], false, false);
		if (tempElement.number() != 0)
		{
			
			// An element is already saved
			if (curElement.number() != 0)
				res += iso.atoms(curElement);
			
			// Save element
			curElement = tempElement;
		}
		
		// Found an integer
		else if (Language::isInteger(function.arguments()[i]))
		{
			
			// Get the atom
			index = atoi(function.arguments()[i].array()) - 1;
			if (curElement.number() == 0)
				atom = iso.atom(index);
			else
				atom = iso.atom(curElement, index);
			
			// Atom is not in the structure
			if (atom == 0)
			{
				Output::newline(ERROR);
				Output::print("Atom ");
				if (curElement.number() != 0)
				{
					Output::print(curElement.symbol());
					Output::print(" ");
				}
				Output::print(index + 1);
				Output::print(" is not in the structure");
				Output::quit();
			}
			
			// Add current atom
			res += atom;
			
			// Clear data
			index = -1;
			curElement.clear();
		}
		
		// Found a float
		else if ((Language::isDecimal(function.arguments()[i])) && (allowPositions))
		{
			
			// Save value
			pos[posIndex++] = atof(function.arguments()[i].array());
			
			// Reached end of position
			if (posIndex == 3)
			{
				
				// Get atom
				atom = iso.atom(pos, Settings::value<double>(TOLERANCE));
				if (atom == 0)
				{
					Output::newline(ERROR);
					Output::print("No atom found at ");
					Output::print(pos, 8, true);
					Output::quit();
				}
				
				// Save atom and reset
				res += atom;
				posIndex = 0;
			}
		}
	}
	
	// Position was not finished
	if (posIndex != 0)
	{
		Output::newline(ERROR);
		Output::print("A partially declared position was found");
		Output::quit();
	}
	
	// An element is saved
	if (curElement.number() != 0)
		res += iso.atoms(curElement);
	
	// Save all atoms if none were set
	if ((!res.length()) && (autoPopAll == true))
	{
		int j;
		res.length(iso.numAtoms());
		for (i = 0, index = 0; i < iso.atoms().length(); ++i)
		{
			for (j = 0; j < iso.atoms()[i].length(); ++j)
				res[index++] = &iso.atoms()[i][j];
		}
	}
	
	// Sort the atoms by index
	sortAtoms(res, 0, res.length() - 1);
	
	// Return list
	return res;
}



/* void Launcher::sortAtoms(List<Atom*>& atoms, int left, int right)
 *
 * Sort atoms by index
 */

void Launcher::sortAtoms(List<Atom*>& atoms, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	int pivot = atoms[pivotIndex]->atomNumber();
	
	// Move pivot to end
	atoms.swap(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (atoms[i]->atomNumber() < pivot)
		{
			atoms.swap(i, newPivotIndex);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	atoms.swap(newPivotIndex, right);
	
	// Recursive calls to next sorts
	sortAtoms(atoms, left, newPivotIndex - 1);
	sortAtoms(atoms, newPivotIndex + 1, right);
}
