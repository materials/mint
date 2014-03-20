/* Copyright 2011-2014 Kyle Michel, Logan Ward
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



#ifndef LAUNCHER_H
#define LAUNCHER_H



#include "iso.h"
#include "structureIO.h"
#include "symmetry.h"
#include "potential.h"
#include "phonons.h"
#include "kmc.h"
#include "diffraction.h"
#include "random.h"
#include "text.h"
#include "list.h"



// Types of runs
enum Keyword {KEY_NONE, KEY_SETTINGS, KEY_HELP, KEY_NUMPROCS, KEY_OUTPUT, KEY_TIME, KEY_TOLERANCE, KEY_PRINT, \
	KEY_NAME, KEY_FIX, KEY_REMOVE, KEY_NEIGHBORS, KEY_SHELLS, KEY_COORDINATION, KEY_REDUCED, KEY_PRIMITIVE, \
	KEY_CONVENTIONAL, KEY_IDEAL, KEY_SHIFT, KEY_TRANSFORM, KEY_ROTATE, KEY_SYMMETRY, KEY_UNIQUE, KEY_EQUIVALENT, \
	KEY_ABOUT, KEY_POINTGROUP, KEY_SPACEGROUP, KEY_REFINE, KEY_ENERGY, KEY_FORCES, KEY_PHONONS, KEY_KMC, \
	KEY_DIFFRACTION, KEY_PERTURB, KEY_RELAX, KEY_OPT, KEY_INTERSTITIAL, KEY_COMPARE};



// Class to store data during run
class Storage
{

	// Structure variables
	int _curID;
	List<int> _id;
	OList<ISO> _iso;
	List<bool> _needsGeneration;
	
	// Symmetry variables
	OList<Symmetry> _symmetry;
	List<bool> _updateSymmetry;
	
	// Write variables
	List<StructureFormat> _format;
	OList<Word> _baseName;
	OList<Word> _history;
	
	// Other variables
	Potential _potential;
	Phonons _phonons;
	Text _kmc;
	Diffraction _diffraction;
	Random _random;
	
public:
	
	// Constructor
	Storage()	{ _curID = 0; }
	
	// Add/remove structure
	void addISO();
	void removeISO(int index);
	
	// Print functions
	void printLabel(int index, int length = -1);
	
	// Structure access functions
	List<int>& id()					{ return _id; }
	OList<ISO>& iso()				{ return _iso; }
	List<bool>& needsGeneration()	{ return _needsGeneration; }
	
	// Symmetry access functions
	OList<Symmetry>& symmetry()		{ return _symmetry; }
	List<bool>& updateSymmetry()	{ return _updateSymmetry; }
	
	// Write access functions
	List<StructureFormat>& format()	{ return _format; }
	OList<Word>& baseName()			{ return _baseName; }
	OList<Word>& history()			{ return _history; }
	
	// Other access functions
	Potential& potential()		{ return _potential; }
	Phonons& phonons()			{ return _phonons; }
	Text& kmc()					{ return _kmc; }
	Diffraction& diffraction()	{ return _diffraction; }
	Random& random()			{ return _random; }
};



// Class to store a single function and its arguments
class Function
{
	
	// Variables
	Keyword _keyword;
	Words _arguments;
	
public:
	
	// Setup functions
	void keyword(Keyword input)			{ _keyword = input; }
	void addArgument(const Word& input)	{ _arguments += input; }
	void removeArgument(int index)		{ _arguments.remove(index); }
	
	// Access functions
	Keyword keyword() const				{ return _keyword; }
	const Words& arguments() const		{ return _arguments; }
};



// Class to store functions for Mint
class Launcher
{
	
	// Type definitions
	typedef OList<Function> Functions;
	
	// Startup functions
	static Words getArguments(int argc, char** argv);
	static Functions getFunctions(const Words& arguments);
	static Keyword getKeyword(const Word& argument);
	static void runSetup(Storage& data, Functions& functions, bool& keepFree);
	static void getCommandLineSettings(Functions& functions, bool& keepFree);
	static void setupJobSize(Functions& functions);
	static void setupOutput(Functions& functions);
	static void setupTolerance(Functions& functions);
	static void setupTime(Functions& functions);
	static void setupPrint(Functions& functions, bool& keepFree);
	static void fixCellParams(Storage& data, Functions& functions);
	
	// General functions
	static void output(const Function& function);
	static void tolerance(const Function& function);
	static void generateStructure(Storage& data, bool keepFree = false);
	static void printStructure(Storage& data, const Functions& functions);
	
	// Change cell
	static void removeAtoms(Storage& data, const Function& function);
	static void reduced(Storage& data);
	static void primitive(Storage& data);
	static void conventional(Storage& data);
	static void ideal(Storage& data, const Function& function);
	static void shift(Storage& data, const Function& function);
	static void transform(Storage& data, const Function& function);
	static void rotate(Storage& data, const Function& function);
	
	// Symmetry functions
	static void updateSymmetry(Storage& data);
	static void symmetry(Storage& data, const Function& function);
	static void unique(Storage& data, const Function& function);
	static void about(Storage& data, const Function& function);
	static void aboutStructure(Storage& data);
	static void aboutAtoms(Storage& data, const Function& function);
	static void pointGroup(Storage& data, const Function& function);
	static void spaceGroup(Storage& data, const Function& function);
	static void refine(Storage& data, const Function& function);
	
	// Potential evaluation functions
	static void energy(Storage& data);
	static void forces(Storage& data);
	static void phonons(Storage& data, const Function& function);
	
	// Diffusion functions
	static void kmc(Storage& data, const Function& function);
	
	// Characterization functions
	static void diffraction(Storage& data, const Function& function);
	
	// Random structure and perturbation
	static void perturb(Storage& data, const Function& function);
	
	// Optimize functions
	static void relax(Storage& data, const Function& function);
	static void optimize(Storage& data, const Function& function);
	
	// Other functions
	static void neighbors(Storage& data, const Function& function);
	static void shells(Storage& data, const Function& function);
	static void coordination(Storage& data, const Function& function);
	static void interstitial(Storage& data, const Function& function);
	
	// Compare two structures
	static void compare(Storage& data, const Function& function);
	
	// Alternate runs
	static bool runSettings(const Functions& functions);
	static bool runHelp(const Functions& functions);
	
	// Helper functions
	static bool fileNameUsed(const OList<Word>& fileNames, const Word& curName, int numToCheck);
	static List<Atom*> getAtoms(const ISO& iso, const Function& function, bool allowPositions = true, \
		bool autoPopAll = true);
	static void sortAtoms(List<Atom*>& atoms, int left, int right);
	
public:
	
	// Main function
	static void start(int argc, char** argv);	
};



#endif
