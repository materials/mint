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



#include "potential.h"
#include "locPotential.h"
#include "extPotential.h"
#include "language.h"



/* void Potential::set(const Text& text)
 *
 * Set potential information
 */

void Potential::set(const Text& text)
{
	
	// Output
	Output::newline();
	Output::print("Setting potential");
	Output::increase();
	
	// Get parsed input
	OList<Text> data = parseInput(text);
	
	// Return if no potential was found
	if (!data.length())
	{
		Output::newline();
		Output::print("File contains no potential information");
		Output::decrease();
		return;
	}
	
	// Loop over parsed input to check whether reference states will be used
	int i;
	for (i = data.length() - 1; i >= 0; --i)
	{
		if (data[i][0][0].equal("references", false, 3))
		{
			_useReferences = true;
			setReferenceData(data[i]);
			data.remove(i);
		}
	}
	
	// Loop over parsed input to figure out what kind of potential to use
	bool useVasp = false;
	bool useLocal = false;
	bool useEspresso = false;
	PotentialType firstLocal = PT_UNKNOWN;
	List<PotentialType> types(data.length());
	for (i = 0; i < data.length(); ++i)
	{
		types[i] = potentialType(data[i][0][0]);
		if (types[i] == PT_VASP)
			useVasp = true;
		else if (types[i] == PT_QE)
			useEspresso = true;
		else if (types[i] != PT_UNKNOWN)
		{
			useLocal = true;
			firstLocal = (firstLocal == PT_UNKNOWN) ? types[i] : firstLocal;
		}
		else
		{
			Output::newline(ERROR);
			Output::print("Internal error in setting potential");
			Output::quit();
		}
	}
	
	// Set local potential
	if (useLocal)
	{
		
		// Make sure Vasp and QE potentials were not set
		if (useVasp)
			error(firstLocal, PT_VASP);
		if (useEspresso)
			error(firstLocal, PT_QE);
		
		// Set potential
		_ipo = new LocalPotential;
	}
	
	// Set Vasp potential
	if (useVasp)
	{
		
		// Make sure local and QE potentials were not set
		if (useLocal)
			error(PT_VASP, firstLocal);
		if (useEspresso)
			error(PT_QE, PT_VASP);
		
		// Set potential
		_ipo = new VaspPot;
	}
	
	// Set quantum espresso file
	if (useEspresso)
	{
		
		// Make sure local and Vasp potentials were not set
		if (useLocal)
			error(PT_QE, firstLocal);
		if (useVasp)
			error(PT_QE, PT_VASP);
		
		// Set potential
		_ipo = new QEPot;
	}
	
	// Add first potential only if vasp or QE
	if ((useVasp) || (useEspresso))
	{
		
		// Add potential
		_ipo->add(data[0], types[0]);
		
		// Print warning if more potentials were passed
		if (data.length() > 1)
		{
			Output::print(WARNING);
			Output::print("Multiple ");
			Output::print(potentialType(types[0]));
			Output::print(" potentials were defined but only first is used");
		}
	}
	
	// Add all local potentials
	else
	{
		for (i = 0; i < data.length(); ++i)
			_ipo->add(data[i], types[i]);
	}
	
	// Print information about references
	if (!_useReferences)
	{
		Output::newline();
		Output::print("Reference states will not be used in energy calculations");
	}
	else
	{
		Output::newline();
		Output::print("Reference states will be used in energy calculations");
		Output::increase();
		for (i = 0; i < _references.length(); ++i)
		{
			Output::newline();
			Output::print("Reference state for ");
			Output::print(_references[i].element().symbol());
			Output::print(" has a set value of ");
			Output::print(_references[i].energyPerAtom());
			Output::print(" eV/atom");
		}
		Output::decrease();
	}
	
	// Output
	Output::decrease();
}



/* void Potential::setReferenceData(const Text& input)
 *
 * Parse reference data text to look for values
 */

void Potential::setReferenceData(const Text& input)
{
	
	// Check for unphysical input setting
	if (input[0].length() > 1)
	{
		if (Language::isNumber(input[0][1], true))
			_unphysicalCutoff = Language::fractionToNumber(input[0][1]);
		else
		{
			Output::newline(ERROR);
			Output::print("Did not recognize value as a number in reference setting (");
			Output::print(input[0][1]);
			Output::print(")");
			Output::quit();
		}
	}
	
	// Loop over contents and check for predefined values
	Element element;
	for (int i = 1; i < input.length(); ++i)
	{
		
		// Skip blank lines and comments
		if (!input[i].length())
			continue;
		if (Language::isComment(input[i][0]))
			continue;
		
		// Check for errors
		element = Element::find(input[i][0], true, true);
		if (input[i].length() < 2)
		{
			Output::newline(ERROR);
			Output::print("Did not find an energy for reference state of ");
			Output::print(element.symbol());
			Output::quit();
		}
		if (!Language::isNumber(input[i][1], true))
		{
			Output::newline(ERROR);
			Output::print("Did not recognize \"");
			Output::print(input[i][1]);
			Output::print("\" as an energy for reference state of ");
			Output::print(element.symbol());
			Output::quit();
		}
		
		// Save referece
		_references.add();
		_references.last().set(element, false);
		_references.last().set(Language::fractionToNumber(input[i][1]));
	}
}



/* OList<Text> Potential::parseInput(const Text& input)
 *
 * Convert input into different potentials
 */

OList<Text> Potential::parseInput(const Text& input)
{
	
	// Copy text and remove = signs
	Text modInput = input;
	modInput.split('=');
	
	// Loop until a potential or references is found
	int i;
	for (i = 0; i < modInput.length(); ++i)
	{
		if (!modInput[i].length())
			continue;
		if (modInput[i][0].equal("references", false, 3))
			break;
		if (potentialType(modInput[i][0]) != PT_UNKNOWN)
			break;
	}
	
	// Get potentials
	OList<Text> res;
	for (; i < modInput.length(); ++i)
	{
		
		// Line has no data
		if (!modInput[i].length())
			continue;
		
		// Found start of a new potential or references
		if ((potentialType(modInput[i][0]) != PT_UNKNOWN) || (modInput[i][0].equal("references", false, 3)))
		{
			res.add();
			res.last().addLine();
			res.last().addWords(modInput[i]);
		}
		
		// Continuing current potential
		else
		{
			res.last().addLine();
			res.last().addWords(modInput[i]);
		}
	}
	
	// Return parsed data
	return res;
}



/* double Potential::getReferenceEnergy(const ISO& iso) const
 *
 * Get reference energy for current structure
 */

double Potential::getReferenceEnergy(const ISO& iso) const
{
	
	// Return if not using reference energy
	if (!_useReferences)
		return 0;
	
	// Loop over elements in the structure
	int i, j;
	bool found;
	double res = 0;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Loop over reference states that have been set
		found = false;
		for (j = 0; j < _references.length(); ++j)
		{
			
			// Found current element so save
			if (iso.atoms()[i][0].element() == _references[j].element())
			{
				res += iso.atoms()[i].length() * _references[j].energyPerAtom();
				found = true;
				break;
			}
		}
		
		// Found a new element
		if (!found)
		{
			
			// Set reference state and get energy
			_references.add();
			_references.last().set(iso.atoms()[i][0].element(), true);
			_references.last().set(*this);
			res += iso.atoms()[i].length() * _references.last().energyPerAtom();
		}
	}
	
	// Return reference energy
	return res;
}



/* void Potential::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
 *
 * Perform nudged elastic band calculation
 */

void Potential::neb(OList<ISO>& isos, double* tsEnergy, ISO* tsISO) const
{
	
	// Set run
	errorIfNotSet();
	if (tsEnergy)
		tsEnergy = 0;
	
	// Perform calculation
	_ipo->neb(isos, tsEnergy, tsISO);
}



/* Word Potential::potentialType(PotentialType type)
 *
 * Return word for a potential
 */

Word Potential::potentialType(PotentialType type)
{
	
	// Found vasp
	if (type == PT_VASP)
		return Word("VASP");
	
	// Found quantum espresso
	if (type == PT_QE)
		return Word("Quantum Espresso");
	
	// Found ewald
	if (type == PT_EWALD)
		return Word("Ewald");
	
	// Found lennard jones
	if (type == PT_LENNARDJONES)
		return Word("Lennard-Jones");
	
	// Found buckingham
	if (type == PT_BUCKINGHAM)
		return Word("Buckingham");
	
	// Found power
	if (type == PT_POWER)
		return Word("power");
	
	// Found exponential
	if (type == PT_EXPONENTIAL)
		return Word("exponential");
	
	// Found covalent
	if (type == PT_COVALENT)
		return Word("covalent");
	
	// Return that word is not a potential type if at this point
	return Word("unknown");
}



/* PotentialType Potential::potentialType(const Word& word)
 *
 * Return type of potential
 */

PotentialType Potential::potentialType(const Word& word)
{
	
	// Found vasp
	if (word.equal("vasp", false, 4))
		return PT_VASP;
	
	// Found quantum espresso
	if ((word.equal("quantum", false, 4)) || (word.equal("espresso", false, 3)) || (word.equal("qe", false)))
		return PT_QE;
	
	// Found ewald
	if (word.equal("ewald", false, 5))
		return PT_EWALD;
	
	// Found electrostatic
	if (word.equal("electrostatic", false, 7))
		return PT_ELECTROSTATIC;
	
	// Found lennard jones
	if (word.equal("lennard", false, 4))
		return PT_LENNARDJONES;
	
	// Found buckingham
	if (word.equal("buckingham", false, 4))
		return PT_BUCKINGHAM;
	
	// Found power
	if (word.equal("power", false, 3))
		return PT_POWER;
	
	// Found exponential
	if (word.equal("exponential", false, 3))
		return PT_EXPONENTIAL;
	
	// Found covalent
	if (word.equal("covalent", false, 3))
		return PT_COVALENT;
	
	// Return that word is not a potential type if at this point
	return PT_UNKNOWN;
}



/* void IPO::readError(const OList<Word>& line)
 *
 * Print error when reading file
 */

void IPO::readError(const OList<Word>& line)
{
	Output::newline(ERROR);
	Output::print("Did not recognize potential setting on line \"");
	for (int i = 0; i < line.length(); ++i)
	{
		Output::print(line[i]);
		if (i != line.length() - 1)
			Output::print(" ");
	}
	Output::print("\"");
	Output::quit();
}



/* void Reference::set(const Potential& potential)
 *
 * Set the reference potential
 */

void Reference::set(const Potential& potential)
{
	
	// Make sure that potential has been set
	potential.errorIfNotSet();
	
	// Make sure that structure has been set
	if (!_iso.numAtoms())
	{
		Output::newline(ERROR);
		Output::print("Internal error: structure has not been set for reference state");
		Output::quit();
	}
	
	// Output
	Output::newline();
	Output::print("Calculating reference state energy for ");
	Output::print(_iso.atoms()[0][0].element().symbol());
	Output::increase();
	
	// Evaluate the potential
	_energyPerAtom = 0;
	potential._ipo->relax(_iso, &_energyPerAtom);
	_energyPerAtom /= _iso.atoms()[0].length();
	
	// Output
	Output::newline();
	Output::print("Energy of ");
	Output::print(_energyPerAtom);
	Output::print(" eV/atom");
	
	// Output
	Output::decrease();
}



/* void Reference::set(const Element& element, bool setISO)
 *
 * Set the structure for the reference element
 */

void Reference::set(const Element& element, bool setISO)
{
	
	// Save element
	_element = element;
	if (!setISO)
		return;
	
	// Set by symbol
	Atom* atom;
	if (element.symbol() == "Ac")
	{
		_iso.basis(Vector3D(3.75544411, 3.75544411, 3.75544411), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ag")
	{
		_iso.basis(Vector3D(2.92999962, 2.92999962, 4.79000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Al")
	{
		_iso.basis(Vector3D(2.87032320, 2.87032320, 2.87032320), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ar")
	{
		_iso.basis(Vector3D(3.80000040, 3.80000040, 6.20000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "As")
	{
		_iso.basis(Vector3D(3.75970033, 3.75970033, 4.10182476), \
			Vector3D(1.09471919194, 1.09471919194, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.22764000019155, 0.22764000019155, 0.31707999942535);
		atom = _iso.addAtom(element);
		atom->fractional(0.77235999980845, 0.77235999980845, 0.68292000057465);
	}
	else if (element.symbol() == "Au")
	{
		_iso.basis(Vector3D(2.88471282, 2.88471282, 2.88471282), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "B")
	{
		_iso.basis(Vector3D(4.90799604, 4.90799604, 5.05700021), \
			Vector3D(1.06412721126, 1.06412721126, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.65729996594790, 0.01040001166617, 0.32190001071976);
		atom = _iso.addAtom(element);
		atom->fractional(0.34270003405210, 0.98959998833383, 0.67809998928024);
		atom = _iso.addAtom(element);
		atom->fractional(0.98959998833383, 0.98959998833383, 0.67809998928024);
		atom = _iso.addAtom(element);
		atom->fractional(0.01040001166617, 0.01040001166617, 0.32190001071976);
		atom = _iso.addAtom(element);
		atom->fractional(0.98959998833383, 0.34270003405210, 0.67809998928024);
		atom = _iso.addAtom(element);
		atom->fractional(0.01040001166617, 0.65729996594790, 0.32190001071976);
		atom = _iso.addAtom(element);
		atom->fractional(0.63229997941912, 0.22060000772391, 0.92650000513307);
		atom = _iso.addAtom(element);
		atom->fractional(0.36770002058088, 0.77939999227609, 0.07349999486693);
		atom = _iso.addAtom(element);
		atom->fractional(0.77939999227609, 0.77939999227609, 0.07349999486693);
		atom = _iso.addAtom(element);
		atom->fractional(0.22060000772391, 0.22060000772391, 0.92650000513307);
		atom = _iso.addAtom(element);
		atom->fractional(0.77939999227609, 0.36770002058088, 0.07349999486693);
		atom = _iso.addAtom(element);
		atom->fractional(0.22060000772391, 0.63229997941912, 0.92650000513307);
	}
	else if (element.symbol() == "Ba")
	{
		_iso.basis(Vector3D(4.35437573, 4.35437573, 4.35437573), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Be")
	{
		_iso.basis(Vector3D(2.28580011, 2.28580011, 3.58430000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.75000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.25000000000000);
	}
	else if (element.symbol() == "Bi")
	{
		_iso.basis(Vector3D(3.30399968, 6.11700000, 6.33534114), \
			Vector3D(1.57079632679, 1.72695032987, 1.57079632679), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.87499997072440, 0.00000000000000, 0.75000002236756);
		atom = _iso.addAtom(element);
		atom->fractional(0.12500002927560, 0.00000000000000, 0.24999997763244);
	}
	else if (element.symbol() == "Br")
	{
		_iso.basis(Vector3D(4.01744011, 4.01744011, 8.72000000), \
			Vector3D(1.57079632679, 1.57079632679, 1.95868378044), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.13500000000000, 0.13500000000000, 0.89000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.86500000000000, 0.86500000000000, 0.11000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.36500000000000, 0.36500000000000, 0.39000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.63500000000000, 0.63500000000000, 0.61000000000000);
	}
	else if (element.symbol() == "C")
	{
		_iso.basis(Vector3D(2.47000022, 2.47000022, 6.79000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.75000000000000);
	}
	else if (element.symbol() == "Ca")
	{
		_iso.basis(Vector3D(3.98999969, 3.98999969, 6.53000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Cd")
	{
		_iso.basis(Vector3D(2.97912035, 2.97912035, 5.61827000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Ce")
	{
		_iso.basis(Vector3D(3.62272017, 3.62272017, 3.62272017), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Cl")
	{
		_iso.basis(Vector3D(3.86697879, 3.86697879, 8.21000000), \
			Vector3D(1.57079632679, 1.57079632679, 1.89959060065), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.11730000000000, 0.11730000000000, 0.89840000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.88270000000000, 0.88270000000000, 0.10160000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.38270000000000, 0.38270000000000, 0.39840000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.61730000000000, 0.61730000000000, 0.60160000000000);
	}
	else if (element.symbol() == "Co")
	{
		_iso.basis(Vector3D(2.50539996, 2.50539996, 4.08930000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Cr")
	{
		_iso.basis(Vector3D(2.49843133, 2.49843133, 2.49843133), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Cs")
	{
		_iso.basis(Vector3D(5.25417612, 5.25417612, 5.25417612), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Cu")
	{
		_iso.basis(Vector3D(2.55622637, 2.55622637, 2.55622637), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Dy")
	{
		_iso.basis(Vector3D(3.58000006, 3.58000006, 8.78332608), \
			Vector3D(1.36556336479, 1.36556336479, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.22200001085514, 0.22200001085514, 0.33399996743459);
		atom = _iso.addAtom(element);
		atom->fractional(0.77799998914487, 0.77799998914486, 0.66600003256541);
	}
	else if (element.symbol() == "Er")
	{
		_iso.basis(Vector3D(3.55899964, 3.55899964, 5.59200000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Eu")
	{
		_iso.basis(Vector3D(3.96466430, 3.96466430, 3.96466430), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "F")
	{
		_iso.basis(Vector3D(3.20189007, 3.20189007, 7.28000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.06610982669), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.10200000000000, 0.03200000000000, 0.09970000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.89800000000000, 0.96800000000000, 0.90030000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.53200000000000, 0.60200000000000, 0.40030000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.46800000000000, 0.39800000000000, 0.59970000000000);
	}
	else if (element.symbol() == "Fe")
	{
		_iso.basis(Vector3D(2.53875347, 2.53875347, 2.53875347), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ga")
	{
		_iso.basis(Vector3D(4.44827073, 4.44827073, 4.52400000), \
			Vector3D(1.57079632679, 1.57079632679, 2.07494164038), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.84510000000000, 0.15490000000000, 0.08100000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.15490000000000, 0.84510000000000, 0.91900000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.65490000000000, 0.34510000000000, 0.58100000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.34510000000000, 0.65490000000000, 0.41900000000000);
	}
	else if (element.symbol() == "Gd")
	{
		_iso.basis(Vector3D(3.89615836, 3.89615836, 3.89615836), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ge")
	{
		_iso.basis(Vector3D(3.97627356, 3.97627356, 3.97627356), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.25000000000000, 0.25000000000000, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.50000000000000, 0.50000000000000, 0.50000000000000);
	}
	else if (element.symbol() == "H")
	{
		_iso.basis(Vector3D(3.60000000, 3.65581728, 3.65581728), \
			Vector3D(1.32593247261, 1.05599037162, 1.05599037162), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.78000000000000, 0.00000000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.22000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "He")
	{
		_iso.basis(Vector3D(3.55936441, 3.55936441, 3.55936441), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Hf")
	{
		_iso.basis(Vector3D(3.19799979, 3.19799979, 5.06100000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Hg")
	{
		_iso.basis(Vector3D(2.75660000, 3.26620884, 8.83152015), \
			Vector3D(1.39431183429, 1.41408974393, 1.13515963467), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ho")
	{
		_iso.basis(Vector3D(3.57609961, 3.57609961, 5.61740000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "I")
	{
		_iso.basis(Vector3D(4.26852117, 4.26852117, 9.78400000), \
			Vector3D(1.57079632679, 1.57079632679, 1.9794923557), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.15434000000000, 0.15434000000000, 0.88259000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.84566000000000, 0.84566000000000, 0.11741000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.34566000000000, 0.34566000000000, 0.38259000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.65434000000000, 0.65434000000000, 0.61741000000000);
	}
	else if (element.symbol() == "In")
	{
		_iso.basis(Vector3D(3.24800000, 3.24800000, 3.37571148), \
			Vector3D(2.07268684691, 2.07268684691, 1.57079632679), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ir")
	{
		_iso.basis(Vector3D(2.75771645, 2.75771645, 2.75771645), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "K")
	{
		_iso.basis(Vector3D(4.61418335, 4.61418335, 4.61418335), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Kr")
	{
		_iso.basis(Vector3D(4.00000033, 4.00000033, 6.53000000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "La")
	{
		_iso.basis(Vector3D(3.74130198, 3.74130198, 3.74130198), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Li")
	{
		_iso.basis(Vector3D(3.03888314, 3.03888314, 3.03888314), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Lu")
	{
		_iso.basis(Vector3D(3.50499997, 3.50499997, 5.54860000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Mg")
	{
		_iso.basis(Vector3D(3.20850043, 3.20850043, 5.21060000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Mn")
	{
		_iso.basis(Vector3D(2.46780267, 2.46780267, 2.46780267), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Mo")
	{
		_iso.basis(Vector3D(2.72570238, 2.72570238, 2.72570238), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "N")
	{
		_iso.basis(Vector3D(3.95700000, 3.95700000, 5.10900000), \
			Vector3D(1.57079632679, 1.57079632679, 1.57079632679), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.90500000000000, 0.09500000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.09500000000000, 0.90500000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.59500000000000, 0.59500000000000, 0.50000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.40500000000000, 0.40500000000000, 0.50000000000000);
	}
	else if (element.symbol() == "Na")
	{
		_iso.basis(Vector3D(3.76700026, 3.76700026, 6.15400000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Nb")
	{
		_iso.basis(Vector3D(2.85277428, 2.85277428, 2.85277428), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Nd")
	{
		_iso.basis(Vector3D(3.39411255, 3.39411255, 3.39411255), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ne")
	{
		_iso.basis(Vector3D(3.19612265, 3.19612265, 3.19612265), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ni")
	{
		_iso.basis(Vector3D(2.43951840, 2.43951840, 2.43951840), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Np")
	{
		_iso.basis(Vector3D(3.38800000, 4.89700000, 4.89700000), \
			Vector3D(1.57079632679, 1.57079632679, 1.57079632679), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.50000000000000, 0.50000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.37500000000000, 0.00000000000000, 0.50000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.62500000000000, 0.50000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "O")
	{
		_iso.basis(Vector3D(3.30699970, 3.30699970, 4.20985974), \
			Vector3D(1.1671563518, 1.1671563518, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.38693333333333, 0.38693333333333, 0.83920000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.61306666666667, 0.61306666666667, 0.16080000000000);
	}
	else if (element.symbol() == "Os")
	{
		_iso.basis(Vector3D(2.72399983, 2.72399983, 4.29500000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "P")
	{
		_iso.basis(Vector3D(3.37700024, 3.37700024, 3.52385537), \
			Vector3D(1.07109578387, 1.07109578387, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.78500003236431, 0.78500003236431, 0.64499990290708);
		atom = _iso.addAtom(element);
		atom->fractional(0.21499996763570, 0.21499996763570, 0.35500009709291);
	}
	else if (element.symbol() == "Pa")
	{
		_iso.basis(Vector3D(3.55745422, 3.55745422, 3.55745422), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pb")
	{
		_iso.basis(Vector3D(3.50074425, 3.50074425, 3.50074425), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pd")
	{
		_iso.basis(Vector3D(2.75771645, 2.75771645, 2.75771645), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pm")
	{
		_iso.basis(Vector3D(3.51965400, 3.51965400, 3.65468380), \
			Vector3D(2.07319189807, 2.07319189807, 1.57079632679), false);
			atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pr")
	{
		_iso.basis(Vector3D(3.66705577, 3.66705577, 3.66705577), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pt")
	{
		_iso.basis(Vector3D(2.80721392, 2.80721392, 2.80721392), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Pu")
	{
		_iso.basis(Vector3D(3.15870000, 3.28821825, 5.32079951), \
			Vector3D(1.4277418254, 1.26942951781, 1.06979310153), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.75000000000000, 0.50000000000000, 0.50000000000000);
	}
	else if (element.symbol() == "Rb")
	{
		_iso.basis(Vector3D(4.85407239, 4.85407239, 4.85407239), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Re")
	{
		_iso.basis(Vector3D(2.75999990, 2.75999990, 4.45800000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Rh")
	{
		_iso.basis(Vector3D(2.70821897, 2.70821897, 2.70821897), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ru")
	{
		_iso.basis(Vector3D(2.75080028, 2.75080028, 4.28190000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "S")
	{
		_iso.basis(Vector3D(4.27999881, 6.40664375, 6.40664375), \
			Vector3D(2.01043971561, 1.79536454018, 1.79536454018), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.08269998955882, 0.23100000731161, 0.33360000657511);
		atom = _iso.addAtom(element);
		atom->fractional(0.91730001044118, 0.76899999268839, 0.66639999342489);
		atom = _iso.addAtom(element);
		atom->fractional(0.74909998298370, 0.66639999342489, 0.89740000073650);
		atom = _iso.addAtom(element);
		atom->fractional(0.25090001701630, 0.33360000657511, 0.10259999926350);
		atom = _iso.addAtom(element);
		atom->fractional(0.14830001775279, 0.89740000073650, 0.23100000731161);
		atom = _iso.addAtom(element);
		atom->fractional(0.85169998224721, 0.10259999926350, 0.76899999268839);
	}
	else if (element.symbol() == "Sb")
	{
		_iso.basis(Vector3D(4.26817705, 4.26817705, 4.49999957), \
			Vector3D(1.07669354155, 1.07669354155, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.56633330004700, 0.56633330004700, 0.30100009985900);
		atom = _iso.addAtom(element);
		atom->fractional(0.43366669995300, 0.43366669995300, 0.69899990014100);
	}
	else if (element.symbol() == "Sc")
	{
		_iso.basis(Vector3D(3.56764099, 3.56764099, 3.56764099), \
			Vector3D(1.68057790658, 2.03221559807, 2.03221559807), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Se")
	{
		_iso.basis(Vector3D(4.36619990, 4.36619980, 4.95360000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439508895), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.77459999841322, 0.00000000515704, 0.66666666666667);
		atom = _iso.addAtom(element);
		atom->fractional(0.22540000158678, 0.22540000674382, 0.33333333333334);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.27460002261163, 0.00000000000000);
	}
	else if (element.symbol() == "Si")
	{
		_iso.basis(Vector3D(3.84009187, 3.84009187, 3.84009187), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.25000000000000, 0.25000000000000, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.50000000000000, 0.50000000000000, 0.50000000000000);
	}
	else if (element.symbol() == "Sm")
	{
		_iso.basis(Vector3D(3.64200042, 3.64200042, 5.77200000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Sn")
	{
		_iso.basis(Vector3D(4.58855732, 4.58855732, 4.58855732), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.25000000000000, 0.25000000000000, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.50000000000000, 0.50000000000000, 0.50000000000000);
	}
	else if (element.symbol() == "Sr")
	{
		_iso.basis(Vector3D(4.29638080, 4.29638080, 4.29638080), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ta")
	{
		_iso.basis(Vector3D(2.84056332, 2.84056332, 2.84056332), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Tb")
	{
		_iso.basis(Vector3D(3.67695526, 3.67695526, 3.67695526), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Tc")
	{
		_iso.basis(Vector3D(2.75771645, 2.75771645, 2.75771645), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Te")
	{
		_iso.basis(Vector3D(2.94713603, 2.94713603, 2.94713603), \
			Vector3D(1.79229902888, 1.79229902888, 1.79229902888), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Th")
	{
		_iso.basis(Vector3D(3.59641580, 3.59641580, 3.59641580), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Ti")
	{
		_iso.basis(Vector3D(2.87085353, 2.87085353, 2.87085353), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Tl")
	{
		_iso.basis(Vector3D(3.36000000, 3.36000000, 3.60025346), \
			Vector3D(2.05627730247, 2.05627730247, 1.57079632679), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Tm")
	{
		_iso.basis(Vector3D(3.57796031, 3.57796031, 3.57796031), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "U")
	{
		_iso.basis(Vector3D(2.85400000, 3.26351865, 4.95500000), \
			Vector3D(1.57079632679, 1.57079632679, 2.02334393938), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.10250000000000, 0.20500000000000, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.89750000000000, 0.79500000000000, 0.75000000000000);
	}
	else if (element.symbol() == "V")
	{
		_iso.basis(Vector3D(2.62041967, 2.62041967, 2.62041967), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "W")
	{
		_iso.basis(Vector3D(2.74075390, 2.74075390, 2.74075390), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Xe")
	{
		_iso.basis(Vector3D(4.33576665, 4.33576665, 4.33576665), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Y")
	{
		_iso.basis(Vector3D(4.12243253, 4.12243253, 4.12243253), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Yb")
	{
		_iso.basis(Vector3D(3.84515279, 3.84515279, 3.84515279), \
			Vector3D(1.91063323617, 1.91063323617, 1.91063323617), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else if (element.symbol() == "Zn")
	{
		_iso.basis(Vector3D(2.66450027, 2.66450027, 4.94530000), \
			Vector3D(1.57079632679, 1.57079632679, 2.09439510239), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.33333333333333, 0.66666666666667, 0.25000000000000);
		atom = _iso.addAtom(element);
		atom->fractional(0.66666666666667, 0.33333333333333, 0.75000000000000);
	}
	else if (element.symbol() == "Zr")
	{
		_iso.basis(Vector3D(3.18905158, 3.18905158, 3.18905158), \
			Vector3D(1.0471975512, 1.0471975512, 1.0471975512), false);
		atom = _iso.addAtom(element);
		atom->fractional(0.00000000000000, 0.00000000000000, 0.00000000000000);
	}
	else
	{
		
		// Element does not have a structure yet
		Output::newline(WARNING);
		Output::print("No internal structure for ");
		Output::print(element.symbol());
		Output::print(": assuming simple cubic");
		
		// Set structure
		const double length = 2*element.radius();
		const double angle = 1.57079632679;
		_iso.basis(Vector3D(length, length, length), Vector3D(angle, angle, angle), false);
		atom = _iso.addAtom(element);
		atom->fractional(0, 0, 0);
	}
}

