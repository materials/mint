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



#include "gaPredict.h"
#include "structureIO.h"
#include "randomStructure.h"
#include "fileSystem.h"
#include "language.h"
#include "output.h"
#include "num.h"



/* ISOSymmetryPair& ISOSymmetryPair::operator= (const ISOSymmetryPair& rhs)
 *
 * Assignment operator for ISOSymmetryPair object
 */

ISOSymmetryPair& ISOSymmetryPair::operator= (const ISOSymmetryPair& rhs)
{
	if (this != &rhs)
	{
		
		// Copy structure and symmetry
		_iso = rhs._iso;
		_symmetry = rhs._symmetry;
		
		// Assign atoms in symmetry orbits
		int i, j, k, m;
		int curAtomNum;
		for (i = 0; i < _symmetry.orbits().length(); ++i)
		{
			
			// Loop over elements in the structure
			for (j = 0; j < _iso.atoms().length(); ++j)
			{
				
				// Elements are the same
				if (_symmetry.orbits()[i].atoms()[0]->element() == _iso.atoms()[j][0].element())
				{
					
					// Loop over atoms in orbit
					for (k = 0; k < _symmetry.orbits()[i].atoms().length(); ++k)
					{
						
						// Loop over atoms of current element
						curAtomNum = _symmetry.orbits()[i].atoms()[k]->atomNumber();
						for (m = 0; m < _iso.atoms()[j].length(); ++m)
						{
							
							// Found atom
							if (_iso.atoms()[j][m].atomNumber() == curAtomNum)
							{
								_symmetry.orbits()[i].set(k, &(_iso.atoms()[j][m]));
								break;
							}
						}
					}
					
					// Break since element was found
					break;
				}
			}
		}
	}
	return *this;
}



/* void GAPredict::run(ISO& iso, Random& random, Potential* potential, Diffraction* diffraction)
 *
 * Run structure prediction using genetic algorithm
 */

void GAPredict::run(ISO& iso, Random& random, Potential* potential, Diffraction* diffraction)
{
	
	// Figure out which metrics to run
	_metrics.length(0);
	if (potential)
	{
		if (potential->isSet())
		{
			_metrics += GAPM_POTENTIAL;
			_potential = potential;
		}
	}
	if (diffraction)
	{
		if (diffraction->isSet())
		{
			_metrics += GAPM_DIFFRACTION;
			_refDiffraction = diffraction;
		}
	}
	
	// No metrics were passed
	if (!_metrics.length())
	{
		Output::newline(ERROR);
		Output::print("Attempting to predict structure with genetic algorithm when no metrics have been defined");
		Output::quit();
	}
	
	// Set the metric to optimize
	GAPredictMetric metToOpt = (_optMetric == GAPM_UNKNOWN) ? _metrics[0] : _optMetric;
	
	// Look for metric to optimize
	int i;
	int metToOptNumber = -1;
	for (i = 0; i < _metrics.length(); ++i)
	{
		if (_metrics[i] == metToOpt)
		{
			metToOptNumber = i;
			break;
		}
	}
	
	// Metric to optimize has not been defined
	if (metToOptNumber < 0)
	{
		Output::newline(ERROR);
		Output::print("Attempting to optimize structure ");
		Output::print(metric(metToOpt).tolower());
		Output::print(" in genetic algorithm when it has not been defined");
		Output::quit();
	}
	
	// Set metric tolerance
	switch (metToOpt)
	{
		case GAPM_POTENTIAL:
			_ga.metricTolerance(_energyTolerance * iso.numAtoms());
			break;
		case GAPM_DIFFRACTION:
			_ga.metricTolerance(_diffractionTolerance);
			break;
		default:
			Output::newline(ERROR);
			Output::print("Internal error: no tolerance for current metric to optimize");
			Output::quit();
	}
	
	// Setup GA metrics
	_ga.metricToOptimize(metToOptNumber);
	_ga.numMetrics(_metrics.length());
	for (i = 0; i < _metrics.length(); ++i)
	{
		if (_metrics[i] == GAPM_POTENTIAL)
		{
			_ga.fitness().optLowest()[i] = true;
			_ga.fitness().names()[i] = metric(_metrics[i]);
			_ga.fitness().units()[i] = "eV";
		}
		else if (_metrics[i] == GAPM_DIFFRACTION)
		{
			_ga.fitness().optLowest()[i] = true;
			_ga.fitness().names()[i] = metric(_metrics[i]);
		}
	}
	
	// Setup crossover
	_ga.crossover().set(this, &GAPredict::crossover);
	
	// Setup mutations
	_ga.mutations().length(3);
	_ga.mutations()[0].set(this, &GAPredict::mutateBasis, _cellMutationProb);
	_ga.mutations()[1].set(this, &GAPredict::mutatePositions, _posMutationProb);
	_ga.mutations()[2].set(this, &GAPredict::mutateWyckoff, _wyckMutationProb);
	
	// Setup screen if needed
	if ((_screenMetric != GAPM_UNKNOWN) && (_numScreens != 0))
	{
		
		// Loop over metrics and check that current is available
		for (i = 0; i < _metrics.length(); ++i)
		{
			if (_metrics[i] == _screenMetric)
			{
				_ga.screen().set(this, &GAPredict::screen, _ga.fitness().optLowest()[i], _numScreens);
				break;
			}
		}
		
		// Did not find metric
		if (i >= _metrics.length())
		{
			Output::newline(WARNING);
			Output::print("Ignoring screen in GA using ");
			Output::print(metric(_screenMetric).tolower());
			Output::print(" because it has not been supplied");
		}
	}
	
	// Check if Gaussian comparisons should be used for diffraction
	bool uesGaussianMatch = true;
	for (i = 0; i < 3; ++i)
	{
		if ((iso.basis().angleFixed()[i]) || (iso.basis().lengthFixed()[i]))
			uesGaussianMatch = false;
	}
	
	// Get current directory
	Word origDir;
	Directory::current(origDir);
	
	// File to write during optimization
	int logID;
	int startID;
	int fitnessID;
	int origID = Output::streamID();
	Word structureFile = "best";
	StructureFormat structureFormat = SF_MINT;
	StructureIO::addExtension(structureFile, structureFormat);
	
	// Create diffraction object
	CalculatedPattern curDiffraction;
	_diffraction = &curDiffraction;
	
	// Loop over number of simulations to run
	int j, k, m;
	double bestMetric;
	Word curDir;
	ISO bestStructure = iso;
	for (i = 0; i < _numSimulations; ++i)
	{
		
		// Reset ga
		_ga.initNewRun();
		
		// Create directory for current run
		curDir = Directory::makePath(origDir, Word("OPT_"));
		curDir += Language::numberToWord(i+1);
		Directory::create(curDir, true);
		Directory::change(curDir);
		
		// Open startup stream
		startID = Output::addStream("start.out");
		Output::setStream(startID);
		Output::method(RESTRICTED);
		
		// Output
		Output::newline();
		Output::print("Initializing structure");
		if (_ga.population().length() != 1)
			Output::print("s");
		Output::increase();
		
		// Generate structures
		for (j = 0; j < _ga.population().length(); ++j)
		{
			_ga.population()[j].iso() = iso;
			RandomStructure::generate(_ga.population()[j].iso(), random, _wyckoffBias, \
				&_ga.population()[j].symmetry());
		}
		
		// Output
		Output::decrease();
		
		// Open log steeam
		logID = Output::addStream("log.out");
		Output::setStream(logID);
		Output::removeStream(startID);
		Output::method(RESTRICTED);
		
		// Open fitness stream
		fitnessID = Output::addStream("fitness.out");
		
		// Output
		Output::newline();
		Output::print("Running simulation");
		Output::increase();
		
		// Loop until max generations is reached
		for (j = 1; j <= _ga.maxGenerations(); ++j)
		{
			
			// Advance the GA
			_ga.advance(random);
			
			// Output
			Output::newline();
			Output::print("Calculating properties of structures in generation ");
			Output::print(j);
			Output::increase();
			
			// Calculate the fitness of each structure
			for (k = 0; k < _ga.population().length(); ++k)
			{
				
				// Relax against the metric to optimize
				Output::quietOn();
				if (metToOpt == GAPM_POTENTIAL)
					potential->relax(_ga.population()[k].iso(), _ga.population()[k].symmetry(), \
						&_ga.fitness().values()[metToOptNumber][k], 0, false, true);
				else if (metToOpt == GAPM_DIFFRACTION)
					_ga.fitness().values()[metToOptNumber][k] = curDiffraction.refine(_ga.population()[k].iso(), \
						_ga.population()[k].symmetry(), *diffraction, false);
				Output::quietOff();
				
				// Print metric value
				Output::newline();
				Output::print("Structure ");
				Output::print(k+1);
				Output::print(": ");
				if (_ga.fitness().names()[metToOptNumber].length())
				{
					Output::print(_ga.fitness().names()[metToOptNumber]);
					Output::print(" of ");
				}
				Output::print(_ga.fitness().values()[metToOptNumber][k]);
				Output::print(" ");
				Output::print(_ga.fitness().units()[metToOptNumber]);
				
				// Loop over other metrics
				for (m = 0; m < _metrics.length(); ++m)
				{
					
					// Skip if the optimized metric
					if (m == metToOptNumber)
						continue;
					
					// Get metric value
					Output::quietOn();
					if (_metrics[m] == GAPM_POTENTIAL)
						potential->single(_ga.population()[k].iso(), _ga.population()[k].symmetry(), \
							&_ga.fitness().values()[m][k], 0, false, true);
					else if (_metrics[m] == GAPM_DIFFRACTION)
						_ga.fitness().values()[m][k] = curDiffraction.set(_ga.population()[k].iso(), \
							_ga.population()[k].symmetry(), diffraction, true);
					Output::quietOff();
					
					// Print metric value
					Output::newline();
					Output::print("Structure ");
					Output::print(k+1);
					Output::print(": ");
					if (_ga.fitness().names()[m].length())
					{
						Output::print(_ga.fitness().names()[m]);
						Output::print(" of ");
					}
					Output::print(_ga.fitness().values()[m][k]);
					Output::print(" ");
					Output::print(_ga.fitness().units()[m]);
				}
			}
			Output::quietOff();
			
			// Output
			Output::decrease();
			
			// Set the fitness
			_ga.setFitness();
	
			// Print the fitness
			Output::setStream(fitnessID);
			printFitness(j);
			Output::setStream(logID);
			
			// Output
			Output::newline();
			Output::print("Best structure after ");
			Output::print(j);
			Output::print(" generation");
			if (j != 1)
				Output::print("s");
			Output::print(": ");
			if (_ga.fitness().names()[metToOptNumber].length())
			{
				Output::print(_ga.fitness().names()[metToOptNumber]);
				Output::print(" of ");
			}
			Output::print(_ga.bestFitness());
			Output::print(" ");
			Output::print(_ga.fitness().units()[metToOptNumber]);
			
			// Found a new best
			if (_ga.gensSinceLastBest() == 0)
				StructureIO::write(structureFile, _ga.bestIndividual().iso(), structureFormat);
			
			// Break if converged
			if (_ga.isConverged())
				break;
		}
		
		// Output
		Output::decrease();
		
		// Print convergence
		if (j >= _ga.maxGenerations())
		{
			Output::newline();
			Output::print("Reached the maximum number of generations (");
			Output::print(_ga.maxGenerations());
			Output::print(")");
		}
		else
		{
			Output::newline();
			Output::print("Reached convergence criterion after ");
			Output::print(j);
			Output::print(" generation");
			if (j != 1)
				Output::print("s");
		}
		
		// Print fitness
		Output::newline();
		if (_ga.fitness().names()[metToOptNumber].length() > 0)
			Output::print(_ga.fitness().names()[metToOptNumber]);
		else
			Output::print("Fitness");
		Output::print(" of optimal structure: ");
		Output::print(_ga.bestFitness());
		Output::print(" ");
		Output::print(_ga.fitness().units()[metToOptNumber]);
		
		// Found a new best
		if ((!i) || ((_ga.fitness().optLowest()[metToOptNumber]) && (_ga.bestFitness() < bestMetric)) || \
			((!_ga.fitness().optLowest()[metToOptNumber]) && (_ga.bestFitness() > bestMetric)))
		{
			bestStructure = _ga.bestIndividual().iso();
			bestMetric = _ga.bestFitness();
		}
		
		// Reset stream
		Output::removeStream(logID);
		Output::removeStream(fitnessID);
		Output::setStream(origID);
		
		// Change to original directory
		Directory::change(origDir);
	}
	
	// Print best result
	Output::newline();
	if (_ga.fitness().names()[metToOptNumber].length() > 0)
		Output::print(_ga.fitness().names()[metToOptNumber]);
	else
		Output::print("Fitness");
	Output::print(" of optimal structure: ");
	Output::print(bestMetric);
	Output::print(" ");
	Output::print(_ga.fitness().units()[metToOptNumber]);
	
	// Save the result
	iso = bestStructure;
}



/* void GAPredict::mutateBasis(ISOSymmetryPair& pair, Random& random)
 *
 * Basis mutation
 */

void GAPredict::mutateBasis(ISOSymmetryPair& pair, Random& random)
{
	Output::quietOn();
	RandomStructure::perturbBasis(pair.iso(), pair.symmetry(), .1, .5, random);
	Output::quietOff();
}



/* void GAPredict::mutatePositions(ISOSymmetryPair& pair, Random& random)
 *
 * Position mutation
 */

void GAPredict::mutatePositions(ISOSymmetryPair& pair, Random& random)
{
	Output::quietOn();
	RandomStructure::perturbAtoms(pair.iso(), pair.symmetry(), .1, .5, random);
	Output::quietOff();
}



/* void GAPredict::mutateWyckoff(ISOSymmetryPair& pair, Random& random)
 *
 * Mutate which Wyckoff positions are occupied
 */

void GAPredict::mutateWyckoff(ISOSymmetryPair& pair, Random& random)
{
	
	// Turn on quiet
	Output::quietOn();
	
	// Loop until at least one position is mutated
	int i, j;
	int numSelected = 0;
	double selectionProb = 0.5;
	for (int trial = 0; trial < 10; ++trial)
	{
		
		// Loop over orbits
		for (i = 0; i < pair.symmetry().orbits().length(); ++i)
		{
			
			// Check whether any atoms in the orbit are fixed
			if (pair.symmetry().orbits()[i].anyAtomsFixed())
				continue;
			
			// Pick whether orbit will be changed
			if (random.decimal(0, 1) < selectionProb)
			{
				
				// Set that all atoms in orbit will be assigned new positions
				++numSelected;
				for (j = 0; j < pair.symmetry().orbits()[i].atoms().length(); ++j)
					pair.symmetry().orbits()[i].atoms()[j]->assigned(false);
			}
		}
		
		// Break if any positions were selected
		if (numSelected > 0)
			break;
	}
	
	// Generate new positions
	RandomStructure::generate(pair.iso(), random, _wyckoffBias, &pair.symmetry());
	
	// Turn off quiet
	Output::quietOff();
}



/* void GAPredict::crossover(ISOSymmetryPair& pair, int mother, int father, Random& random)
 *
 * Crossover operator for prediction
 */

void GAPredict::crossover(ISOSymmetryPair& pair, int mother, int father, Random& random)
{
	
	// Clear space
	pair.iso().clear();
	pair.symmetry().clear();
	
	// Save space group of parents
	pair.iso().spaceGroup(_ga.population()[mother].iso().spaceGroup());
	pair.iso().basis().latticeSystem(_ga.population()[mother].iso().basis().latticeSystem());
	
	// Pick random number to adjust between mother and father basis
	int i;
	double scale = random.decimal(0, 1);
	Vector3D lengths = _ga.population()[mother].iso().basis().lengths() * scale + \
		_ga.population()[father].iso().basis().lengths() * (1 - scale);
	Vector3D angles = _ga.population()[mother].iso().basis().angles() * scale + \
		_ga.population()[father].iso().basis().angles() * (1 - scale);
	pair.iso().basis(lengths, angles, false);
	for (i = 0; i < 3; ++i)
	{
		if ((_ga.population()[mother].iso().basis().lengthFixed()[i]) || \
			(_ga.population()[father].iso().basis().lengthFixed()[i]))
			pair.iso().basis().lengthFixed(i, true);
		if ((_ga.population()[mother].iso().basis().angleFixed()[i]) || \
			(_ga.population()[father].iso().basis().angleFixed()[i]))
			pair.iso().basis().angleFixed(i, true);
	}
	
	// Get orbits that have atoms with fixed positions
	List<Orbit*>::D2 fixed = getFixedOrbits(mother, father);
	
	// Group together orbits in both structures based on multiplicity
	List<Orbit*>::D3 orbitGroups = groupOrbits(mother, father, fixed);
	
	// Get number of atoms to add for each element
	int j;
	List<int> numToAdd(orbitGroups.length());
	for (i = 0; i < orbitGroups.length(); ++i)
	{
		for (j = 0; j < _ga.population()[mother].iso().atoms().length(); ++j)
		{
			if (orbitGroups[i][0][0]->atoms()[0]->element() == _ga.population()[mother].iso().atoms()[j][0].element())
			{
				numToAdd[i] = _ga.population()[mother].iso().atoms()[j].length();
				break;
			}
		}
		for (j = 0; j < fixed[0].length(); ++j)
		{
			if (orbitGroups[i][0][0]->atoms()[0]->element() == fixed[0][j]->atoms()[0]->element())
				numToAdd[i] -= fixed[0][j]->atoms().length();
		}
	}
	
	// Get combinations of orbits that are allowed
	OList<Linked<List<int> > > allowedGroups(orbitGroups.length());
	for (i = 0; i < orbitGroups.length(); ++i)
		getOrbitCombinations(allowedGroups[i], orbitGroups[i], numToAdd[i]);
	
	// Loop until a good set of positions is found
	int k, m, n;
	bool found = true;
	int randGroup;
	double minBondLength;
	List<Orbit*> curGroups;
	List<Orbit*>::D2 groupsRemaining;
	Linked<List<int> >::iterator set;
	for (i = 0; i < _maxCrossoverLoops; ++i)
	{
		
		// Loop over elements
		curGroups.length(0);
		for (j = 0; j < allowedGroups.length(); ++j)
		{
			
			// Skip if empty
			if (!allowedGroups[j].length())
				continue;
			
			// Pick set of sites at random
			set = allowedGroups[j].begin() + random.integer(0, allowedGroups[j].length() - 1);
			
			// Save list of sites
			groupsRemaining = orbitGroups[j];
			
			// Loop over multiplicities
			for (k = 0; k < (*set).length(); ++k)
			{
				
				// Loop over number to add
				for (m = 0; m < (*set)[k]; ++m)
				{
					
					// Pick group at random
					randGroup = random.integer(0, groupsRemaining[k].length() - 1);
					
					// Add group and remove
					curGroups += groupsRemaining[k][randGroup];
					groupsRemaining[k].remove(randGroup);
				}
			}
		}
		
		// Loop over groups that are being saved
		found = false;
		for (j = 0; j < curGroups.length(); ++j)
		{
			
			// Loop over atoms in group
			for (k = 0; k < curGroups[j]->atoms().length(); ++k)
			{
				
				// Loop over fixed groups to get distances
				for (m = 0; m < fixed[0].length(); ++m)
				{
					minBondLength = _minBondFraction * (curGroups[j]->atoms()[0]->element().radius() + \
						fixed[0][m]->atoms()[0]->element().radius());
					for (n = 0; n < fixed[0][m]->atoms().length(); ++n)
					{
						
						// Atoms within cutoff distance
						if (pair.iso().basis().distance(curGroups[j]->atoms()[k]->fractional(), FRACTIONAL, \
							fixed[0][m]->atoms()[n]->fractional(), FRACTIONAL) < minBondLength)
						{
							found = true;
							break;
						}
					}
					if (found)
						break;
				}
				if (found)
					break;
				
				// Loop over all other groups to get distance
				for (m = j; m < curGroups.length(); ++m)
				{
					minBondLength = _minBondFraction * (curGroups[j]->atoms()[0]->element().radius() + \
						curGroups[m]->atoms()[0]->element().radius());
					for (n = (m == j) ? k + 1 : 0; n < curGroups[m]->atoms().length(); ++n)
					{
						
						// Atoms within cutoff distance
						if (pair.iso().basis().distance(curGroups[j]->atoms()[k]->fractional(), FRACTIONAL, \
							curGroups[m]->atoms()[n]->fractional(), FRACTIONAL) < minBondLength)
						{
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
			if (found)
				break;
		}
		
		// Break if a good combination was found
		if (!found)
			break;
	}
	
	// If a good structure was not found, then just save mother or father
	if (found)
	{
		if (random.integer(0, 1) == 0)
			pair = _ga.population()[mother];
		else
			pair = _ga.population()[father];
		return;
	}
	
	// Save symmetry properties
	pair.symmetry().operations(_ga.population()[mother].symmetry().operations());
	
// NEED TO SET SPECIAL POSITION OPERATORS AND POINT ROTATIONS/TRANSLATIONS
	
	// Loop over fixed groups to add
	Atom* curAtom;
	int orbitIndex = 0;
	for (i = 0; i < fixed[0].length(); ++i)
	{
		
		// Add orbit
		pair.symmetry().addOrbit(*(fixed[0][i]));
		
		// Loop over atoms in orbit
		for (j = 0; j < fixed[0][i]->atoms().length(); ++j)
		{
			
			// Add atom to structure
			curAtom = pair.iso().addAtom(*(fixed[0][i]->atoms()[j]));
			
			// Set pointer to atom in orbit
			pair.symmetry().orbit(orbitIndex).set(j, curAtom);
		}
		
		// Go to next orbit
		++orbitIndex;
	}
	
	// Loop over mixed groups to add
	for (i = 0; i < curGroups.length(); ++i)
	{
		
		// Add orbit
		pair.symmetry().addOrbit(*(curGroups[i]));
		
		// Loop over atoms in orbit
		for (j = 0; j < curGroups[i]->atoms().length(); ++j)
		{
			
			// Add atom to structure
			curAtom = pair.iso().addAtom(*(curGroups[i]->atoms()[j]));
			
			// Set pointer to atom in orbit
			pair.symmetry().orbit(orbitIndex).set(j, curAtom);
		}
		
		// Go to next orbit
		++orbitIndex;
	}
}



/* List<Orbit*>::D2 GAPredict::getFixedOrbits(int mother, int father)
 *
 * Get fixed orbits
 */

List<Orbit*>::D2 GAPredict::getFixedOrbits(int mother, int father)
{
	
	// Loop over orbits in mother
	int i, j;
	bool found;
	List<Orbit*>::D2 res (2);
	for (i = 0; i < _ga.population()[mother].symmetry().orbits().length(); ++i)
	{
		
		// Loop over atoms in orbit and check if any coordinates are fixed
		found = false;
		for (j = 0; j < _ga.population()[mother].symmetry().orbits()[i].atoms().length(); ++j)
		{
			if (_ga.population()[mother].symmetry().orbits()[i].atoms()[j]->anyFixed())
			{
				found = true;
				break;
			}
		}
		
		// Orbit not fixed
		if (!found)
			continue;
		
		// Loop over orbits in father to find equivalent orbit
		found = false;
		for (j = 0; j < _ga.population()[father].symmetry().orbits().length(); ++j)
		{
			if (_ga.population()[mother].symmetry().orbits()[i].atoms().length() == \
				_ga.population()[father].symmetry().orbits()[j].atoms().length())
			{
				if (_ga.population()[mother].symmetry().orbits()[i].atoms()[0]->equal(\
					*_ga.population()[father].symmetry().orbits()[j].atoms()[0], 1e-4))
				{
					
					// Found match so save
					res[0] += &(_ga.population()[mother].symmetry().orbits()[i]);
					res[1] += &(_ga.population()[father].symmetry().orbits()[j]);
					found = true;
					break;
				}
			}
		}
		
		// Did not find orbit in father
		if (!found)
		{
			Output::newline(ERROR);
			Output::print("Internal error: fixed orbit in mother not found in father");
			Output::quit();
		}
	}
	
	// Return result
	return res;
}



/* List<Orbit*>::D3 GAPredict::groupOrbits(int mother, int father, const List<Orbit*>::D2& fixedOrbits)
 *
 * Group orbits from mother and father by multiplicity
 */

List<Orbit*>::D3 GAPredict::groupOrbits(int mother, int father, const List<Orbit*>::D2& fixedOrbits)
{
	
	// Loop over parents
	int i, j, k;
	bool fixed;
	bool foundElem;
	bool foundMult;
	int parentIndex = 0;
	List<Orbit*>::D3 res;
	for (int parent = mother; ((parent == mother) || (parent == father)); parent += father - mother, ++parentIndex)
	{
		
		// Loop over orbits of current parent
		for (i = 0; i < _ga.population()[parent].symmetry().orbits().length(); ++i)
		{
			
			// Check if orbit is fixed
			fixed = false;
			for (j = 0; j < fixedOrbits[parentIndex].length(); ++j)
			{
				if (fixedOrbits[parentIndex][j] == &(_ga.population()[parent].symmetry().orbits()[i]))
				{
					fixed = true;
					break;
				}
			}
			if (fixed)
				continue;
			
			// Loop over elements
			foundElem = false;
			for (j = 0; j < res.length(); ++j)
			{

				// Elements are the same
				if (_ga.population()[parent].symmetry().orbits()[i].atoms()[0]->element() == \
					res[j][0][0]->atoms()[0]->element())
				{

					// Loop over multiplicities of current element
					foundMult = false;
					for (k = 0; k < res[j].length(); ++k)
					{

						// Orbit multiplicities are the same
						if (_ga.population()[parent].symmetry().orbits()[i].atoms().length() == \
							res[j][k][0]->atoms().length())
						{
							res[j][k] += &(_ga.population()[parent].symmetry().orbits()[i]);
							foundMult = true;
							break;
						}
					}

					// New multiplicity
					if (!foundMult)
					{
						res[j].add();
						res[j].last() += &(_ga.population()[parent].symmetry().orbits()[i]);
					}

					// Break since element was found
					foundElem = true;
					break;
				}
			}
			
			// Found new element
			if (!foundElem)
			{
				res.add();
				res.last().add();
				res.last().last() += &(_ga.population()[parent].symmetry().orbits()[i]);
			}
		}
	}
	
	// Return groups
	return res;
}



/* void GAPredict::getOrbitCombinations(Linked<List<int> >& allowedGroups, const List<Orbit*>::D2 orbitGroups,
 *		int numToAdd)
 *
 * Get combinations of allowed orbits
 */

void GAPredict::getOrbitCombinations(Linked<List<int> >& allowedGroups, const List<Orbit*>::D2 orbitGroups, \
	int numToAdd)
{
	
	// Return if no atoms to add
	if (numToAdd < 1)
		return;
	
	// Make list of multiplicities
	List<int> lengths(orbitGroups.length());
	List<int> mults(orbitGroups.length());
	for (int i = 0; i < orbitGroups.length(); ++i)
	{
		lengths[i] = orbitGroups[i].length();
		mults[i] = orbitGroups[i][0]->atoms().length();
	}
	
	// Make allowed combinations
	List<int> startList(orbitGroups.length());
	startList.fill(0);
	recurseOrbits(allowedGroups, lengths, mults, startList, 0, numToAdd);
}



/* void GAPredict::recurseOrbits(Linked<List<int> >& allowedGroups, const List<int>& lengths, const List<int>& mults,
 *		const List<int>& curList, int curGroup, int numToAdd)
 *
 * Recursively generate list of allowed orbit combinations
 */

void GAPredict::recurseOrbits(Linked<List<int> >& allowedGroups, const List<int>& lengths, const List<int>& mults, \
	const List<int>& curList, int curGroup, int numToAdd)
{
	
	// Count current number of atoms
	int i;
	int total = 0;
	for (i = 0; i < mults.length(); ++i)
		total += curList[i] * mults[i];
	
	// Return if out of range
	if (total > numToAdd)
		return;
	
	// Have the correct number of atoms
	if (total == numToAdd)
	{
		allowedGroups.add(curList);
		return;
	}
	
	// Return if on last element
	if (curGroup >= lengths.length())
		return;
	
	// Figure out max number of instances to add
	int max = (int)Num<double>::floor((double)(numToAdd - total) / mults[curGroup]);
	if (max > lengths[curGroup])
		max = lengths[curGroup];

	// Loop over multiplicities of current group
	List<int> newList = curList;
	for (i = 0; i <= max; ++i)
	{
		newList[curGroup] = i;
		recurseOrbits(allowedGroups, lengths, mults, newList, curGroup + 1, numToAdd);
	}
}



/* Return the static value of the current screening function
 */

double GAPredict::screen(ISOSymmetryPair& pair)
{
	Output::quietOn(true);
	double value = 0;
	if (_screenMetric == GAPM_POTENTIAL)
		_potential->single(pair.iso(), pair.symmetry(), &value, 0, false, true);
	else if (_screenMetric == GAPM_DIFFRACTION)
		value = _diffraction->set(pair.iso(), pair.symmetry(), _refDiffraction, _useReitveld, true);
	Output::quietOff();
	return value;
}



/* void GAPredict::printFitness(int generation)
 *
 * Print the current fitness
 */

void GAPredict::printFitness(int generation)
{
	
	// Print generation
	Output::newline();
	if (generation > 1)
		Output::newline();
	Output::print("Generation ");
	Output::print(generation);
	
	// Initialize message
	int i;
	Output message;
	message.addLines(_ga.population().length() + 1);
	message.addLine();
	message.add("Rank");
	for (i = 0; i < _ga.fitness().names().length(); ++i)
	{
		message.add("  ");
		if (_ga.fitness().names()[i].length() > 0)
		{
			if (_ga.fitness().units()[i].length() > 0)
				message.add(_ga.fitness().names()[i] + Word(" (") + _ga.fitness().units()[i] + ')');
			else
				message.add(_ga.fitness().names()[i]);
		}
		else
			message.add(Word("Metric ") + Language::numberToWord(i+1));
	}
	
	// Loop over population
	int j, k;
	for (i = 0; i < _ga.population().length(); ++i)
	{
		
		// Print number
		message.addLine();
		message.add(i + 1);
		
		// Loop over order to find current
		for (j = 0; j < _ga.fitness().totalOrder().length(); ++j)
		{
			
			// Found current structure in order
			if (_ga.fitness().totalOrder()[j] == i)
			{
				
				// Loop over metrics
				for (k = 0; k < _ga.fitness().values().length(); ++k)
				{
					message.add("  ");
					message.add(_ga.fitness().values()[k][j]);
				}
				
				// Break since found
				break;
			}
		}
	}
	
	// Print fitness
	Output::newline();
	Output::print(message, RIGHT);
}
