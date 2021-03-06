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



#ifndef GAPREDICT_H
#define GAPREDICT_H



#include "ga.h"
#include "iso.h"
#include "symmetry.h"
#include "potential.h"
#include "diffraction.h"
#include "random.h"
#include "text.h"
#include "list.h"



// Fitness metrics
enum GAPredictMetric {GAPM_UNKNOWN, GAPM_POTENTIAL, GAPM_DIFFRACTION};



// Class to store structure and symmetry pair
class ISOSymmetryPair
{
	
	// Variables
	ISO _iso;
	Symmetry _symmetry;
	
public:
	
	// Functions
	ISOSymmetryPair& operator= (const ISOSymmetryPair& rhs);
	
	// Access functions
	ISO& iso()				{ return _iso; }
	Symmetry& symmetry()	{ return _symmetry; }
};



// Class to run a genetic algorithm structure prediction
class GAPredict
{

	// Settings
	int _numScreens;
	int _numSimulations;
	int _maxCrossoverLoops;
	double _minBondFraction;
	double _wyckoffBias;
	double _cellMutationProb;
	double _posMutationProb;
	double _wyckMutationProb;
	double _energyTolerance;
	double _diffractionTolerance;
	bool _userietveld;
	bool _saveAllResults; // Whether to print out each candidate structure
	GAPredictMetric _optMetric;
	GAPredictMetric _screenMetric;
	
	// Variables
	GeneticAlgorithm<ISOSymmetryPair, GAPredict> _ga;
	List<GAPredictMetric> _metrics;
	
	// Storage
	Potential* _potential;
	CalculatedPattern* _diffraction;
	Diffraction* _refDiffraction;
	
	// Functions
	void printFitness(int generation);
	
	// Mutation functions
	void mutateBasis(ISOSymmetryPair& pair, Random& random);
	void mutatePositions(ISOSymmetryPair& pair, Random& random);
	void mutateWyckoff(ISOSymmetryPair& pair, Random& random);
	
	// Crossover functions
	void crossover(ISOSymmetryPair& pair, int mother, int father, Random& random);
	List<Orbit*>::D2 getFixedOrbits(int mother, int father);
	List<Orbit*>::D3 groupOrbits(int mother, int father, const List<Orbit*>::D2& fixedOrbits);
	void getOrbitCombinations(Linked<List<int> >& allowedGroups, const List<Orbit*>::D2 orbitGroups, int numToAdd);
	void recurseOrbits(Linked<List<int> >& allowedGroups, const List<int>& lengths, const List<int>& mults, \
		const List<int>& curList, int curGroup, int numToAdd);
	
	// Screen functions
	double screen(ISOSymmetryPair& pair);
	
	void saveResult(int entry);
	
public:
	
	// Constructor
	GAPredict();
	
	// General settings functions
	void numScreens(int input)						{ _numScreens = input; }
	void numSimulations(int input)					{ _numSimulations = input; }
	void maxCrossoverLoops(int input)				{ _maxCrossoverLoops = input; }
	void minBondFraction(double input)				{ _minBondFraction = input; }
	void wyckoffBias(double input)					{ _wyckoffBias = input; }
	void cellMutationProbability(double input)		{ _cellMutationProb = input; }
	void positionMutationProbability(double input)	{ _posMutationProb = input; }
	void wyckoffMutationProbability(double input)	{ _wyckMutationProb = input; }
	void metricToOptimize(GAPredictMetric input)	{ _optMetric = input; }
	void metricToScreen(GAPredictMetric input)		{ _screenMetric = input; }
	void energyTolerance(double input)				{ _energyTolerance = input; }
	void diffractionTolerance(double input)			{ _diffractionTolerance = input; }
	void useRietveld(bool input)					{ _userietveld = input; }
	/**
	 * Set whether to print each candidate structure
     * @param input Desired setting
     */
	void setSaveAllResults(bool input)				{ _saveAllResults = input; }
	
	// GA settings functions
	void populationSize(int input)	{ _ga.populationSize(input); }
	void convergeOver(int input)	{ _ga.convergeOver(input); }
	void maxGenerations(int input)	{ _ga.maxGenerations(input); }
	void numToKeep(int input)		{ _ga.numToKeep(input); }
	
	// Selection settings functions
	void selectionMethod(GASelectionMethod input)	{ _ga.selection().set(input); }
	
	// Run functions
	void run(ISO& iso, Random& random, bool restartable = false,
			Potential* potential = 0, Diffraction* diffraction = 0);
	
	// Used for restarting old runs
	void writeRestartInformation();
	void readRestartInformation(ISO& inputStrc);
	
	// Static member functions
	static GAPredictMetric metric(const Word& input);
	static Word metric(GAPredictMetric input);
};



/* inline GAPredict::GAPredict()
 *
 * Constructor for GAPredict object
 */
inline GAPredict::GAPredict() {
	_numScreens = 0;
	_numSimulations = 1;
	_maxCrossoverLoops = 100;
	_minBondFraction = 0.5;
	_wyckoffBias = 0.5;
	_cellMutationProb = 0.1;
	_posMutationProb = 0.1;
	_wyckMutationProb = 0.1;
	_optMetric = GAPM_UNKNOWN;
	_screenMetric = GAPM_UNKNOWN;
	_energyTolerance = 1e-3;
	_diffractionTolerance = 1e-4;
	_saveAllResults = false;
}



/**
 * Return metric from word
 * @param input [in] Name of metric to be optimized
 * @return Corresponding GAPredictMetric
 */
inline GAPredictMetric GAPredict::metric(const Word& input) {
	if ((input.equal("potential", false, 3)) || (input.equal("energy", false, 4)))
		return GAPM_POTENTIAL;
	if ((input.equal("diffraction", false, 4)) || (input.equal("xray", false)))
		return GAPM_DIFFRACTION;
	return GAPM_UNKNOWN;
}



/**
 * Return word from metric
 * @param input [in] GAPredictMetric value
 * @return Name of that metric
 */
inline Word GAPredict::metric(GAPredictMetric input) {
	switch (input) {
		case GAPM_POTENTIAL:
			return Word("Energy");
		case GAPM_DIFFRACTION:
			return Word("Diffraction match");
		default:
			return Word("Unknown");
	}
}



#endif
