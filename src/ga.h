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



#ifndef GA_H
#define GA_H



#include <cmath>
#include "random.h"
#include "text.h"
#include "list.h"
#include "num.h"
#include "random.h"
#include "output.h"
#include "constants.h"



// Different crossover selection methods
enum GASelectionMethod {GA_UNKNOWN, GA_TOURNAMENT, GA_ROULETTE};



// Class to run a genetic algorithm
template <class T, class FunClass>
class GeneticAlgorithm
{
	
public:
	
	// Member classes
	class Fitness;
	class Selection;
	class Crossover;
	class Mutation;
	class Screen;
	
private:

	// Variables
	OList<T> _population;
	Selection _selection;
	Fitness _fitness;
	Crossover _crossover;
	OList<Mutation> _mutations;
	Screen _screen;

	// Settings
	int _convergeOver;
	int _populationSize;
	int _maxGenerations;
	int _numToKeepBetweenGens;
	double _metricTol;
	
	// Convergence variables
	bool _isConverged;
	int _numGenerations;
	int _gensSinceLastBest;
	
	// Storage variables
	int _metricToOpt;
	double _bestFitness;
	T _bestIndividual;
	
public:
	
	// Constructor
	GeneticAlgorithm();
	
	// Setup variables
	void numMetrics(int input);
	void populationSize(int input);
	void convergeOver(int input)		{ _convergeOver = input; }
	void maxGenerations(int input)		{ _maxGenerations = input; }
	void numToKeep(int input)			{ _numToKeepBetweenGens = input; }
	void metricToOptimize(int input)	{ _metricToOpt = input; }
	void metricTolerance(double input)	{ _metricTol = input; }
	
	// Functions
	void initNewRun();
	void advance(Random& random);
	void setFitness();
	
	// Access functions
	bool isConverged() const		{ return _gensSinceLastBest >= _convergeOver; }
	OList<T>& population()			{ return _population; }
	Fitness& fitness()				{ return _fitness; }
	Selection& selection()			{ return _selection; }
	Crossover& crossover()			{ return _crossover; }
	OList<Mutation>& mutations()	{ return _mutations; }
	Screen& screen()				{ return _screen; }
	double bestFitness() const		{ return _bestFitness; }
	T& bestIndividual()				{ return _bestIndividual; }
	
	// Convergence access functions
	/**
     * @return Number of generations since GA was started
     */
	int numGenerations() const		{ return _numGenerations; } 
	int gensSinceLastBest() const	{ return _gensSinceLastBest; }
	int maxGenerations() const		{ return _maxGenerations; }
	
	// Operations used when restarting operation
	void setGensSinceLastBest(int input) { _gensSinceLastBest = input; }
	void setGenNumber(int input) { _numGenerations = input; }
	void setBestIndividual(T individual, double fitness) {
		_bestIndividual = individual;
		_bestFitness = fitness;
	}
	
	// Static member functions
	static GASelectionMethod selection(const Word& input);
	static Word selection(GASelectionMethod input);
};



// Class to store fitness functions
template <class T, class FunClass>
class GeneticAlgorithm<T, FunClass>::Fitness
{
	
	// Variables
	List<bool> _optLowest;
	OList<Word> _names;
	OList<Word> _units;
	List<double>::D2 _values;
	List<int>::D2 _order;
	List<int> _totalOrder;

	// Functions
	void sortByFitness(int metric);
	void sortByTotalRank(List<int>& totals, int metricToOpt);

public:
	
	// Functions
	void sort(int metricToOpt);
	
	// Access functions
	List<bool>& optLowest()		{ return _optLowest; }
	OList<Word>& names()		{ return _names; }
	OList<Word>& units()		{ return _units; }
	/**
     * @return Fitness values for each entry as 2D, where [m][k] is the value
	 * of entry k for metric m.
     */
	List<double>::D2& values()	{ return _values; }
	List<int>::D2& order()		{ return _order; }
	List<int>& totalOrder()		{ return _totalOrder; }
};



// Class to store selection schemes
template <class T, class FunClass>
class GeneticAlgorithm<T, FunClass>::Selection
{
	
	// Selection methods
	void roulette(int& mother, int& father, int numInPopulation, Random& random);
	void tournament(int& mother, int& father, int numInPopulation, Random& random);

	// Pointer to method
	void (Selection::*_function)(int& mother, int& father, int numInPopulation, Random& random);
	
	// Functions
	int genParent(int numChoices, Random& random)
		{ return random.integerOnNormal(0, numChoices - 1, 0, numChoices / 2); }

public:
	
	// Constructor
	Selection()	{ set(GA_TOURNAMENT); }
	
	// Setup function
	void set(GASelectionMethod method);
	
	// Functions
	void generate(int& mother, int& father, int numInPopulation, Random& random)
		{ (this->*_function)(mother, father, numInPopulation, random); }
};



// Class to store crossover
template <class T, class FunClass>
class GeneticAlgorithm<T, FunClass>::Crossover
{
	
	// Variables
	void (FunClass::*_function)(T& child, int mother, int father, Random& random);
	FunClass* _obj;
	
public:
	
	// Setup
	void set(FunClass* obj, void (FunClass::*function)(T& child, int mother, int father, Random& random))
		{ _obj = obj; _function = function; }
	
	// Evaluate
	void operator() (T& child, int mother, int father, Random& random)
		{ (*_obj.*_function)(child, mother, father, random); }
};



// Class to store a mutation
template <class T, class FunClass>
class GeneticAlgorithm<T, FunClass>::Mutation
{
	
	// Variables
	double _probability;
	void (FunClass::*_function)(T& individual, Random& random);
	FunClass* _obj;
	
public:
	
	// Setup
	void set(FunClass* obj, void (FunClass::*function)(T& individual, Random& random), double probability)
		{ _obj = obj; _function = function; _probability = probability; }
	
	// Evaluation
	void operator() (T& individual, Random& random)
		{ (*_obj.*_function)(individual, random); }
	
	// Access
	double probability() const	{ return _probability; }
};



// Class to store a screening functions
template <class T, class FunClass>
class GeneticAlgorithm<T, FunClass>::Screen
{

	// Variables
	bool _optLowest;
	int _numTrials;
	double (FunClass::*_function)(T& individual);
	FunClass* _obj;
	
public:
	
	// Constructor
	Screen()	{ _numTrials = 0; }
	
	// Setup
	void set(FunClass* obj, double (FunClass::*function)(T& individual), bool optLowest, int numTrials)
		{ _obj = obj; _function = function; _optLowest = optLowest; _numTrials = numTrials; }
	
	// Evaluation
	double operator() (T& individual)
		{ return (*_obj.*_function)(individual); }
	
	// Access functions
	bool optLowest() const	{ return _optLowest; }
	int numTrials() const	{ return _numTrials; }
};



// =====================================================================================================================
// GeneticAlgorithm
// =====================================================================================================================

/* inline GeneticAlgorithm::GeneticAlgorithm()
 *
 * Constructor for GeneticAlgorithm object
 */

template <class T, class FunClass>
inline GeneticAlgorithm<T, FunClass>::GeneticAlgorithm()
{
	_convergeOver = 100;
	_populationSize = 0;
	_maxGenerations = 1000;
	_numToKeepBetweenGens = 0;
	_isConverged = false;
	_numGenerations = 0;
	_gensSinceLastBest = 0;
	_metricToOpt = 0;
	_metricTol = 0;
}



/* inline void GeneticAlgorithm::initNewRun()
 *
 * Set variables for new run
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::initNewRun()
{
	_isConverged = false;
	_numGenerations = 0;
	_gensSinceLastBest = 0;
}



/* inline void GeneticAlgorithm::numMetrics(int input)
 *
 * Set the number of metrics that will be optimized
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::numMetrics(int input)
{
	_fitness.optLowest().length(input);
	_fitness.names().length(input);
	_fitness.units().length(input);
	_fitness.values().length(input);
	_fitness.order().length(input);
	for (int i = 0; i < input; ++i)
	{
		_fitness.values()[i].length(_populationSize);
		_fitness.order()[i].length(_populationSize);
	}
}



/* inline void GeneticAlgorithm::populationSize(int input)
 *
 * Set the genetic algorithm population size
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::populationSize(int input)
{
	_populationSize = input;
	_population.length(input);
	for (int i = 0; i < _fitness.values().length(); ++i)
	{
		_fitness.values()[i].length(input);
		_fitness.order()[i].length(input);
	}
	_fitness.totalOrder().length(input);
}



/* inline GASelectionMethod GeneticAlgorithm::selection(const Word& input)
 *
 * Get selection method from word
 */

template <class T, class FunClass>
inline GASelectionMethod GeneticAlgorithm<T, FunClass>::selection(const Word& input)
{
	if (input.equal("tournament", false, 5))
		return GA_TOURNAMENT;
	if (input.equal("roulette", false, 4))
		return GA_ROULETTE;
	return GA_UNKNOWN;
}


/* inline Word GeneticAlgorithm::selection(GASelectionMethod input)
 *
 * Get word from selection method
 */

template <class T, class FunClass>
inline Word GeneticAlgorithm<T, FunClass>::selection(GASelectionMethod input)
{
	switch (input)
	{
		case GA_TOURNAMENT:
			return Word("Tournament");
		case GA_ROULETTE:
			return Word("Roulette");
		default:
			return Word("Unknown");
	}
}



/* void GeneticAlgorithm::advance(Random& random)
 *
 * Create new generation in genetic algorithm
 */

template <class T, class FunClass>
void GeneticAlgorithm<T, FunClass>::advance(Random& random)
{
	
	// Return if this is the first loop
	if (++_numGenerations == 1)
		return;
	
	// Decide how many children to make
	int numToMake = (_screen.numTrials() > 1) ? _screen.numTrials() : 1;
	OList<T> trials(numToMake);
	List<double> trialValues(numToMake);
	
	// Loop over children to add
	int i, j, k;
	int best;
	int mother;
	int father;
	OList<T> children(_populationSize - _numToKeepBetweenGens);
	for (i = 0; i < children.length(); ++i)
	{
		
		// Loop over trials
		best = 0;
		for (j = 0; j < numToMake; ++j)
		{
			
			// Perform crossover
			_selection.generate(mother, father, _populationSize, random);
			_crossover(trials[j], mother, father, random);
			
			// Make mutations
			for (k = 0; k < _mutations.length(); ++k)
			{
				if (random.decimal(0, 1) < _mutations[k].probability())
					_mutations[k](trials[j], random);
			}
			
			// Save best if needed
			if (numToMake > 1)
			{
				trialValues[j] = _screen(trials[j]);
				if (_screen.optLowest())
				{
					if (trialValues[j] < trialValues[best])
						best = j;
				}
				else
				{
					if (trialValues[j] > trialValues[best])
						best = j;
				}
			}
		}
		
		// Save child
		children[i] = trials[best];
	}
	
	// Replace individuals with children
	for (i = _populationSize - 1; i >= _numToKeepBetweenGens; --i)
	{
		for (j = 0; j < _populationSize; ++j)
		{
			if (_fitness.totalOrder()[j] == i)
			{
				_population[j] = children[_populationSize - i - 1];
				break;
			}
		}
	}
}



/**
 * Sort fitness values and check for convergence
 */
template <class T, class FunClass>
void GeneticAlgorithm<T, FunClass>::setFitness()
{
	
	// Sort the fitness
	_fitness.sort(_metricToOpt);
	
	// Loop over population for fitness metric to optimize
	++_gensSinceLastBest;
	bool newBest = false;
	for (int i = 0; i < _fitness.order()[_metricToOpt].length(); ++i)
	{
		
		// Found best
		if (_fitness.order()[_metricToOpt][i] == 0)
		{
			
			// Check if a new best
			newBest = false;
			if (_fitness.optLowest()[_metricToOpt])
			{
				if (_fitness.values()[_metricToOpt][i] < _bestFitness - _metricTol)
					newBest = true;
			}
			else
			{
				if (_fitness.values()[_metricToOpt][i] > _bestFitness + _metricTol)
					newBest = true;
			}
			
			// This is the first loop or found a new best
			if ((_numGenerations == 1) || (newBest))
			{
				_bestFitness = _fitness.values()[_metricToOpt][i];
				_bestIndividual = _population[i];
				_gensSinceLastBest = 0;
			}
			
			// Break since best was found
			break;
		}
	}
	
	// Check for convergence
	if (_gensSinceLastBest >= _convergeOver)
		_isConverged = true;
}



/* void GeneticAlgorithm::Fitness::sort(int metricToOpt)
 *
 * Sort fitness variables
 */

template <class T, class FunClass>
void GeneticAlgorithm<T, FunClass>::Fitness::sort(int metricToOpt)
{
	
	// Rank each individual under each metric
	int i;
	for (i = 0; i < _values.length(); ++i)
		sortByFitness(i);
	
	// Get total for each individual
	int j;
	List<int> totals;
	if (_order.length())
	{
		totals.length(_order[0].length());
		totals.fill(0);
		for (i = 0; i < _order[0].length(); ++i)
		{
			for (j = 0; j < _order.length(); ++j)
				totals[i] += _order[j][i];
		}
	}
	
	// Sort by total
	sortByTotalRank(totals, metricToOpt);
}



/* void GeneticAlgorithm::Fitness::sortByFitness(int metric)
 *
 * Sort fitness metrics
 */

template <class T, class FunClass>
void GeneticAlgorithm<T, FunClass>::Fitness::sortByFitness(int metric)
{
	
	// Fill with -1
	_order[metric].fill(-1);
	
	// Loop over ranks
	int i, j;
	int curIndex;
	double curValue;
	for (i = 0; i < _values[metric].length(); ++i)
	{
		
		// Loop until the first non-set value is found
		for (j = 0; j < _order[metric].length(); ++j)
		{
			if (_order[metric][j] == -1)
			{
				curIndex = j;
				curValue = _values[metric][j];
				break;
			}
		}
		
		// Loop over remaining values
		for (; j < _order[metric].length(); ++j)
		{
			if (_order[metric][j] != -1)
				continue;
			if (_optLowest[metric])
			{
				if (_values[metric][j] < curValue)
				{
					curIndex = j;
					curValue = _values[metric][j];
				}
			}
			else
			{
				if (_values[metric][j] > curValue)
				{
					curIndex = j;
					curValue = _values[metric][j];
				}
			}
		}
		
		// Save order
		_order[metric][curIndex] = i;
	}
}



/* void GeneticAlgorithm::Fitness::sortByTotalRank(List<int>& totals, int metricToOpt)
 *
 * Quicksort fitness values by total rank
 */

template <class T, class FunClass>
void GeneticAlgorithm<T, FunClass>::Fitness::sortByTotalRank(List<int>& totals, int metricToOpt)
{
	
	// Fill with -1
	_totalOrder.fill(-1);
	
	// Loop over ranks
	int i, j;
	int curIndex;
	int curValue;
	for (i = 0; i < totals.length(); ++i)
	{
		
		// Loop until the first non-set value is found
		for (j = 0; j < _totalOrder.length(); ++j)
		{
			if (_totalOrder[j] == -1)
			{
				curIndex = j;
				curValue = totals[j];
				break;
			}
		}
		
		// Loop over remaining values
		for (; j < _totalOrder.length(); ++j)
		{
			if (_totalOrder[j] != -1)
				continue;
			if (totals[j] < curValue)
			{
				curIndex = j;
				curValue = totals[j];
			}
			else if (totals[j] == curValue)
			{
				if (_optLowest[metricToOpt])
				{
					if (_values[metricToOpt][j] < _values[metricToOpt][curIndex])
						curIndex = j;
				}
				else
				{
					if (_values[metricToOpt][j] > _values[metricToOpt][curIndex])
						curIndex = j;
				}
			}
		}
		
		// Save order
		_totalOrder[curIndex] = i;
	}
}



// =====================================================================================================================
// Selection
// =====================================================================================================================

/* inline void GeneticAlgorithm::Selection::set(GASelectionMethod method)
 *
 * Set selection method
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::Selection::set(GASelectionMethod method)
{
	switch (method)
	{
		case GA_TOURNAMENT:
			_function = &Selection::tournament;
			break;
		case GA_ROULETTE:
			_function = &Selection::roulette;
			break;
		default:
			Output::newline(ERROR);
			Output::print("Unknown genetic algorithm selection method");
			Output::quit();
	}
}



/* inline void GeneticAlgorithm::Selection::roulette(int& mother, int& father, int numInPopulation, Random& random)
 *
 * Roulette selection in genetic algorithm
 * Select on the distribution Exp{[-x / (numInPopulation / 2)]^2} where x is the rank
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::Selection::roulette(int& mother, int& father, int numInPopulation, \
	Random& random)
{
	
	// Pick the mother
	mother = genParent(numInPopulation, random);
	
	// Pick the father
	do
	{
		father = genParent(numInPopulation, random);
	} while (mother == father);
}



/* inline void GeneticAlgorithm::Selection::tournament(int& mother, int& father, int numInPopulation, Random& random)
 *
 * Tournament selection in genetic algorithm
 */

template <class T, class FunClass>
inline void GeneticAlgorithm<T, FunClass>::Selection::tournament(int& mother, int& father, int numInPopulation, \
	Random& random)
{
	
	// Determine tournament size
	int tournamentSize = (int)Num<double>::ceil(sqrt(numInPopulation) * 0.6);
	int tournamentPop[tournamentSize];
	
	// Make tournament for mother
	int i;
	for (i = 0; i < tournamentSize; ++i)
		tournamentPop[i] = random.integer(0, numInPopulation-1);
	
	// Pick the mother
	mother = (tournamentSize > 0) ? tournamentPop[0] : 0;
	for (i = 1; i < tournamentSize; ++i)
	{
		if (tournamentPop[i] < mother)
			mother = tournamentPop[i];
	}
	
	// Loop until a father is found that is not the same as the mother
	do
	{
		
		// Make tournament for father
		for (i = 0; i < tournamentSize; ++i)
			tournamentPop[i] = random.integer(0, numInPopulation-1);

		// Pick the father
		father = (tournamentSize > 0) ? tournamentPop[0] : 0;
		for (i = 1; i < tournamentSize; ++i)
		{
			if (tournamentPop[i] < father)
				father = tournamentPop[i];
		}
	} while (mother == father);
}



#endif
