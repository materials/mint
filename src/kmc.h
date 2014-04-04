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



#ifndef KMC_H
#define KMC_H



#include "num.h"
#include "iso.h"
#include "elements.h"
#include "structureIO.h"
#include "text.h"
#include "random.h"
#include "fileSystem.h"
#include "constants.h"
#include "output.h"
#include "list.h"



// Class to deal with KMC simulations
class KMC
{
	
	// Classes to store data
	class StorageNode;
	class Node;
	class Attempt;
	class Tracer;
	
	// Variables to store results
	double _activationEnergy;
	Vector3D _activationEnergyVector;
	double _D0;
	Vector3D _D0vector;
	
	// Variables to use during simulation
	bool _runComplement;
	Basis _basis;
	Matrix3D _idealToUnit;
	OList<StorageNode> _nodes;
	OList<Attempt>::D2 _attempts;
	
	// Settings variables
	double _convergence;
	double _jumpsPerAtom;
	
	// Functions
	void setJumps(const Text& input);
	void setSimulation(const Text& input);
	
	// Run functions
	void runPrimary(Random& random);
	void runComplement(Random& random);
	void runSimulations(Random& random, const List<double>& sizes, const List<double>& invTemps);
	Vector3D runSingleSimulation(Random& random, int numNodes, Node* curNode, const OList<Tracer>& tracers);
	void createSimulationNodes(OList<Node>& simNodes, double cellSize, Basis& simCellBasis);
	void initializeTracers(OList<Node>& simNodes, OList<Tracer>& tracers, int start);
	void analyzeResults(Basis& unitBasis, const List<double>& invTemps, const OList<Vector3D>& Dfrac);
	
	// Helper functions
	void getTransformations(const Text& input);
	double concentrationPerCell(double invTemp);
	static int pickStartNode(const OList<Node>& simNodes, double invTemp, Random& random);
	
	// Setup helper functions
	void makeStructures(ISO& start, ISO& final, ISO startInitial, const ISO& startRelaxed, ISO finalInitial, \
		const ISO& finalRelaxed, const Matrix3D& rotStart, const Vector3D& transStart, const Matrix3D& rotFinal, \
		const Vector3D& transFinal);
	void transformStructure(ISO& iso, const Matrix3D& rotation, const Vector3D& translation);
	
public:
	
	// Constructors
	KMC();
	
	// Settings functions
	void convergence(double input)		{ _convergence = input; }
	void jumpsPerAtom(double input)		{ _jumpsPerAtom = input; }
	
	// Setup functions
	void clear()	{ _nodes.clear(); _attempts.clear(); _runComplement = false; }
	static void generateJumps(const ISO& iso, const Element& element, bool useInterstitials, double minImageDis, \
		double tol, StructureFormat format, double maxJumpDistance = 7.0);
	
	// Run functions
	bool set(const Text& input);
	bool set(const Word& file)	{ return set(Read::text(file)); }
	void run(Random& random);
	
	// Print functions
	void print(const Word& file) const;
	
	// Access functions
	double D0() const								{ return _D0; }
	double activationEnergy() const					{ return _activationEnergy; }
	const Vector3D& activationEnergyVector() const	{ return _activationEnergyVector; }
	const Vector3D& D0vector() const				{ return _D0vector; }
	
	// Other functions
	static bool isKMCFile(const Text& content);
	static bool isKMCFile(const Word& file)		{ return isKMCFile(Read::text(file)); }
};



// Class to store information about a node during setup of KMC simulation
class KMC::StorageNode
{
	
	// Variables
	double _formationEnergy;
	Vector3D _position;
	OList<Vector3D> _vectors;
	List<Attempt*> _attempts;
	
public:
	
	// Setup functions
	void formationEnergy(double input)		{ _formationEnergy = input; }
	void position(const Vector3D& input)	{ _position = input; }
	void add(const Vector3D& vector, Attempt* attempt);
	
	// Access functions
	double formationEnergy() const			{ return _formationEnergy; }
	const Vector3D& position() const		{ return _position; }
	const OList<Vector3D>& vectors() const	{ return _vectors; }
	const List<Attempt*>& attempts() const	{ return _attempts; }
};



// Class to store information about a node during a KMC simulation
class KMC::Node
{
	
	// Variables
	double _formationEnergy;
	Tracer* _tracer;
	Vector3D _position;
	OList<Vector3D> _vectors;
	List<Attempt*> _attempts;
	List<Node*> _endNodes;
	List<double> _rates;
	
public:
	
	// Setup functions
	void formationEnergy(double input)		{ _formationEnergy = input; }
	void position(const Vector3D& input)	{ _position = input; }
	void addVector(const Vector3D& input)	{ _vectors += input; }
	void addAttempt(Attempt* input)			{ _attempts += input; }
	void addEndNode(Node* input)			{ _endNodes += input; }
	void setRates();
	
	// Access functions
	Tracer*& tracer()						{ return _tracer; }
	double formationEnergy() const			{ return _formationEnergy; }
	const Vector3D& position() const		{ return _position; }
	const OList<Vector3D>& vectors() const	{ return _vectors; }
	const List<Node*>& endNodes() const		{ return _endNodes; }
	const List<double>& rates() const		{ return _rates; }
};



// Class to store information about a transition attempt
class KMC::Attempt
{
	
	// Variables
	double _rate;
	double _frequency;
	double _barrier;
	List<double> _gsModes;
	List<double> _tsModes;
	
public:
	
	// Constructor
	Attempt();
	
	// Setup
	void set(double barrier);
	void set(double barrier, double frequency);
	void set(double barrier, const Word& gsModesFile, const Word& tsModesFile);
	
	// Set the temperature
	double setRate(double temperature);
	
	// Access functions
	double rate() const	{ return _rate; }
};



// Class to store information about a tracer during KMC simulation
class KMC::Tracer
{
	int _numJumps;
	Vector3D _vector;	
public:
	Tracer()								{ _vector = 0.0; _numJumps = 0; }
	void add(const Vector3D& input)			{ _vector += input; ++_numJumps; }
	void subtract(const Vector3D& input)	{ _vector -= input; ++_numJumps; }
	int numJumps() const					{ return _numJumps; }
	const Vector3D& vector() const			{ return _vector; }
};



// =====================================================================================================================
// KMC
// =====================================================================================================================

/* inline KMC::KMC()
 *
 * Constructor for KMC simulation object
 */

inline KMC::KMC()
{
	_runComplement = false;
	_convergence = 0.5;
	_jumpsPerAtom = 100.0;
}



/* inline double KMC::concentrationPerCell(double invTemp)
 *
 * Return the concentration of occupied nodes at a given temperature
 */

inline double KMC::concentrationPerCell(double invTemp)
{
	double res = 0;
	for (int i = 0; i < _nodes.length(); ++i)
		res += exp(-_nodes[i].formationEnergy() * invTemp / Constants::kb);
	return res;
}



/* inline int KMC::pickStartNode(const OList<Node>& simNodes, double invTemp, Random& random)
 *
 * Pick starting node based on site formation energies
 */

inline int KMC::pickStartNode(const OList<Node>& simNodes, double invTemp, Random& random)
{
	int i;
	double total;
	for (i = 0, total = 0; i < simNodes.length(); ++i)
		total += exp(-simNodes[i].formationEnergy() * invTemp / Constants::kb);
	double ran = random.decimal(0, total);
	for (i = 0, total = 0; i < simNodes.length(); ++i)
	{
		total += exp(-simNodes[i].formationEnergy() * invTemp / Constants::kb);
		if (ran <= total)
			return i;
	}
	return simNodes.length() - 1;
}



/* inline bool KMC::isKMCFile(const Text& content)
 *
 * Return whether a file contains kmc simulation data
 */

inline bool KMC::isKMCFile(const Text& content)
{
	
	// Loop over lines in file and look for at least one line that starts with SITE and one with JUMP
	bool foundSite = false;
	bool foundJump = false;
	for (int i = 0; i < content.length(); ++i)
	{
		if (content[i].length() > 0)
		{
			if (!foundSite)
			{
				if (content[i][0].equal("site", false))
					foundSite = true;
			}
			if (!foundJump)
			{
				if (content[i][0].equal("jump", false))
					foundJump = true;
			}
			if ((foundSite) && (foundJump))
				return true;
		}
	}
	
	// Return that not a match if at this point
	return false;
}



// =====================================================================================================================
// KMC::StorageNode
// =====================================================================================================================

/* inline void KMC::StorageNode::add(const Vector3D& vector, Attempt* attempt)
 *
 * Add a transition to a storage node
 */

inline void KMC::StorageNode::add(const Vector3D& vector, Attempt* attempt)
{
	_vectors += vector;
	_attempts += attempt;
}



// =====================================================================================================================
// KMC::Node
// =====================================================================================================================

/* inline void KMC::Node::setRates()
 *
 * Set cumulative rates for a node
 */

inline void KMC::Node::setRates()
{
	_rates.length(_attempts.length());
	if (_attempts.length())
		_rates[0] = _attempts[0]->rate();
	for (int i = 0, j = 1; j < _attempts.length(); ++i, ++j)
		_rates[j] = _rates[i] + _attempts[j]->rate();
}



// =====================================================================================================================
// KMC::Attempt
// =====================================================================================================================

/* inline KMC::Attempt::Attempt()
 *
 * Constructor for KMC attempt
 */

inline KMC::Attempt::Attempt()
{
	_frequency = 1e13;
	_barrier = 1e10;
}



/* inline void KMC::Attempt::set(double barrier)
 *
 * Set barrier for kmc transition
 */

inline void KMC::Attempt::set(double barrier)
{
	_barrier = barrier;
	if (_barrier <= 0)
	{
		Output::newline(ERROR);
		Output::print("Found a barrier of ");
		Output::print(barrier);
		Output::print(" eV when expecting a value greater than 0");
		Output::quit();
	}
}



/* inline void KMC::Attempt::set(double barrier, double frequency)
 *
 * Set barrier and attempt frequency for kmc transition
 */

inline void KMC::Attempt::set(double barrier, double frequency)
{
	set(barrier);
	_frequency = frequency;
}



#endif
