/* randomStructure.h -- Generate a random structure
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef RANDOMSTRUCTURE_H
#define RANDOMSTRUCTURE_H



#include "num.h"
#include "iso.h"
#include "spaceGroup.h"
#include "symmetry.h"
#include "random.h"
#include "list.h"



// Class to generate a random structure
class RandomStructure
{
	
	// Static member variables
	static int _maxTrialLoops;
	static double _minBondFraction;
	static double _weights[143];
	
	// General helper functions
	static Matrix3D generateBasis(const Matrix3D origMat, double min, double max, Random& random, \
		const Symmetry* symmetry = 0);
	static Matrix3D generateBasis(double targetVolume, Random& random, \
		LatticeSystem latticeSystem = LS_TRICLINIC);
	static List<Atom*>::D2 getAtomsToAssign(ISO& iso);
	static int weightBin(double ratio)	{ return (int)((ratio + 1e-6 - 0.95)/0.05); }
	static double weight(int index)		{ return ((index >= 143) || (index < 0)) ? 0 : _weights[index]; }
	static double weight(double ratio)	{ return weight(weightBin(ratio)); }
	
	// Generate structure without symmetry
	static void generateWithoutSymmetry(ISO& iso, double targetVolume, Random& random, Symmetry* symmetry);
	static void generatePositions(ISO& iso, Random& random);
	
	// Generate structure with symmetry
	static void generateWithSymmetry(ISO& iso, double targetVolume, Random& random, double bias, Symmetry* symmetry);
	static void generatePositions(ISO& iso, Random& random, const SpaceGroup& spaceGroup, double bias, \
		Symmetry* symmetry);
	static List<Wyckoff*>::D2 getWyckoffGroups(const ISO& iso, const SpaceGroup& spaceGroup, double minRadius);
	static void sortWyckoffGroups(List<Wyckoff*>::D2& groups, int left, int right);
	static void getAllowedWyckoff(Linked<List<int>::D2>& allowedGroups, const List<Atom*>::D2& atomsToAssign, \
		const List<Wyckoff*>::D2 wyckoffGroups);
	static void recurseAllowedWyckoff(Linked<List<int>::D2>& res, const List<Atom*>::D2& atomsToAssign, \
		const List<Wyckoff*>::D2 wyckoffGroups, const List<bool>& limited, const List<int>& numUsed, \
		const List<int>::D2& curList, int curAtom, int curGroup);
	static void setSiteSymmetry(const ISO& iso, const SpaceGroup& spaceGroup, Symmetry* symmetry);
	
public:
	
	// Setup functions
	static void maxTrialLoops(int input)		{ _maxTrialLoops = input; }
	static void minBondFraction(double input)	{ _minBondFraction = input; }
	
	// Perturb functions
	static void perturbBasis(ISO& iso, double min, double max, Random& random);
	static void perturbAtoms(ISO& iso, double min, double max, Random& random);
	
	// Perturb using symmetry
	static void perturbBasis(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random);
	static void perturbAtoms(ISO& iso, const Symmetry& symmetry, double min, double max, Random& random);
	
	// Generate structure
	static void generate(ISO& iso, Random& random)
		{ generate(iso, random, 1, 0); }
	static void generate(ISO& iso, Random& random, double bias)
		{ generate(iso, random, bias, 0); }
	static void generate(ISO& iso, Random& random, Symmetry* symmetry)
		{ generate(iso, random, 1, symmetry); }
	static void generate(ISO& iso, Random& random, double bias, Symmetry* symmetry);
};



#endif
