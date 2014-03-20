/* kpoints.h -- Create kpoints
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef KPOINTS_H
#define KPOINTS_H



#include "text.h"
#include "num.h"



// Class to store data about kpoints
class KPoints
{
	
public:
	
	// Centering types
	enum Centering {KP_GAMMA, KP_MONKHORST};
	
	// File formats
	enum Format {KP_VASP};
	
private:
  
	// Variables
	int _mesh[3];
	Centering _centering;

	// Functions
	void writeVasp() const;
	
	// Static member variables
	static double _kppraMetalHigh;
	static double _kppraMetalNormal;
	static double _kppraMetalLow;
	static double _kppraInsulatorHigh;
	static double _kppraInsulatorNormal;
	static double _kppraInsulatorLow;
	
public:

	// Constructor
	KPoints()	{ _centering = KP_GAMMA; }

	// Functions
	void mesh(int numAtoms, const Vector3D& lengths, double kppra);
	void mesh(const int* mesh)				{ _mesh[0] = mesh[0]; _mesh[1] = mesh[1]; _mesh[2] = mesh[2]; }
	void centering(Centering centering)		{ _centering = centering; }

	// Access functions
	const int* mesh() const					{ return _mesh; }
	Centering centering() const				{ return _centering; }
	static double kppraMetalHigh()			{ return _kppraMetalHigh; }
	static double kppraMetalNormal()		{ return _kppraMetalNormal; }
	static double kppraMetalLow()			{ return _kppraMetalLow; }
	static double kppraInsulatorHigh()		{ return _kppraInsulatorHigh; }
	static double kppraInsulatorNormal()	{ return _kppraInsulatorNormal; }
	static double kppraInsulatorLow()		{ return _kppraInsulatorLow; }

	// Write functions
	void write(const Word& file, Format format) const;
};



#endif
