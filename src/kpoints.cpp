/* kpoints.cpp -- Create kpoints
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include <cmath>
#include "kpoints.h"
#include "output.h"



// Kpoint densities as points per reciprocal atom
double KPoints::_kppraMetalHigh = 5000;
double KPoints::_kppraMetalNormal = 1000;
double KPoints::_kppraMetalLow = 500;
double KPoints::_kppraInsulatorHigh = 500;
double KPoints::_kppraInsulatorNormal = 100;
double KPoints::_kppraInsulatorLow = 50;



/* void KPoints::mesh(int numAtoms, const Vector3D& lengths, double kppra)
 *
 * Create mesh
 */

void KPoints::mesh(int numAtoms, const Vector3D& lengths, double kppra)
{
	
	// Set the total number of kpoints as # points per reciprocal atom / number of atoms
	// Set the average number of kpoints along each direction as cube root of total number of kpoints
	double avgKPoints = Num<double>::ceil(pow(kppra / numAtoms, 1.0/3.0));
	
	// Get the average lattice vector length
	int i;
	double avgLength = 0;
	for (i = 0; i < 3; i++)
		avgLength += lengths[i];
	avgLength /= 3.0;
	
	// Set the number of kpoints along each direction as average number of kpoints along each direction multiplied
	//		by the fraction of the average length in reciprocal space
	for (i = 0; i < 3; i++)
		_mesh[i] = (int) Num<double>::ceil(Num<double>::round(avgKPoints * avgLength / lengths[i], 0.01));
    
    // Make sure that none of the kpoints are zero
    for (i = 0; i < 3; i++)
		_mesh[i] = (_mesh[i] < 1) ? 1 : _mesh[i];
}



/* void KPoints::write(const Word& file, Format format) const
 *
 * Write points file
 */

void KPoints::write(const Word& file, Format format) const
{
	
	// Open file for writing if needed
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file != "stdout")
		Output::setStream(Output::addStream(file));
	
	// Set output method
	Output::method(STANDARD);
	
	// Write to correct format
	switch (format)
	{
		case KP_VASP:
			writeVasp();
			break;
		default:
			break;
	}
	
	// Reset output
	if (file != "stdout")
		Output::removeStream(Output::streamID());
	Output::setStream(origStream);
	Output::method(origMethod);
}



/* void KPoints::writeVasp() const
 *
 * Write kpoints to vasp file
 */

void KPoints::writeVasp() const
{
	
	// Print header
	Output::newline();
	Output::print("Automatically generated KPOINTS file");
	
	// Write automatic
	Output::newline();
	Output::print("  0");
	
	// Write centering line
	Output::newline();
	if (_centering == KP_MONKHORST)
		Output::print("Monkhorst-Pack");
	else
		Output::print("Gamma");

	// Write mesh
	Output::newline();
	Output::print("  ");
	for (int i = 0; i < 3; ++i)
	{
		Output::print(_mesh[i]);
		Output::print(" ");
	}
}
