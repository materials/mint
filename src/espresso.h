/* espresso.h -- Quantum Espresso functionality
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef ESPRESSO_H
#define ESPRESSO_H



#include "iso.h"
#include "elements.h"
#include "kpoints.h"
#include "text.h"
#include "list.h"
#include "num.h"
#include "fileSystem.h"



// Namespace to store quantum espresso functionality
namespace Espresso
{

	// Classes
	class Job;
	class Files;
	
	// Enumerated types
	enum Calculation {EC_NONE, EC_SINGLE, EC_RELAX_ALL, EC_RELAX_POSITIONS};
	enum Accuracy {EA_UNKNOWN, EA_HIGH, EA_NORMAL, EA_LOW};
	enum Material {EM_UNKNOWN, EM_METAL, EM_INSULATOR};

	// Functions
	Accuracy accuracy(const Word& word);
	Word accuracy(Accuracy accuracy);
}



// Class for running a Quantum Espresso job
class Espresso::Job
{
	
	// Control variables
	bool _saveFiles;
	Word _executable;
	Word _potDirectory;
	Calculation _calculation;
	Accuracy _accuracy;
	Material _material;
	OList<Element> _elements;
	OList<Word> _potentials;
	
	// Variables to store results
	bool _completed;
	List<double> _energy;
	OList<Vector3D >::D2 _forces;
	ISO _structure;
	
	// Static member variables
	static Word _dirName;
	static Word _inputFileName;
	static Word _outputFileName;
	
	// Functions
	void gatherResults(const Basis* basis);
	
public:
	
	// Constructor
	Job();
	
	// Setup functions
	void saveFiles(bool saveFiles)					{ _saveFiles = saveFiles; }
	void executable(const Word& exe)				{ _executable = exe; }
	void potDirectory(const Word& dir)				{ _potDirectory = dir; }
	void calculation(Calculation calculation)		{ _calculation = calculation; }
	void accuracy(Accuracy accuracy)				{ _accuracy = accuracy; }
	void material(Material material)				{ _material = material; }
	void elements(const OList<Element>& elements)	{ _elements = elements; }
	void potentials(const OList<Word>& potentials)	{ _potentials = potentials; }
	
	// Access functions
	const Word& executable() const			{ return _executable; }
	const Word& potDirectory() const		{ return _potDirectory; }
	Calculation calculation() const			{ return _calculation; }
	Accuracy accuracy() const				{ return _accuracy; }
	Material material() const				{ return _material; }
	const OList<Element>& elements() const	{ return _elements; }
	const OList<Word>& potentials() const	{ return _potentials; }
	
	// Functions
	bool submit(const ISO& iso, bool quitIfError = true);
	
	// Access functions
	bool completed() const						{ return _completed; }
	const List<double>& energy() const			{ return _energy; }
	const OList<Vector3D>::D2& forces() const	{ return _forces; }
	void updateStructure(ISO& iso) const;
};



// Class for reading and writing Quantum Espresso files
class Espresso::Files
{
	
	// Keywords
	enum Keyword {EK_UNKNOWN, EK_CONTROL, EK_SYSTEM, EK_ELECTRONS, EK_IONS, EK_CELL, EK_SPECIES, EK_POSITIONS, \
		EK_KPOINTS, EK_PARAMETERS, EK_CONSTRAINTS, EK_OCCUPATIONS};
	
	// Lattice description types
	enum ValueType {VT_ALAT, VT_BOHR, VT_ANGSTROM, VT_CRYSTAL};
	
	// Functions
	static Keyword keyword(const Word& word);
	
	// Helper functions
	static bool isComment(const Word& word) { return (word[0] == '!'); }
	static void setBasis(ISO& iso, double A, double B, double C, double cosAB, double cosAC, double cosBC, \
		bool showOutput);
	static void setBasis(ISO& iso, double m11, double m12, double m13, double m21, double m22, double m23, \
		double m31, double m32, double m33, bool showOutput);
	
public:
	
	// Read structure
	static ISO read(Text content, bool readAsInput = true, const Basis* basis = 0);
	static ISO read(const Word& file, bool readAsInput = true, const Basis* basis = 0)
		{ return read(Read::text(file), readAsInput, basis); }
	
	// Write structure
    static void write(const Word& file, const ISO& iso, CoordinateType coordinates, Job* job = 0, \
		KPoints* kpoints = 0);

	// Check if file is in correct format
	static bool isFormat(const Word& file)		{ return isFormat(Read::text(file)); }
    static bool isFormat(const Text& content);
};



// =====================================================================================================================
// Job
// =====================================================================================================================

/* inline Espresso::Job::Job()
 *
 * Constructor for quantum espresso job
 */

inline Espresso::Job::Job()
{
	_saveFiles = false;
	_calculation = EC_SINGLE;
	_accuracy = EA_NORMAL;
	_material = EM_METAL;
	_completed = false;
}



// =====================================================================================================================
// Files
// =====================================================================================================================

/* inline void Espresso::Files::setBasis(ISO& iso, double A, double B, double C, double cosAB, double cosAC,
 *		double cosBC, bool showOutput)
 *
 * Set metric values for basis
 */

inline void Espresso::Files::setBasis(ISO& iso, double A, double B, double C, double cosAB, double cosAC, \
	double cosBC, bool showOutput)
{
	iso.basis(Vector3D(A, B, C), Vector3D(acos(cosBC), acos(cosAC), acos(cosAB)), showOutput);
}



/* inline void Espresso::Files::setBasis(ISO& iso, double m11, double m12, double m13, double m21, double m22,
 *		double m23, double m31, double m32, double m33, bool showOutput)
 *
 * Set vector values for basis
 */

inline void Espresso::Files::setBasis(ISO& iso, double m11, double m12, double m13, double m21, double m22, \
	double m23, double m31, double m32, double m33, bool showOutput)
{
	iso.basis(Matrix3D(m11, m12, m13, m21, m22, m23, m31, m32, m33), showOutput);
}



#endif
