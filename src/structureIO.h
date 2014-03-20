/* structureIO.h -- Input and output of structure files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef STRUCTUREIO_H
#define STRUCTUREIO_H



#include "iso.h"
#include "text.h"
#include "fileSystem.h"



// Types of structure files
enum StructureFormat {SF_UNKNOWN, SF_AUTO, SF_VASP4, SF_VASP5, SF_MINT, SF_CIF, SF_CRYSTALMAKER, SF_FINDSYM, \
	SF_ESPRESSO, SF_JSON};



// Namespace for storing function to make ISO from file
namespace StructureIO
{
	
	// Test if a file is a structure
	StructureFormat getFormat(const Text& content);
	StructureFormat getFormat(const Word& file);
	
	// Read and write files
	ISO read(const Text& content, StructureFormat format, double tol = 1e-4, double clusterTol = 0.2);
	ISO read(const Text& content, double tol = 1e-4, double clusterTol = 0.2);
	ISO read(const Word& file, double tol = 1e-4, double clusterTol = 0.2);
	ISO read(const Word& file, StructureFormat format, double tol = 1e-4, double clusterTol = 0.2);
	void write(const Word& file, const ISO& iso, StructureFormat format, CoordinateType coordinates = FRACTIONAL, \
		double tol = 1e-4);
	
	// Helper functions
	Word structureFormat(StructureFormat format);
	StructureFormat structureFormat(const Words& words);
	void addExtension(Word& word, StructureFormat format);
}



/* inline StructureFormat StructureIO::getFormat(const Word& file)
 *
 * Get file format
 */

inline StructureFormat StructureIO::getFormat(const Word& file)
{
	return getFormat(Read::text(file));
}



/* inline ISO StructureIO::read(const Text& content, double tol, double clusterTol)
 *
 * Get structure from file
 */

inline ISO StructureIO::read(const Text& content, double tol, double clusterTol)
{
	return read(content, getFormat(content), tol, clusterTol);
}



/* inline ISO StructureIO::read(const Word& file, StructureFormat format, double tol, double clusterTol)
 *
 * Get structure from file
 */

inline ISO StructureIO::read(const Word& file, StructureFormat format, double tol, double clusterTol)
{
	return read(Read::text(file), format, tol, clusterTol);
}



#endif
