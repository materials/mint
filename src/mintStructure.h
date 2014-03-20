/* mintStructure.h -- Structure for mint files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef MINTSTRUCTURE_H
#define MINTSTRUCTURE_H



#include "iso.h"
#include "text.h"
#include "fileSystem.h"



// Mint structure functions
class MintStructure
{
	
	// Labels
	enum Label {UNKNOWN, SPACEGROUP, BASIS, ELEMENT, FRACPOS, CARTPOS, FIXED, MULTIPLICITY, EXPAND, INTERSTITIAL};
	
	// Functions
	static Label getLabel(const Word& word);
	static void addAtom(ISO& iso, const Element& element, const Vector3D& position, bool useFractional, \
		bool setFixed, const bool* fixed, bool interstitial);
	
public:
	
	// Functions
	static ISO read(const Text& content, double tol, double clusterTol);
	static ISO read(const Word& file, double tol, double clusterTol)
		{ return read(Read::text(file), tol, clusterTol); }
    static void write(const Word& file, const ISO& iso, CoordinateType coordinates);
	static bool isFormat(const Word& file)	{ return isFormat(Read::text(file)); }
    static bool isFormat(const Text& content);
};



#endif
