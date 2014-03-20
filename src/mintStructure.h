/* Copyright 2011-2014 Kyle Michel, Logan Ward
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
