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



#ifndef PDF_H
#define PDF_H



#include "iso.h"
#include "diffraction.h"
#include "text.h"
#include "fileSystem.h"



// Class to interpret powder diffraction file
class PDF
{
public:
	
	// Read file
	static void read(ExperimentalPattern& diffraction, Text content, ISO* iso = 0);
	static void read(ExperimentalPattern& diffraction, const Word& file, ISO* iso = 0)
		{ read(diffraction, Read::text(file), iso); }

	// Check if file is in correct format
	static bool isFormat(const Text& content);
	static bool isFormat(const Word& file)		{ return isFormat(Read::text(file)); }
};



#endif
