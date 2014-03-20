/* pdf.h -- Interpret powder diffraction file
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
	static void read(Diffraction& diffraction, Text content, ISO* iso = 0);
	static void read(Diffraction& diffraction, const Word& file, ISO* iso = 0)
		{ read(diffraction, Read::text(file), iso); }

	// Check if file is in correct format
	static bool isFormat(const Text& content);
	static bool isFormat(const Word& file)		{ return isFormat(Read::text(file)); }
};



#endif
