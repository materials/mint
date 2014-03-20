/* structureIO.cpp -- Input and output of structure files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "structureIO.h"
#include "mintStructure.h"
#include "crystalMaker.h"
#include "vasp.h"
#include "findsym.h"
#include "espresso.h"
#include "json.h"
#include "cif.h"
#include "language.h"
#include "output.h"



/* StructureFormat StructureIO::getFormat(const Text& content)
 *
 * Get the format of a structure file
 */

StructureFormat StructureIO::getFormat(const Text& content)
{
	
	// Check if mint
	if (MintStructure::isFormat(content))
		return SF_MINT;
	
	// Check if cif
	if (CIF::isFormat(content))
		return SF_CIF;
	
	// Check if vasp 4
	if (Vasp::Structure::isVersion4(content))
		return SF_VASP4;
	
	// Check if vasp 5
	if (Vasp::Structure::isVersion5(content))
		return SF_VASP5;
	
	// Check if quantum espresso
	if (Espresso::Files::isFormat(content))
		return SF_ESPRESSO;
	
	// Return that type is not known
	return SF_UNKNOWN;
}



/* ISO StructureIO::read(const Text& content, StructureFormat format, double tol, double clusterTol)
 *
 * Read contents of structure file
 */

ISO StructureIO::read(const Text& content, StructureFormat format, double tol, double clusterTol)
{
	switch (format)
	{
		case SF_MINT:
			return MintStructure::read(content, tol, clusterTol);
		case SF_CIF:
			return CIF::read(content, clusterTol);
		case SF_VASP4:
			return Vasp::Structure::read(content, false);
		case SF_VASP5:
			return Vasp::Structure::read(content, true);
		case SF_ESPRESSO:
			return Espresso::Files::read(content);
		default:
			Output::newline(ERROR);
			Output::print("Cannot read structure from file with unknown format");
			Output::quit();
			return ISO();
	}
}



/* ISO StructureIO::read(const Word& file, double tol, double clusterTol)
 *
 * Read contents of structure file
 */

ISO StructureIO::read(const Word& file, double tol, double clusterTol)
{
	
	// Get file contents
	Text content = Read::text(file);
	
	// Read structure
	return read(content, getFormat(content), tol, clusterTol);
}



/* void StructureIO::write(const Word& file, const ISO& iso, StructureFormat format, CoordinateType coordinates,
 *		double tol)
 *
 * Write structure to file
 */

void StructureIO::write(const Word& file, const ISO& iso, StructureFormat format, CoordinateType coordinates, \
	double tol)
{
	switch (format)
	{
		case SF_MINT:
			MintStructure::write(file, iso, coordinates);
			break;
		case SF_CIF:
			CIF::write(file, iso, tol);
			break;
		case SF_VASP4:
			Vasp::Structure::write(file, iso, coordinates, false);
			break;
		case SF_VASP5:
			Vasp::Structure::write(file, iso, coordinates, true);
			break;
		case SF_CRYSTALMAKER:
			CrystalMaker::write(file, iso);
			break;
		case SF_FINDSYM:
			FindSym::write(file, iso, tol);
			break;
		case SF_ESPRESSO:
			Espresso::Files::write(file, iso, coordinates);
			break;
		case SF_JSON:
			JSON::write(file, iso);
			break;
		default:
			Output::newline(ERROR);
			Output::print("Cannot print to unknown format");
			Output::quit();
	}
}



/* Word StructureIO::structureFormat(StructureFormat format)
 *
 * Return the name of a format
 */

Word StructureIO::structureFormat(StructureFormat format)
{
	switch (format)
	{
		case SF_AUTO:
			return Word("Automatic");
		case SF_MINT:
			return Word("Mint");
		case SF_CIF:
			return Word("CIF");
		case SF_VASP4:
			return Word("VASP 4");
		case SF_VASP5:
			return Word("VASP 5");
		case SF_CRYSTALMAKER:
			return Word("CrystalMaker");
		case SF_FINDSYM:
			return Word("FindSym");
		case SF_ESPRESSO:
			return Word("Quantum Espresso");
		case SF_JSON:
			return Word("JSON");
		default:
			return Word("unknown");
	}
}



/* StructureFormat StructureIO::structureFormat(const Words& words)
 *
 * Return the format of a file from name
 */

StructureFormat StructureIO::structureFormat(const Words& words)
{
	
	// Loop over words in file
	for (int i = 0; i < words.length(); i++)
	{
		
		// Break if a comment
		if (Language::isComment(words[i]))
			break;
		
		// Found automatic
		if (words[i].equal("automatic", false, 4))
			return SF_AUTO;
		
		// Found mint
		if ((words[i].equal("mint", false, 4)) || (words[i].equal(".mint", false, 5)))
			return SF_MINT;
		
		// Found CIF
		if ((words[i].equal("cif", false)) || (words[i].equal(".cif", false)))
			return SF_CIF;
		
		// Found vasp
		if ((words[i].equal("vasp", false)) || (words[i].equal(".vasp", false)))
		{
			
			// Check for vasp 4
			if (i < words.length() - 1)
			{
				if (words[++i] == "4")
					return SF_VASP4;
			}
			
			// Must be vasp 5
			return SF_VASP5;
		}
		
		// Found vasp 4
		if ((words[i].equal("vasp4", false)) || (words[i].equal(".vasp4", false)))
			return SF_VASP4;
		
		// Found vasp 5
		if ((words[i].equal("vasp5", false)) || (words[i].equal(".vasp5", false)))
			return SF_VASP5;
		
		// Found crystalMaker
		if ((words[i].equal("crystalmaker", false, 5)) || (words[i].equal("cm", false)) || \
			(words[i].equal("cmtx", false)) || (words[i].equal(".cmtx", false)))
			return SF_CRYSTALMAKER;
		
		// Found findsym
		if ((words[i].equal("findsym", false)) || (words[i].equal(".findsym", false)) || \
			(words[i].equal("fs", false)) || (words[i].equal(".fs", false)))
			return SF_FINDSYM;
		
		// Found quantum espresso
		if ((words[i].equal("quantum", false, 4)) || (words[i].equal("espresso", false, 3)) || \
			(words[i].equal("qe", false)) || (words[i].equal(".qe", false)))
			return SF_ESPRESSO;
		
		// Found json
		if ((words[i].equal("json", false)) || (words[i].equal(".json", false)))
			return SF_JSON;
	}
	
	// Did not find a type
	return SF_UNKNOWN;
}



/* void StructureIO::addExtension(Word& word, StructureFormat format)
 *
 * Add extension
 */

void StructureIO::addExtension(Word& word, StructureFormat format)
{
	switch (format)
	{
		case SF_MINT:
			word += ".mint";
			break;
		case SF_CIF:
			word += ".cif";
			break;
		case SF_VASP4:
		case SF_VASP5:
			word += ".vasp";
			break;
		case SF_CRYSTALMAKER:
			word += ".cmtx";
			break;
		case SF_FINDSYM:
			word += ".fs";
			break;
		case SF_ESPRESSO:
			word += ".qe";
			break;
		case SF_JSON:
			word += ".json";
			break;
		default:
			break;
	}
}
