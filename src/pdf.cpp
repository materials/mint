/* pdf.cpp -- Interpret powder diffraction file
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "pdf.h"
#include "elements.h"
#include "num.h"
#include "language.h"
#include "output.h"
#include <cstdlib>



/* void PDF::read(Diffraction& diffraction, Text content, ISO* iso)
 *
 * Read PDF file and save structure information and diffraction pattern
 */

void PDF::read(Diffraction& diffraction, Text content, ISO* iso)
{
	
	// Output
	Output::newline();
	Output::print("Reading PDF file");
	Output::increase();
	
	// Clear space
	if (iso)
		iso->clear();
	diffraction.clear();
	
	// Remove brackets
	content.split('<');
	content.split('>');
	
	// Variables to store cell information
	int totalAtoms = 0;
	Word spaceGroup;
	OList<Element> elements;
	List<double> composition;
	Vector3D lengths(0.0);
	Vector3D angles(0.0);
	
	// Variables to store diffraction data
	double wavelength = 0;
	Linked<double> theta;
	Linked<double> intensity;
	
	// Loop over contents
	int i, j, k;
	bool foundTheta = false;
	bool foundIntensity = false;
	Word tempWord;
	Word tempNumber;
	double curTheta;
	double curIntensity;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is too short
		if (content[i].length() == 0)
			continue;
		
		// Skip if line is commented
		if (Language::isComment(content[i][0]))
			continue;
		
		// Found intensity
		if (content[i][0] == "intensity")
		{
			
			// Start of new peak
			if (content[i].length() == 1)
			{
				foundTheta = false;
				foundIntensity = false;
			}
			
			// Save peak
			else if (content[i].length() > 1)
			{
				content[i][1].strip('m');
				curIntensity = atof(content[i][1].array());
				foundIntensity = true;
			}
		}
		
		// Found theta
		else if ((content[i][0] == "theta") && (content[i].length() > 1))
		{
			curTheta = atof(content[i][1].array());
			foundTheta = true;
		}
		
		// Found end of intensity
		else if (content[i][0] == "/intensity")
		{
			
			// Intensity was not set
			if (!foundIntensity)
			{
				Output::newline(ERROR);
				Output::print("Did not find intensity of current peak");
				Output::quit();
			}
			
			// Theta was not set
			if (!foundTheta)
			{
				Output::newline(ERROR);
				Output::print("Did not find theta of current peak");
				Output::quit();
			}
			
			// Save peak
			if (curIntensity != 0)
			{
				theta += curTheta;
				intensity += curIntensity;
			}
		}
		
		// Found lattice parameter
		else if ((content[i][0] == "xtla") && (content[i].length() > 1))
			lengths[0] = atof(content[i][1].array());
		else if ((content[i][0] == "xtlb") && (content[i].length() > 1))
			lengths[1] = atof(content[i][1].array());
		else if ((content[i][0] == "xtlc") && (content[i].length() > 1))
			lengths[2] = atof(content[i][1].array());
		else if ((content[i][0] == "xtlal") && (content[i].length() > 1))
			angles[0] = Num<double>::toRadians(atof(content[i][1].array()));
		else if ((content[i][0] == "xtlbe") && (content[i].length() > 1))
			angles[1] = Num<double>::toRadians(atof(content[i][1].array()));
		else if ((content[i][0] == "xtlga") && (content[i].length() > 1))
			angles[2] = Num<double>::toRadians(atof(content[i][1].array()));
		
		// Found wavelength
		else if (((content[i][0] == "wave_length") || (content[i][0] == "lambda")) && (content[i].length() > 1))
			wavelength = atof(content[i][1].array());
		
		// Found spacegroup
		else if ((content[i][0] == "spgr") && (content[i].length() > 1))
			spaceGroup = content[i][1];
		
		// Found Pearson symbol
		else if ((content[i][0] == "pearson") && (content[i].length() > 1))
		{
			tempWord.length(0);
			for (j = 2; j < content[i][1].length(); ++j)
				tempWord += content[i][1][j];
			if (Language::isNumber(tempWord))
				totalAtoms = atoi(tempWord.array());
		}
		
		// Found atomic percent
		else if (content[i][0] == "atomic_percent")
		{
			
			for (j = 1; j < content[i].length(); ++j)
			{
				if (content[i][j][0] == '/')
					break;
				tempWord.length(0);
				for (k = 0; k < content[i][j].length(); ++k)
				{
					if ((content[i][j][k] == '0') || (content[i][j][k] == '1') || (content[i][j][k] == '2') || \
						(content[i][j][k] == '3') || (content[i][j][k] == '4') || (content[i][j][k] == '5') || \
						(content[i][j][k] == '6') || (content[i][j][k] == '7') || (content[i][j][k] == '8') || \
						(content[i][j][k] == '9'))
						break;
					tempWord += content[i][j][k];
				}
				tempNumber.length(0);
				for (; k < content[i][j].length(); ++k)
					tempNumber += content[i][j][k];
				elements += Element::find(tempWord, false, true);
				composition += atof(tempNumber.array()) / 100.0;
			}
		}
	}
	
	// Lattice parameters were not set
	if (iso)
	{
		for (i = 0; i < 3; ++i)
		{
			if ((lengths[i] == 0) || (angles[i] == 0))
			{
				Output::newline(ERROR);
				Output::print("Basis not properly set using xt(la)(b)(c) and xt(al)(ba)(ga)");
				Output::quit();
			}
		}
	}
	
	// Total number of atoms was not correct
	if ((iso) && (!totalAtoms))
	{
		Output::newline(ERROR);
		Output::print("Could not determine total number of atoms");
		Output::quit();
	}
	
	// Get the number of each element to add
	List<int> numToAdd(composition.length());
	if (iso)
	{
		for (i = 0; i < composition.length(); ++i)
		{
			if (Num<double>::abs(Num<double>::round(composition[i] * totalAtoms, 1) - composition[i] * totalAtoms) > 1e-2)
			{
				Output::newline(ERROR);
				Output::print("There is not an integral number of ");
				Output::print(elements[i].symbol());
				Output::print(" atoms in the structure (");
				Output::print(composition[i] * totalAtoms);
				Output::print(")");
				Output::quit();
			}
			numToAdd[i] = (int) Num<double>::round(composition[i] * totalAtoms, 1);
		}
	}
	
	// Output
	if (iso)
	{
		Output::newline();
		Output::print("Saving structure properties");
		Output::increase();
	}
	
	// Set the space group
	if (iso)
	{
		if (spaceGroup.length())
		{
			Output::newline();
			Output::print("Setting space group to ");
			Output::print(spaceGroup);
			iso->spaceGroup(spaceGroup);
		}
		else
		{
			Output::newline(WARNING);
			Output::print("Space group was not found - assuming P1");
		}
	}
	
	// Set the basis
	if (iso)
		iso->basis(lengths, angles, true);
	
	// Add the atoms
	if (iso)
	{
		for (i = 0; i < elements.length(); ++i)
		{

			// Output
			Output::newline();
			Output::print("Adding ");
			Output::print(numToAdd[i]);
			Output::print(" ");
			Output::print(elements[i].symbol());
			Output::print(" atom");
			if (numToAdd[i] != 1)
				Output::print("s");
			Output::print(" that will be assigned random positions");

			// Add atoms
			for (j = 0; j < numToAdd[i]; ++j)
				iso->addAtom(elements[i]);
		}
	}
	
	// Output
	if (iso)
		Output::decrease();
	Output::newline();
	Output::print("Saving x-ray properties");
	Output::increase();
	
	// Save wavelength
	if (wavelength > 0)
	{
		Output::newline();
		Output::print("Setting wavelength to ");
		Output::print(wavelength);
		Output::print(" Ang");
		diffraction.wavelength(wavelength);
	}
	else
	{
		Output::newline(WARNING);
		Output::print("Wavelength was not found - assuming ");
		Output::print(diffraction.wavelength());
		Output::print(" Ang");
	}
	
	// Print number of peaks
	Output::newline();
	Output::print("Found ");
	Output::print(theta.length());
	Output::print(" peak");
	if (theta.length() != 1)
		Output::print("s");
	Output::increase();
	
	// Print peaks
	Linked<double>::iterator itTheta = theta.begin();
	Linked<double>::iterator itIntens = intensity.begin();
	for (; itTheta != theta.end(); ++itTheta, ++itIntens)
	{
		Output::newline();
		Output::print("Two-theta and intensity of ");
		Output::print(*itTheta, 4);
		Output::print(" ");
		Output::print(*itIntens, 4);
	}
	
	// Save peaks
	diffraction.set(theta, intensity);
	
	// Output
	Output::decrease();
	
	// Output
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* bool PDF::isFormat(const Text& content)
 *
 * Return whether a file is in PDF format
 */

bool PDF::isFormat(const Text& content)
{
	
	// Loop over lines in file
	for (int i = 0; i < content.length(); ++i)
	{
		
		// Skip if line is empty
		if (!content[i].length())
			continue;
		
		// Found <pdfcard>
		if (content[i][0] == "<pdfcard>")
			return true;
	}
	
	// Did not find <pdfcard>
	return false;
}
