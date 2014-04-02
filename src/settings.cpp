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



#include "settings.h"
#include "language.h"
#include "fileSystem.h"
#include "output.h"
#include <cmath>
#include <cstdlib>



/* bool Setting::set(const OList<Word>& input)
 *
 * Get setting from file
 */

bool Setting::set(const OList<Word>& input)
{
	
	// Check if first word is the tag
	if (!input.length())
		return false;
	if (!isTag(input[0]))
		return false;
	
	// Value is a bool
	if (_valueType == VT_BOOL)
	{
		char temp;
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			temp = tolower(input[i][0]);
			if ((temp == '1') || (temp == 't'))
			{
				_valueBool = true;
				return true;
			}
			if ((temp == '0') || (temp == 'f'))
			{
				_valueBool = false;
				return true;
			}
		}
	}
	
	// Value is an integer
	else if (_valueType == VT_INT)
	{
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			if (Language::isInteger(input[i]))
			{
				_valueInt = atoi(input[i].array());
				return true;
			}
		}
	}
	
	// Value is a double
	else if (_valueType == VT_DOUBLE)
	{
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			if (Language::isNumber(input[i]))
			{
				_valueDouble = atof(input[i].array());
				return true;
			}
		}
	}
	
	// Value is a coordinate type
	else if (_valueType == VT_COORDINATES)
	{
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			if ((input[i].equal("fractional", false, 4)) || (input[i].equal("direct", false, 3)))
			{
				_valueCoordinates = FRACTIONAL;
				return true;
			}
			if (input[i].equal("cartesian", false, 4))
			{	
				_valueCoordinates = CARTESIAN;
				return true;
			}
		}
	}
	
	// Value is a structure format
	else if (_valueType == VT_STRFORMAT)
	{
		_valueStrFormat = StructureIO::structureFormat(input);
		if (_valueStrFormat != SF_UNKNOWN)
			return true;
	}
	
	// Value is a GA structure prediction metric
	else if (_valueType == VT_GAOPTMETRIC)
	{
		GAPredictMetric curMetric;
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			curMetric = GAPredict::metric(input[i]);
			if (curMetric != GAPM_UNKNOWN)
			{
				_valueGAPredictMetric = curMetric;
				return true;
			}
		}
	}
	
	// Value is a GA selection method
	else if (_valueType == VT_GASELECTION)
	{
		GASelectionMethod curMethod;
		for (int i = 1; i < input.length(); ++i)
		{
			if (Language::isComment(input[i]))
				break;
			curMethod = GeneticAlgorithm<int, Setting>::selection(input[i]);
			if (curMethod != GA_UNKNOWN)
			{
				_valueGASelection = curMethod;
				return true;
			}
		}
	}
	
	// Something went wrong
	error(input);
	return false;
}



/* void Setting::printDefault()
 *
 * Print the default value of a setting
 */

void Setting::printDefault()
{
	
	// Print bool
	if (_valueType == VT_BOOL)
	{
		if (_defaultBool == true)
			Output::print("true");
		else
			Output::print("false");
	}
	
	// Print integer
	else if (_valueType == VT_INT)
		Output::print(_defaultInt);
	
	// Print double
	else if (_valueType == VT_DOUBLE)
	{
		if ((fabs(_defaultDouble) < 1e-2) || (fabs(_defaultDouble) > 1e3))
			Output::printSci(_defaultDouble);
		else
			Output::print(_defaultDouble);
	}
	
	// Print coordinate type
	else if (_valueType == VT_COORDINATES)
	{
		if (_defaultCoordinates == FRACTIONAL)
			Output::print("Fractional");
		else
			Output::print("Cartesian");
	}
	
	// Print structure format
	else if (_valueType == VT_STRFORMAT)
		Output::print(StructureIO::structureFormat(_defaultStrFormat));
	
	// Print GA structure prediction metric
	else if (_valueType == VT_GAOPTMETRIC)
		Output::print(GAPredict::metric(_defaultGAPredictMetric));
	
	// Print GA selection method
	else if (_valueType == VT_GASELECTION)
		Output::print(GeneticAlgorithm<int, Setting>::selection(_defaultGASelection));
}



/* void Setting::printValue()
 *
 * Print the current value of a setting
 */

void Setting::printValue()
{
	
	// Print bool
	if (_valueType == VT_BOOL)
	{
		if (_valueBool == true)
			Output::print("true");
		else
			Output::print("false");
	}
	
	// Print integer
	else if (_valueType == VT_INT)
		Output::print(_valueInt);
	
	// Print double
	else if (_valueType == VT_DOUBLE)
	{
		if ((fabs(_valueDouble) < 1e-2) || (fabs(_valueDouble) > 1e3))
			Output::printSci(_valueDouble);
		else
			Output::print(_valueDouble);
	}
	
	// Print coordinate type
	else if (_valueType == VT_COORDINATES)
	{
		if (_valueCoordinates == FRACTIONAL)
			Output::print("Fractional");
		else
			Output::print("Cartesian");
	}
	
	// Print structure format
	else if (_valueType == VT_STRFORMAT)
		Output::print(StructureIO::structureFormat(_valueStrFormat));
	
	// Print GA structure prediction metric
	else if (_valueType == VT_GAOPTMETRIC)
		Output::print(GAPredict::metric(_valueGAPredictMetric));
	
	// Print GA selection method
	else if (_valueType == VT_GASELECTION)
		Output::print(GeneticAlgorithm<int, Setting>::selection(_valueGASelection));
}



/* void Setting::error(const OList<Word>& input)
 *
 * Error message when reading a setting
 */

void Setting::error(const OList<Word>& input)
{
	Output::newline(ERROR);
	Output::print("Could not recognize setting in \"");
	Output::print(input, false);
	Output::print("\"");
	Output::quit();
}



// Static member values of Settings
Word Settings::_globalFile;
const int Settings::_numSettings = 34;
Setting Settings::_settings[Settings::_numSettings];



// Settings list
Settings::Helper Settings::_helper;
Settings::Helper::Helper()
{
	
	// Set the global file name
	_globalFile = Directory::makePath(Directory::home(), ".mint_settings");
	
	// NUMPROCS
	Settings::_settings[(int)NUMPROCS].setup(1, "numprocs");
	
	// OUTPUT_LEVEL
	Settings::_settings[(int)OUTPUT_LEVEL].setup(0, "displaylevel");
		
	// OUTPUT_TAB
	Settings::_settings[(int)OUTPUT_TAB].setup(4, "displaytab");
	
	// TIME_SHOW
	Settings::_settings[(int)TIME_SHOW].setup(false, "timeshow");
	
	// TIME_PRECISION
	Settings::_settings[(int)TIME_PRECISION].setup(4, "timeprec");
	
	// TIME_FORMAT
	Settings::_settings[(int)TIME_FORMAT].setup(true, "formattime");
	
	// TOLERANCE
	Settings::_settings[(int)TOLERANCE].setup(1.0e-4, "tolerance");
	
	// CLUSTERTOL
	Settings::_settings[(int)CLUSTERTOL].setup(0.2, "clustertol");
	
	// USE_STDOUT
	Settings::_settings[(int)USE_STDOUT].setup(false, "usestdout");
	
	// COORDINATES
	Settings::_settings[(int)COORDINATES].setup(FRACTIONAL, "coordinates");
	
	// STRUCTURE_FORMAT
	Settings::_settings[(int)STRUCTURE_FORMAT].setup(SF_AUTO, "strformat");
	
	// OVERWRITE
	Settings::_settings[(int)OVERWRITE].setup(false, "overwrite");
	
	// ADD_EXTENSION
	Settings::_settings[(int)ADD_EXTENSION].setup(true, "addextension");
	
	// RANDSTR_MAXLOOPS
	Settings::_settings[(int)RANDSTR_MAXLOOPS].setup(100, "randstrmaxloops");
	
	// RANDSTR_MINBOND
	Settings::_settings[(int)RANDSTR_MINBOND].setup(0.5, "randstrminbond");
	
	// GAOPT_NUMSIM
	Settings::_settings[(int)GAOPT_NUMSIM].setup(1, "gaoptnumsim");
	
	// GAOPT_POPSIZE
	Settings::_settings[(int)GAOPT_POPSIZE].setup(10, "gaoptpopsize");
	
	// GAOPT_CELLMUTPROB
	Settings::_settings[(int)GAOPT_CELLMUTPROB].setup(0.1, "gaoptcellmutprob");
	
	// GAOPT_POSMUTPROB
	Settings::_settings[(int)GAOPT_POSMUTPROB].setup(0.1, "gaoptposmutprob");
	
	// GAOPT_WYCKMUTPROB
	Settings::_settings[(int)GAOPT_WYCKMUTPROB].setup(0.1, "gaoptwyckmutprob");
	
	// GAOPT_METRICTOOPT
	Settings::_settings[(int)GAOPT_METRICTOOPT].setup(GAPM_UNKNOWN, "gaoptmetric");
	
	// GAOPT_CONVERGEOVER
	Settings::_settings[(int)GAOPT_CONVERGEOVER].setup(10, "gaoptconverge");
	
	// GAOPT_MAXGENS
	Settings::_settings[(int)GAOPT_MAXGENS].setup(1000, "gaoptmaxgens");
	
	// GAOPT_NUMTOKEEP
	Settings::_settings[(int)GAOPT_NUMTOKEEP].setup(0, "gaoptnumtokeep");
	
	// GAOPT_SELECTION
	Settings::_settings[(int)GAOPT_SELECTION].setup(GA_TOURNAMENT, "gaoptselection");
	
	// GAOPT_ENERGYTOL
	Settings::_settings[(int)GAOPT_ENERGYTOL].setup(1e-3, "gaoptenergytol");
	
	// GAOPT_DIFFRACTIONTOL
	Settings::_settings[(int)GAOPT_DIFFRACTIONTOL].setup(1e-4, "gaoptdifftol");
	
	// GAOPT_SCREENMETHOD
	Settings::_settings[(int)GAOPT_SCREENMETHOD].setup(GAPM_UNKNOWN, "gaoptscreenmethod");
	
	// GAOPT_SCREENNUM
	Settings::_settings[(int)GAOPT_SCREENNUM].setup(0, "gaoptscreennum");
	
	// WYCKOFFBIAS
	Settings::_settings[(int)WYCKOFFBIAS].setup(0.5, "wyckoffbias");
	
	// MINIMAGEDISTANCE
	Settings::_settings[(int)MINIMAGEDISTANCE].setup(8.0, "minimagedis");
	
	// MAXJUMPDISTANCE
	Settings::_settings[(int)MAXJUMPDISTANCE].setup(7.0, "maxjumpdistance");
	
	// KMC_MINJUMPS
	Settings::_settings[(int)KMC_JUMPSPERATOM].setup(100.0, "kmcjumpsperatom");
	
	// KMC_CONVERGENCE
	Settings::_settings[(int)KMC_CONVERGENCE].setup(0.01, "kmcconverge");
}



/* void Settings::set(const Text& content)
 *
 * Set settings from file contents
 */

void Settings::set(const Text& content)
{
	
	// Loop over file contents
	int i, j;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Loop over settings and save
		for (j = 0; j < _numSettings; ++j)
		{
			if (_settings[j].set(content[i]))
				break;
		}
	}
}



/* bool Settings::isFormat(const Text& content, double target)
 *
 * Return whether contents are a settings file
 */

bool Settings::isFormat(const Text& content, double target)
{
	
	// Loop over lines in file
	int i, j;
	double tagCount = 0;
	double lineCount = 0;
	for (i = 0; i < content.length(); ++i)
	{
		
		// Current line is empty
		if (!content[i].length())
			continue;
		
		// Current line starts with a comment
		if (Language::isComment(content[i][0]))
			continue;
		
		// Count current line and check if a tag
		++lineCount;
		for (j = 0; j < _numSettings; ++j)
		{
			if (_settings[j].isTag(content[i][0]))
			{
				++tagCount;
				break;
			}
		}
	}
	
	// Return false if file does not contain any content
	if (!lineCount)
		return false;
	
	// Return false if number of tags vs. number of lines is less than target
	if (tagCount / lineCount < target)
		return false;
	
	// Return that a settings file
	return true;
}
