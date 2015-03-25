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

/**
 * =========================
 *  Developer's guide 
 * =========================
 * 
 * --> How to add a new setting
 * 
 * 1. Add a new SettingsLabel
 * 2. Increase number of slots in which to store settings
 *		a. Look for "const int Settings::_numSettings =" in settings.cpp
 * 3. Assign a name for this setting and a default value in the constructor
 *    of Settings::Help (look in settings.cpp)
 * 4. Add a description of the new setting into README.md
 * 
 * --> How to use settings
 * 
 * Easy enough, just call "Settings::value([SettingsLabel], [variable in which to store value])"
 *   or "Settings::value&lt;[type of setting (i.e. double)]&gt;([SettingLabel])"
 */


#ifndef SETTINGS_H
#define SETTINGS_H



#include "iso.h"
#include "structureIO.h"
#include "ga.h"
#include "gaPredict.h"
#include "text.h"
#include "list.h"
#include <cmath>



// Types of settings
enum SettingsLabel {NUMPROCS, OUTPUT_LEVEL, OUTPUT_TAB, TIME_SHOW, TIME_PRECISION, TIME_FORMAT, TOLERANCE, CLUSTERTOL,\
	USE_STDOUT, COORDINATES, STRUCTURE_FORMAT, OVERWRITE, ADD_EXTENSION, RANDSTR_MAXLOOPS, RANDSTR_MINBOND, \
	GAOPT_NUMSIM, GAOPT_POPSIZE, GAOPT_CELLMUTPROB, GAOPT_POSMUTPROB, GAOPT_WYCKMUTPROB, GAOPT_METRICTOOPT, \
	GAOPT_CONVERGEOVER, GAOPT_MAXGENS, GAOPT_NUMTOKEEP, GAOPT_SELECTION, GAOPT_ENERGYTOL, GAOPT_DIFFRACTIONTOL, \
	GAOPT_USERIETVELD, GAOPT_SCREENMETHOD, GAOPT_SCREENNUM, GAOPT_ALLOWRESTART, GAOPT_SAVEALLRESULTS, \
	WYCKOFFBIAS, MINIMAGEDISTANCE, MAXJUMPDISTANCE, KMC_JUMPSPERATOM, KMC_CONVERGENCE, \
	XRD_BACKGROUNDCOUNT, XRD_LATPARAM};



// Class to store a single setting
class Setting
{
	
public:
	
	// Possible types
	enum ValueType {VT_BOOL, VT_INT, VT_DOUBLE, VT_COORDINATES, VT_STRFORMAT, VT_GAOPTMETRIC, VT_GASELECTION};
	
private:

	// General variables
	ValueType _valueType;
	Word _tag;
	
	// Store default values
	bool _defaultBool;
	int _defaultInt;
	double _defaultDouble;
	CoordinateType _defaultCoordinates;
	StructureFormat _defaultStrFormat;
	GAPredictMetric _defaultGAPredictMetric;
	GASelectionMethod _defaultGASelection;
	
	// Store current values
	bool _valueBool;
	int _valueInt;
	double _valueDouble;
	CoordinateType _valueCoordinates;
	StructureFormat _valueStrFormat;
	GAPredictMetric _valueGAPredictMetric;
	GASelectionMethod _valueGASelection;
	
	// Functions
	void commonSetup(const Word& tag);
	void error(const OList<Word>& input);
	
public:
	
	// Setup functions
	void setup(             bool defValue, const Word& tag);
	void setup(              int defValue, const Word& tag);
	void setup(           double defValue, const Word& tag);
	void setup(   CoordinateType defValue, const Word& tag);
	void setup(  StructureFormat defValue, const Word& tag);
	void setup(  GAPredictMetric defValue, const Word& tag);
	void setup(GASelectionMethod defValue, const Word& tag);
	
	// Print functions
	void printDefault();
	void printValue();
	
	// Set value type of setting
	void valueType(ValueType input)	{ _valueType = input; }
	
	// Get value from input
	bool set(const OList<Word>& input);
	bool isTag(const Word& input) const	{ return _tag.equal(input, false, _tag.length(), '='); }
	
	// Set value of setting
	void value(bool input)				{ _valueBool            = input; }
	void value(int input)				{ _valueInt             = input; }
	void value(double input)			{ _valueDouble          = input; }
	void value(CoordinateType input)	{ _valueCoordinates     = input; }
	void value(StructureFormat input)	{ _valueStrFormat       = input; }
	void value(GAPredictMetric input)	{ _valueGAPredictMetric = input; }
	void value(GASelectionMethod input)	{ _valueGASelection     = input; }
	
	// Access functions
	bool valueBool() const							{ return _valueBool; }
	int valueInt() const							{ return _valueInt; }
	double valueDouble() const						{ return _valueDouble; }
	CoordinateType valueCoordinates() const			{ return _valueCoordinates; }
	StructureFormat valueStrFormat() const			{ return _valueStrFormat; }
	GAPredictMetric valueGAPredictMetric() const	{ return _valueGAPredictMetric; }
	GASelectionMethod valueGASelection() const		{ return _valueGASelection; }
};



// Class to store all settings
class Settings
{
	
	// Settings objects
	static const int _numSettings;
	static Setting * _settings;
	
	// Local of global settings file
	static Word _globalFile;
	
	// Initialize variable
	class Helper { public: Helper(); };
	static Helper _helper; 

public:
	
	// Setup functions
	static void set(const Text& content);
	template <class T> static void value(SettingsLabel setting, const T& value)
		{ _settings[(int)setting].value(value); }
	
	// Access functions
	template <class T> static T value(SettingsLabel setting);
	static const Word globalFile()	{ return _globalFile; }
	
	// Helper functions
	static bool isFormat(const Text& content, double target = 0.5);
};



// =====================================================================================================================
// Setting
// =====================================================================================================================

/* inline void Setting::setup(bool defValue, const Word& tag)
 *
 * Setup Setting object
 */

inline void Setting::setup(bool defValue, const Word& tag)
{
	_valueType = VT_BOOL;
	_defaultBool = _valueBool = defValue;
	commonSetup(tag);
}

inline void Setting::setup(int defValue, const Word& tag)
{
	_valueType = VT_INT;
	_defaultInt = _valueInt = defValue;
	commonSetup(tag);
}

inline void Setting::setup(double defValue, const Word& tag)
{
	_valueType = VT_DOUBLE;
	_defaultDouble = _valueDouble = defValue;
	commonSetup(tag);
}

inline void Setting::setup(CoordinateType defValue, const Word& tag)
{
	_valueType = VT_COORDINATES;
	_defaultCoordinates = _valueCoordinates = defValue;
	commonSetup(tag);
}

inline void Setting::setup(StructureFormat defValue, const Word& tag)
{
	_valueType = VT_STRFORMAT;
	_defaultStrFormat = _valueStrFormat = defValue;
	commonSetup(tag);
}

inline void Setting::setup(GAPredictMetric defValue, const Word& tag)
{
	_valueType = VT_GAOPTMETRIC;
	_defaultGAPredictMetric = _valueGAPredictMetric = defValue;
	commonSetup(tag);
}

inline void Setting::setup(GASelectionMethod defValue, const Word& tag)
{
	_valueType = VT_GASELECTION;
	_defaultGASelection = _valueGASelection = defValue;
	commonSetup(tag);
}



/* inline void Setting::commonSetup(const Word& tag)
 *
 * Finish setup of Settings object
 */

inline void Setting::commonSetup(const Word& tag)
{
	_tag = tag;
}



// =====================================================================================================================
// Settings
// =====================================================================================================================

/* inline T Settings::value<bool>(SettingsLabel setting)
 *
 * Return the value of a setting
 */

template <>
inline bool Settings::value<bool>(SettingsLabel setting)
{ return _settings[(int)setting].valueBool(); }

template <>
inline int Settings::value<int>(SettingsLabel setting)
{ return _settings[(int)setting].valueInt(); }

template <>
inline double Settings::value<double>(SettingsLabel setting)
{ return _settings[(int)setting].valueDouble(); }

template <>
inline CoordinateType Settings::value<CoordinateType>(SettingsLabel setting)
{ return _settings[(int)setting].valueCoordinates(); }

template <>
inline StructureFormat Settings::value<StructureFormat>(SettingsLabel setting)
{ return _settings[(int)setting].valueStrFormat(); }

template <>
inline GAPredictMetric Settings::value<GAPredictMetric>(SettingsLabel setting)
{ return _settings[(int)setting].valueGAPredictMetric(); }

template <>
inline GASelectionMethod Settings::value<GASelectionMethod>(SettingsLabel setting)
{ return _settings[(int)setting].valueGASelection(); }



#endif
