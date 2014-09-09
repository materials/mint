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



#include "multi.h"
#include "output.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>



/* StreamInfo::StreamInfo()
 *
 * Constructor for StreamInfo object
 */

StreamInfo::StreamInfo()
{
	_streamIsSet = false;
	_stream = 0;
	_havePrinted = 0;
	_numLinesBefore = 0;
	_numLinesAfter = 0;
	_hold = false;
	_method = STANDARD;
	_type = ORDINARY;
	_showTermination = true;
	_ordinaryOn = true;
	_warningsOn = true;
	_errorsOn = true;
	_maxLevel = 10000;
	_spacesPerLevel = 4;
	_currentLevel = 0;
	_print = true;
}



/* StreamInfo::~StreamInfo()
 * 
 * Destructor for StreamInfo object
 */

StreamInfo::~StreamInfo()
{
	if (!_streamIsSet)
		return;
	if (_havePrinted)
	{
		for (int i = 0; i < _numLinesAfter; ++i)
			*_stream << endl;
	}
	if (_outfile.is_open())
		_outfile.close();
}



/* void StreamInfo::setup(ostream* input, int numBlankLinesAtStart, int numBlankLinesAtEnd)
 *
 * Constructor for StreamInfo object
 */

void StreamInfo::setup(ostream* input, int numBlankLinesAtStart, int numBlankLinesAtEnd)
{
	if (!Multi::rank())
	{
		_stream = input;
		_streamIsSet = true;
		_numLinesBefore = numBlankLinesAtStart;
		_numLinesAfter = numBlankLinesAtEnd;
	}
}



/* void StreamInfo::setup(const Word& file, int numBlankLinesAtStart, int numBlankLinesAtEnd)
 *
 * Constructor for StreamInfo object when writing to file
 */

void StreamInfo::setup(const Word& file, int numBlankLinesAtStart, int numBlankLinesAtEnd)
{
	
	// Only run on root
	if (Multi::rank() == 0)
	{
		
		// Open file and save properties
		_outfile.open(file.array());
		if (!_outfile.is_open())
		{
			Output::newline(ERROR);
			Output::print("Could not open file for writing: ");
			Output::print(file);
			Output::quit();
		}
		_outfileName = file;
		_stream = &_outfile;
		_streamIsSet = true;
		_numLinesBefore = numBlankLinesAtStart;
		_numLinesAfter = numBlankLinesAtEnd;
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void StreamInfo::quietOn(bool warningsOff, bool errorsOff)
 *
 * Turn on quiet
 */

void StreamInfo::quietOn(bool warningsOff, bool errorsOff)
{
	
	// Save current states
	_prevOrdinaryOn += _ordinaryOn;
	_prevWarningsOn += _warningsOn;
	_prevErrorsOn += _errorsOn;
	
	// Turn off messages
	_ordinaryOn = false;
	_warningsOn = (!warningsOff);
	_errorsOn = (!errorsOff);
}



/* void StreamInfo::quietOff()
 *
 * Turn quiet off
 */

void StreamInfo::quietOff()
{
	
	// Quiet is already off
	if (!_prevOrdinaryOn.length())
		return;
	
	// Save previous states
	_ordinaryOn = _prevOrdinaryOn.last();
	_warningsOn = _prevWarningsOn.last();
	_errorsOn = _prevErrorsOn.last();
	
	// Remove previous states
	_prevOrdinaryOn.remove(_prevOrdinaryOn.length() - 1);
	_prevWarningsOn.remove(_prevWarningsOn.length() - 1);
	_prevErrorsOn.remove(_prevErrorsOn.length() - 1);
}



/* void StreamInfo::updatePrint()
 *
 * Decide whether output should be printed based on current settings
 */

void StreamInfo::updatePrint()
{
	
	// Stream is not on
	if (!_streamIsSet)
	{
		_print = false;
		return;
	}
	
	// Current type is an error
	if (_type == ERROR)
	{
		
		// Error printing is on
		if (_errorsOn)
			_print = true;
		
		// Error printing is off
		else
			_print = false;
	}
	
	// Current type is a warning
	else if (_type == WARNING)
	{
		
		// Warning printing is on
		if (_warningsOn)
			_print = true;
		
		// Warning printing is off
		else
			_print = false;
	}
	
	// Current type is ordinary
	else
	{
		
		// Ordinary printing is off
		if (!_ordinaryOn)
			_print = false;
		
		// Ordinary printing is on
		else
		{
			
			// Printing in standard mode
			if (_method == STANDARD)
				_print = true;
			
			// Printing in restricted mode
			else
			{
				
				// Inside of range
				if (_currentLevel < _maxLevel)
					_print = true;
				
				// Outside of range
				else
					_print = false;
			}
		}
	}
}



/* void Output::initialize()
 *
 * Start output framework
 */

void Output::initialize()
{
	Output::_ids += 1;
	Output::_streams.add();
	Output::_streams.last().setup(&cout);
}



// Initialize static variables
List<int> Output::_ids;
OList<StreamInfo> Output::_streams;
int Output::_primary = 0;
int Output::_curStream = 0;
char Output::_arg[10];
char Output::_buffer[50];



/* double Output::checkZero(double in, int prec)
 *
 * Return zero without number or original number
 */

double Output::checkZero(double in, int prec)
{
	double comp = in;
	comp *= pow(10, (double)(prec+1));
	if ((comp < 1) && (comp > -1))
		return 0;
	return in;
}



/* int Output::roundToOne(double in)
 *
 * Round value to nearest integer
 */

int Output::roundToOne(double in)
{
	if (in > 0)
    {
        if (in - floor(in) >= 0.5)
            return (int) ceil(in);
        else
            return (int) floor(in);
    }
    else
    {
        if (in - ceil(in) <= -0.5)
            return (int) floor(in);
        else
            return (int) ceil(in);
    }
}



/* void Output::add(int input)
 *
 * Add integer to output object
 */

void Output::add(int input)
{
	sprintf(_buffer, "%d", input);
	add(_buffer);
}



/* void Output::add(double input, int prec)
 * 
 * Add double to output object
 */

void Output::add(double input, int prec)
{
	sprintf(_arg, "%s%d%s", "%.", prec, "f");
	sprintf(_buffer, _arg, checkZero(input, prec));
	add(_buffer);
}



/* void Output::addSci(double input, int prec)
 * 
 * Add double to output object in scientific notation
 */

void Output::addSci(double input, int prec)
{
	sprintf(_arg, "%s%d%s", "%.", prec, "e");
	sprintf(_buffer, _arg, checkZero(input, abs((int)floor(log10(fabs(input)))) + prec));
	add(_buffer);
}



/* void Output::addTab()
 *
 * Add tab to output object
 */

void Output::addTab()
{
	for (int i = 0; i < _streams[_curStream].spacesPerLevel() - 1; ++i)
		_buffer[i] = ' ';
	if (_addSpaces)
		_buffer[_streams[_curStream].spacesPerLevel() - 1] = '\0';
	else
	{
		_buffer[_streams[_curStream].spacesPerLevel() - 1] = ' ';
		_buffer[_streams[_curStream].spacesPerLevel()] = '\0';
	}
	add(_buffer);
}



/** 
 * Add a new stream to available output streams
 * 
 * @param input Pointer to an output stream
 * @param numBlankLinesAtStart Number of blank lines to insert at head of file
 * @param numBlankLinesAtEnd Number of blank lines to insert when closing file
 * @return ID number of new stream
 */
int Output::addStream(ostream* input, int numBlankLinesAtStart, int numBlankLinesAtEnd)
{

	// Save stream
	_ids += _ids.last() + 1;
	_streams.add();
	_streams.last().setup(input, numBlankLinesAtStart, numBlankLinesAtEnd);
	
	// Return new stream id
	return _ids.last();
}


/** 
 * Add a new file to available output streams
 * 
 * @param input Path to file 
 * @param numBlankLinesAtStart Number of blank lines to insert at head of file
 * @param numBlankLinesAtEnd Number of blank lines to insert when closing file
 * @return ID number of new stream
 */
int Output::addStream(const Word& input, int numBlankLinesAtStart, int numBlankLinesAtEnd)
{

	// Save stream
	_ids += _ids.last() + 1;
	_streams.add();
	_streams.last().setup(input, numBlankLinesAtStart, numBlankLinesAtEnd);
	
	// Return new stream id
	return _ids.last();
}



/* void Output::removeStream(int inID)
 *
 * Delete a stream
 */

void Output::removeStream(int inID)
{
	int num = 0;
	for (; num < _ids.length(); ++num)
	{
		if (_ids[num] == inID)
			break;
	}
	if (num < _ids.length())
	{
		_ids.remove(num);
		_streams.remove(num);
		if (_curStream >= num)
			_curStream--;
		if (_primary >= num)
			_primary--;
	}
}



/* void Output::setPrimary(int inID)
 *
 * Set the primary stream
 */

void Output::setPrimary(int inID)
{
	for (int i = 0; i < _ids.length(); ++i)
	{
		if (_ids[i] == inID)
		{
			_primary = i;
			return;
		}
	}
}



/* void Output::setStream(int inID)
 *
 * Set the current stream
 */

void Output::setStream(int inID)
{
	for (int i = 0; i < _ids.length(); ++i)
	{
		if (_ids[i] == inID)
		{
			_curStream = i;
			return;
		}
	}
}



/* List<bool> Output::onOff()
 *
 * Return a list with whether ordinary messages, warnings, and errors will be printed
 */

List<bool> Output::onOff()
{
	List<bool> res(3);
	res[0] = _streams[_curStream].ordinaryOn();
	res[1] = _streams[_curStream].warningsOn();
	res[2] = _streams[_curStream].errorsOn();
	return res;
}



/* void Output::method(PrintMethod input)
 *
 * Set the print method
 */

void Output::method(PrintMethod input)
{
	if (((_streams[_curStream].havePrinted()) && (input == STANDARD)) || \
		((_streams[_curStream].havePrinted() == 2) && (_streams[_curStream].method() == STANDARD)))
		_streams[_curStream].hold(true);
	if (_streams[_curStream].havePrinted())
		_streams[_curStream].havePrinted(1);
	_streams[_curStream].method(input);
}



/* void Output::newline(PrintType input)
 *
 * Print a new line
 */

void Output::newline(PrintType input)
{
	
	// Set current type
	_streams[_curStream].type(input);
	
	// Check if printing
	if (!_streams[_curStream].print())
		return;
	
	// This is the first print
	int i;
	if (!_streams[_curStream].havePrinted())
	{
		for (i = 0; i < _streams[_curStream].numLinesBefore() - 1; ++i)
			stream() << endl;
	}
	
	// Print blank line if holding
	if (_streams[_curStream].hold())
	{
		stream() << endl;
		_streams[_curStream].hold(false);
	}
	
	// Print new line in standard setting
	if (_streams[_curStream].method() == STANDARD)
	{
		if ((_streams[_curStream].havePrinted()) || (_streams[_curStream].numLinesBefore()))
			stream() << endl;
		_streams[_curStream].havePrinted(2);
		return;
	}
	
	// Print new line in restricted setting
	if ((_streams[_curStream].havePrinted()) || (_streams[_curStream].numLinesBefore()))
		stream() << endl;
	
	// Print ordinary line
	if (_streams[_curStream].type() == ORDINARY)
	{
		stream() << "O: ";
		for (i = 0; i < _streams[_curStream].currentLevel() * _streams[_curStream].spacesPerLevel(); ++i)
			stream() << " ";
	}
	
	// Print warning line or error line
	else
	{
		
		// Print warning line
		if (_streams[_curStream].type() == WARNING)
			stream() << "W: ";
		
		// Print error line
		else
			stream() << "E: ";
		
		// Print spaces
		if (_streams[_curStream].currentLevel() < _streams[_curStream].maxLevel())
		{
			for (i = 0; i < _streams[_curStream].currentLevel() * _streams[_curStream].spacesPerLevel(); ++i)
				stream() << " ";
		}
		else
		{
			for (i = 0; i < _streams[_curStream].maxLevel() * _streams[_curStream].spacesPerLevel(); ++i)
				stream() << " ";
		}
	}
	
	// Print message if needed
	if (_streams[_curStream].type() == WARNING)
		stream() << "WARNING: ";
	else if (_streams[_curStream].type() == ERROR)
		stream() << "ERROR: ";
	stream() << flush;
	
	// Save that print was made
	_streams[_curStream].havePrinted(2);
}



/* void Output::tab()
 * 
 * Print a tab
 */

void Output::tab()
{
	
	// Check if printing
	if (!_streams[_curStream].print())
		return;
	
	// Print tab
	for (int i = 0; i < _streams[_curStream].spacesPerLevel(); ++i)
		stream() << " ";
	stream() << flush;
}



/* void Output::print(char message)
 *
 * Print a single character
 */

void Output::print(char message)
{
	if (!_streams[_curStream].print())
		return;
	stream() << message << flush;
}



/* void Output::print(const char* message)
 *
 * Print a char message
 */

void Output::print(const char* message)
{
	if (!_streams[_curStream].print())
		return;
	stream() << message << flush;
}



/* void Output::print(int message)
 *
 * Print an integer message
 */

void Output::print(int message)
{
	if (!_streams[_curStream].print())
		return;
	stream() << message << flush;
}



/* void Output::print(unsigned long int message)
 *
 * Print an unsigned long integer message
 */

void Output::print(unsigned long int message)
{
	if (!_streams[_curStream].print())
		return;
	stream() << message << flush;
}



/* void Output::print(double message, int prec)
 *
 * Print a double message
 */

void Output::print(double message, int prec)
{
	if (!_streams[_curStream].print())
		return;
	int places = prec;
	if (places == -1)
	{
		double temp = message;
		for (places = 0; places <= 8; ++places)
		{
			if ((temp - roundToOne(temp) > -1e-9) && (temp - roundToOne(temp) < 1e-9))
				break;
			temp *= 10;
		}
	}
	streamsize origPrec = stream().precision();
	stream() << setiosflags(ios::fixed) << setprecision(places) << checkZero(message, places);
	stream() << resetiosflags(ios::fixed) << setprecision(origPrec) << flush;
}



/* void Output::printSci(double message, int prec)
 *
 * Print double in scientific notation
 */

void Output::printSci(double message, int prec)
{
	if (!_streams[_curStream].print())
		return;
	int order = (int)floor(log10(fabs(message)));
	int places = prec;
	if (places == -1)
	{
		double temp = message / pow(10.0, order);
		for (places = 0; places <= 8; places++)
		{
			if ((temp - roundToOne(temp) > -1e-9) && (temp - roundToOne(temp) < 1e-9))
				break;
			temp *= 10;
		}
	}
	streamsize origPrec = stream().precision();
	stream() << setiosflags(ios::scientific) << setprecision(places);
	stream() << checkZero(message, abs(order) + places);
	stream() << resetiosflags(ios::scientific) << setprecision(origPrec) << flush;
}



/* void Output::print(const Words& message, bool useComma, bool useAnd)
 *
 * Print a list of char arrays
 */

void Output::print(const Words& message, bool useComma, bool useAnd)
{
	if (!_streams[_curStream].print())
		return;
	for (int i = 0; i < message.length(); ++i)
	{
		stream() << message[i];
		if (useComma)
		{
			if ((message.length() > 2) && (i != message.length() - 1))
				stream() << ",";
			if ((i == message.length() - 2) && (useAnd))
				stream() << " and";
		}
		if (i != message.length() - 1)
			stream() << " ";
	}
	stream() << flush;
}



/* void Output::print(const List<int>& message, bool useComma, bool useAnd)
 *
 * Print a list of integers
 */

void Output::print(const List<int>& message, bool useComma, bool useAnd)
{
	if (!_streams[_curStream].print())
		return;
	for (int i = 0; i < message.length(); ++i)
	{
		stream() << message[i];
		if (useComma)
		{
			if ((message.length() > 2) && (i != message.length() - 1))
				stream() << ",";
			if ((i == message.length() - 2) && (useAnd))
				stream() << " and";
		}
		if (i != message.length() - 1)
			stream() << " ";
	}
	stream() << flush;
}



/* void Output::print(const List<double>& message, int prec, bool useComma, bool useAnd)
 *
 * Print a list of doubles
 */

void Output::print(const List<double>& message, int prec, bool useComma, bool useAnd)
{
	if (!_streams[_curStream].print())
		return;
	streamsize origPrec = stream().precision();
	for (int i = 0; i < message.length(); ++i)
	{
		stream() << setiosflags(ios::fixed) << setprecision(prec) << checkZero(message[i], prec);
		stream() << resetiosflags(ios::fixed) << setprecision(origPrec);
		if (useComma)
		{
			if ((message.length() > 2) && (i != message.length() - 1))
				stream() << ",";
			if ((i == message.length() - 2) && (useAnd))
				stream() << " and";
		}
		if (i != message.length() - 1)
			stream() << " ";
	}
	stream() << flush;
}



/* void Output::print(const Vector3D& message, int prec, bool useComma)
 *
 * Print vector
 */

void Output::print(const Vector3D& message, int prec, bool useComma)
{
	if (!_streams[_curStream].print())
		return;
	streamsize origPrec = stream().precision();
	for (int i = 0; i < 3; ++i)
	{
		stream() << setiosflags(ios::fixed) << setprecision(prec) << checkZero(message[i], prec);
		stream() << resetiosflags(ios::fixed) << setprecision(origPrec);
		if (i != 2)
		{
			if (useComma)
				stream() << ",";
			stream() << " ";
		}
	}
	stream() << flush;
}



/* void Output::print(const Output& message, PrintAlign align, PrintType input)
 *
 * Print Output object
 */

void Output::print(const Output& message, PrintAlign align, PrintType input)
{
	
	// Set the current type
	_streams[_curStream].type(input);
	
	// Check if printing
	if (!_streams[_curStream].print())
		return;
	
	// Figure out the maximum number of columns
	int i;
	int numCols = 0;
	for (i = 0; i < message._data.length(); ++i)
	{
		if (message._data[i].length() > numCols)
			numCols = message._data[i].length();
	}
	
	// Make alignment vector
	List<PrintAlign> alignList (numCols);
	for (i = 0; i < numCols; ++i)
		alignList[i] = align;
	
	// Print
	print(message, alignList, input);
}



/* void Output::print(const Output& message, const List<PrintAlign>& align, PrintType input)
 *
 * Print Output object
 */

void Output::print(const Output& message, const List<PrintAlign>& align, PrintType input)
{
	
	// Set the current type
	_streams[_curStream].type(input);
	
	// Check if printing
	if (!_streams[_curStream].print())
		return;
	
	// Figure out the maximum number of columns
	int i;
	int numCols = 0;
	for (i = 0; i < message._data.length(); ++i)
	{
		if (message._data[i].length() > numCols)
			numCols = message._data[i].length();
	}
	
	// Set the width for each column
	int j;
	int width[numCols];
	for (i = 0; i < numCols; ++i)
		width[i] = 0;
	for (i = 0; i < message._data.length(); ++i)
	{
		for (j = 0; j < message._data[i].length(); ++j)
		{
			if (message._data[i][j].length() > width[j])
				width[j] = message._data[i][j].length();
		}
	}
	
	// Print all values
	for (i = 0; i < message._data.length(); ++i)
	{
		
		// New line
		if (i)
			newline(_streams[_curStream].type());
		
		// Print words
		for (j = 0; j < message._data[i].length(); ++j)
		{
			
			// Set left alignment if needed
			if (align[j] == LEFT)
				stream() << resetiosflags(ios::right) << setiosflags(ios::left);
			
			// Print value
			stream() << setw(width[j]) << message._data[i][j];
			if (j != message._data[i].length() - 1)
			{
				if (message._data[i][j].length() == 1)
				{
					if ((message._data[i][j][0] != '(') && (message._addSpaces))
						stream() << " ";
				}
				else if (message._addSpaces)
					stream() << " ";
			}
			
			// Unset left alignment if needed
			if (align[j] == LEFT)
				stream() << resetiosflags(ios::left) << setiosflags(ios::right);
		}
	}	
	
	// Flush stream
	stream() << flush;
}



/* void Output::printPadded(int message, int width, PrintAlign align)
 *
 * Print padded integer
 */

void Output::printPadded(int message, int width, PrintAlign align)
{
	if (!_streams[_curStream].print())
		return;
	if (align == LEFT)
		stream() << resetiosflags(ios::right) << setiosflags(ios::left);
	stream() << setw(width) << message << flush;
	if (align == LEFT)
		stream() << resetiosflags(ios::left) << setiosflags(ios::right);
}



/* void Output::printPadded(double message, int width, PrintAlign align, int prec)
 *
 * Print padded double
 */

void Output::printPadded(double message, int width, PrintAlign align, int prec)
{
	if (!_streams[_curStream].print())
		return;
	int places = prec;
	if (places == -1)
	{
		double temp = message;
		for (places = 0; places <= 8; places++)
		{
			if ((temp - roundToOne(temp) > -1e-9) && (temp - roundToOne(temp) < 1e-9))
				break;
			temp *= 10;
		}
	}
	streamsize origPrec = stream().precision();
	if (align == LEFT)
		stream() << resetiosflags(ios::right) << setiosflags(ios::left);
	stream() << setiosflags(ios::fixed) << setprecision(places) << setw(width) << checkZero(message, places);
	stream() << resetiosflags(ios::fixed) << setprecision(origPrec) << flush;
	if (align == LEFT)
		stream() << resetiosflags(ios::left) << setiosflags(ios::right);
}



/* void Output::printPaddedSci(double message, int width, PrintAlign align, int prec)
 *
 * Print padded double in scientific format
 */

void Output::printPaddedSci(double message, int width, PrintAlign align, int prec)
{
	if (!_streams[_curStream].print())
		return;
	int order = (int)floor(log10(fabs(message)));
	int places = prec;
	if (places == -1)
	{
		double temp = message / pow(10.0, order);
		for (places = 0; places <= 8; places++)
		{
			if ((temp - roundToOne(temp) > -1e-9) && (temp - roundToOne(temp) < 1e-9))
				break;
			temp *= 10;
		}
	}
	streamsize origPrec = stream().precision();
	if (align == LEFT)
		stream() << resetiosflags(ios::right) << setiosflags(ios::left);
	stream() << setiosflags(ios::scientific) << setprecision(places);
	stream() << setw(width) << checkZero(message, abs(order) + places);
	stream() << resetiosflags(ios::scientific) << setprecision(origPrec) << flush;
	if (align == LEFT)
		stream() << resetiosflags(ios::left) << setiosflags(ios::right);
}



/* void Output::quit()
 *
 * Exit program
 */

void Output::quit()
{
	for (int i = 0; i < _streams.length(); ++i)
	{
		if ((_streams[i].showTermination()) && (_streams[i].streamIsSet()) && (!Multi::rank()))
			(_streams[i].stream()) << endl << "Terminating before completion" << flush;
	}
	exit(1);
}
