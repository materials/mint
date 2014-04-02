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



#ifndef TIMER_H
#define TIMER_H



#include "text.h"



// Class for timer functions
class Timer
{
	
	// Variables
	double _start;
	
public:
	
	// Constructor
	Timer(bool startTimer = true) { if (startTimer) start(); }
	
	// Functions
	void start();
	
	// Access functions
	double currentNumber();
	Word current(bool formatTime, int precision);
	Word current(bool formatTime)					{ return current(formatTime, 4); }
	Word current(int precision)						{ return current(true, precision); }
	Word current()									{ return current(true, 4); }
};



#endif
