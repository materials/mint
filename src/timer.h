/* timer.h -- Timer functions
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
