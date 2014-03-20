/* help.h -- Print help information
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef MINTHELP_H
#define MINTHELP_H



#include "text.h"



// Class for help functions
class Help
{

	// Functions
	static void general();
	static void interface();
	static void functions();
	static void settings();

public:
	
	// Functions
	static void run(const Words& args);
};



#endif
