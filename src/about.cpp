/* about.cpp -- About message
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "about.h"
#include "output.h"



/* void About::print()
 *
 * Print information about the program
 */

void About::print()
{
	PrintMethod origMethod = Output::method();
	Output::method(STANDARD);
	
	Output::newline();
	Output::print("Materials Interface (mint) ");
	#if defined(VERSION)
		Output::print("v");
		Output::print(VERSION);
	#else
		Output::print("unknown version");
	#endif
	
	Output::newline(); Output::print("\nCopyright (C) 2011-2014 by Northwestern University, All Rights Reserved");
	Output::newline(); Output::print("Contact: Kyle Michel (kylemichel@gmail.com)");
     
	Output::newline();
	Output::newline();
	Output::print("Compiled on ");
	#ifdef COMP
		Output::print(COMP);
	#else
		Output::print("unknown date at unknown time");
	#endif
	Output::newline();
	Output::print("    MPI: ");
	#ifdef MINT_MPI
		Output::print("on");
	#else
		Output::print("off");
	#endif
		
	Output::method(origMethod);
 }
