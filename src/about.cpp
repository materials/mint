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
	
	Output::newline();
	Output::newline();
	Output::print("Compiled on ");
	#ifdef COMP
		Output::print(COMP);
	#else
		Output::print("unknown date at unknown time");
	#endif
	Output::newline();
	Output::print("MPI support: ");
	#ifdef MINT_MPI
		Output::print("on");
	#else
		Output::print("off");
	#endif
	
	Output::newline();
	Output::newline();
	Output::newline(); Output::print("Copyright 2011-2014 Kyle Michel, Logan Ward, Christopher Wolverton, Christopher Wolverton");
	Output::newline();
	Output::newline(); Output::print("Mint is free software: you can redistribute it and/or modify it under the");
	Output::newline(); Output::print("terms of the GNU Lesser General Public License as published by the Free");
	Output::newline(); Output::print("Software Foundation, either version 3 of the License, or (at your option)");
	Output::newline(); Output::print("any later version.");
	Output::newline();
	Output::newline(); Output::print("Mint is distributed in the hope that it will be useful, but WITHOUT ANY");
	Output::newline(); Output::print("WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS");
	Output::newline(); Output::print("FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for");
	Output::newline(); Output::print("more details.");
	Output::newline();
	Output::newline(); Output::print("You should have received a copy of the GNU Lesser General Public License");
	Output::newline(); Output::print("along with Mint.  If not, see <http://www.gnu.org/licenses/>.");
	Output::newline();
	Output::newline(); Output::print("Contact: Kyle Michel (kylemichel@gmail.com)");
	Output::newline(); Output::print("         Logan Ward (LoganWard2012@u.northwestern.edu)");
		
	Output::method(origMethod);
 }

