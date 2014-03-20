/* Copyright 2011-2014 Kyle Michel, Logan Ward
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
#include "launcher.h"



// Main function
int main(int argc, char** argv)
{
	
	// Start MPI
	Multi::initialize(argc, argv);
	
	// Start output
	Output::initialize();
    
    // Call launcher to process arguments and call functions
    Launcher::start(argc, argv);

	// Close MPI
	Multi::finalize();
    
	// Finished
    return 0;
}
