/* mint.cpp -- Main function for mint program
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
