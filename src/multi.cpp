/* multi.cpp -- MPI functions
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifdef MINT_MPI
	#include <mpi.h>
#endif
#include "multi.h"
#include "language.h"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
using namespace std;



// Initialize static variables
int Multi::_rank = 0;
int Multi::_worldSize = 1;
#ifdef MINT_MPI
	bool Multi::_mpiOn = true;
#else
	bool Multi::_mpiOn = false;
#endif
bool Multi::_setExitFun = false;
bool Multi::_runningFun = false;
int Multi::_jobSize = 1;



/* void Multi::initialize(int argc, char** argv)
 *
 * Initialize MPI functionality
 */

void Multi::initialize(int argc, char** argv)
{
	#ifdef MINT_MPI
	
		// Start MPI
    	MPI_Init(&argc, &argv);
	
		// Get the rank and number of processors
    	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &_worldSize);
		
	#endif
}



/* void Multi::finalize()
 *
 * End MPI functionality
 */

void Multi::finalize()
{
	#ifdef MINT_MPI
	
		// Close MPI
    	MPI_Finalize();
		
	#endif
}



/* bool Multi::safeCall(void (&function)(), const char* stdoutFile, const char* stderrFile)
 *
 * Call a function and avoid premature exit
 */

bool Multi::safeCall(void (&function)(), const char* stdoutFile, const char* stderrFile)
{
	
	// Set atexit function if needed
	if (!_setExitFun)
	{
		atexit(exitfun);
		_setExitFun = true;
	}
	
	// Redirect stdout if needed
	int out;
	if (stdoutFile)
	{
		out = dup(fileno(stdout));
		freopen(stdoutFile, "w", stdout);
	}
	
	// Redirect stderr if needed
	int err;
	if (stderrFile)
	{
		err = dup(fileno(stderr));
		freopen(stderrFile, "w", stderr);
	}
	
	// Call function
	bool ranWithoutError = true;
	_runningFun = true;
	try { function(); }
	catch (...)	{ ranWithoutError = false; }
	_runningFun = false;
	
	// Reset stdout if needed
	if (stdoutFile)
	{
		fflush(stdout);
		dup2(out, fileno(stdout));
		close(out);
	}
	
	// Reset stderr if needed
	if (stderrFile)
	{
		fflush(stderr);
		dup2(err, fileno(stderr));
		close(err);
	}
	
	// Return result
	return ranWithoutError;
}



/* void Multi::external(const Word& exe, const char* flags, const char* stdoutFile)
 *
 * Run external program
 */

void Multi::external(const Word& exe, const char* flags, const char* stdoutFile)
{
	
	// Set mpirun command
	Word mpirun;
	
	// Only run if there is more than one processor
	if (_jobSize > 1)
	{
		
		// Save mpirun call
		#ifdef MPIRUN
			mpirun = MPIRUN;
		#endif
		
		// Look for -np
		int i;
		if (mpirun.contains("-np", false))
		{
			for (i = 0; i < mpirun.length() - 2; ++i)
			{
				if ((mpirun[i] == '-') && ((mpirun[i+1] == 'n') || (mpirun[i+1] == 'N')) && \
					((mpirun[i+2] == 'p') || (mpirun[i+2] == 'P')))
				{
					mpirun.insert(" ", i+3);
					mpirun.insert(Language::numberToWord(_jobSize), i+4);
				}
			}
		}
		
		// Look for -n
		else if (mpirun.contains("-n", false))
		{
			for (i = 0; i < mpirun.length() - 1; ++i)
			{
				if ((mpirun[i] == '-') && ((mpirun[i+1] == 'n') || (mpirun[i+1] == 'N')))
				{
					mpirun.insert(" ", i+2);
					mpirun.insert(Language::numberToWord(_jobSize), i+3);
				}
			}
		}
	}
	
	// Add executable name
	if (mpirun.length() > 0)
		mpirun += " ";
	mpirun += exe;
	
	// Add flags
	if (flags)
	{
		mpirun += " ";
		mpirun += flags;
	}
	
	// Add output redirect
	if (stdoutFile)
	{
		mpirun += " &> ";
		mpirun += stdoutFile;
	}
	
	// Add wait
	mpirun += "; wait";
	
	// Call job for root only
	if (!_rank)
		system(mpirun.array());
	barrier();
}
