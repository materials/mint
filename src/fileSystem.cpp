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
#include "fileSystem.h"
#include "output.h"
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>



/* Text Read::text(const Word& file)
 *
 * Read text file into memory
 */

Text Read::text(const Word& file)
{
	
	// Variable to store result
	Text content;
	
	// Open file for reading
    ifstream input (file.array());
    if (!input.is_open())
    {
        Output::newline(ERROR);
        Output::print("Error opening file: ");
        Output::print(file);
        Output::quit();
    }
	
	// Count the number of lines
	int numLines = 0;
	int bufferSize = 10000;
	char buffer[bufferSize];
	while (input.getline(buffer, bufferSize))
		++numLines;
	
	// Allocate space for all lines
	content.addLines(numLines);
	
	// Return to beginning of file
	input.clear();
	input.seekg(0, ios::beg);

	// Loop over lines and save
	while (input.getline(buffer, bufferSize))
	{
		buffer[bufferSize-1] = '\0';
		content.addLine(buffer);
	}
	
	// Close file
	input.close();
	
	// Return result
	return content;
}



/* void File::copy(const Word& from, const Word& to)
 *
 * Copy one file to another
 */

void File::copy(const Word& from, const Word& to)
{
	
	// Only run on root
	if (Multi::rank() == 0)
	{
		
		// Open input file
		ifstream infile (from.array());
	    if (!infile.is_open())
	    {
	        Output::newline(ERROR);
	        Output::print("Could not open file ");
	        Output::print(from);
	        Output::print(" for reading");
	        Output::quit();
	    }

		// Open the output file
	    ofstream outfile (to.array());
	    if (!outfile.is_open())
	    {
	        Output::newline(ERROR);
	        Output::print("Could not open file ");
	        Output::print(to);
	        Output::print(" for writing");
	        Output::quit();
	    }

		// Copy contents
		outfile << infile.rdbuf();

	    // Close files
	    infile.close();
	    outfile.close();
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* bool File::exists(const Word& file, bool exitIfNotFound)
 *
 * Return whether a file exists
 */

bool File::exists(const Word& file, bool exitIfNotFound)
{
    
    // Try to open the file
    ifstream infile;
    infile.open(file.array());
    if (infile.is_open())
    {
        infile.close();
        return true;
    }
    
    // File does not exist and exit
    if (exitIfNotFound)
    {
        Output::newline(ERROR);
        Output::print("The file ");
        Output::print(file);
        Output::print(" does not exist");
        Output::quit();
    }
    
    // File does not exist and do not exit
    return false;
}



/* void File::create(const Word& file, bool clearIfExists)
 *
 * Create a file
 */

void File::create(const Word& file, bool clearIfExists)
{
	
	// Run on root processor only
	if (Multi::rank() == 0)
	{
		
		// Create the file if it does not exist or should be cleared
		if ((!exists(file)) || (clearIfExists))
		{
			ofstream outfile (file.array());
			outfile << "";
			outfile.close();
		}
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void File::remove(const Word& file)
 *
 * Remove a file
 */

void File::remove(const Word& file)
{
	
	// Only run on root
	if (Multi::rank() == 0)
	{
		
		// Remove file if it exists
		if (exists(file))
			std::remove(file.array());
	}
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* char* Directory::home()
 *
 * Return the name of the home directory
 */

char* Directory::home()
{
	char* res = getenv("HOME");
	if (res == 0)
	{
		Output::newline(ERROR);
		Output::print("Could not determine home directory");
		Output::quit();
	}
	return res;
}



/* bool Directory::exists(const Word& path)
 *
 * Return whether directory exists
 */
 
bool Directory::exists(const Word& path)
{
	
	// Open directory
    DIR* pdir = opendir(path.array());
    
    // Directory exists
    if (pdir)
    {
        closedir(pdir);
        return true;
    }
    
    // Directory does not exist
    return false;
}



/* void Directory::create(const Word& path, bool deleteIfExists)
 *
 * Create a directory at path
 * If the directory already exists and deleteIfExists == true, then delete first
 */
 
void Directory::create(const Word& path, bool deleteIfExists)
{
	
	// Run on root only
	if (Multi::rank() == 0)
	{
		
		// Removing directory
		if (deleteIfExists)
			recurseRemove(path);
		
		// Not sure why this works but is needed when using multiple nodes
		DIR *pdir = opendir(".");
		closedir(pdir);
		
		// Create the directory if it does not exist
		if (!exists(path))
		{
			if (mkdir(path.array(), 0777))
			{
				Output::newline(ERROR);
				Output::print("Failed to create directory: ");
				Output::print(path);
				Output::quit();
			}
		}
	}

	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void Directory::remove(const Word& path)
 *
 * Delete a directory tree and all files within
 */
 
void Directory::remove(const Word& path)
{
	
	// Remove file if it exists
	if (Multi::rank() == 0)
		recurseRemove(path);
	
	// Wait for all processors to reach this point
	Multi::barrier();
}



/* void Directory::recurseRemove(const Word& path)
 *
 * Recursively delete a directory tree and all files within
 */
 
void Directory::recurseRemove(const Word& path)
{
	
	// Variables
    DIR* pdir;
    struct dirent* pfiles;

    // Directory exists
    if (Directory::exists(path))
    {

        // Open directory
        pdir = opendir(path.array());

        // Loop over all file/folders
        while ((pfiles = readdir(pdir)))
        {
            if ((!Text::equal(pfiles->d_name, ".")) && (!Text::equal(pfiles->d_name, "..")))
				recurseRemove(makePath(path, pfiles->d_name));
        }

        // Close directory
        closedir(pdir);
    }

    // Remove file
	std::remove(path.array());
}



/* bool Directory::current(Word& path, int length)
 *
 * Return the path of the current directory
 */

bool Directory::current(Word& path)
{
	int len = 250;
	char temp[len];
	if (!getcwd(temp, len))
		return false;
	path = temp;
	return true;
}



/* Word Directory::makePath(const Word& base, const Word& append)
 *
 * Create path name
 */

Word Directory::makePath(const Word& base, const Word& append)
{
	if (!base.length())
		return append;
	if (!append.length())
		return base;
	if ((base.last() == '/') || (append[0] == '/'))
		return base + append;
	return base + '/' + append;
}



/* void Directory::change(const Word& newDir)
 *
 * Change the current directory
 */

void Directory::change(const Word& newDir)
{
	if (chdir(newDir.array()) != 0)
	{
		Output::newline(ERROR);
		Output::print("Failed to change directory to ");
		Output::print(newDir);
		Output::quit();
	}
}
