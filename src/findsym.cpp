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



#include "findsym.h"
#include "num.h"
#include "output.h"



/* void FindSym::write(const Word& file, const ISO& iso, double tol)
 *
 * Write structure to findsym file
 */

void FindSym::write(const Word& file, const ISO& iso, double tol)
{
	
	// Precision for printing numbers
	int prec = 14;
	
	// Setup output if file was set
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file.length() > 0)
	{
		
		// Open file for writing if needed
		if (file != "stdout")
			Output::setStream(Output::addStream(file));

		// Set output method
		Output::method(STANDARD);
	}
	
	// Write elements on comment line
	int i;
	Output::newline();
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		Output::print(iso.atoms()[i][0].element().symbol());
		Output::print(" ");
	}
	
	// Write tolerance
	Output::newline();
	Output::printSci(tol);
	
	// Print lattice parameters
	Output::newline();
	Output::print("2");
	Output::newline();
	for (i = 0; i < 3; ++i)
	{
		Output::print(iso.basis().lengths()[i]);
		Output::print(" ");
	}
	for (i = 0; i < 3; ++i)
	{
		Output::print(Num<double>::toDegrees(iso.basis().angles()[i]));
		Output::print(" ");
	}
	
	// Write form for lattice vectors
	Output::newline();
	Output::print("2");
	Output::newline();
	Output::print("P");
	
	// Print number of atoms
	Output::newline();
	Output::print(iso.numAtoms());
	
	// Print type of each atom
	int j;
	Output::newline();
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			Output::print(i + 1);
			Output::print(" ");
		}
	}
	
	// Make list of atom positions
	int k;
	Output message;
	message.addLines(iso.numAtoms());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			message.addLine();
			for (k = 0; k < 3; ++k)
				message.add(iso.atoms()[i][j].fractional()[k], prec);
		}
	}
	
	// Print coordinates
	Output::newline();
	Output::print(message, RIGHT);
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}
