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



#include "json.h"
#include "bonds.h"
#include "language.h"
#include "output.h"
#include "list.h"



/* void JSON::write(const Word& file, const ISO& iso)
 *
 * Write structure to json file
 */

void JSON::write(const Word& file, const ISO& iso)
{
	
	// Precision for printing numbers
	int prec = 8;
	
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
	
	// Print opening bracket
	Output::newline();
	Output::print("{");
	
	// Print basis line
	Output::newline();
	Output::tab();
	Output::print("\"basis\": [");
	
	// Loop over basis vectors to write
	int i, j, k;
	Words endStrings;
	Words startStrings;
	Vector3D mults;
	Vector3D start;
	Output message;
	message.addLines(12);
	for (i = 0; i < 4; ++i)
	{
		
		// Set starting point for current vertex and directions for each lattice vector
		if (i == 0)
		{
			start = 0.0;
			mults.set(1, 1, 1);
		}
		else if (i == 1)
		{
			start = iso.basis().vectors()[0];
			start += iso.basis().vectors()[1];
			mults.set(-1, -1, 1);
		}
		else if (i == 2)
		{
			start = iso.basis().vectors()[0];
			start += iso.basis().vectors()[2];
			mults.set(-1, 1, -1);
		}
		else if (i == 3)
		{
			start = iso.basis().vectors()[1];
			start += iso.basis().vectors()[2];
			mults.set(1, -1, -1);
		}
		
		// Add vectors
		startStrings = makeString(start, prec);
		for (j = 0; j < 3; ++j)
		{
			message.addLine();
			message.addTab();
			message.addTab();
			message.add("{");
			message.add("\"start\":");	
			message.add("[");
			for (k = 0; k < 3; ++k)
				message.add(startStrings[k]);
			message.add("],");
			message.add("\"end\":");	
			message.add("[");
			endStrings = makeString(start + Vector3D(iso.basis().vectors()[j]) * mults[j], prec);
			for (k = 0; k < 3; ++k)
				message.add(endStrings[k]);
			message.add("]");
			if ((i == 3) && (j == 2))
				message.add("}");
			else
				message.add("},");
		}
	}
	
	// Alignment
	List<PrintAlign> align(16);
	align.fill(RIGHT);
	align.last() = LEFT;
	
	// Print vectors
	Output::newline();
	Output::print(message, align);
	
	// Close basis line
	Output::newline();
	Output::tab();
	Output::print("],");
	
	// Print atoms line
	Output::newline();
	Output::tab();
	Output::print("\"atoms\": [");
	
	// Set atoms
	int m;
	Words posStrings;
	OList<Vector3D > vecsToAdd;
	vecsToAdd += Vector3D(0.0);
	message.clear();
	message.addLines(iso.numAtoms());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Get images to add
			vecsToAdd.length(1);
			if (iso.atoms()[i][j].fractional()[0] < 1e-2)
			{
				vecsToAdd += Vector3D(1, 0, 0);
				if (iso.atoms()[i][j].fractional()[1] < 1e-2)
				{
					vecsToAdd += Vector3D(1, 1, 0);
					if (iso.atoms()[i][j].fractional()[2] < 1e-2)
						vecsToAdd += Vector3D(1, 1, 1);
				}
				if (iso.atoms()[i][j].fractional()[2] < 1e-2)
					vecsToAdd += Vector3D(1, 0, 1);
			}
			if (iso.atoms()[i][j].fractional()[1] < 1e-2)
			{
				vecsToAdd += Vector3D(0, 1, 0);
				if (iso.atoms()[i][j].fractional()[2] < 1e-2)
					vecsToAdd += Vector3D(0, 1, 1);
			}
			if (iso.atoms()[i][j].fractional()[2] < 1e-2)
				vecsToAdd += Vector3D(0, 0, 1);
			
			// Add atoms
			for (k = 0; k < vecsToAdd.length(); ++k)
			{
				message.addLine();
				message.addTab();
				message.addTab();
				message.add("{");
				message.add("\"element\":");
				message.add(Word("\"") + iso.atoms()[i][j].element().symbol() + "\",");
				message.add("\"location\":");
				message.add("[");
				posStrings = makeString(iso.basis().getCartesian(iso.atoms()[i][j].fractional() + vecsToAdd[k]), prec);
				for (m = 0; m < 3; ++m)
					message.add(posStrings[m]);
				message.add("]");
				if ((i == iso.atoms().length() - 1) && (j == iso.atoms()[i].length() - 1) && \
					(k == vecsToAdd.length() - 1))
					message.add("}");
				else
					message.add("},");
			}
		}
	}
	
	// Alignment
	align.length(12);
	align.fill(RIGHT);
	align.last() = LEFT;
	
	// Print atoms
	Output::newline();
	Output::print(message, align);
	
	// Print end of atoms section
	Output::newline();
	Output::tab();
	Output::print("],");
	
	// Print bonds line
	Output::newline();
	Output::tab();
	Output::print("\"bonds\": [");
	
	
	
	// Print end of bonds
	Output::newline();
	Output::tab();
	Output::print("]");
	
	// Print closing bracket
	Output::newline();
	Output::print("}");
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}
