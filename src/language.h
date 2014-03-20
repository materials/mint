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



#ifndef LANGUAGE_H
#define LANGUAGE_H



#include "text.h"



// Namespace with language functions
namespace Language
{
    
    // Tests
    bool isNumber(const Word& input, bool allowFraction = false);
    bool isInteger(const Word& input);
    bool isDecimal(const Word& input);
    bool isComment(const Word& input);
    
    // Functions
	Word numberToWord(int number);
	Word numberToFraction(double number, double tol = 1e-4, bool addDot = false);
	double fractionToNumber(const Word& input);
}



#endif
