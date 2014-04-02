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



#include "language.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>



/* bool Language::isNumber(const Word& input, bool allowFraction)
 *
 * Return if word is any number
 */

bool Language::isNumber(const Word& input, bool allowFraction)
{
	
	// String is empty
	if (!input.length())
		return false;
		
	// Check word
	bool foundNumber = false;
	for (int i = 0; i < input.length(); ++i)
	{
		
		// Found a number
		if ((input[i] == '0') || (input[i] == '1') || (input[i] == '2') || (input[i] == '3') || (input[i] == '4') || \
			(input[i] == '5') || (input[i] == '6') || (input[i] == '7') || (input[i] == '8') || (input[i] == '9'))
			foundNumber = true;
		
		// Found something else that could be in a number
		else if ((input[i] == '.') || (input[i] == 'e') || (input[i] == 'E') || (input[i] == 'd') || \
			(input[i] == 'D') || (input[i] == '+') || (input[i] == '-') || ((input[i] == '/') && (allowFraction)))
			continue;
		
		// Found something else
		else
			return false;
	}
	
	// Did whether a number was found
	return foundNumber;
}



/* bool Language::isInteger(const Word& input)
 *
 * Return if word is an integer
 */

bool Language::isInteger(const Word& input)
{
	
	// Check if a number
	if (!isNumber(input))
		return false;
	
	// Check for a .
	for (int i = 0; i < input.length(); ++i)
	{
		if (input[i] == '.')
			return false;
	}
	
	// Found an integer if at this point
	return true;
}



/* bool Language::isDecimal(const Word& input)
 *
 * Return if word is a decimal number
 */

bool Language::isDecimal(const Word& input)
{

	// Check if a number
	if (!isNumber(input))
		return false;
	
	// Check for a .
	for (int i = 0; i < input.length(); ++i)
	{
		if (input[i] == '.')
			return true;
	}
	
	// Did not find a decimal if at this point
	return false;
}



/* bool Language::isComment(const Word& input)
 * 
 * Return TRUE if a comment and FALSE otherwise
 */

bool Language::isComment(const Word& input)
{
	if (!input.length())
		return false;
	if (input[0] == '#')
		return true;
	return false;
}



/* Word Language::numberToWord(int number)
 *
 * Convert an integer number to Word object
 */

Word Language::numberToWord(int number)
{
	char buffer[20];
	sprintf(buffer, "%d", number);
	return Word(buffer);
}



/* Word Language::numberToFraction(double number, double tol, bool addDot)
 *
 * Convert a number to fractional representation
 */

Word Language::numberToFraction(double number, double tol, bool addDot)
{

    // Result
	Word res;
    
    // Add negative if needed
    if (number < 0)
        res += "-";
    number = fabs(number);
    
    // Number is very small
    if (number < 1e-12)
    {
		if (addDot)
			return "0.0";
		else
			return "0";
        return res;
    }
    
    // Number is very large
    if (number > 1e14)
    {
		if (addDot)
			res += "99999999999999.0";
		else
        	res += "99999999999999";
        return res;
    }

	// Number is 1
	if ((number < 1 + tol) && (number > 1 - tol))
	{
		if (addDot)
			res += "1.0";
		else
			res += "1";
		return res;
	}
 
    // Loop until fraction is found
	double temp;
    double Z = number;
    double denom = 1;
    double num = number;
	double prevDenom = 0;
	if (fabs(Z - (int)Z) > tol/10)
	{
		do
	    {
			Z = 1.0/(Z - (int)Z);
	        temp = denom;
	        denom = denom * (int)Z + prevDenom;
	        prevDenom = temp;
	        num = (int)(number * denom + 0.5);
	    }
		while ((fabs(number - num/denom) > tol) && (Z != (int)Z));
	}

    // Save result
	char buffer[20];
	if (fabs(fabs(denom) - 1) < tol)
	{
		if (addDot)
			sprintf(buffer, "%d%s", (int)num, ".0");
		else
			sprintf(buffer, "%d", (int)num);
	}
	else
	{
		if (addDot)
			sprintf(buffer, "%d%s%d%s", (int)num, ".0/", (int)denom, ".0");
		else
			sprintf(buffer, "%d%s%d", (int)num, "/", (int)denom);
	}
	res += buffer;
	
	// Return result
	return res;
}



/* double Language::fractionToNumber(const Word& input)
 * 
 * Turn a fraction into number
 */

double Language::fractionToNumber(const Word& input)
{

	// Get the numerator
	int i;
	char temp[25];
	for (i = 0; ((input[i] != '/') && (input[i] != '\0')); ++i)
		temp[i] = input[i];
	temp[i] = '\0';
	double res = atof(temp);
	
	// Get the denominator
	if (input[i] == '/')
	{
		int j = 0;
		for (++i; input[i] != '\0'; ++i, ++j)
			temp[j] = input[i];
		temp[j] = '\0';
		if (j)
			res /= atof(temp);
	}
	
	// Return result
	return res;
}
