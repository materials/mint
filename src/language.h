/* language.h -- Deal with user input
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
