/* findsym.h -- Structure for FindSym files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef FINDSYM_H
#define FINDSYM_H



#include "iso.h"
#include "text.h"



// FindSym structure functions
class FindSym
{	
public:
    static void write(const Word& file, const ISO& iso, double tol);
};



#endif
