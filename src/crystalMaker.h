/* crystalMaker.h -- Structure for CrystalMaker files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef CRYSTALMAKER_H
#define CRYSTALMAKER_H



#include "iso.h"
#include "text.h"



// CrystalMaker structure functions
class CrystalMaker
{	
public:
    static void write(const Word& file, const ISO& iso);
};



#endif
