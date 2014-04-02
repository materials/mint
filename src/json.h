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



#ifndef JSON_H
#define JSON_H



#include "iso.h"
#include "text.h"
#include "num.h"
#include <cstdio>
#include <cstdlib>



// JSON structure functions
class JSON
{	
	static Words makeString(const Vector3D& vector, int prec);
public:
    static void write(const Word& file, const ISO& iso);
};



/* inline Words JSON::makeString(const Vector3D& vector, int prec)
 *
 * Convert number to string
 */

inline Words JSON::makeString(const Vector3D& vector, int prec)
{
	Words res;
	char arg[10];
	char temp[25];
	for (int i = 0; i < 3; ++i)
	{
		if (i != 2)
			sprintf(arg, "%s%d%s", "%.", prec, "f,");
		else
			sprintf(arg, "%s%d%s", "%.", prec, "f");
		sprintf(temp, arg, vector[i]);
		res += temp;
	}
	return res;
}



#endif
