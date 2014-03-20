/* json.h -- Structure for json files
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
