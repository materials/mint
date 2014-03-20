/* interstitial.h -- Search for interstitial sites
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef INTERSTITIAL_H
#define INTERSTITIAL_H



#include "iso.h"
#include "symmetry.h"
#include "num.h"
#include "list.h"



// Class to compute interstitial sites
class Interstitial
{
	
	// Variables to store results
	OList<Vector3D> _sites;
	
	// Functions
	bool minimizePoint(Vector3D& point, const ISO& iso, double scale);
	Vector3D getStep(const ISO& iso, const Vector3D& point, Vector3D& deriv, Matrix3D& hessian, double scale);
	
public:
	
	// Functions
	void clear()	{ _sites.length(0); }
	void evaluate(const ISO& iso, const Symmetry& symmetry, int numPointsPerAtom, double tol, double scale = 1.0);
	void voronoi(const ISO& iso, const Symmetry& symmetry, double tol);
	
	// Access functions
	const OList<Vector3D>& sites() const	{ return _sites; }
};



#endif
