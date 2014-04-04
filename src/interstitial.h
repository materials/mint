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
