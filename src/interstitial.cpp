/* interstitial.cpp -- Search for interstitial sites
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "multi.h"
#include "interstitial.h"
#include "constants.h"
#include "output.h"
#include <cmath>



/* void Interstitial::evaluate(const ISO& iso, const Symmetry& symmetry, int numPointsPerAtom, double tol, double scale)
 *
 * Find interstitial sites in a structure
 */

void Interstitial::evaluate(const ISO& iso, const Symmetry& symmetry, int numPointsPerAtom, double tol, double scale)
{
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Searching for interstitial sites using ");
	Output::print(numPointsPerAtom);
	Output::print(" starting point");
	if (numPointsPerAtom != 1)
		Output::print("s");
	Output::print(" per atom and a scale of ");
	Output::print(scale);
	Output::increase();
	
	// Constants used in generating points on sphere around each atom
	double phiScale = Constants::pi * (3 - sqrt(5));
	double yScale = 2.0 / numPointsPerAtom;
	
	// Loop over unique atoms in the structure
	int i, j, k;
	int count = 0;
	double y;
	double r;
	double phi;
	double curDistance;
	double nearDistance;
	double startDistance;
	Vector3D curPoint;
	Linked<Vector3D > points;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		
		// Get the distance to nearest atom in the structure
		nearDistance = -1;
		for (j = 0; j < iso.atoms().length(); ++j)
		{
			for (k = 0; k < iso.atoms()[j].length(); ++k)
			{
				curDistance = iso.basis().distance(symmetry.orbits()[i].atoms()[0]->fractional(), FRACTIONAL, \
					iso.atoms()[j][k].fractional(), FRACTIONAL);
				if (curDistance > 0)
				{
					if ((nearDistance == -1) || (curDistance < nearDistance))
						nearDistance = curDistance;
				}
			}
		}
		
		// Set the starting distance away from atom
		startDistance = nearDistance / 2;
		
		// Loop over starting points
		for (j = 0; j < numPointsPerAtom; ++j)
		{
			
			// Check if running current point on current processor
			if ((++count + Multi::rank()) % Multi::worldSize() == 0)
			{
				
				// Get current starting point
				y = j * yScale - 1 + (yScale / 2);
				r = sqrt(1 - y*y);
				phi = j * phiScale; 
				curPoint.set(symmetry.orbits()[i].atoms()[0]->cartesian()[0] + startDistance*r*cos(phi), \
					symmetry.orbits()[i].atoms()[0]->cartesian()[1] + startDistance*y, \
					symmetry.orbits()[i].atoms()[0]->cartesian()[2] + startDistance*r*sin(phi));
				
				// Minimize the current point
				if (!minimizePoint(curPoint, iso, scale))
					continue;
				
				// Save current point in fractional coordinates
				points += iso.basis().getFractional(curPoint);
				ISO::moveIntoCell(*points.last());
			}
		}
	}
	
	// Reduce list of points to unique ones
	int m;
	bool found;
	int numLoops;
	Vector3D rotPoint;
	Vector3D equivPoint;
	Vector3D origin(0.0);
	Linked<double> distances;
	Linked<double>::iterator itDist;
	Linked<Vector3D> uniquePoints;
	Linked<Vector3D>::iterator it;
	Linked<Vector3D>::iterator itUnique;
	for (i = 0; i < Multi::worldSize(); ++i)
	{
		
		// Send number of points in list on current processor
		numLoops = points.length();
		Multi::broadcast(numLoops, i);
		
		// Loop over points
		if (i == Multi::rank())
			it = points.begin();
		for (j = 0; j < numLoops; ++j)
		{
			
			// Send out current point
			if (i == Multi::rank())
			{
				curPoint = *it;
				++it;
			}
			Multi::broadcast(curPoint, i);
			
			// Get current distance to origin
			curDistance = iso.basis().distance(curPoint, FRACTIONAL, origin, FRACTIONAL);
			
			// Loop over points that were already saved
			found = false;
			itDist = distances.begin();
			itUnique = uniquePoints.begin();
			for (; itDist != distances.end(); ++itDist, ++itUnique)
			{
				
				// Current points are not the same
				if (Num<double>::abs(curDistance - *itDist) <= tol)
				{
					if (iso.basis().distance(curPoint, FRACTIONAL, *itUnique, FRACTIONAL) <= tol)
					{
						found = true;
						break;
					}
				}
				
				// Loop over symmetry operations
				for (k = 0; k < symmetry.operations().length(); ++k)
				{
					
					// Loop over translations
					rotPoint = symmetry.operations()[k].rotation() * curPoint;
					for (m = 0; m < symmetry.operations()[k].translations().length(); ++m)
					{
						
						// Check if points are the same
						equivPoint = rotPoint;
						equivPoint += symmetry.operations()[k].translations()[m];
						if (iso.basis().distance(equivPoint, FRACTIONAL, *itUnique, FRACTIONAL) <= tol)
						{
							found = true;
							break;
						}
					}
					if (found)
						break;
				}
				if (found)
					break;
			}
			
			// Found a new point
			if (!found)
			{
				distances += curDistance;
				uniquePoints += curPoint;
			}
		}
	}
	
	// Save unique points
	_sites.length(uniquePoints.length());
	for (i = 0, it = uniquePoints.begin(); it != uniquePoints.end(); ++i, ++it)
		_sites[i] = *it;
	
	// Output
	Output::newline();
	Output::print("Found ");
	Output::print(_sites.length());
	Output::print(" possible interstitial site");
	if (_sites.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _sites.length(); ++i)
	{
		Output::newline();
		Output::print("Site ");
		Output::print(i+1);
		Output::print(":");
		for (j = 0; j < 3; ++j)
		{
			Output::print(" ");
			Output::print(_sites[i][j], 8);
		}
	}
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* bool Interstitial::minimizePoint(Vector3D& point, const ISO& iso, double scale)
 *
 * Minimize point
 */

bool Interstitial::minimizePoint(Vector3D& point, const ISO& iso, double scale)
{
	
	// Loop until converged
	int i;
	int minIndex;
	int loopNum = 0;
	int maxLoops = 50;
	int curNegCount;
	int origNegCount;
	double mag;
	double tol = 1e-8;
	double stepSize = 1;
	double minEigenvalue;
	Vector3D step;
	Vector3D deriv;
	Vector3D eigenvalues;
	Vector3D minEigenvector;
	Matrix3D eigenvectors;
	Matrix3D hessian;
	for (; loopNum < maxLoops; ++loopNum)
	{
				
		// Get current step, derivatives, and second derivatives
		step = getStep(iso, point, deriv, hessian, scale);
		
		// Check if at a stationary point
		if (deriv.magnitude() < tol)
		{
			
			// Count the number of negative eigenvalues
			origNegCount = 0;
			eigenvalues = hessian.eigenvalues(&eigenvectors);
			for (i = 0; i < 3; ++i)
			{
				if (eigenvalues[i] < -tol)
					++origNegCount;
			}
			
			// Found a minimum so break
			if (origNegCount == 0)
				break;
			
			// Get the most negative eigenvalue
			minIndex = 0;
			if (eigenvalues[1] < eigenvalues[minIndex])
				minIndex = 1;
			if (eigenvalues[2] < eigenvalues[minIndex])
				minIndex = 2;
			minEigenvalue = eigenvalues[minIndex];
			minEigenvector[0] = eigenvectors(0, minIndex);
			minEigenvector[1] = eigenvectors(1, minIndex);
			minEigenvector[2] = eigenvectors(2, minIndex);
			
			// Make step along lowest eigenvector
			step = minEigenvector * stepSize;
			if (step.magnitude() < tol)
			{
				loopNum = maxLoops;
				break;
			}
			point -= step;
			
			// Perform gradient based minimization until there is one less negative eigenvalue
			for (++loopNum; loopNum < maxLoops; ++loopNum)
			{
				
				// Get current derivative
				getStep(iso, point, deriv, hessian, scale);
				
				// Count the number of negative eigenvalues
				curNegCount = 0;
				eigenvalues = hessian.eigenvalues();
				for (i = 0; i < 3; ++i)
				{
					if (eigenvalues[i] < -tol)
						++curNegCount;
				}
				
				// Break if a negative eigenvalue was removed
				if (curNegCount < origNegCount)
					break;
				
				// Make step
				step = deriv;
				step *= stepSize / deriv.magnitude();
				point -= step;
			}
		}
		
		// Make sure that step is not too large and update position
		else
		{
			mag = step.magnitude();
			if (mag > 0.25)
				step *= 0.25 / mag;
			point -= step;
		}
	}
	
	// Return that point was not found if not converged
	if (loopNum >= maxLoops)
		return false;
	
	// Return that point was found
	return true;
}



/* Vector3D Interstitial::getStep(const ISO& iso, const Vector3D& point, Vector3D& deriv, Matrix3D& H, double scale)
 *
 * Get derivatives of density function at current point
 */

Vector3D Interstitial::getStep(const ISO& iso, const Vector3D& point, Vector3D& deriv, Matrix3D& H, double scale)
{
	
	// Clear data
	deriv = 0.0;
	H = 0.0;
	
	// Get fractional coordinates of current point
	Vector3D fracPoint = iso.basis().getFractional(point);
	ISO::moveIntoCell(fracPoint);
	
	// Loop over atoms
	int i, j, k;
	double dis;
	double norm;
	double denom;
	double cutoff;
	double exponent;
	Vector3D dif;
	ImageIterator images;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		
		// Set cutoff for current element
		norm = scale * iso.atoms()[i][0].element().radius();
		cutoff = -log(1e-8) * norm;
		
		// Set image iterator
		images.setCell(iso.basis(), cutoff);
		
		// Loop over atoms of current element
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			
			// Reset image iterator for current atom
			images.reset(fracPoint, iso.atoms()[i][j].fractional());
			
			// Loop over images
			while (!images.finished())
			{
				
				// Get current value
				dis = ++images;
				if (dis < 1e-8)
					continue;
				exponent = exp(-dis / norm);
			
				// Get derivative vector
				dif = images.cartVector();
				dif *= -1;
				for (k = 0; k < 3; ++k)
					deriv[k] += -exponent * dif[k] / (norm * dis);
				
				// Get second derivative vector
				denom = norm * norm * pow(dis, 3);
				H(0, 0) += exponent * (dif[0]*dif[0]*dis - norm * (dif[1]*dif[1] + dif[2]*dif[2])) / denom;
			 	H(1, 1) += exponent * (dif[1]*dif[1]*dis - norm * (dif[0]*dif[0] + dif[2]*dif[2])) / denom;
				H(2, 2) += exponent * (dif[2]*dif[2]*dis - norm * (dif[0]*dif[0] + dif[1]*dif[1])) / denom;
				
				// Add to Hessian
				H(0, 1) += exponent * dif[0]*dif[1] * (norm + dis) / denom;
				H(0, 2) += exponent * dif[0]*dif[2] * (norm + dis) / denom;
				H(1, 2) += exponent * dif[1]*dif[2] * (norm + dis) / denom;
			}
		}
	}

	// Finish building Hessian
	H(1, 0) = H(0, 1);
	H(2, 0) = H(0, 2);
	H(2, 1) = H(1, 2);
	
	// Return Newton step
	return H.inverse() * deriv;
}



/* void Interstitial::voronoi(const ISO& iso, const Symmetry& symmetry, double tol)
 *
 * Search for interstitial sites using vertices of the Voronoi volumes around each atom 
 */

void Interstitial::voronoi(const ISO& iso, const Symmetry& symmetry, double tol)
{	
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Searching for interstitial sites using Voronoi method");
	Output::increase();
	
	// Set up image iterator
	ImageIterator images;
	images.setCell(iso.basis(), 12);
	
	// Loop over unique atoms in the structure
	int i, j, k;
	List<double> weights;
	OList<Vector3D> points;
	OList<Vector3D> vertices;
	Linked<Vector3D> intPoints;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		
		// Reset variables
		weights.length(0);
		points.length(0);
		vertices.length(0);
		
		// Loop over atoms in the structure
		for (j = 0; j < iso.atoms().length(); ++j)
		{
			for (k = 0; k < iso.atoms()[j].length(); ++k)
			{
				
				// Loop over images
				images.reset(symmetry.orbits()[i].atoms()[0]->fractional(), iso.atoms()[j][k].fractional());
				while (!images.finished())
				{
					
					// Skip if atoms are the same
					if (++images < 1e-8)
						continue;
					
					// Save current point
					points += symmetry.orbits()[i].atoms()[0]->cartesian() + images.cartVector();
					weights += 0.5;
				}
			}
		}
		
		// Calculate Voronoi volume
		symmetry.orbits()[i].atoms()[0]->cartesian().voronoi(points, weights, tol, &vertices);
		
		// Save points
		for (j = 0; j < vertices.length(); ++j)
		{
			intPoints += iso.basis().getFractional(vertices[j]);
			ISO::moveIntoCell(*intPoints.last());
		}
	}
	
	// Reduce points to unique ones
	bool found;
	double curDistance;
	Vector3D rotPoint;
	Vector3D equivPoint;
	Vector3D origin(0.0);
	Linked<double> distances;
	Linked<double>::iterator itDist;
	Linked<Vector3D>::iterator it;
	Linked<Vector3D> uniquePoints;
	Linked<Vector3D>::iterator itUnique;
	for (it = intPoints.begin(); it != intPoints.end(); ++it)
	{
		
		// Get current distance to origin
		curDistance = iso.basis().distance(*it, FRACTIONAL, origin, FRACTIONAL);
		
		// Loop over points that were already saved
		found = false;
		itDist = distances.begin();
		itUnique = uniquePoints.begin();
		for (; itDist != distances.end(); ++itDist, ++itUnique)
		{
			
			// Current points are not the same
			if (Num<double>::abs(curDistance - *itDist) <= tol)
			{
				if (iso.basis().distance(*it, FRACTIONAL, *itUnique, FRACTIONAL) <= tol)
				{
					found = true;
					break;
				}
			}
			
			// Loop over symmetry operations
			for (i = 0; i < symmetry.operations().length(); ++i)
			{
				
				// Loop over translations
				rotPoint = symmetry.operations()[i].rotation() * *it;
				for (j = 0; j < symmetry.operations()[i].translations().length(); ++j)
				{
					
					// Check if points are the same
					equivPoint = rotPoint;
					equivPoint += symmetry.operations()[i].translations()[j];
					if (iso.basis().distance(equivPoint, FRACTIONAL, *itUnique, FRACTIONAL) <= tol)
					{
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			if (found)
				break;
		}
		
		// Found a new point
		if (!found)
		{
			distances += curDistance;
			uniquePoints += *it;
		}
	}
	
	// Save unique points
	_sites.length(uniquePoints.length());
	for (i = 0, it = uniquePoints.begin(); it != uniquePoints.end(); ++i, ++it)
		_sites[i] = *it;
	
	// Output
	Output::newline();
	Output::print("Found ");
	Output::print(_sites.length());
	Output::print(" possible interstitial site");
	if (_sites.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _sites.length(); ++i)
	{
		Output::newline();
		Output::print("Site ");
		Output::print(i+1);
		Output::print(":");
		for (j = 0; j < 3; ++j)
		{
			Output::print(" ");
			Output::print(_sites[i][j], 8);
		}
	}
	Output::decrease();
	
	// Output
	Output::decrease();
}
