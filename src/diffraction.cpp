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



#include "multi.h"
#include "diffraction.h"
#include "language.h"
#include "output.h"
#include <cstdlib>



// Static member variables
double Diffraction::_resolution = 1e-2;



/* double Diffraction::set(const ISO& iso, const Symmetry& symmetry, const Diffraction* ref, bool fitBfactors)
 *
 * Set diffraction pattern given a structure
 */

double Diffraction::set(const ISO& iso, const Symmetry& symmetry, const Diffraction* ref, bool fitBfactors)
{
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Calculating peak intensities for the structure");
	Output::increase();
	
	// Pull data from reference pattern if available
	if (ref)
	{
		_method = ref->_method;
		_wavelength = ref->_wavelength;
		_minTwoTheta = ref->_minTwoTheta;
		_maxTwoTheta = ref->_maxTwoTheta;
	}
	
	// Set the peak locations
	setATFParams(symmetry);
	setPeakLocations(iso, symmetry);
	
	// If a reference pattern was passed then optimize max intensity and B factors
	double rFactor = 0;
	Vector Bfactors(symmetry.orbits().length(), 1.0);
	if (ref)
	{
		
		// Output
		Output::newline();
		Output::print("Optimizing against reference pattern");
		Output::increase();
		
		// Match peaks to reference pattern
		matchPeaksToReference(*ref);
		
		// Optimize the pattern
		optStructure(*ref, symmetry, fitBfactors, false, true);
		rFactor = curRFactor(*ref, DR_ABS);
		
		// Output
		Output::newline();
		Output::print("Optimal R factor: ");
		Output::print(rFactor);
		Output::decrease();
	}
	
	// Set the peak intensities if no reference was passed
	else
	{
		setPeakIntensities(symmetry, Bfactors, false, false);
		scaleIntensity(1000);
	}
	
	// Extract the intensities
	extractIntensities();
	
	// Print intensities
	Output::newline();
	Output::print("Generated ");
	Output::print(_twoTheta.length());
	Output::print(" peak");
	if (_twoTheta.length() != 1)
		Output::print("s");
	Output::increase();
	for (int i = 0; i < _twoTheta.length(); ++i)
	{
		Output::newline();
		Output::print("Two-theta and intensity of ");
		Output::print(_twoTheta[i]);
		Output::print(" ");
		Output::print(_intensity[i]);
	}
	Output::decrease();
	
	// Output
	Output::decrease();
	
	// Return the R factor
	return rFactor;
}



/* double Diffraction::refine(ISO& iso, Symmetry& symmetry, const Diffraction& reference, bool showWarnings)
 *
 * Refine a structure against reference pattern
 */

double Diffraction::refine(ISO& iso, Symmetry& symmetry, const Diffraction& reference, bool showWarnings)
{
	
	// Clear space
	clear();
	
	// Output
	Output::newline();
	Output::print("Refining structure against reference pattern");
	Output::increase();
	
	// Pull data from reference pattern if available
	_method = reference._method;
	_wavelength = reference._wavelength;
	_minTwoTheta = reference._minTwoTheta;
	_maxTwoTheta = reference._maxTwoTheta;
	
	// Set the peak locations
	setATFParams(symmetry);
	setPeakLocations(iso, symmetry);
	matchPeaksToReference(reference);
	
	// Optimize max intensity, B factors, and positions
	optStructure(reference, symmetry, true, true, showWarnings, &iso.basis().lengths());
	double rFactor = curRFactor(reference, DR_ABS);
	
	// Output
	Output::newline();
	Output::print("Optimal R factor: ");
	Output::print(rFactor);
	
	// Output
	Output::decrease();
	
	// Return R factor
	return rFactor;
}



/* void Diffraction::setPeakLocations(const ISO& iso, const Symmetry& symmetry)
 *
 * Set the peak locations for a pattern
 */

void Diffraction::setPeakLocations(const ISO& iso, const Symmetry& symmetry)
{
	
	// Clear space
	_peaks.length(0);
	
	// Calculate the hkl range
	int i, j;
	double range[3];
	double maxMag = 2 * sin(Num<double>::toRadians(_maxTwoTheta / 2)) / _wavelength;
	Vector3D vec;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			vec[j] = iso.basis().reducedInverse()(j, i);
		range[i] = Num<double>::abs(Num<double>::ceil(maxMag / vec.magnitude()));
	}
	
	// Conversion matrix to take reduced cell reciprocal lattice vector to unit cell reciprocal lattice vector
	Matrix3D convHKL = iso.basis().unitPointToReduced().transpose();
	
	// Generate symmetry operations for reduced cell
	Matrix3D P = iso.basis().unitToReduced().transpose();
	Matrix3D Q = P.inverse();
	OList<Matrix3D > operations(symmetry.operations().length());
	for (i = 0; i < symmetry.operations().length(); ++i)
	{
		operations[i] = P;
		operations[i] *= symmetry.operations()[i].rotation();
		operations[i] *= Q;
		operations[i] = operations[i].transpose();
	}
	
	// Remove identity operation
	Matrix3D identity;
	identity.makeIdentity();
	for (i = 0; i < operations.length(); ++i)
	{
		if (operations[i] == identity)
		{
			operations.remove(i);
			break;
		}
	}
	
	// Get the intrinsic part of all symmetry operations
	OList<Vector3D >::D2 translations(symmetry.operations().length());
	for (i = 0; i < symmetry.operations().length(); ++i)
	{
		translations[i].length(symmetry.operations()[i].translations().length());
		for (j = 0; j < symmetry.operations()[i].translations().length(); ++j)
			translations[i][j] = Symmetry::intrinsicTranslation(symmetry.operations()[i].rotation(), \
				symmetry.operations()[i].translations()[j]);
	}
	
	// Get the peak intensities
	bool found;
	int mult;
	double product;
	double twoTheta;
	Vector3D hkl;
	Vector3D redHKL;
	Vector3D symHKL;
	Linked<double>::iterator itTT;
	Linked<double>::iterator itI;
	Linked<Vector3D > equivPoints;
	Linked<Vector3D >::iterator itEquiv;
	for (redHKL[0] = -range[0]; redHKL[0] <= range[0]; ++redHKL[0])
	{
		for (redHKL[1] = -range[1]; redHKL[1] <= range[1]; ++redHKL[1])
		{
			for (redHKL[2] = -range[2]; redHKL[2] <= range[2]; ++redHKL[2])
			{
				
				// Loop over operations to generate equivalent points
				mult = 1;
				equivPoints.clear();
				equivPoints += redHKL;
				for (i = 0; i < operations.length(); ++i)
				{
					
					// Get new point
					symHKL = operations[i] * redHKL;
					for (j = 0; j < 3; ++j)
						symHKL[j] = Num<double>::round(symHKL[j], 1);
					
					// Generating an earlier point
					if (symHKL[0] < redHKL[0] - 1e-4)
						mult = 0;
					else if (Num<double>::abs(symHKL[0] - redHKL[0]) < 1e-4)
					{
						if (symHKL[1] < redHKL[1] - 1e-4)
							mult = 0;
						else if (Num<double>::abs(symHKL[1] - redHKL[1]) < 1e-4)
						{
							if (symHKL[2] < redHKL[2] - 1e-4)
								mult = 0;
						}
					}
					if (mult == 0)
						break;
					
					// Check if point was already used
					found = false;
					for (itEquiv = equivPoints.begin(); itEquiv != equivPoints.end(); ++itEquiv)
					{
						if ((Num<double>::abs((*itEquiv)[0] - symHKL[0]) < 1e-4) && \
							(Num<double>::abs((*itEquiv)[1] - symHKL[1]) < 1e-4) && \
							(Num<double>::abs((*itEquiv)[2] - symHKL[2]) < 1e-4))
						{
							found = true;
							break;
						}
					}
					if (!found)
					{
						++mult;
						equivPoints += symHKL;
					}
				}
				
				// Multiplier is zero so skip
				if (mult == 0)
					continue;
			
				// Convert current reduced basis hkl to unit cell
				hkl = convHKL * redHKL;
				
				// Check if direction will be a systematic absence
				found = false;
				for (i = 0; i < symmetry.operations().length(); ++i)
				{
					
					// Check if R*hkl = hkl
					symHKL = symmetry.operations()[i].rotation() * hkl;
					if ((Num<double>::abs(symHKL[0] - hkl[0]) > 1e-4) || \
						(Num<double>::abs(symHKL[1] - hkl[1]) > 1e-4) || \
						(Num<double>::abs(symHKL[2] - hkl[2]) > 1e-4))
						continue;
					
					// Loop over intrinsic translations and check if ti*hkl = integer for all
					for (j = 0; j < translations[i].length(); ++j)
					{
						product = translations[i][j] * hkl;
						if (Num<double>::abs(Num<double>::round(product, 1) - product) > 1e-4)
						{
							found = true;
							break;
						}
					}
					
					// Break if hkl is a system absence
					if (found)
						break;
				}
				
				// Current hkl will be a systematic absence
			//	if (found)
			//		continue;
			
				// Get the current angle
				twoTheta = 2 * Num<double>::toDegrees(angle(iso.basis(), hkl));
				
				// Skip if too small or too large
				if ((twoTheta < _minTwoTheta) || (twoTheta > _maxTwoTheta))
					continue;
				
				// Loop over known points to look for accidental overlap
				found = false;
				for (i = 0; i < _peaks.length(); ++i)
				{
					
					// Peaks overlap
					if (Num<double>::abs(twoTheta - _peaks[i][0].twoThetaDeg) < _resolution/2)
					{
						_peaks[i].add();
						_peaks[i].last().twoThetaDeg = twoTheta;
						_peaks[i].last().twoThetaRad = Num<double>::toRadians(twoTheta);
						_peaks[i].last().lpFactor = lpFactor(_peaks[i].last().twoThetaRad / 2.0);
						_peaks[i].last().multiplicity = mult;
						_peaks[i].last().hkl = hkl;
						found = true;
						break;
					}
				}
				
				// Found a new peak
				if (!found)
				{
					_peaks.add();
					_peaks.last().add();
					_peaks.last()[0].twoThetaDeg = twoTheta;
					_peaks.last()[0].twoThetaRad = Num<double>::toRadians(twoTheta);
					_peaks.last()[0].lpFactor = lpFactor(_peaks.last()[0].twoThetaRad / 2.0);
					_peaks.last()[0].multiplicity = mult;
					_peaks.last()[0].hkl = hkl;
				}
			}
		}
	}
}



/* void Diffraction::setPeakIntensities(const Symmetry& symmetry, const Vector& Bfactors, bool getBderivs,
 *		bool getPosDerivs, double scale)
 *
 * Calculate the peak intensities
 */

void Diffraction::setPeakIntensities(const Symmetry& symmetry, const Vector& Bfactors, bool getBderivs, \
	bool getPosDerivs, double scale)
{
	
	// Loop over peaks
	int i, j;
	int count = 0;
	for (i = 0; i < _peaks.length(); ++i)
	{
		for (j = 0; j < _peaks[i].length(); ++j)
		{
			
			// Get peak intensity
			if ((++count + Multi::rank()) % Multi::worldSize() == 0)
			{
				
				// Make sure derivs are allocated
				_peaks[i][j].derivBfactors.length(getBderivs == true ? Bfactors.length() : 0);
				_peaks[i][j].derivPositions.length(getPosDerivs == true ? 3 * symmetry.orbits().length() : 0); 
				_peaks[i][j].derivBfactors.fill(0.0);
				_peaks[i][j].derivPositions.fill(0.0);
				
				// Get intensity
				_peaks[i][j].intensity = structureFactorSquared(symmetry, _peaks[i][j].twoThetaRad/2, \
					_peaks[i][j].hkl, Bfactors, _peaks[i][j].derivBfactors, _peaks[i][j].derivPositions);
				_peaks[i][j].intensity *= scale * _peaks[i][j].multiplicity * _peaks[i][j].lpFactor;

				// Save B factor derivatives
				_peaks[i][j].derivBfactors *= scale * _peaks[i][j].multiplicity * _peaks[i][j].lpFactor;
				
				// Save position derivatives
				_peaks[i][j].derivPositions *= scale * _peaks[i][j].multiplicity * _peaks[i][j].lpFactor;
			}
		}
	}
	
	// Send peak intensities between processors
	int root = Multi::worldSize() - 1;
	for (i = 0; i < _peaks.length(); ++i)
	{
		for (j = 0; j < _peaks[i].length(); ++j)
		{
			Multi::broadcast(_peaks[i][j].intensity, root);
			Multi::broadcast(_peaks[i][j].derivBfactors, root);
			Multi::broadcast(_peaks[i][j].derivPositions, root);
			root = Num<int>::prev(root, Multi::worldSize() - 1);
		}
	}
}



/* void Diffraction::scaleIntensity(double max)
 *
 * Scale intensities to max value
 */

void Diffraction::scaleIntensity(double max)
{

	// Find the current maximum intensity
	int i, j;
	double curTot;
	double curMax = 0;
	for (i = 0; i < _peaks.length(); ++i)
	{
		curTot = _peaks[i][0].intensity;
		for (j = 1; j < _peaks[i].length(); ++j)
			curTot += _peaks[i][j].intensity;
		if (curTot > curMax)
			curMax = curTot;
	}
	
	// Scale all intensities
	double scale = (curMax == 0) ? 1 : max / curMax;
	for (i = 0; i < _peaks.length(); ++i)
	{
		for (j = 0; j < _peaks[i].length(); ++j)
		{
			_peaks[i][j].intensity *= scale;
			_peaks[i][j].derivBfactors *= scale;
			_peaks[i][j].derivPositions *= scale;
		}
	}
}



/* void Diffraction::extractIntensities()
 *
 * Pull the intensities from working peaks
 */

void Diffraction::extractIntensities()
{
	
	// Build a list of intensities
	int i, j;
	double maxIntensity = 0;
	List<double> intensity (_peaks.length());
	for (i = 0; i < _peaks.length(); ++i)
	{
		
		// Get current intensity
		intensity[i] = 0;
		for (j = 0; j < _peaks[i].length(); ++j)
			intensity[i] += _peaks[i][j].intensity;
		
		// Check if a new max
		if (intensity[i] > maxIntensity)
			maxIntensity = intensity[i];
	}
	
	// Save peaks
	_twoTheta.length(0);
	_intensity.length(0);
	for (i = 0; i < intensity.length(); ++i)
	{
		if (intensity[i] / maxIntensity > _minRelativeIntensity)
		{
			_twoTheta += _peaks[i][0].twoThetaDeg;
			_intensity += intensity[i];
		}
	}
	
	// Sort peaks
	sortByTwoTheta(0, _twoTheta.length() - 1);
}



/* void Diffraction::matchPeaksToReference(const Diffraction& reference)
 *
 * Match working peaks to reference pattern peaks
 */

void Diffraction::matchPeaksToReference(const Diffraction& reference)
{
	
	// Return if there are no reference peaks
	if (reference._twoTheta.length() == 0)
		return;
	
	// Tolerance for peaks to be aligned
	double tol = 0.15;
	
	// Loop until no changes are made to working peaks
	int i, j, k;
	bool changeMade = true;
	int localStart = 0;
	int nearIndex;
	double curDif;
	double nearDif;
	while (changeMade)
	{
		
		// Loop over working peaks
		changeMade = false;
		for (i = localStart; i < _peaks.length(); ++i)
		{
			
			// Loop over reference peaks to find nearest
			nearIndex = 0;
			nearDif = Num<double>::abs(_peaks[i][0].twoThetaDeg - reference._twoTheta[0]);
			for (j = 1; j < reference._twoTheta.length(); ++j)
			{
				curDif = Num<double>::abs(_peaks[i][0].twoThetaDeg - reference._twoTheta[j]);
				if (curDif < nearDif)
				{
					nearIndex = j;
					nearDif = curDif;
				}
			}
			
			// Skip if nearest peak is not within range
			if (nearDif > tol)
			{
				for (j = 0; j < _peaks[i].length(); ++j)
					_peaks[i][j].patternIndex = -1;
				++localStart;
				continue;
			}
			
			// Save pattern index for peak
			for (j = 0; j < _peaks[i].length(); ++j)
				_peaks[i][j].patternIndex = nearIndex;
			
			// Loop over working peaks and check for peak that matches same reference
			for (j = 0; j < i; ++j)
			{
				
				// Current peaks match to same reference
				if (nearIndex == _peaks[j][0].patternIndex)
				{
					for (k = 0; k < _peaks[i].length(); ++k)
						_peaks[j] += _peaks[i][k];
					_peaks.remove(i);
					changeMade = true;
					break;
				}
			}
			
			// Restart if change was made
			if (changeMade)
				break;
		}
	}
}



/* double Diffraction::structureFactorSquared(const Symmetry& symmetry, double angle, const Vector3D& hkl,
 *		const Vector& Bfactors, Vector& Bderivs, Vector& posDerivs)
 *
 * Get the squared structure factor
 */

double Diffraction::structureFactorSquared(const Symmetry& symmetry, double angle, const Vector3D& hkl, \
	const Vector& Bfactors, Vector& Bderivs, Vector& posDerivs)
{
	
	// Variables to store real and imaginary parts of derivatives
	Vector realBderivs (Bderivs.length(), 0.0);
	Vector imagBderivs (Bderivs.length(), 0.0);
	Vector realPosDerivs (posDerivs.length(), 0.0);
	Vector imagPosDerivs (posDerivs.length(), 0.0);
	
	// Loop over all atoms and calculate magnitude squared
	int i, j, k;
	double dot;
	double pre;
	double real = 0;
	double imag = 0;
	double sinTerm;
	double cosTerm;
	double preBderiv;
	double prePosDeriv;
	double thermFactor;
	double scatteringFactor;
	double thermFactorDeriv = 0;
	Atom* curAtom;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{

		// Get scaling factors for current atom
		curAtom = symmetry.orbits()[i].atoms()[0];
		scatteringFactor = atomicScatteringFactor(i, angle);
		thermFactor = (_method == DM_SIMPLE) ? 1 : thermalFactor(angle, Bfactors[i]);
		if (Bderivs.length() > 0)
			thermFactorDeriv = (_method == DM_SIMPLE) ? 0 : thermalFactorDeriv(angle, Bfactors[i]);
		
		// Loop over equivalent atoms
		for (j = 0; j < symmetry.orbits()[i].atoms().length(); ++j)
		{
	
			// Calculate contribution from current atom
			curAtom = symmetry.orbits()[i].atoms()[j];
			dot = 2 * Constants::pi * (hkl * curAtom->fractional());
			sinTerm = sin(dot);
			cosTerm = cos(dot);
			
			// Add contribution
			pre = scatteringFactor * thermFactor * curAtom->occupancy();
			real += pre * cosTerm;
			imag += pre * sinTerm;
			
			// Add the B factor derivative if passed
			if (Bderivs.length() > 0)
			{
				preBderiv = scatteringFactor * thermFactorDeriv * curAtom->occupancy();
				realBderivs[i] += preBderiv * cosTerm;
				imagBderivs[i] += preBderiv * sinTerm;
			}
			
			// Add position derivatives if passed
			if ((j == 0) && (posDerivs.length() > 0))
			{
				for (k = 0; k < 3; ++k)
				{
					prePosDeriv = pre * 2 * Constants::pi * hkl[k] * symmetry.orbits()[i].atoms().length();
					realPosDerivs[3*i + k] += prePosDeriv * sinTerm;
					imagPosDerivs[3*i + k] += prePosDeriv * cosTerm;
				}
			}
		}
	}
	
	// Save B factor derivatives
	for (i = 0; i < Bderivs.length(); ++i)
		Bderivs[i] = 2 * (real * realBderivs[i] + imag * imagBderivs[i]);
	
	// Save position derivatives
	if (posDerivs.length() > 0)
	{
		Vector3D tempDeriv;
		for (i = 0; i < posDerivs.length(); ++i)
			posDerivs[i] = 2 * (-real * realPosDerivs[i] + imag * imagPosDerivs[i]);
		for (i = 0; i < symmetry.orbits().length(); ++i)
		{
			for (j = 0; j < 3; ++j)
				tempDeriv[j] = posDerivs[3*i + j];
			tempDeriv *= symmetry.orbits()[i].specialPositions()[0].rotation();
			if (symmetry.orbits()[i].anyAtomsFixed())
				tempDeriv = 0.0;
			for (j = 0; j < 3; ++j)
				posDerivs[3*i + j] = tempDeriv[j];
		}
	}
	
	// Return the square of the magnitude
	return real*real + imag*imag;
}



/* void Diffraction::set(const Linked<double>& twoTheta, const Linked<double>& intensity)
 *
 * Set pattern from data
 */

void Diffraction::set(const Linked<double>& twoTheta, const Linked<double>& intensity)
{
	
	// Clear space
	clear();
	
	// Loop over two theta values and get max and min distances between them
	double curDif;
	double minDif = 0;
	double maxDif = 0;
	Linked<double>::iterator prev = twoTheta.begin();
	Linked<double>::iterator cur = prev + 1;
	if (twoTheta.length() >= 2)
		minDif = maxDif = *cur - *prev;
	for (++cur, ++prev; cur != twoTheta.end(); ++cur, ++prev)
	{
		curDif = *cur - *prev;
		if (curDif < minDif)
			minDif = curDif;
		else if (curDif > maxDif)
			maxDif = curDif;
	}
	
	// Save peaks if data is already processed
	if ((maxDif > 1.1*minDif) || (maxDif == 0))
	{
		_twoTheta.length(twoTheta.length());
		_intensity.length(intensity.length());
		Linked<double>::iterator itTheta = twoTheta.begin();
		Linked<double>::iterator itIntens = intensity.begin();
		for (int i = 0; itTheta != twoTheta.end(); ++itTheta, ++itIntens, ++i)
		{
			_twoTheta[i] = *itTheta;
			_intensity[i] = *itIntens;
		}
	}
	
	// Save peaks after processing
	else
	{
		
		// Output
		Output::newline();
		Output::print("Processing raw data");
		Output::increase();
		
		// Make a copy of the data to process
		Linked<double> twoThetaCopy = twoTheta;
		Linked<double> intensityCopy = intensity;
		
		// Process data
		smoothData(twoThetaCopy, intensityCopy);
		removeBackground(twoThetaCopy, intensityCopy);

		// Get peaks
		List<double>::D2 peakTwoTheta;
		List<double>::D2 peakIntensity;
		getPeaks(peakTwoTheta, peakIntensity, twoThetaCopy, intensityCopy);
		fitPeaks(peakTwoTheta, peakIntensity);
		
		// Output
		Output::decrease();
	}
	
	// Sort results
	sortByTwoTheta(0, _twoTheta.length() - 1);	
	
	// Save min and max two theta
	_minTwoTheta = _twoTheta[0] - _resolution/2;
	_maxTwoTheta = _twoTheta.last() + _resolution/2;
}



/* double Diffraction::rFactor(const Diffraction& reference)
 *
 * Return the R factor between two patterns
 */

double Diffraction::rFactor(const Diffraction& reference)
{
	
	// Output
	Output::newline();
	Output::print("Calculating R factor compared to reference pattern");
	Output::increase();
	
	// Save peaks as working peaks
	int i;
	_peaks.length(_twoTheta.length());
	for (i = 0; i < _twoTheta.length(); ++i)
	{
		_peaks[i].add();
		_peaks[i][0].twoThetaDeg = _twoTheta[i];
		_peaks[i][0].twoThetaRad = Num<double>::toRadians(_twoTheta[i]);
		_peaks[i][0].intensity = _intensity[i];
	}
	
	// Match peaks to reference pattern
	matchPeaksToReference(reference);
	
	// Get the R factor
	double optMaxIntensity;
	double rFactor = optIntensity(reference, &optMaxIntensity);
	
	// Output
	Output::newline();
	Output::print("Optimal R factor of ");
	Output::print(rFactor);
	Output::print(" with a max intensity of ");
	Output::print(optMaxIntensity);
	
	// Output
	Output::decrease();
	
	// Return result
	return rFactor;
}



/* double Diffraction::optIntensity(const Diffraction& reference, double* maxIntensity)
 *
 * Find the max intensity that optimizes match against reference pattern
 */

double Diffraction::optIntensity(const Diffraction& reference, double* maxIntensity)
{

	// Loop over working peaks and save scaling to reference peak
	int i, j;
	double workingMax = 0;
	double curIntensity;
	Linked<double> trialScales;
	for (i = 0; i < _peaks.length(); ++i)
	{
		
		// Get current intensity
		curIntensity = _peaks[i][0].intensity;
		for (j = 1; j < _peaks[i].length(); ++j)
			curIntensity += _peaks[i][j].intensity;
		
		// Found a new max
		if (curIntensity > workingMax)
			workingMax = curIntensity;
		
		// Current working peak does not match to any reference peaks
		if (_peaks[i][0].patternIndex < 0)
			continue;
		
		// Save scaling factor
		trialScales += reference._intensity[_peaks[i][0].patternIndex] / curIntensity;
	}
	
	// Patterns do not overlap
	if (trialScales.length() == 0)
		trialScales += 0;
	
	// Initialize best fit
	double prevScale = *trialScales.begin();
	double bestScale = *trialScales.begin();
	setNewIntensities(workingMax, bestScale*workingMax);
	double bestRfactor = curRFactor(reference, DR_ABS);
	
	// Loop over trials
	double curRfactor;
	for (Linked<double>::iterator itScale = trialScales.begin() + 1; itScale != trialScales.end(); ++itScale)
	{
		
		// Get current R factor
		setNewIntensities(prevScale*workingMax, *itScale * workingMax);
		prevScale = *itScale;
		curRfactor = curRFactor(reference, DR_ABS);
		
		// Found a new best
		if (curRfactor < bestRfactor)
		{
			bestRfactor = curRfactor;
			bestScale = *itScale;
		}
	}
	
	// Set peaks to best
	setNewIntensities(prevScale*workingMax, bestScale*workingMax);
	if (maxIntensity)
		*maxIntensity = bestScale * workingMax;
	
	// Return best fit
	return bestRfactor;
}



/* double Diffraction::optStructure(const Diffraction& reference, const Symmetry& symmetry, bool optBfactors,
 *		bool optPositions, bool showWarnings, const Vector3D* lengths)
 *
 * Optimize R factor by varying B factors and intensity scaling
 */

double Diffraction::optStructure(const Diffraction& reference, const Symmetry& symmetry, bool optBfactors, \
	bool optPositions, bool showWarnings, const Vector3D* lengths)
{
	
	// Set the initial intensities
	Vector Bfactors(symmetry.orbits().length(), 1.0);
	setPeakIntensities(symmetry, Bfactors, optBfactors, optPositions);
	
	// Get current max intensity
	int i, j;
	double curIntensity;
	double origMaxIntensity = 0;
	for (i = 0; i < _peaks.length(); ++i)
	{
		curIntensity = 0;
		for (j = 0; j < _peaks[i].length(); ++j)
			curIntensity += _peaks[i][j].intensity;
		if (curIntensity > origMaxIntensity)
			origMaxIntensity = curIntensity;
	}
	
	// Optimize the intensity
	double optMaxIntensity;
	optIntensity(reference, &optMaxIntensity);
	
	// Initialize scaling factor parameter
	double IScale = optMaxIntensity / origMaxIntensity;
	double nextIScale;
	double IScaleDeriv;
	double nextIScaleDeriv;
	double IScaleStep;
	double IScaleStepScale = IScale / 20;
	
	// Initialize B factor paramaters
	double BfactorScale = 0.1;
	Vector nextBfactors;
	Vector BfactorDerivs;
	Vector prevBfactorDerivs;
	Vector nextBfactorDerivs;
	Vector BfactorStep;
	Vector BfactorMaxStep;
	Vector cgVectorBfactors;
	Vector cgNormBfactors;
	if (optBfactors)
	{
		nextBfactors.length(symmetry.orbits().length());
		BfactorDerivs.length(symmetry.orbits().length());
		prevBfactorDerivs.length(symmetry.orbits().length());
		nextBfactorDerivs.length(symmetry.orbits().length());
		BfactorStep.length(symmetry.orbits().length());
		BfactorMaxStep.length(symmetry.orbits().length());
		BfactorMaxStep.fill(0.1);
		cgVectorBfactors.length(symmetry.orbits().length());
		cgNormBfactors.length(symmetry.orbits().length());
	}
	
	// Initialize position parameters
	Vector3D posScale;
	Vector positions;
	Vector nextPositions;
	Vector positionDerivs;
	Vector prevPositionDerivs;
	Vector nextPositionDerivs;
	Vector positionStep;
	Vector positionMaxStep;
	Vector cgVectorPositions;
	Vector cgNormPositions;
	if (optPositions)
	{
		positions.length(3 * symmetry.orbits().length());
		nextPositions.length(3 * symmetry.orbits().length());
		positionDerivs.length(3 * symmetry.orbits().length());
		prevPositionDerivs.length(3 * symmetry.orbits().length());
		nextPositionDerivs.length(3 * symmetry.orbits().length());
		positionStep.length(3 * symmetry.orbits().length());
		positionMaxStep.length(3 * symmetry.orbits().length());
		for (i = 0; i < 3; ++i)
		{
			positionMaxStep[i] = 0.1 / (*lengths)[i];
			for (j = 1; j < symmetry.orbits().length(); ++j)
				positionMaxStep[i + 3*j] = positionMaxStep[i];
			posScale[i] = positionMaxStep[i];
		}
		cgVectorPositions.length(3 * symmetry.orbits().length());
		cgNormPositions.length(3 * symmetry.orbits().length());
	}
	
	// Variable to store best R factor
	double rFactor = -1;
	
	// Make conjugate gradient moves
	bool scaleConverged;
	bool BfactorsConverged;
	bool positionsConverged;
	int curMove = -1;
	int prevCGMove = 0;
	int moveLoop;
	int loopNum = 0;
	double a;
	double b;
	double c;
	double norm;
	double scale;
	double projDeriv0;
	double projDeriv1;
	double projDeriv2;
	Vector3D scaleVec;
	int numVars = 1 + BfactorDerivs.length() + positionDerivs.length();
	int maxLoops = 50 * numVars;
	for (loopNum = 0; loopNum < maxLoops; ++loopNum)
	{
		
		// Set the current pattern and check for convergence
		setPeakIntensities(symmetry, Bfactors, optBfactors, optPositions, IScale);
		if ((optBfactors) && (optPositions))
			rFactor = curRFactor(reference, DR_SQUARED, &IScale, &IScaleDeriv, &BfactorDerivs, &positionDerivs);
		else if (optBfactors)
			rFactor = curRFactor(reference, DR_SQUARED, &IScale, &IScaleDeriv, &BfactorDerivs, 0);
		else if (optPositions)
			rFactor = curRFactor(reference, DR_SQUARED, &IScale, &IScaleDeriv, 0, &positionDerivs);
		else
			rFactor = curRFactor(reference, DR_SQUARED, &IScale, &IScaleDeriv, 0, 0);
		if (optPositions)
			symDerivatives(symmetry, positionDerivs);
		
		// Break if on last loop
		if (loopNum >= maxLoops - 1)
			break;
		
		// Check for convergence
		scaleConverged = isConverged(IScaleDeriv);
		BfactorsConverged = isConverged(BfactorDerivs, &Bfactors);
		positionsConverged = isConverged(positionDerivs);
		if ((scaleConverged) && (BfactorsConverged) && (positionsConverged))
			break;
		
		// Set current move
		if (!scaleConverged)
			curMove = 0;
		else if ((!BfactorsConverged) && ((prevCGMove != 1) || (positionsConverged)))
			curMove = prevCGMove = 1;
		else if (!positionsConverged)
			curMove = prevCGMove = 2;
		
		// Move B factors or positions
		if ((curMove == 1) || (curMove == 2))
		{
			
			// Set references to variables
			Vector& params     = (curMove == 1) ? Bfactors          : positions;
			Vector& nextParams = (curMove == 1) ? nextBfactors      : nextPositions;
			Vector& derivs     = (curMove == 1) ? BfactorDerivs     : positionDerivs;
			Vector& prevDerivs = (curMove == 1) ? prevBfactorDerivs : prevPositionDerivs;
			Vector& nextDerivs = (curMove == 1) ? nextBfactorDerivs : nextPositionDerivs;
			Vector& step       = (curMove == 1) ? BfactorStep       : positionStep;
			Vector& maxStep    = (curMove == 1) ? BfactorMaxStep    : positionMaxStep;
			Vector& cgVector   = (curMove == 1) ? cgVectorBfactors  : cgVectorPositions;
			Vector& cgNorm     = (curMove == 1) ? cgNormBfactors    : cgNormPositions;
			
			// Initialize variables
			nextParams.fill(0);
			prevDerivs.fill(0);
			nextDerivs.fill(0);
			step.fill(0);
			cgVector.fill(0);
			
			// Save positions if needed
			if (curMove == 2)
			{
				for (i = 0; i < symmetry.orbits().length(); ++i)
				{
					for (j = 0; j < 3; ++j)
						positions[3*i + j] = symmetry.orbits()[i].atoms()[0]->fractional()[j];
				}
			}
	
			// Perform conjugate gradient loop
			for (moveLoop = 0; (loopNum < maxLoops) && (moveLoop < Num<int>::max(5, params.length())); \
				++loopNum, ++moveLoop)
			{
				
				// Set the current scale
				scale = prevDerivs*prevDerivs;
				if (Num<double>::neq(scale, 0, 1e-12))
					scale = (derivs - prevDerivs)*derivs / scale;
				if ((scale > 0.5) || (scale < 0))
					scale = 0;
				
				// Set current conjugate gradient vector
				cgVector *= scale;
				cgVector -= derivs;
				cgNorm = cgVector;
				norm = cgNorm.magnitude();
				if (Num<double>::neq(norm, 0, 1e-12))
					cgNorm /= norm;
				
				// Set the step
				step = cgNorm;
				if (curMove == 1)
				{
					scale = BfactorScale;
					for (i = 0; i < step.length(); ++i)
					{
						if (Num<double>::abs(scale*step[i]) > maxStep[i])
							scale = maxStep[i] / Num<double>::abs(step[i]);
					}
					step *= scale;
				}
				else
				{
					scaleVec = posScale;
					for (i = 0; i < step.length(); ++i)
					{
						if (Num<double>::abs(scaleVec[(i+3)%3]*step[i]) > maxStep[i])
						{
							scale = maxStep[i] / Num<double>::abs(step[i]) / scaleVec[(i+3)%3];
							scaleVec *= scale;
						}
					}
					for (i = 0; i < step.length(); ++i)
						step[i] *= scaleVec[(i+3)%3];
				}
				
				// Update parameters
				nextParams = params;
				nextParams += step;
				
				// Set next pattern
				if (curMove == 1)
				{
					constrainBFactors(nextParams);
					setPeakIntensities(symmetry, nextParams, true, false, IScale);
					curRFactor(reference, DR_SQUARED, 0, 0, &nextDerivs, 0);
				}
				else
				{
					setPositions(symmetry, nextParams);
					setPeakIntensities(symmetry, Bfactors, false, true, IScale);
					curRFactor(reference, DR_SQUARED, 0, 0, 0, &nextDerivs);
					symDerivatives(symmetry, nextDerivs);
				}
				
				// Get projected derivatives
				projDeriv0 = derivs * cgNorm;
				projDeriv2 = nextDerivs * cgNorm;
				
				// Set next pattern
				step /= 2;
				nextParams = params;
				nextParams += step;
				if (curMove == 1)
				{
					constrainBFactors(nextParams);
					setPeakIntensities(symmetry, nextParams, true, false, IScale);
					curRFactor(reference, DR_SQUARED, 0, 0, &nextDerivs, 0);
				}
				else
				{
					setPositions(symmetry, nextParams);
					setPeakIntensities(symmetry, Bfactors, false, true, IScale);
					curRFactor(reference, DR_SQUARED, 0, 0, 0, &nextDerivs);
					symDerivatives(symmetry, nextDerivs);
				}

				// Get projected derivative
				projDeriv1 = nextDerivs * cgNorm;
				
				// Get optimal scaling
				a = 0.5*(projDeriv0 + projDeriv2) - projDeriv1;
				if (Num<double>::neq(a, 0, 1e-12))
				{
					b = -1.5*projDeriv0 + 2*projDeriv1 - 0.5*projDeriv2;
					c = projDeriv0;
					scale = b*b - 4*a*c;
					if (scale < 0)
						scale = 0;
					scale = (-b + sqrt(scale)) / (2*a) / 2;
				}
				else
				{
					scale = projDeriv0 - projDeriv2;
					if (Num<double>::eq(scale, 0, 1e-12))
						scale = 0.5;
					else
						scale = projDeriv0 / scale;
				}
				
				// Adjust trial step if needed
				if (Num<double>::abs(scale) < 0.25)
				{
					if (curMove == 1)
						BfactorScale /= 2;
					else
						posScale /= 2;
				}
				if (Num<double>::abs(scale) > 2)
				{
					if (curMove == 1)
						BfactorScale *= 2;
					else
						posScale *= 2;
				}
				
				// Make step
				step *= 2;
				if (curMove == 1)
				{
					for (i = 0; i < step.length(); ++i)
					{
						if (Num<double>::abs(scale*step[i]) > maxStep[i])
							scale = Num<double>::sign(scale) * maxStep[i] / Num<double>::abs(step[i]);
					}
					step *= scale;
				}
				else
				{
					scaleVec = scale;
					for (i = 0; i < step.length(); ++i)
					{
						if (Num<double>::abs(scaleVec[(i+3)%3]*step[i]) > maxStep[i])
						{
							scale = Num<double>::abs(maxStep[i] / step[i] / scaleVec[(i+3)%3]);
							scaleVec *= scale;
						}
					}
					for (i = 0; i < step.length(); ++i)
						step[i] *= scaleVec[(i+3)%3];
				}
				params += step;
				
				// Save current values
				prevDerivs = derivs;
				
				// Get current values
				if (curMove == 1)
				{
					constrainBFactors(params);
					setPeakIntensities(symmetry, Bfactors, true, false, IScale);
					rFactor = curRFactor(reference, DR_SQUARED, 0, 0, &derivs, 0);
				}
				else
				{
					setPositions(symmetry, positions);
					setPeakIntensities(symmetry, Bfactors, false, true, IScale);
					rFactor = curRFactor(reference, DR_SQUARED, 0, 0, 0, &derivs);
					symDerivatives(symmetry, derivs);
				}
				
				// Check for convergence
				if (isConverged(derivs, curMove == 1 ? &Bfactors : 0))
				{
					--loopNum;
					break;
				}
			}
		}
		
		// Move scaling parameter
	//	if ((curMove == 0) || (loopNum >= maxLoops))
		{
			
			// Make test move
			nextIScale = IScale;
			IScaleStep = -Num<double>::sign(IScaleDeriv) * IScaleStepScale;
			nextIScale += IScaleStep;
			
			// Set current pattern
			setPeakIntensities(symmetry, Bfactors, false, false, nextIScale);
			curRFactor(reference, DR_SQUARED, &nextIScale, &nextIScaleDeriv, 0, 0);
			
			// Adjust trial step if needed
			scale = IScaleDeriv - nextIScaleDeriv;
			if (Num<double>::neq(scale, 0, 1e-12))
			{
				if (Num<double>::abs(IScaleDeriv / scale) < 0.5)
					IScaleStepScale /= 2;
			}
			
			// Make move
			if (Num<double>::eq(scale, 0, 1e-12))
				IScaleStep /= 2;
			else
				IScaleStep *= IScaleDeriv / scale;
			IScale += IScaleStep;
		}
	}
	
	// Did not converge
	if ((loopNum >= maxLoops - 1) && (showWarnings))
	{
		Output::newline(WARNING);
		Output::print("Failed to converge during conjugate gradient loop in diffraction routine");
	}
	
	// Output
	if (optBfactors)
	{
		for (i = 0; i < Bfactors.length(); ++i)
		{
			Output::newline();
			Output::print("Optimized B factor for atom ");
			Output::print(symmetry.orbits()[i].atoms()[0]->atomNumber() + 1);
			Output::print(" (");
			Output::print(symmetry.orbits()[i].atoms()[0]->element().symbol());
			Output::print("): ");
			Output::print(Bfactors[i]);
		}
	}
	if (optPositions)
	{
		for (i = 0; i < symmetry.orbits().length(); ++i)
		{
			Output::newline();
			Output::print("Optimized position for atom ");
			Output::print(symmetry.orbits()[i].atoms()[0]->atomNumber() + 1);
			Output::print(" (");
			Output::print(symmetry.orbits()[i].atoms()[0]->element().symbol());
			Output::print("): ");
			for (j = 0; j < 3; ++j)
			{
				Output::print(symmetry.orbits()[i].atoms()[0]->fractional()[j]);
				if (j != 2)
					Output::print(", ");
			}
		}
	}
	
	// Return optimized R factor
	return rFactor;
}



/* void Diffraction::setPositions(const Symmetry& symmetry, const Vector& positions)
 *
 * Set positions during optimization
 */

void Diffraction::setPositions(const Symmetry& symmetry, const Vector& positions)
{
	int i, j, k;
	Vector3D newPos;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		for (j = 0; j < symmetry.orbits()[i].atoms().length(); ++j)
		{
			for (k = 0; k < 3; ++k)
				newPos[k] = positions[3*i + k];
			newPos *= symmetry.orbits()[i].generators()[j].rotation();
			newPos += symmetry.orbits()[i].generators()[j].translations()[0];
			ISO::moveIntoCell(newPos);
			symmetry.orbits()[i].atoms()[j]->fractional(newPos);
		}
	}
}



/* void Diffraction::symDerivatives(const Symmetry& symmetry, Vector& derivs)
 *
 * Make position derivatives obey symmetry
 */

void Diffraction::symDerivatives(const Symmetry& symmetry, Vector& derivs)
{
	int i, j;
	Vector3D tempDeriv;
	for (i = 0; i < symmetry.orbits().length(); ++i)
	{
		for (j = 0; j < 3; ++j)
			tempDeriv[j] = derivs[3*i + j];
		tempDeriv *= symmetry.orbits()[i].specialPositions()[0].rotation();
		if (symmetry.orbits()[i].anyAtomsFixed())
			tempDeriv = 0.0;
		for (j = 0; j < 3; ++j)
			derivs[3*i + j] = tempDeriv[j];
	}
}



/* double Diffraction::curRFactor(const Diffraction& reference, Rmethod method, double* scale, double* scaleDeriv,
 *		Vector* Bderivs, Vector* posDerivs)
 *
 * Set the current R factor
 */

double Diffraction::curRFactor(const Diffraction& reference, Rmethod method, double* scale, double* scaleDeriv, \
	Vector* Bderivs, Vector* posDerivs)
{
	
	// Variable to store result
	double rFactor = 0;
	Vector curBderivs;
	Vector curPosDerivs;
	if (Bderivs)
	{
		Bderivs->fill(0);
		curBderivs = *Bderivs;
	}
	if (posDerivs)
	{
		posDerivs->fill(0);
		curPosDerivs = *posDerivs;
	}
	
	// Clear scale derivative if passed
	if (scaleDeriv)
		*scaleDeriv = 0;
	
	// Variable to store which peaks in reference pattern are not used
	int i;
	bool refPeaksUsed[reference._twoTheta.length()];
	for (i = 0; i < reference._twoTheta.length(); ++i)
		refPeaksUsed[i] = false;
	
	// Loop over working peaks
	int j;
	double curIntensity;
	double intensityDif;
	for (i = 0; i < _peaks.length(); ++i)
	{
		
		// Get total intensity of current peak
		curIntensity = _peaks[i][0].intensity;
		for (j = 1; j < _peaks[i].length(); ++j)
			curIntensity += _peaks[i][j].intensity;
		
		// Current peak does not align with one in reference
		if (_peaks[i][0].patternIndex < 0)
			intensityDif = -curIntensity;
		
		// Get difference in peak intensities
		else
		{
			refPeaksUsed[_peaks[i][0].patternIndex] = true;
			intensityDif = reference._intensity[_peaks[i][0].patternIndex] - curIntensity;
		}
		
		// Update result
		if (method == DR_ABS)
			rFactor += Num<double>::abs(intensityDif);
		else if (method == DR_SQUARED)
			rFactor += intensityDif*intensityDif;
		
		// Save scale derivative
		if (scaleDeriv)
		{
			if (method == DR_SQUARED)
				*scaleDeriv += -2 * curIntensity * intensityDif / *scale;
		}
		
		// Save B factor derivative
		if (Bderivs)
		{			
			curBderivs = _peaks[i][0].derivBfactors;
			for (j = 1; j < _peaks[i].length(); ++j)
				curBderivs += _peaks[i][j].derivBfactors;
			for (j = 0; j < Bderivs->length(); ++j)
			{
				if (method == DR_SQUARED)
					(*Bderivs)[j] += -2 * intensityDif * curBderivs[j];
			}
		}
		
		// Save position derivatives
		if (posDerivs)
		{
			curPosDerivs = _peaks[i][0].derivPositions;
			for (j = 1; j < _peaks[i].length(); ++j)
				curPosDerivs += _peaks[i][j].derivPositions;
			for (j = 0; j < posDerivs->length(); ++j)
			{
				if (method == DR_SQUARED)
					(*posDerivs)[j] += -2 * intensityDif * curPosDerivs[j];
			}
		}
	}
	
	// Loop over peaks in reference pattern
	double norm = 0;
	for (i = 0; i < reference._twoTheta.length(); ++i)
	{
		
		// Peak was not matched to one in current pattern
		if (refPeaksUsed[i] == false)
		{
			if (method == DR_ABS)
				rFactor += Num<double>::abs(reference._intensity[i]);
			else if (method == DR_SQUARED)
				rFactor += reference._intensity[i] * reference._intensity[i];
		}
		
		// Update norm
		if (method == DR_ABS)
			norm += Num<double>::abs(reference._intensity[i]);
		else
			norm += reference._intensity[i] * reference._intensity[i];
	}
	
	// Scale derivs if needed
	if (Num<double>::neq(norm, 0, 1e-8))
	{
		if (scaleDeriv)
			*scaleDeriv /= norm;
		if (Bderivs)
			*Bderivs /= norm;
		if (posDerivs)
			*posDerivs /= norm;
	}
	
	// Return R factor
	return (Num<double>::eq(norm, 0, 1e-8)) ? rFactor : rFactor / norm;
}



/* void Diffraction::setNewIntensities(double oldMax, double newMax)
 *
 * Loop over intensities and set new scales
 */

void Diffraction::setNewIntensities(double oldMax, double newMax)
{
	int i, j;
	double scale = newMax / oldMax;
	for (i = 0; i < _peaks.length(); ++i)
	{
		for (j = 0; j < _peaks[i].length(); ++j)
		{
			_peaks[i][j].intensity *= scale;
			_peaks[i][j].derivBfactors *= scale;
			_peaks[i][j].derivPositions *= scale;
		}
	}
}



/* bool Diffraction::isFormat(const Text& text)
 *
 * Return whether file contains diffraction data
 */

bool Diffraction::isFormat(const Text& text)
{
	int pairCount = 0;
	int lineCount = 0;
	for (int i = 0; i < text.length(); ++i)
	{
		if (!text[i].length())
			continue;
		if (Language::isComment(text[i][0]))
			continue;
		++lineCount;
		if (text[i].length() < 2)
			continue;
		if ((Language::isNumber(text[i][0])) && (Language::isNumber(text[i][1])))
			++pairCount;
	}
	if (!lineCount)
		return false;
	if ((double)pairCount / lineCount < 0.5)
		return false;
	return true;
}



/* void Diffraction::set(const Text& text)
 *
 * Set diffraction pattern from file
 */

void Diffraction::set(const Text& text)
{
	
	// Output
	Output::newline();
	Output::print("Reading diffraction data from file");
	Output::increase();
	
	// Clear space
	clear();
	
	// Loop over lines in file
	int i;
	Linked<double> rawTwoTheta;
	Linked<double> rawIntensity;
	for (i = 0; i < text.length(); ++i)
	{
		
		// Skip if line is blank
		if (!text[i].length())
			continue;
		
		// Skip if line is too short
		if (text[i].length() < 2)
			continue;
		
		// Found wavelength
		if (text[i][0].equal("wavelength", false, 4))
		{
			if (Language::isNumber(text[i][1]))
				_wavelength = atof(text[i][1].array());
			else
			{
				Output::newline(ERROR);
				Output::print("Did not recognize wavelength value in diffraction file (");
				Output::print(text[i][1]);
				Output::print(")");
				Output::quit();
			}
		}
		
		// Found fwhm
		else if (text[i][0].equal("fwhm", false, 4))
		{
			if (Language::isNumber(text[i][1]))
				fwhm(atof(text[i][1].array()));
			else
			{
				Output::newline(ERROR);
				Output::print("Did not recognize FWHM value in diffraction file (");
				Output::print(text[i][1]);
				Output::print(")");
				Output::quit();
			}
		}
		
		// Found variance
		else if (text[i][0].equal("variance", false, 3))
		{
			if (Language::isNumber(text[i][1]))
				variance(atof(text[i][1].array()));
			else
			{
				Output::newline(ERROR);
				Output::print("Did not recognize variance value in diffraction file (");
				Output::print(text[i][1]);
				Output::print(")");
				Output::quit();
			}
		}
		
		// Found a data line
		else if ((Language::isNumber(text[i][0])) && (Language::isNumber(text[i][1])))
		{
			rawTwoTheta += atof(text[i][0].array());
			rawIntensity += atof(text[i][1].array());
		}
	}
	
	// Process data
	set(rawTwoTheta, rawIntensity);
	
	// Output
	Output::newline();
	Output::print("Found ");
	Output::print(_twoTheta.length());
	Output::print(" peak");
	if (_twoTheta.length() != 1)
		Output::print("s");
	Output::increase();
	for (i = 0; i < _twoTheta.length(); ++i)
	{
		Output::newline();
		Output::print("Two-theta and intensity of ");
		Output::print(_twoTheta[i]);
		Output::print(" ");
		Output::print(_intensity[i]);
	}
	Output::decrease();
	
	// Output
	Output::decrease();
}



/* void Diffraction::smoothData(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity)
 *
 * Apply smoothing function to intensities
 */

void Diffraction::smoothData(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity)
{
	
	// Get average spacing
	double spacing = 0;
	Linked<double>::iterator prevTwoTheta = rawTwoTheta.begin();
	Linked<double>::iterator itTwoTheta = prevTwoTheta + 1;
	Linked<double>::iterator prevIntensity = rawIntensity.begin();
	Linked<double>::iterator itIntensity = prevIntensity + 1;
	for (; itTwoTheta != rawTwoTheta.end(); ++itTwoTheta, ++prevTwoTheta, ++itIntensity, ++prevIntensity)
		spacing += (*itTwoTheta - *prevTwoTheta) / (rawTwoTheta.length() - 1);
	
	// Number of points to smooth over (should be an odd number)
	int numSmoothPoints = 2*(int)Num<double>::ceil(0.5 / spacing) + 1;
	
	// Weight for point at max distance away
	double farWeight = 0.05;
	
	// Set scale for Gaussian used in smoothing
	double scale = -0.5*0.5 / log(farWeight);
	
	// Variables to store two-theta and intensity values
	int i;
	double twoTheta[numSmoothPoints];
	double intensity[numSmoothPoints];
	for (i = 0; i < numSmoothPoints; ++i)
		twoTheta[i] = intensity[i] = -1;
	
	// Initialize fitting points
	int fitIndex = numSmoothPoints / 2;
	i = fitIndex;
	itTwoTheta = rawTwoTheta.begin();
	itIntensity = rawIntensity.begin();
	for (; i < numSmoothPoints; ++i, ++itTwoTheta, ++itIntensity)
	{
		twoTheta[i] = *itTwoTheta;
		intensity[i] = *itIntensity;
	}
	
	// Variable to store new points
	double newIntensity[rawIntensity.length()];
	newIntensity[0] = smoothPoint(twoTheta, intensity, scale, numSmoothPoints);
	
	// Set new points
	int intensityIndex = 1;
	for (; itIntensity != rawIntensity.end(); ++itIntensity, ++itTwoTheta)
	{
		
		// Move current points
		for (i = 0; i < numSmoothPoints - 1; ++i)
		{
			twoTheta[i] = twoTheta[i+1];
			intensity[i] = intensity[i+1];
		}
		
		// Add new point
		twoTheta[numSmoothPoints - 1] = *itTwoTheta;
		intensity[numSmoothPoints - 1] = *itIntensity;
		
		// Smooth point
		newIntensity[intensityIndex++] = smoothPoint(twoTheta, intensity, scale, numSmoothPoints);
	}
	
	// Set final points
	int j;
	for (i = numSmoothPoints - 1; i > fitIndex; --i)
	{
		
		// Move current points
		for (j = 0; j < numSmoothPoints - 1; ++j)
		{
			twoTheta[j] = twoTheta[j+1];
			intensity[j] = intensity[j+1];
		}
		
		// Smooth point
		newIntensity[intensityIndex++] = smoothPoint(twoTheta, intensity, scale, numSmoothPoints);
	}
	
	// Save smoothed points
	for (i = 0, itIntensity = rawIntensity.begin(); itIntensity != rawIntensity.end(); ++i, ++itIntensity)
		*itIntensity = newIntensity[i];
}



/* double Diffraction::smoothPoint(const double* twoTheta, const double* intensity, double scale, int numPoints)
 *
 * Smooth single point
 */

double Diffraction::smoothPoint(const double* twoTheta, const double* intensity, double scale, int numPoints)
{
	
	// Get the peak index
	int i;
	int peakIndex = numPoints / 2;
	for (i = numPoints - 1; i > peakIndex; --i)
	{
		if (intensity[i] < 0)
			--peakIndex;
	}
	
	// Loop over points
	double weight;
	double totalWeight = 0;
	double totalIntensity = 0;
	for (i = 0; i < numPoints; ++i)
	{
		
		// Skip if negative intensity
		if (intensity[i] < 0)
			continue;
		
		// Calculate value
		weight = exp(-pow(twoTheta[i] - twoTheta[peakIndex], 2.0) / scale);
		totalIntensity += weight * intensity[i];
		totalWeight += weight;
	}
	
	// Return smoothed value
	return totalIntensity / totalWeight;
}



/* void Diffraction::removeBackground(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity)
 *
 * Remove background data
 */

void Diffraction::removeBackground(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity)
{
	
	// Get max deviation between two points
	double curDev;
	double maxDev = 0;
	Linked<double>::iterator prevIntensity = rawIntensity.begin();
	Linked<double>::iterator itIntensity = prevIntensity + 1;
	for (; itIntensity != rawIntensity.end(); ++itIntensity, ++prevIntensity)
	{
		curDev = Num<double>::abs(*itIntensity - *prevIntensity);
		if (curDev > maxDev)
			maxDev = curDev;
	}
	double devTol = 0.03 * maxDev;
	
	// Make copy of points to use in fit
	Linked<double> fitTwoTheta(rawTwoTheta);
	Linked<double> fitIntensity(rawIntensity);
	
	// Loop over points and look for those that will be used in fit
	int i;
	double prevValue;
	itIntensity = fitIntensity.begin() + 1;
	Linked<double>::iterator itTwoTheta = fitTwoTheta.begin() + 1;
	Linked<double>::iterator itRemoveTwoTheta;
	Linked<double>::iterator itRemoveIntensity;
	for (; itTwoTheta != fitTwoTheta.end(); ++itTwoTheta, ++itIntensity)
	{
		
		// Check if current deviation is outside fit range
		if (Num<double>::abs(*itIntensity - *(itIntensity - 1)) > devTol)
		{
			
			// Remove points
			for (i = -(int)Num<double>::sign(*itIntensity - *(itIntensity - 1)); i <= 1; i += 2)
			{
				while (1)
				{		
					prevValue = *itIntensity;
					itRemoveTwoTheta = itTwoTheta;
					itRemoveIntensity = itIntensity;
					++itTwoTheta;
					++itIntensity;
					fitTwoTheta.remove(itRemoveTwoTheta);
					fitIntensity.remove(itRemoveIntensity);
					if (itIntensity == fitIntensity.end())
						break;
					if (i*(*itIntensity - prevValue) > -devTol)
						break;
				}
				if (itIntensity == fitIntensity.end())
					break;
			}
		}
		if (itIntensity == fitIntensity.end())
			break;
	}
	
	// Save fitting points
	itTwoTheta = fitTwoTheta.begin();
	itIntensity = fitIntensity.begin();
	List<double>::D2 fitPoints(fitTwoTheta.length());
	for (i = 0; itTwoTheta != fitTwoTheta.end(); ++itTwoTheta, ++itIntensity, ++i)
	{
		fitPoints[i].length(2);
		fitPoints[i][0] = *itTwoTheta;
		fitPoints[i][1] = *itIntensity;
	}
	
	// Get fitting parameters
	int start = 0;
	int numTerms = (fitTwoTheta.length() > 16) ? 4 : fitTwoTheta.length() / 4;
	Vector fit = Fit::polynomial(fitPoints, start, numTerms);
	
	// Subtract background
	double tot;
	double curPower;
	itTwoTheta = rawTwoTheta.begin();
	itIntensity = rawIntensity.begin();
	for (; itTwoTheta != rawTwoTheta.end(); ++itTwoTheta, ++itIntensity)
	{
		tot = 0;
		if (numTerms != 0)
		{
			curPower = pow(*itTwoTheta, start);
			tot += fit[0] * curPower;
		}
		for (i = 1; i < numTerms; ++i)
		{
			curPower *= *itTwoTheta;
			tot += fit[i] * curPower;
		}		
		*itIntensity -= tot;
	}
}



/* void Diffraction::getPeaks(List<double>::D2& peakTwoTheta, List<double>::D2& peakIntensity,
 *		const Linked<double>& rawTwoTheta, const Linked<double>& rawIntensity)
 *
 * Get peaks from raw data
 */

void Diffraction::getPeaks(List<double>::D2& peakTwoTheta, List<double>::D2& peakIntensity, \
	const Linked<double>& rawTwoTheta, const Linked<double>& rawIntensity)
{
	
	// Tolerance for point being on a peak
	double peakTol = 0.03;
	
	// Clear space for results
	peakTwoTheta.clear();
	peakIntensity.clear();
	
	// Loop over points to get maximum
	double max = 0;
	Linked<double>::iterator itIntensity = rawIntensity.begin();
	for (; itIntensity != rawIntensity.end(); ++itIntensity)
	{
		if (*itIntensity > max)
			max = *itIntensity;
	}
	
	// Set absolute peak tolerance
	peakTol *= max;
	
	// Loop over points and find peaks
	double prev;
	itIntensity = rawIntensity.begin() + 1;
	Linked<double>::iterator itTwoTheta = rawTwoTheta.begin() + 1;
	for (; itTwoTheta != rawTwoTheta.end(); ++itTwoTheta, ++itIntensity)
	{
		
		// Current point is on a peak
		if (*itIntensity > peakTol)
		{
			
			// Loop backwards until negative slope is found
			prev = *itIntensity;
			for (--itTwoTheta, --itIntensity; itTwoTheta != rawTwoTheta.begin(); --itTwoTheta, --itIntensity)
			{
				if ((prev - *itIntensity < 0) || (*itIntensity < 0))
				{
					++itTwoTheta;
					++itIntensity;
					break;
				}
				prev = *itIntensity;
			}
			
			// Add peak and point
			peakTwoTheta.add();
			peakIntensity.add();
			peakTwoTheta.last() += *itTwoTheta;
			peakIntensity.last() += *itIntensity;
			
			// Add points while slope is positive
			for (++itTwoTheta, ++itIntensity; itTwoTheta != rawTwoTheta.end(); ++itTwoTheta, ++itIntensity)
			{
				if (*itIntensity - peakIntensity.last().last() < 0)
					break;
				peakTwoTheta.last() += *itTwoTheta;
				peakIntensity.last() += *itIntensity;
			}
			
			// Remove current peak if there are no negative slope points
			if (itTwoTheta == rawTwoTheta.end())
			{
				peakTwoTheta.remove(peakTwoTheta.length() - 1);
				peakIntensity.remove(peakIntensity.length() - 1);
				break;
			}
			
			// Add points while slope is negative
			for (; itTwoTheta != rawTwoTheta.end(); ++itTwoTheta, ++itIntensity)
			{
				if ((*itIntensity - peakIntensity.last().last() > 0) || (*itIntensity < 0))
				{
					--itTwoTheta;
					--itIntensity;
					break;
				}
				peakTwoTheta.last() += *itTwoTheta;
				peakIntensity.last() += *itIntensity;
			}
		}
		
		// Break if needed
		if (itTwoTheta == rawTwoTheta.end())
			break;
	}
}



/* void Diffraction::fitPeaks(const List<double>::D2& peakTwoTheta, const List<double>::D2& peakIntensity)
 *
 * Fit functions to each peak to get intensities and locations
 */

void Diffraction::fitPeaks(const List<double>::D2& peakTwoTheta, const List<double>::D2& peakIntensity)
{
	
	// Allocate space
	_twoTheta.length(peakTwoTheta.length());
	_intensity.length(peakTwoTheta.length());
	
	// Functors used in fitting
	Functor<Diffraction> gaussFun(this, &Diffraction::gaussian);
	Functor<Diffraction> psFun(this, &Diffraction::PS);
	Functor<Diffraction> psTT(this, &Diffraction::PStwoTheta);
	VectorFunctor<Diffraction> gaussDeriv(this, &Diffraction::gaussianDerivs);
	VectorFunctor<Diffraction> psDeriv(this, &Diffraction::PSderivs);
	
	// Loop over peaks
	int i, j;
	double temp;
	double initialTwoTheta;
	double twoThetaStep;
	List<double>::D2 points;
	Vector gaussianParams;
	Vector initialGaussian(3);
	Vector initialPS(8);
	Functor<Diffraction> ttFun (this, &Diffraction::PStwoTheta);
	for (i = 0; i < peakTwoTheta.length(); ++i)
	{
		
		// Save two-theta/intensity pairs
		points.length(peakTwoTheta[i].length());
		for (j = 0; j < peakTwoTheta[i].length(); ++j)
		{
			points[j].length(2);
			points[j][0] = peakTwoTheta[i][j];
			points[j][1] = peakIntensity[i][j];
		}
		
		// Initialize fitting parameters for Gaussian function
		initialGaussian[0] = 0.25;
		initialGaussian[1] = points[0][0];
		initialGaussian[2] = points[0][1];
		for (j = 1; j < points.length(); ++j)
		{
			if (points[j][1] > initialGaussian[2])
			{
				initialGaussian[1] = points[j][0];
				initialGaussian[2] = points[j][1];
			}
		}
		
		// Fit Gaussian function
		gaussianParams = Fit::LM<Diffraction>(points, gaussFun, gaussDeriv, initialGaussian, 1e-5);
		
		// Initialize pseudo-voigt fitting parameters
		initialPS[0] = initialPS[1] = initialPS[2] = 0.01;
		initialPS[3] = gaussianParams[1];
		temp = tan(Num<double>::toRadians(initialPS[3]/2));
		initialPS[4] = initialPS[5] = initialPS[6] = gaussianParams[0] / (1 + temp + temp*temp);
		initialPS[7] = gaussianParams[2];
		
		// Fit data
		_PSparams = Fit::LM<Diffraction>(points, psFun, psDeriv, initialPS, 1e-5);

		// Save peak
		initialTwoTheta = _PSparams[3];
		twoThetaStep = 1e-3;	
		_intensity[i] = Solve<Diffraction>::maximize(psTT, 1e-8, initialTwoTheta, twoThetaStep, _twoTheta[i]);		
		_intensity[i] = Num<double>::integrate<Diffraction>(ttFun, _twoTheta[i] - 0.25, _twoTheta[i] + 0.25);
	}
}



/* void Diffraction::print(const Word& file, bool broaden) const
 *
 * Print diffraction data
 */

void Diffraction::print(const Word& file, bool broaden) const
{
	
	// Open file for writing if needed
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file != "stdout")
		Output::setStream(Output::addStream(file));
	
	// Set output method
	Output::method(STANDARD);
	
	// If printing to file, then print settings
	Output message;
	if (file != "stdout")
	{
		Output::newline();
		Output::print("Wavelength ");
		Output::print(_wavelength);
		Output::newline();
		Output::print("FWHM ");
		Output::printSci(_fwhm);
		Output::newline();
		Output::print("Resolution ");
		Output::print(_resolution);
	}
	
	// If printing to screen then add header
	else
	{
		message.addLine();
		message.add("Two-theta");
		message.add("Intensity");
		message.addLine();
		message.add("---------");
		message.add("---------");
	}
	
	// Add peaks
	if (!broaden)
	{
		message.addLines(_twoTheta.length());
		for (int i = 0; i < _twoTheta.length(); ++i)
		{
			message.addLine();
			message.add(_twoTheta[i], 10);
			message.add(_intensity[i], 10);
		}
	}
	
	// Add broadened data
	else
	{
		int i;
		double intensity;
		double twoTheta = _minTwoTheta - 3 * _fwhm;
		message.addLines((int)(_maxTwoTheta - _minTwoTheta)/_resolution + 1);
		for (twoTheta = (twoTheta < 3) ? 3 : twoTheta; twoTheta <= _maxTwoTheta + 3 * _fwhm; twoTheta += _resolution)
		{
			intensity = 0;
			for (i = 0; i < _twoTheta.length(); ++i)
				intensity += _intensity[i] * exp(-pow(twoTheta - _twoTheta[i], 2.0) / (2*_variance));
			message.addLine();
			message.add(twoTheta, 10);
			message.add(intensity, 10);
		}
	}
	
	// Print peaks
	Output::newline();
	Output::print(message, RIGHT);
	
	// Reset output
	if (file != "stdout")
		Output::removeStream(Output::streamID());
	Output::setStream(origStream);
	Output::method(origMethod);
}



/* void Diffraction::sortByTwoTheta(int left, int right)
 *
 * Sort diffraction data by two-theta values
 */

void Diffraction::sortByTwoTheta(int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	double pivot = _twoTheta[pivotIndex];
	
	// Move pivot to end
	Num<double>::swap(_twoTheta[pivotIndex], _twoTheta[right]);
	Num<double>::swap(_intensity[pivotIndex], _intensity[right]);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (_twoTheta[i] < pivot)
		{
			Num<double>::swap(_twoTheta[i], _twoTheta[newPivotIndex]);
			Num<double>::swap(_intensity[i], _intensity[newPivotIndex]);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	Num<double>::swap(_twoTheta[newPivotIndex], _twoTheta[right]);
	Num<double>::swap(_intensity[newPivotIndex], _intensity[right]);
	
	// Recursive calls to next sorts
	sortByTwoTheta(left, newPivotIndex - 1);
	sortByTwoTheta(newPivotIndex + 1, right);
}



/* void Diffraction::setATFParams(const Symmetry& symmetry)
 *
 * Set the parameters used for atomic form factor
 */

void Diffraction::setATFParams(const Symmetry& symmetry)
{
	
	// Allocate space
	_atfParams.length(symmetry.orbits().length());
	
	// Get parameters
	double a1, a2, a3, a4, b1, b2, b3, b4, c;
	for (int i = 0; i < symmetry.orbits().length(); ++i)
	{
		const Element& element = symmetry.orbits()[i].atoms()[0]->element();
		if (element.number() == 1)
		{
			a1 = 0.489918; b1 = 20.659300; a2 = 0.262003; b2 = 7.740390; a3 = 0.196767; b3 = 49.551899; 
			a4 = 0.049879; b4 = 2.201590; c = 0.001305;
		}
		else if (element.number() == 2)
		{
			a1 = 0.873400; b1 = 9.103700; a2 = 0.630900; b2 = 3.356800; a3 = 0.311200; b3 = 22.927601;
			a4 = 0.178000; b4 = 0.982100; c = 0.006400;
		}
		else if (element.number() == 3)
		{
			a1 = 1.128200; b1 = 3.954600; a2 = 0.750800; b2 = 1.052400; a3 = 0.617500; b3 = 85.390503;
			a4 = 0.465300; b4 = 168.261002; c = 0.037700;
		}
		else if (element.number() == 4)
		{
			a1 = 1.591900; b1 = 43.642700; a2 = 1.127800; b2 = 1.862300; a3 = 0.539100; b3 = 103.483002;
			a4 = 0.702900; b4 = 0.542000; c = 0.038500;
		}
		else if (element.number() == 5)
		{
			a1 = 2.054500; b1 = 23.218500; a2 = 1.332600; b2 = 1.021000; a3 = 1.097900; b3 = 60.349800;
			a4 = 0.706800; b4 = 0.140300; c = -0.193200;
		}
		else if (element.number() == 6)
		{
			a1 = 2.310000; b1 = 20.843901; a2 = 1.020000; b2 = 10.207500; a3 = 1.588600; b3 = 0.568700;
			a4 = 0.865000; b4 = 51.651199; c = 0.215600;
		}
		else if (element.number() == 7)
		{
			a1 = 12.212600; b1 = 0.005700; a2 = 3.132200; b2 = 9.893300; a3 = 2.012500; b3 = 28.997499;
			a4 = 1.166300; b4 = 0.582600; c = -11.529000;
		}
		else if (element.number() == 8)
		{
			a1 = 3.048500; b1 = 13.277100; a2 = 2.286800; b2 = 5.701100; a3 = 1.546300; b3 = 0.323900;
			a4 = 0.867000; b4 = 32.908901; c = 0.250800;
		}
		else if (element.number() == 9)
		{
			a1 = 3.539200; b1 = 10.282500; a2 = 2.641200; b2 = 4.294400; a3 = 1.517000; b3 = 0.261500;
			a4 = 1.024300; b4 = 26.147600; c = 0.277600;
		}
		else if (element.number() == 10)
		{
			a1 = 3.955300; b1 = 8.404200; a2 = 3.112500; b2 = 3.426200; a3 = 1.454600; b3 = 0.230600;
			a4 = 1.125100; b4 = 21.718399; c = 0.351500;
		}
		else if (element.number() == 11)
		{
			a1 = 4.762600; b1 = 3.285000; a2 = 3.173600; b2 = 8.842200; a3 = 1.267400; b3 = 0.313600;
			a4 = 1.112800; b4 = 129.423996; c = 0.676000;
		}
		else if (element.number() == 12)
		{
			a1 = 5.420400; b1 = 2.827500; a2 = 2.173500; b2 = 79.261101; a3 = 1.226900; b3 = 0.380800;
			a4 = 2.307300; b4 = 7.193700; c = 0.858400;
		}
		else if (element.number() == 13)
		{
			a1 = 6.420200; b1 = 3.038700; a2 = 1.900200; b2 = 0.742600; a3 = 1.593600; b3 = 31.547199;
			a4 = 1.964600; b4 = 85.088600; c = 1.115100;
		}
		else if (element.number() == 14)
		{
			a1 = 6.291500; b1 = 2.438600; a2 = 3.035300; b2 = 32.333698; a3 = 1.989100; b3 = 0.678500; 
			a4 = 1.541000; b4 = 81.693703; c = 1.140700;
		}
		else if (element.number() == 15)
		{
			a1 = 6.434500; b1 = 1.906700; a2 = 4.179100; b2 = 27.157000; a3 = 1.780000; b3 = 0.526000;
			a4 = 1.490800; b4 = 68.164497; c = 1.114900;
		}
		else if (element.number() == 16)
		{
			a1 = 6.905300; b1 = 1.467900; a2 = 5.203400; b2 = 22.215099; a3 = 1.437900; b3 = 0.253600;
			a4 = 1.586300; b4 = 56.172001; c = 0.866900;
		}
		else if (element.number() == 17)
		{
			a1 = 11.460400; b1 = 0.010400; a2 = 7.196400; b2 = 1.166200; a3 = 6.255600; b3 = 18.519400;
			a4 = 1.645500; b4 = 47.778400; c = -9.557400;
		}
		else if (element.number() == 18)
		{
			a1 = 7.484500; b1 = 0.907200; a2 = 6.772300; b2 = 14.840700; a3 = 0.653900; b3 = 43.898300;
			a4 = 1.644200; b4 = 33.392899; c = 1.444500;
		}
		else if (element.number() == 19)
		{
			a1 = 8.218600; b1 = 12.794900; a2 = 7.439800; b2 = 0.774800; a3 = 1.051900; b3 = 213.186996;
			a4 = 0.865900; b4 = 41.684101; c = 1.422800;
		}
		else if (element.number() == 20)
		{
			a1 = 8.626600; b1 = 10.442100; a2 = 7.387300; b2 = 0.659900; a3 = 1.589900; b3 = 85.748398;
			a4 = 1.021100; b4 = 178.436996; c = 1.375100;
		}
		else if (element.number() == 21)
		{
			a1 = 9.189000; b1 = 9.021300; a2 = 7.367900; b2 = 0.572900; a3 = 1.640900; b3 = 136.108002;
			a4 = 1.468000; b4 = 51.353100; c = 1.332900;
		}
		else if (element.number() == 22)
		{
			a1 = 9.759500; b1 = 7.850800; a2 = 7.355800; b2 = 0.500000; a3 = 1.699100; b3 = 35.633801;
			a4 = 1.902100; b4 = 116.105003; c = 1.280700;
		}
		else if (element.number() == 23)
		{
			a1 = 10.297100; b1 = 6.865700; a2 = 7.351100; b2 = 0.438500; a3 = 2.070300; b3 = 26.893801;
			a4 = 2.057100; b4 = 102.477997; c = 1.219900;
		}
		else if (element.number() == 24)
		{
			a1 = 10.640600; b1 = 6.103800; a2 = 7.353700; b2 = 0.392000; a3 = 3.324000; b3 = 20.262600;
			a4 = 1.492200; b4 = 98.739899; c = 1.183200;
		}
		else if (element.number() == 25)
		{
			a1 = 11.281900; b1 = 5.340900; a2 = 7.357300; b2 = 0.343200; a3 = 3.019300; b3 = 17.867399;
			a4 = 2.244100; b4 = 83.754303; c = 1.089600;
		}
		else if (element.number() == 26)
		{
			a1 = 11.769500; b1 = 4.761100; a2 = 7.357300; b2 = 0.307200; a3 = 3.522200; b3 = 15.353500;
			a4 = 2.304500; b4 = 76.880501; c = 1.036900;
		}
		else if (element.number() == 27)
		{
			a1 = 12.284100; b1 = 4.279100; a2 = 7.340900; b2 = 0.278400; a3 = 4.003400; b3 = 13.535900;
			a4 = 2.348800; b4 = 71.169197; c = 1.011800;
		}
		else if (element.number() == 28)
		{
			a1 = 12.837600; b1 = 3.878500; a2 = 7.292000; b2 = 0.256500; a3 = 4.443800; b3 = 12.176300;
			a4 = 2.380000; b4 = 66.342102; c = 1.034100;
		}
		else if (element.number() == 29)
		{
			a1 = 13.338000; b1 = 3.582800; a2 = 7.167600; b2 = 0.247000; a3 = 5.615800; b3 = 11.396600;
			a4 = 1.673500; b4 = 64.812599; c = 1.191000;
		}
		else if (element.number() == 30)
		{
			a1 = 14.074300; b1 = 3.265500; a2 = 7.031800; b2 = 0.233300; a3 = 5.165200; b3 = 10.316300;
			a4 = 2.410000; b4 = 58.709702; c = 1.304100;
		}
		else if (element.number() == 31)
		{
			a1 = 15.235400; b1 = 3.066900; a2 = 6.700600; b2 = 0.241200; a3 = 4.359100; b3 = 10.780500;
			a4 = 2.962300; b4 = 61.413502; c = 1.718900;
		}
		else if (element.number() == 32)
		{
			a1 = 16.081600; b1 = 2.850900; a2 = 6.374700; b2 = 0.251600; a3 = 3.706800; b3 = 11.446800;
			a4 = 3.683000; b4 = 54.762501; c = 2.131300;
		}
		else if (element.number() == 33)
		{
			a1 = 10.672300; b1 = 2.634500; a2 = 6.070100; b2 = 0.264700; a3 = 3.431300; b3 = 12.947900;
			a4 = 4.277900; b4 = 47.797199; c = 2.531000;
		}
		else if (element.number() == 34)
		{
			a1 = 17.000601; b1 = 2.409800; a2 = 5.819600; b2 = 0.272600; a3 = 3.973100; b3 = 15.237200;
			a4 = 4.354300; b4 = 43.816299; c = 2.840900;
		}
		else if (element.number() == 35)
		{
			a1 = 17.178900; b1 = 2.172300; a2 = 5.235800; b2 = 16.579599; a3 = 5.637700; b3 = 0.260900;
			a4 = 3.985100; b4 = 41.432800; c = 2.955700;
		}
		else if (element.number() == 36)
		{
			a1 = 17.355499; b1 = 1.938400; a2 = 6.728600; b2 = 16.562300; a3 = 5.549300; b3 = 0.226100;
			a4 = 3.537500; b4 = 39.397202; c = 2.825000;
		}
		else if (element.number() == 37)
		{
			a1 = 17.178400; b1 = 1.788800; a2 = 9.643500; b2 = 17.315100; a3 = 5.139900; b3 = 0.274800;
			a4 = 1.529200; b4 = 164.934006; c = 3.487300;
		}
		else if (element.number() == 38)
		{
			a1 = 17.566299; b1 = 1.556400; a2 = 9.818400; b2 = 14.098800; a3 = 5.422000; b3 = 0.166400;
			a4 = 2.669400; b4 = 132.376007; c = 2.506400;
		}
		else if (element.number() == 39)
		{
			a1 = 17.775999; b1 = 1.402900; a2 = 10.294600; b2 = 12.800600; a3 = 5.726290; b3 = 0.125599;
			a4 = 3.265880; b4 = 104.353996; c = 1.912130;
		}
		else if (element.number() == 40)
		{
			a1 = 17.876499; b1 = 1.276180; a2 = 10.948000; b2 = 11.916000; a3 = 5.417320; b3 = 0.117622;
			a4 = 3.657210; b4 = 87.662697; c = 2.069290;
		}
		else if (element.number() == 41)
		{
			a1 = 17.614201; b1 = 1.188650; a2 = 12.014400; b2 = 11.766000; a3 = 4.041830; b3 = 0.204785;
			a4 = 3.533460; b4 = 69.795700; c = 3.755910;
		}
		else if (element.number() == 42)
		{
			a1 = 3.702500; b1 = 0.277200; a2 = 17.235600; b2 = 1.095800; a3 = 12.887600; b3 = 11.004000;
			a4 = 3.742900; b4 = 61.658401; c = 4.387500;
		}
		else if (element.number() == 43)
		{
			a1 = 19.130100; b1 = 0.864132; a2 = 11.094800; b2 = 8.144870; a3 = 4.649010; b3 = 21.570700;
			a4 = 2.712630; b4 = 86.847198; c = 5.404280;
		}
		else if (element.number() == 44)
		{
			a1 = 19.267401; b1 = 0.808520; a2 = 12.918200; b2 = 8.434670; a3 = 4.863370; b3 = 24.799700;
			a4 = 1.567560; b4 = 94.292801; c = 5.378740;
		}
		else if (element.number() == 45)
		{
			a1 = 19.295700; b1 = 0.751536; a2 = 14.350100; b2 = 8.217580; a3 = 4.734250; b3 = 25.874901;
			a4 = 1.289180; b4 = 98.606201; c = 5.328000;
		}
		else if (element.number() == 46)
		{
			a1 = 19.331900; b1 = 0.698655; a2 = 15.501700; b2 = 7.989290; a3 = 5.295370; b3 = 25.205200;
			a4 = 0.605844; b4 = 76.898598; c = 5.265930;
		}
		else if (element.number() == 47)
		{
			a1 = 19.280800; b1 = 0.644600; a2 = 16.688499; b2 = 7.472600; a3 = 4.804500; b3 = 24.660500;
			a4 = 1.046300; b4 = 99.815598; c = 5.179000;
		}
		else if (element.number() == 48)
		{
			a1 = 19.221399; b1 = 0.594600; a2 = 17.644400; b2 = 6.908900; a3 = 4.461000; b3 = 24.700800;
			a4 = 1.602900; b4 = 87.482498; c = 5.069400;
		}
		else if (element.number() == 49)
		{
			a1 = 19.162399; b1 = 0.547600; a2 = 18.559601; b2 = 6.377600; a3 = 4.294800; b3 = 25.849899;
			a4 = 2.039600; b4 = 92.802902; c = 4.939100;
		}
		else if (element.number() == 50)
		{
			a1 = 19.188900; b1 = 5.830300; a2 = 19.100500; b2 = 0.503100; a3 = 4.458500; b3 = 26.890900;
			a4 = 2.466300; b4 = 83.957100; c = 4.782100;
		}
		else if (element.number() == 51)
		{
			a1 = 19.641800; b1 = 5.303400; a2 = 19.045500; b2 = 0.460700; a3 = 5.037100; b3 = 27.907400;
			a4 = 2.682700; b4 = 75.282501; c = 4.590900;
		}
		else if (element.number() == 52)
		{
			a1 = 19.964399; b1 = 4.817420; a2 = 19.013800; b2 = 0.420885; a3 = 6.144870; b3 = 28.528400;
			a4 = 2.523900; b4 = 70.840302; c = 4.352000;
		}
		else if (element.number() == 53)
		{
			a1 = 20.147200; b1 = 4.347000; a2 = 18.994900; b2 = 0.381400; a3 = 7.513800; b3 = 27.766001;
			a4 = 2.273500; b4 = 66.877602; c = 4.071200;
		}
		else if (element.number() == 54)
		{
			a1 = 20.293301; b1 = 3.928200; a2 = 19.029800; b2 = 0.344000; a3 = 8.976700; b3 = 26.465900;
			a4 = 1.990000; b4 = 64.265800; c = 3.711800;
		}
		else if (element.number() == 55)
		{
			a1 = 20.389200; b1 = 3.569000; a2 = 19.106199; b2 = 0.310700; a3 = 10.662000; b3 = 24.387899;
			a4 = 1.495300; b4 = 213.904007; c = 3.335200;
		}
		else if (element.number() == 56)
		{
			a1 = 20.336100; b1 = 3.216000; a2 = 19.297001; b2 = 0.275600; a3 = 10.888000; b3 = 20.207300;
			a4 = 2.695900; b4 = 167.201996; c = 2.773100;
		}
		else if (element.number() == 57)
		{
			a1 = 20.577999; b1 = 2.948170; a2 = 19.599001; b2 = 0.244475; a3 = 11.372700; b3 = 18.772600;
			a4 = 3.287190; b4 = 133.123993; c = 2.146780;
		}
		else if (element.number() == 58)
		{
			a1 = 21.167101; b1 = 2.812190; a2 = 19.769501; b2 = 0.226836; a3 = 11.851300; b3 = 17.608299;
			a4 = 3.330490; b4 = 127.112999; c = 1.862640;
		}
		else if (element.number() == 59)
		{
			a1 = 22.044001; b1 = 2.773930; a2 = 19.669701; b2 = 0.222087; a3 = 12.385600; b3 = 16.766899;
			a4 = 2.824280; b4 = 143.643997; c = 2.058300;
		}
		else if (element.number() == 60)
		{
			a1 = 22.684500; b1 = 2.662480; a2 = 19.684700; b2 = 0.210628; a3 = 12.774000; b3 = 15.885000;
			a4 = 2.851370; b4 = 137.903000; c = 1.984860;
		}
		else if (element.number() == 61)
		{
			a1 = 23.340500; b1 = 2.562700; a2 = 19.609501; b2 = 0.202088; a3 = 13.123500; b3 = 15.100900;
			a4 = 2.875160; b4 = 132.720993; c = 2.028760;
		}
		else if (element.number() == 62)
		{
			a1 = 24.004200; b1 = 2.472740; a2 = 19.425800; b2 = 0.196451; a3 = 13.439600; b3 = 14.399600;
			a4 = 2.896040; b4 = 128.007004; c = 2.209630;
		}
		else if (element.number() == 63)
		{
			a1 = 24.627399; b1 = 2.387900; a2 = 19.088600; b2 = 0.194200; a3 = 13.760300; b3 = 13.754600;
			a4 = 2.922700; b4 = 123.174004; c = 2.574500;
		}
		else if (element.number() == 64)
		{
			a1 = 25.070900; b1 = 2.253410; a2 = 19.079800; b2 = 0.181951; a3 = 13.851800; b3 = 12.933100;
			a4 = 3.545450; b4 = 101.398003; c = 2.419600;
		}
		else if (element.number() == 65)
		{
			a1 = 25.897600; b1 = 2.242560; a2 = 18.218500; b2 = 0.196143; a3 = 14.316700; b3 = 12.664800;
			a4 = 2.953540; b4 = 115.362000; c = 3.589240;
		}
		else if (element.number() == 66)
		{
			a1 = 26.507000; b1 = 2.180200; a2 = 17.638300; b2 = 0.202172; a3 = 14.559600; b3 = 12.189900;
			a4 = 2.965770; b4 = 111.874001; c = 4.297280;
		}
		else if (element.number() == 67)
		{
			a1 = 26.904900; b1 = 2.070510; a2 = 17.294001; b2 = 0.197940; a3 = 14.558300; b3 = 11.440700;
			a4 = 3.638370; b4 = 92.656601; c = 4.567960;
		}
		else if (element.number() == 68)
		{
			a1 = 27.656300; b1 = 2.073560; a2 = 16.428499; b2 = 0.223545; a3 = 14.977900; b3 = 11.360400;
			a4 = 2.982330; b4 = 105.703003; c = 5.920460;
		}
		else if (element.number() == 69)
		{
			a1 = 28.181900; b1 = 2.028590; a2 = 15.885100; b2 = 0.238849; a3 = 15.154200; b3 = 10.997500;
			a4 = 2.987060; b4 = 102.960999; c = 6.756210;
		}
		else if (element.number() == 70)
		{
			a1 = 28.664101; b1 = 1.988900; a2 = 15.434500; b2 = 0.257119; a3 = 15.308700; b3 = 10.664700;
			a4 = 2.989630; b4 = 100.417000; c = 7.566720;
		}
		else if (element.number() == 71)
		{
			a1 = 28.947599; b1 = 1.901820; a2 = 15.220800; b2 = 9.985190; a3 = 15.100000; b3 = 0.261033;
			a4 = 3.716010; b4 = 84.329803; c = 7.976280;
		}
		else if (element.number() == 72)
		{
			a1 = 29.143999; b1 = 1.832620; a2 = 15.172600; b2 = 9.599900; a3 = 14.758600; b3 = 0.275116;
			a4 = 4.300130; b4 = 72.028999; c = 8.581540;
		}
		else if (element.number() == 73)
		{
			a1 = 29.202400; b1 = 1.773330; a2 = 15.229300; b2 = 9.370460; a3 = 14.513500; b3 = 0.295977;
			a4 = 4.764920; b4 = 63.364399; c = 9.243540;
		}
		else if (element.number() == 74)
		{
			a1 = 29.081800; b1 = 1.720290; a2 = 15.430000; b2 = 9.225900; a3 = 14.432700; b3 = 0.321703;
			a4 = 5.119820; b4 = 57.056000; c = 9.887500;
		}
		else if (element.number() == 75)
		{
			a1 = 28.762100; b1 = 1.671910; a2 = 15.718900; b2 = 9.092270; a3 = 14.556400; b3 = 0.350500;
			a4 = 5.441740; b4 = 52.086102; c = 10.472000;
		}
		else if (element.number() == 76)
		{
			a1 = 28.189400; b1 = 1.629030; a2 = 16.155001; b2 = 8.979480; a3 = 14.930500; b3 = 0.382661;
			a4 = 5.675890; b4 = 48.164700; c = 11.000500;
		}
		else if (element.number() == 77)
		{
			a1 = 27.304899; b1 = 1.592790; a2 = 16.729601; b2 = 8.865530; a3 = 15.611500; b3 = 0.417916;
			a4 = 5.833770; b4 = 45.001099; c = 11.472200;
		}
		else if (element.number() == 78)
		{
			a1 = 27.005899; b1 = 1.512930; a2 = 17.763901; b2 = 8.811740; a3 = 15.713100; b3 = 0.424593;
			a4 = 5.783700; b4 = 38.610298; c = 11.688300;
		}
		else if (element.number() == 79)
		{
			a1 = 16.881901; b1 = 0.461100; a2 = 18.591299; b2 = 8.621600; a3 = 25.558201; b3 = 1.482600;
			a4 = 5.860000; b4 = 36.395599; c = 12.065800;
		}
		else if (element.number() == 80)
		{
			a1 = 20.680901; b1 = 0.545000; a2 = 19.041700; b2 = 8.448400; a3 = 21.657499; b3 = 1.572900;
			a4 = 5.967600; b4 = 38.324600; c = 12.608900;
		}
		else if (element.number() == 81)
		{
			a1 = 27.544600; b1 = 0.655150; a2 = 19.158400; b2 = 8.707510; a3 = 15.538000; b3 = 1.963470;
			a4 = 5.525930; b4 = 45.814899; c = 13.174600;
		}
		else if (element.number() == 82)
		{
			a1 = 31.061701; b1 = 0.690200; a2 = 13.063700; b2 = 2.357600; a3 = 18.441999; b3 = 8.618000;
			a4 = 5.969600; b4 = 47.257900; c = 13.411800;
		}
		else if (element.number() == 83)
		{
			a1 = 33.368900; b1 = 0.704000; a2 = 12.951000; b2 = 2.923800; a3 = 16.587700; b3 = 8.793700;
			a4 = 6.469200; b4 = 48.009300; c = 13.578200;
		}
		else if (element.number() == 84)
		{
			a1 = 34.672600; b1 = 0.700999; a2 = 15.473300; b2 = 3.550780; a3 = 13.113800; b3 = 9.556420;
			a4 = 7.025800; b4 = 47.004501; c = 13.677000;
		}
		else if (element.number() == 85)
		{
			a1 = 35.316299; b1 = 0.685870; a2 = 19.021099; b2 = 3.974580; a3 = 9.498870; b3 = 11.382400;
			a4 = 7.425180; b4 = 45.471500; c = 13.710800;
		}
		else if (element.number() == 86)
		{
			a1 = 35.563099; b1 = 0.663100; a2 = 21.281601; b2 = 4.069100; a3 = 8.003700; b3 = 14.042200;
			a4 = 7.443300; b4 = 44.247299; c = 13.690500;
		}
		else if (element.number() == 87)
		{
			a1 = 35.929901; b1 = 0.646453; a2 = 23.054701; b2 = 4.176190; a3 = 12.143900; b3 = 23.105200;
			a4 = 2.112530; b4 = 150.645004; c = 13.724700;
		}
		else if (element.number() == 88)
		{
			a1 = 35.763000; b1 = 0.616341; a2 = 22.906401; b2 = 3.871350; a3 = 12.473900; b3 = 19.988701;
			a4 = 3.210970; b4 = 142.324997; c = 13.621100;
		}
		else if (element.number() == 89)
		{
			a1 = 35.659698; b1 = 0.589092; a2 = 23.103201; b2 = 3.651550; a3 = 12.597700; b3 = 18.599001;
			a4 = 4.086550; b4 = 117.019997; c = 13.526600;
		}
		else if (element.number() == 90)
		{
			a1 = 35.564499; b1 = 0.563359; a2 = 23.421900; b2 = 3.462040; a3 = 12.747300; b3 = 17.830900;
			a4 = 4.807030; b4 = 99.172203; c = 13.431400;
		}
		else if (element.number() == 91)
		{
			a1 = 35.884701; b1 = 0.547751; a2 = 23.294800; b2 = 3.415190; a3 = 14.189100; b3 = 16.923500;
			a4 = 4.172870; b4 = 105.250999; c = 13.428700;
		}
		else if (element.number() == 92)
		{
			a1 = 36.022800; b1 = 0.529300; a2 = 23.412800; b2 = 3.325300; a3 = 14.949100; b3 = 16.092699;
			a4 = 4.188000; b4 = 100.612999; c = 13.396600;
		}
		else if (element.number() == 93)
		{
			a1 = 36.187401; b1 = 0.511929; a2 = 23.596399; b2 = 3.253960; a3 = 15.640200; b3 = 15.362200;
			a4 = 4.185500; b4 = 97.490799; c = 13.357300;
		}
		else if (element.number() == 94)
		{
			a1 = 36.525398; b1 = 0.499384; a2 = 23.808300; b2 = 3.263710; a3 = 16.770700; b3 = 14.945500;
			a4 = 3.479470; b4 = 105.980003; c = 13.381200;
		}
		else if (element.number() == 95)
		{
			a1 = 36.670601; b1 = 0.483629; a2 = 24.099199; b2 = 3.206470; a3 = 17.341499; b3 = 14.313600;
			a4 = 3.493310; b4 = 102.273003; c = 13.359200;
		}
		else if (element.number() == 96)
		{
			a1 = 36.648800; b1 = 0.465154; a2 = 24.409599; b2 = 3.089970; a3 = 17.399000; b3 = 13.434600;
			a4 = 4.216650; b4 = 88.483398; c = 13.288700;
		}
		else if (element.number() == 97)
		{
			a1 = 36.788101; b1 = 0.451018; a2 = 24.773600; b2 = 3.046190; a3 = 17.891899; b3 = 12.894600;
			a4 = 4.232840; b4 = 86.002998; c = 13.275400;
		}
		else if (element.number() == 98)
		{
			a1 = 36.918499; b1 = 0.437533; a2 = 25.199499; b2 = 3.007750; a3 = 18.331699; b3 = 12.404400;
			a4 = 4.243910; b4 = 83.788101; c = 13.267400;
		}
		else
		{
			Output::newline(ERROR);
			Output::print("Atomic scattering factor is not defined for ");
			Output::print(element.symbol());
			Output::quit();
		}
		_atfParams[i].length(9);
		_atfParams[i][0] = a1;
		_atfParams[i][1] = a2;
		_atfParams[i][2] = a3;
		_atfParams[i][3] = a4;
		_atfParams[i][4] = b1;
		_atfParams[i][5] = b2;
		_atfParams[i][6] = b3;
		_atfParams[i][7] = b4;
		_atfParams[i][8] = c;
	}
}



/* double Diffraction::atomicScatteringFactor(int index, double angle)
 *
 * Return the atomic scattering factor
 */

double Diffraction::atomicScatteringFactor(int index, double angle)
{
	
	// Get s value
	double s = sin(angle) / _wavelength;
	double s2 = s * s;
	if (s > 2)
	{
		Output::newline(WARNING);
		Output::print("Atomic scattering factor is not optimized for s greater than 2");
	}
	
	// Return result
	return _atfParams[index][0]*exp(-_atfParams[index][4]*s2) + \
		_atfParams[index][1]*exp(-_atfParams[index][5]*s2) + \
		_atfParams[index][2]*exp(-_atfParams[index][6]*s2) + \
		_atfParams[index][3]*exp(-_atfParams[index][7]*s2) + _atfParams[index][8];
}
